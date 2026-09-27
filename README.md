# ohc_derive

`ohc_derive` builds combined-level ("synthetic-level") ocean-heat-content analysis quantities — OHCA, OHU, their trends, and gridded anomaly maps — from native-level ME4OH submissions and a standard bathymetry, writing one NetCDF per synthetic level. It works from the ME4OH-shaped submissions the ingest `publish` step writes (any quantity), or from any group's bare ME4OH-protocol files with `--contract ME4OH`, plus a standard bathy.

## What it computes

**Inputs** are the native-level submission NetCDFs that make up a synthetic level — its *constituents* — plus a standard bathymetry on the common grid. A submission is an ME4OH-shaped file: `DATA(LONGITUDE, LATITUDE, TIME)` on the common 1° grid, monthly, missing as NaN. `cell_area` is regenerated from the grid.

### What a submission must carry

On top of the ME4OH protocol itself, derive reads four things from a submission:

| what | where | used for |
|---|---|---|
| the **quantity table** | `quantity` attr (the ingest `[quantity]` table as compact JSON: `name`, `kind`, `units`, `scale_terms`, published units) | what the field is, its kind (extensive/intensive), its units; `cp0`/`rho0` for the gcos emitter |
| the **native-level tag** | `mapped_layer` attr, `<top>_<bottom>` | keying the constituent the level plan asks for |
| the **ensemble** | a `<NAME>ENS_…` sibling file, `DATA(MEMBER, LONGITUDE, LATITUDE, TIME)`, found by swapping the `<NAME>_` filename prefix from `quantity.name` | the `_sd` companions |
| the **provenance chain** | `*_run_config` / `*_run_facts` / `*_code_version` attrs | rolled forward into the blob (optional — absent blocks are simply not forwarded) |

**Files from the LocalGP ingest stage carry all of this.** `ohc_ingest`'s `publish.py` stamps the `quantity` table, `mapped_layer` and the provenance blocks, and writes the `<NAME>ENS_` sibling with `--ensemble`; nothing further is needed, and this is the default.

**Files that stop at the ME4OH protocol — other groups' submissions — carry none of it, and don't have to.** Pass `--contract ME4OH` and derive fills the gaps from the protocol, which fixes them: `DATA` is ocean heat content density in TJ/m² computed with the protocol's `cp0 = 3989.244 J/kg/K` and `rho0 = 1030 kg/m³` (so the quantity table is the protocol's OHC table); the layer is the filename's `lev<low>_<high>` token; there is no ensemble, so the run is mean-only (a `DATA_SD`, if present, is not used). Everything inferred is recorded on the blob in `ohc_derive_run_facts.inferred_config`, keyed by constituent, so the record says the table came from the specification and not from the file. Use `contiguous_from_top` (the default) for such files, not `as_published`: the latter verifies a footprint our publish step enforced, which a third party's file never went through. A bare protocol file passed *without* `--contract` is refused with a message saying so — derive never infers silently. The table in [`contracts.py`](contracts.py) is the one place the protocol's assumptions live.

A **synthetic level** is an `n_fac`-weighted sum of native ME4OH levels, shallowest first — for example `0_2000` is `15_20`(×3) + `15_300` + `300_700` + `700_1850` + `1800_1850`(×3). `n_fac` scales a thin measured layer up to the slab it stands in for; each constituent carries its own dbar `top`/`bottom`, used against the bathy in the mask. The level table is a **plan file** passed with `--levels` (TOML, one `[[level]]` per synthetic level with its `name`, `require_top` and `contributors`); LocalGP's is [`levels/localgp.toml`](levels/localgp.toml): `0_300`, `0_700`, `0_1000`, `700_2000`, `0_2000`. A group whose native layers decompose the target layers differently writes its own plan — "`0_2000` is one `0_700` plus one `700_2000`" — and nothing else changes. On load, derive checks each level: contributors shallowest first, no repeated tag, the `n_fac`-weighted thicknesses summing to the level's thickness (`high − low`), and `require_top` within it; a plan that doesn't tile its layer is refused, since it would combine silently into the wrong number. The plan used is recorded in `ohc_derive_run_facts` (`levels_file`, `level_plan`).

The submissions' `quantity` table names the field's **kind**: `extensive` (a per-area density, like OHC in J/m², that sums over area and stacks over layers) or `intensive` (a per-cell value, like a mixed layer depth, that does neither). Each primitive and mask prescription declares the kinds it applies to (`map_transforms.KINDS`, `masks.KINDS`), and a run is refused up front if its plan doesn't suit the field: an intensive field can't take the `integral`-based quantities (that would need an area mean, not written), can't use `fully_wet_nan` (its 0-fill means "adds nothing", an extensive idea), and can't be combined across several constituents (a thickness-weighted mean, not written). The error names each mismatch.

One run builds one synthetic level (`--level`), so levels parallelize across jobs. A level is a plan for combining constituents, and one constituent needs no plan: omit `--level` with exactly one submission in the pool and the run builds that submission's **identity level** — its native tag as the level, one contributor, `n_fac = 1` — which is how a single-layer quantity such as a mixed layer depth goes through the factory (the combine returns the constituent). It proceeds in six steps:

1. **load** — read the constituents' submissions (mean field + members) and the standard bathy.
2. **mask** — apply the cross-layer mask prescription over the constituents.
3. **primitives** — reduce each constituent to its map-level primitives: the area-weighted `integral` (a time series) and the gridded field (`map`).
4. **build** — compose the primitives into the requested deliverables.
5. **collapse** — per constituent, take the central value from the mean field and the 1σ spread across the members.
6. **combine** — `n_fac` sum the constituents into the synthetic level: values linearly, standard deviations worst-case.

The mean field and its members ride together on a leading `realization` axis (index 0 is the mean field, the rest are members). Steps 3–4 transform every realization the same way, and step 5 reads the central value off index 0 and the spread off the members. So an anomaly demean along time hits every realization, and each member is referenced to its own window mean.

### Masking

Cross-layer masking is pluggable (`masks.REGISTRY`); each prescription returns the masked constituents, a footprint (the cells counted for area), and a per-cell column height (metres). `apply` turns those into the level's `area_m2` (footprint cell area) and `volume_m3` (the cell-area-weighted sum of the height), so the volume tapers with the kept column rather than assuming a slab. Alongside, it writes `mask_<tag>_<data>_tw<baseline>_<level>_<mask>_<product_name>_<author>.png` (the footprint) and `coverage_<tag>_<data>_tw<baseline>_<level>_<mask>_<product_name>_<author>.nc` on the product grid (the `<data>_tw<baseline>` token matches the main `.nc`, so different-baseline runs don't overwrite each other's auxiliaries; `<product_name>_<author>` is the publication descriptor since these two files are published) — per cell, the `kept_thickness` (the height above) and the `uncaptured_thickness`, the in-layer water the kept column didn't reach: from the bottom of the kept column down to `min(bathy, layer bottom)`, floored at 0. The coverage `.nc` also carries the full `--citation` sentence as a top-level `citation` attr.

`contiguous_from_top` is the default prescription. `require_top` is a **thickness of the layer's own top** — the metres below `level.low` that must be present — and each level carries its own in the plan file (the top 300 m: every LocalGP level is `300`; `--require-top` overrides it). A cell **survives** only where every constituent reaching into the layer's top `require_top` metres is defined at all times and members. Constituents are atomic (a whole LocalGP layer or nothing), so a `require_top` that lands partway into one requires that whole constituent. Selection is by the constituent's declared bounds; each constituent is a single already-integrated value per cell (present or NaN — no internal depth). With `require_top = 300`: `0_2000` requires `15_20` + `15_300`, whose bounds happen to tile exactly 0–300; `700_2000` requires `700_1850`, because `require_depth = 1000` falls within its 700–1850 bounds — so that one value must be present. The shallowest constituent's `n_fac` reaches the layer top, so any positive `require_top` requires it. Within a surviving cell, walk the constituents from the layer top down and **keep** them until the first undefined one, then discard it and everything below (set to NaN, so they vanish from the nan-aware integral). The kept run's `n_fac`-weighted thickness is the cell's height (a kept `1800_1850` adds 3·50 = 150 m), so the volume is the true tapering ocean.

`as_published` is for an identity level built from our own publish output, and is asked for by name. With one constituent there is nothing to coordinate across layers, and the mask publish applied is final: the data is untouched, the footprint is the cells finite at every time and member, and the prescription **verifies** that this is time-constant — a cell finite at some times or members but not others is a hard error (every stage downstream assumes a time-constant footprint; publish with `incomplete_timeseries` and `ensemble_incomplete` honored makes it so). It reports no column height, so the blob carries `area_m2` but no `volume_m3`, and no coverage `.nc` is written; the footprint png is.

`fully_wet_nan` is the earlier prescription (bathy-driven): each constituent is classified against the standard bathy as fully wet (floor ≥ `bottom`), intersecting, or dry; a cell drops where any fully-wet constituent is undefined, a constituent contributes where not dry and defined (else 0), and the height is the level's full nominal thickness everywhere in-footprint (a slab).

"Defined" is judged across every member and timestep together, so the footprint is identical for all realizations and the ensemble spread reflects real spread rather than footprint jitter.

### Quantities

`--quantities` selects deliverables from `temporal_transforms.REGISTRY`. `--time-window` (`YEAR0:YEAR1`) sets the anomaly baseline and the years the trends are fit over; omitted, it spans all years.

| `--quantities` name | dims | what it is |
|---|---|---|
| `ohca` | (year,) | monthly area-integrated anomaly (window baseline), annual-averaged |
| `ohu` | (year,) | month-to-month tendency of the integral, annual-averaged; NaN-seeded at t0 |
| `ohca_trend` | scalar | OLS slope of the annual integral over the window |
| `ohu_trend` | scalar | OLS slope of the annual tendency over the window |
| `map` | (time, lat, lon) | per-cell monthly anomaly (window baseline) |
| `field` | (time, lat, lon) | the masked field itself, per cell and month — no transform; its `_sd` is the per-cell member spread. Ignores the window |

The integral-based quantities are **extensive** — the per-area submission field integrated over the footprint (e.g. TJ from a TJ/m² field), with the trends and tendency carrying the matching per-year and per-month scaling; `map` stays in the submission's per-area unit. Each variable (and its `_sd`) says so in its attrs: `field_units` is the submission's published units (from the `quantity` table) and `reduction` is `area_integral` (summed over the footprint's cell areas — units are field units × m², and dividing by `area_m2` recovers a per-area density) or `grid` (still per cell). A variable that is a rate also carries `per`, the step it is a rate over: `ohu` is the annual mean of a month-to-month difference, so `per = "month"`; the two trends are slopes against the annual axis, so `per = "year"`. `ohca` and `map` are states and carry no `per`. Turning these into per-area target densities (OHCA in J/m², OHU in W/m², and so on) is the downstream packaging step's job, done from that metadata and the geometry below.

### Output

One NetCDF per level, `derive_<tag>_<data>_tw<baseline>_<level>.nc` — the token carries both year spans: `<data>` (`YYYY_YYYY`) is the data span actually present, and `tw<baseline>` (`twYYYY_YYYY`) is the `--time-window` baseline, defaulting to the whole data span (so a windowless run reads `…_2004_2025_tw2004_2025_…` and a 2005-2024 baseline `…_2004_2025_tw2005_2024_…`, and runs differing only in baseline don't collide): each requested quantity plus its `_sd` companion (omitted under `--no-ensemble`), with header attrs `level`, `area_m2` and `volume_m3` (both from the mask step — the tapered footprint area and column volume; `volume_m3` is absent under `as_published`, which reports no height), the constituents' `quantity` table (the ingest `[quantity]` attr, one compact-JSON attr — the gcos emitter reads `cp0`/`rho0` from its `scale_terms`), and `provenance_tag` / `provenance_link`. The factory emits these extensive quantities and geometry; packaging derives the intensive per-area densities from them.

**Provenance chain.** derive is a fan-in — one synthetic level is built from N constituent submissions — so it rolls each upstream provenance block forward *grouped by constituent*: for every `*_run_config` / `*_run_facts` / `*_code_version` on the inputs (e.g. `localgp_ingest_*`, `localgp_publish_*`), the output carries `<block> = {constituent_tag: block}` (compact JSON, one line). Blocks are opaque — parsed only to nest, never read — so nothing is deduplicated and there's no coupling to the upstream schema. On top of that, derive stamps its own `ohc_derive_run_config` (resolved args), `ohc_derive_run_facts` (level and whether it is an identity level, geometry, constituents, `n_fac`, axis, the `quantity` table), and `ohc_derive_code_version`. Each step namespaces its block by identity, so the chain accretes at every stage; the global `provenance_tag` / `provenance_link` are this derive run's own.

## Usage

### Environment

Described in `Dockerfile` to generate a containerized environment; make a similar env in anaconda on blanca when running on the cluster.

### Test

Tests live in `tests/` and run under pytest — in the containerized environment:

```
docker image build -t ohc_derive:test .
docker container run -v $(pwd):/app ohc_derive:test pytest
```

### Run

See `derive.slurm` for a real example of running this on blanca at CU. With the ensemble on (the default), each constituent's `<NAME>ENS_` sibling **must** sit next to its `<NAME>_` file or the loader exits, and every quantity gets a collapsed `_sd`. `--no-ensemble` is the central-only path (mean field, no `_sd`); central values are identical either way.

#### run.py options

All configuration is on the command line — no env, no config file. The available quantities are `temporal_transforms.REGISTRY`; the mask prescriptions are `masks.REGISTRY`.

| option | default | effect |
|---|---|---|
| `SUBMISSION.nc …` (positional) | *(required)* | the constituent submissions (`<NAME>_…`, or bare ME4OH files with `--contract ME4OH`); the `<NAME>ENS_` member siblings are found automatically. Each level selects the native constituents it needs by tag, so you can pass the whole pool of submissions and let each run pick — but the pool must hold **exactly one file per native level** (a duplicate tag, e.g. a stray window/experiment/rerun, is a hard error, not a silent last-wins). |
| `--level` | *(none)* | the synthetic level to build, by name in the `--levels` plan, e.g. `0_2000`. Omitted with exactly one submission: that submission's identity level. Omitted with any other number of submissions: an error. |
| `--levels` | *(none)* | the level plan (TOML; LocalGP's is `levels/localgp.toml`). Required whenever `--level` is given; an identity level needs no plan. |
| `--bathy` | *(required)* | standard bathymetry NetCDF on the common grid. |
| `--quantities` | *(required)* | comma list from `ohca,ohu,ohca_trend,ohu_trend,map,field`. Unknown names error. |
| `--mask` | `contiguous_from_top` | mask prescription (`masks.REGISTRY`): `contiguous_from_top`, `fully_wet_nan`, or `as_published`. Pass `as_published` by name for an identity level built from our own publish output. |
| `--require-top` | *(the level's own)* | metres of the layer's own top (from `level.low`) that must be defined for a cell to survive; overrides the level's `require_top` (in the plan file). Used by `contiguous_from_top`, ignored by `fully_wet_nan`. |
| `--time-window` | *(all years)* | `YEAR0:YEAR1` — the anomaly baseline and the trend-fit years. Separator `:`, `-`, or `_`, so the filename token `2004_2025` works verbatim. |
| `--no-ensemble` | off (ensemble **on**) | mean field only — skip the `_sd` companions and do not read the `<NAME>ENS_` siblings. |
| `--contract` | *(none)* | `ME4OH`: the submissions are bare protocol files (another group's); infer the quantity table, the layer (from the filename) and mean-only from the protocol, and record it as `inferred_config`. Omit for files from the ingest `publish` step, which carry everything. |
| `--tag` | *(required)* | provenance tag: the **run token** in the filename (`derive_<tag>_<data>_tw<baseline>_<level>.nc`) **and** the `provenance_tag` header attr. Whitespace-stripped, never lowercased — must match the provenance record char-for-char. |
| `--provenance-link` | *(none)* | URL/path to the provenance record; written to the `provenance_link` header attr. |
| `--code-version` | *(required)* | URL to the exact ohc_derive code (commit/release); written to the `ohc_derive_code_version` header attr. |
| `--product-name` | *(required)* | product name; trailing filename token on the **published mask/coverage auxiliaries** only (not the blob). Whitespace-stripped, case preserved. A placeholder is fine if you don't care. |
| `--author` | *(required)* | author; last filename token on the mask/coverage auxiliaries (e.g. `Giglio_etal2026`). |
| `--citation` | *(required)* | citation sentence; written to the coverage `.nc`'s top-level `citation` attr. |
| `--out` | `.` | output directory. |

## Adding a quantity or a mask

A **quantity**: write a recipe `f(primitives, window) -> DataArray(realization, …)` over the helpers in [`temporal_transforms.py`](temporal_transforms.py) and register it in `temporal_transforms.REGISTRY` as `(recipe, primitive)`, naming the step-3 primitive it draws on (`integral` or `map`) — that sets its `reduction` attr. A **mask prescription**: write `f(level, constituents, reference_bathy, require_top) -> (masked, footprint, height)` and register it in `masks.REGISTRY` — `footprint` (lat, lon bool) gives the area, `height` (lat, lon metres) gives the volume. In both cases the runner and the combine do the rest — no other file changes.
