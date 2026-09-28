# ogp_derive

`ogp_derive` is the pipeline's factory. It takes ME4OH-shaped submissions — native-level fields on the
common grid — and builds analysis quantities from them, one NetCDF per level: for ocean heat content,
the combined-layer OHCA and OHU series, their trends, and the gridded anomaly map; for a single-layer
quantity such as a mixed layer depth, the field itself with its ensemble spread. All physics and statistics lives here — cross-layer masking, area integration, anomaly referencing, the ensemble → standard-deviation
collapse, and the cross-layer combine. The emitters downstream only package.

Submissions come from two places: the LocalGP ingest stage
([`localgp_ogp_ingest`](../localgp_ogp_ingest/)), whose `publish` step annotates its files with
everything derive reads, or any group that stopped at the ME4OH protocol, with `--contract ME4OH`.

## Inputs

A submission is `DATA(LONGITUDE, LATITUDE, TIME)` on the common 1° grid, monthly, missing as NaN —
the ME4OH shape. Derive also needs a standard bathymetry on the same grid (`--bathy`; `etopo60.cdf`,
the file ingest pins the grid to); `cell_area` is regenerated from the grid.

### What a submission must carry

On top of the protocol shape, derive reads four things from a submission:

| what | where | used for |
|---|---|---|
| the **quantity table** | `quantity` attr — the ingest `[quantity]` table as compact JSON: `name`, `kind`, `units`, `scale_terms`, published units | what the field is, its kind (extensive/intensive), its units, and scale factors |
| the **native-level tag** | `mapped_layer` attr, `<top>_<bottom>` | keying the constituent a level plan asks for |
| the **ensemble** | a `<NAME>ENS_…` sibling file, `DATA(MEMBER, LONGITUDE, LATITUDE, TIME)`, found by swapping the `<NAME>_` filename prefix (from `quantity.name`) | the `_sd` companions |
| the **provenance chain** | `*_run_config` / `*_run_facts` / `*_code_version` attrs | rolled forward into the blob (optional — absent blocks are simply not forwarded) |

**Files from `localgp_ogp_ingest` carry all of this** — its `publish.py` stamps the `quantity` table,
`mapped_layer` and the provenance blocks, and writes the `<NAME>ENS_` sibling with `--ensemble`. Nothing
further is needed.

**Files that stop at the ME4OH protocol — other groups' submissions — carry none of it, and don't
have to.** Pass `--contract ME4OH` and derive fills the gaps from the protocol, which fixes them:
`DATA` is ocean heat content density in TJ/m² computed with the protocol's `cp0 = 3989.244 J/kg/K`
and `rho0 = 1030 kg/m³`, so the quantity table is the protocol's OHC table; the layer is the
filename's `lev<low>_<high>` token; there is no ensemble, so the run is mean-only (a `DATA_SD`, if
present, is not used). Everything inferred is recorded on the blob in
`ohc_derive_run_facts.inferred_config`, keyed by constituent, so the record says the table came from
the specification and not from the file. A bare protocol file passed *without* `--contract` is
refused with a message pointing at it — derive never infers silently. [`contracts.py`](contracts.py)
is the one place the protocol's assumptions live.

## The six steps

One run builds one synthetic level. The mean field and its members ride together on a leading `realization`
axis (index 0 is the mean field, the rest are members); steps 3–4 transform every realization the same
way, and step 5 reads the central value off index 0 and the spread off the members — so an anomaly
demean along time hits every realization, and each member is referenced to its own window mean.

1. **load** — read the constituents' submissions (mean field + members) and the standard bathy.
2. **mask** — apply the cross-layer mask prescription over the constituents (below).
3. **primitives** — reduce each constituent to its map-level primitives: the area-weighted `integral`
   (a time series) and the gridded field (`map`).
4. **build** — compose the primitives into the requested quantities (below).
5. **collapse** — per constituent, take the central value from the mean field and the 1σ spread
   across the members.
6. **combine** — `n_fac`-sum the constituents into the level: values linearly.

## Levels: the plan file and the identity level

A **synthetic level** is an `n_fac`-weighted sum of native levels, shallowest first — for example
LocalGP's `0_2000` is `15_20`(×3) + `15_300` + `300_700` + `700_1850` + `1800_1850`(×3). `n_fac`
scales a thin measured layer up to the slab it stands in for; each constituent carries its own dbar
`top`/`bottom`, used against the bathy in the mask. The table of levels is a **plan file** passed with
`--levels`: TOML, one `[[level]]` per synthetic level with its `name`, `require_top` and
`contributors`. LocalGP's is [`levels/localgp.toml`](levels/localgp.toml) (`0_300`, `0_700`,
`0_1000`, `700_2000`, `0_2000`). A group whose native layers decompose the target layers differently
writes its own plan — "`0_2000` is one `0_700` plus one `700_2000`," for example — and nothing else changes. On
load, derive checks each level: contributors shallowest first, no repeated tag, the `n_fac`-weighted
thicknesses summing to the level's thickness (`high − low`), and `require_top` within it; a plan that
doesn't tile its layer is refused, since it would combine silently into the wrong number. The plan used
is recorded in `ohc_derive_run_facts` (`levels_file`, `level_plan`).

A level is a plan for combining constituents, and one constituent needs no plan: omit `--level` with
exactly one submission in the pool and the run builds that submission's **identity level** — its
native tag as the level, one contributor, `n_fac = 1`. This is how a single-layer quantity (a mixed
layer depth for example) goes through the factory, and how another group's directly-submitted target layer does:
mask, primitives, collapse, and a combine that returns the constituent.

## Kinds

The `quantity` table names the field's **kind**: `extensive` (a per-area density, like OHC in J/m²,
that sums over area and stacks over layers) or `intensive` (a per-cell value, like a mixed layer
depth, that does neither). Each primitive and mask prescription declares the kinds it applies to
(`map_transforms.KINDS`, `masks.KINDS`), and a run is refused up front if its plan doesn't suit the
field: an intensive field can't take the `integral`-based quantities, can't use `fully_wet_nan` (its 0-fill means "adds nothing", an extensive idea), and
can't be combined across several constituents. The error
names each mismatch.

## Cross-layer masking

Cross-layer masking (dropping cells based on logic that considers the full column rather than just individual native layers as is donw in the ingest shims) is pluggable (`masks.REGISTRY`). A prescription returns the masked constituents, a footprint
(the cells counted for area), and a per-cell column height (metres); `apply` turns those into the
level's `area_m2` (footprint cell area) and, when there is a height, `volume_m3` (the cell-area-weighted
sum of the height, so the volume tapers with the kept column rather than assuming a slab). Alongside it
writes the footprint as `mask_<tag>_<data>_tw<baseline>_<level>_<mask>_<product_name>_<author>.png`
and, when there is a height, `coverage_…_<product_name>_<author>.nc` on the product grid — per cell,
the `kept_thickness` and the `uncaptured_thickness` (in-layer water the kept column didn't reach, from
the bottom of the kept column down to `min(bathy, layer bottom)`, floored at 0), with the `--citation`
sentence as a top-level attr. The `<data>_tw<baseline>` token matches the main `.nc`, so runs that
differ only in baseline don't overwrite each other's auxiliaries; `<product_name>_<author>` is the
publication descriptor, since these two files are published.

"Defined" is judged across every member and timestep together, so the footprint is identical for all
realizations and the ensemble spread reflects real spread rather than footprint jitter. Every stage
downstream assumes the footprint is time-constant.

**`contiguous_from_top`** (the default when `--mask` is omitted). `require_top` is a thickness of the
layer's own top — the metres below `level.low` that must be present — and each level carries its own in
the plan file (every LocalGP level: `300`; `--require-top` overrides it). A cell **survives** only
where every constituent reaching into the layer's top `require_top` metres is defined at all times and
members. Constituents are atomic (a whole native layer or nothing), so a `require_top` that lands
partway into one requires that whole constituent; selection is by the constituent's declared bounds.
With `require_top = 300`: `0_2000` requires `15_20` + `15_300`, whose bounds tile exactly 0–300;
`700_2000` requires `700_1850`, because depth 1000 falls within its bounds. The shallowest
constituent's `n_fac` reaches the layer top, so any positive `require_top` requires it. Within a
surviving cell, walk the constituents from the layer top down and **keep** them until the first
undefined one, then discard it and everything below (NaN, so they vanish from the nan-aware integral).
The kept run's `n_fac`-weighted thickness is the cell's height (a kept `1800_1850` adds 3·50 = 150 m),
so the volume is the true tapering ocean. With one constituent this reduces to "keep a cell where the
field is defined at every time step" — the right choice for a bare protocol file, whose footprint
derive has to coordinate itself.

**`as_published`** — asked for by name, for an identity level built from `localgp_ogp_ingest` output.
The mask publish applied is final: the data is untouched, the footprint is the cells finite at every
time and member, and the prescription **verifies** that this is time-constant — a cell finite at some
times or members but not others is a hard error (publish with `incomplete_timeseries` and
`ensemble_incomplete` honored makes it so). It reports no column height, so the blob carries `area_m2`
but no `volume_m3`, and no coverage `.nc` is written; the footprint png is.

**`fully_wet_nan`** — the earlier, bathy-driven prescription: each constituent is classified against
the standard bathy as fully wet (floor ≥ `bottom`), intersecting, or dry; a cell drops where any
fully-wet constituent is undefined, a constituent contributes where not dry and defined (else 0), and
the height is the level's full nominal thickness everywhere in-footprint (a slab). Extensive only.

## Quantities

`--quantities` selects from `temporal_transforms.REGISTRY`. `--time-window` (`YEAR0:YEAR1`) sets the
anomaly baseline and the years the trends are fit over; omitted, it spans all years.

| name | dims | primitive | what it is |
|---|---|---|---|
| `ohca` | (year,) | integral | monthly area-integrated anomaly (window baseline), annual-averaged |
| `ohu` | (year,) | integral | month-to-month tendency of the integral, annual-averaged; NaN-seeded at t0 |
| `ohca_trend` | scalar | integral | OLS slope of the annual integral over the window |
| `ohu_trend` | scalar | integral | OLS slope of the annual tendency over the window |
| `map` | (time, lat, lon) | map | per-cell monthly anomaly (window baseline) |
| `field` | (time, lat, lon) | map | the masked field itself, per cell and month — no transform; its `_sd` is the per-cell member spread. Ignores the window |

The integral-based quantities are extensive — the per-area field integrated over the footprint (TJ
from a TJ/m² field), with the trends and tendency carrying the matching per-year and per-month
scaling; `map` and `field` stay in the field's per-area unit. Each variable (and its `_sd`) says so in
its attrs: `field_units` is the submission's published units (from the `quantity` table) and
`reduction` is `area_integral` (summed over the footprint's cell areas — units are field units × m²,
and dividing by `area_m2` recovers a per-area density) or `grid` (still per cell). A rate also carries
`per`, the step it is a rate over: `ohu` is the annual mean of a month-to-month difference, so
`per = "month"`; the two trends are slopes against the annual axis, so `per = "year"`. `ohca`, `map`
and `field` are states and carry no `per`. Converting to the deliverables' units (J/m², W/m², …) is
the emitters' job, from this metadata and the geometry attrs.

## Output

One NetCDF per level, `derive_<tag>_<data>_tw<baseline>_<level>.nc`. The token carries both year
spans: `<data>` (`YYYY_YYYY`) is the data span present, and `tw<baseline>` is the `--time-window`
baseline, defaulting to the whole data span — so a windowless run reads `…_2004_2025_tw2004_2025_…`, a
2005–2024 baseline `…_2004_2025_tw2005_2024_…`, and runs differing only in baseline don't collide. It
holds each requested quantity plus its `_sd` companion (omitted under `--no-ensemble` or `--contract
ME4OH`), and the attrs `level`, `area_m2`, `volume_m3` (absent under `as_published`), `time_window`,
`quantity` (the table, one compact-JSON attr — the gcos emitter reads `cp0`/`rho0` from its
`scale_terms`), `provenance_tag`, `provenance_link`, and the provenance chain.

**Provenance chain.** derive is a fan-in — one level from N constituents — so it rolls each upstream
block forward *grouped by constituent*: for every `*_run_config` / `*_run_facts` / `*_code_version` on
the inputs (`localgp_ingest_*`, `localgp_publish_*`), the output carries `<block> = {constituent_tag:
block}` (compact JSON, one line). Blocks are opaque — parsed only to nest, never read. On top of that
derive stamps its own `ohc_derive_run_config` (the resolved arguments, including the mask actually
used), `ohc_derive_run_facts` (level, `identity_level`, `levels_file` and `level_plan`, mask,
`require_top`, window, ensemble, geometry, `n_fac`, the `quantity` table, and `inferred_config` under a
contract), and `ohc_derive_code_version`. The stage key is `ohc_derive` — the repository's former name,
kept because downstream files key on it. The global `provenance_tag` / `provenance_link` are this run's
own.

## Usage

### Environment

`Dockerfile` builds the test/run environment (`numpy`, `xarray`, `netCDF4`, `matplotlib`, `pytest`;
Python ≥ 3.11 for `tomllib`, or `tomli`); make the equivalent conda env on the cluster.

### Test

```
docker image build -t ogp_derive:test .
docker container run -v $(pwd):/app ogp_derive:test pytest
```

### Run

Three happy paths, each with a launcher in [`examples/`](examples/):

- **LocalGP OHC** — [`examples/derive_ohc.slurm`](examples/derive_ohc.slurm): the pool of `OHC_`
  submissions, `--levels levels/localgp.toml --level <level>`, `--time-window <baseline>`,
  `--quantities ohca,ohu,ohca_trend,ohu_trend,map`, `--mask contiguous_from_top`. `run.sh` loops the
  windows and levels over it.
- **LocalGP single-layer quantity (MLD)** — [`examples/derive_mld.slurm`](examples/derive_mld.slurm):
  the one `MLD_` submission, no `--level`, `--quantities field`, `--mask as_published`.
- **Another group's ME4OH files** — as the OHC path, plus `--contract ME4OH`; no `--level` if they
  submit the target layers directly, or `--levels <their plan>` if they decompose them. Mean-only by
  construction.

With the ensemble on (the default), each constituent's `<NAME>ENS_` sibling **must** sit next to its
`<NAME>_` file or the loader exits, and every quantity gets a collapsed `_sd`; central values are the
same either way.

#### run.py options

All configuration is on the command line — no env, no config file. Every resolved option lands in
`ohc_derive_run_config`.

| option | required | default | what it does |
|---|:--:|---|---|
| `SUBMISSION.nc …` (positional) | **yes** | | the constituent submissions (`<NAME>_…`, or bare ME4OH files with `--contract ME4OH`); `<NAME>ENS_` siblings are found automatically. A level selects the constituents it needs by tag, so the whole pool can be passed — but it must hold **exactly one file per native level** (a duplicate tag is a hard error) |
| `--bathy` | **yes** | | standard bathymetry NetCDF on the common grid |
| `--quantities` | **yes** | | comma list from `ohca,ohu,ohca_trend,ohu_trend,map,field`; unknown names error |
| `--tag` | **yes** | | provenance tag: the run token in the filename and the `provenance_tag` attr. Whitespace-stripped, never lowercased — must match the provenance record char-for-char |
| `--code-version` | **yes** | | URL to the exact `ogp_derive` commit/release → `ohc_derive_code_version` |
| `--product-name` | **yes** | | trailing filename token on the published mask/coverage auxiliaries (not the blob); whitespace-stripped, case preserved. A placeholder is fine |
| `--author` | **yes** | | last filename token on the mask/coverage auxiliaries (e.g. `Giglio_etal2026`) |
| `--citation` | **yes** | | citation sentence → the coverage `.nc`'s `citation` attr |
| `--level` | | *(identity level)* | the synthetic level to build, by name in the `--levels` plan. Omitted: the identity level of the one submission in the pool (any other pool size is an error) |
| `--levels` | with `--level` | | the level plan (TOML); LocalGP's is `levels/localgp.toml`. Required whenever `--level` is given; an identity level needs no plan |
| `--mask` | | `contiguous_from_top` | mask prescription: `contiguous_from_top`, `as_published`, or `fully_wet_nan` |
| `--require-top` | | *(the level's own)* | metres of the layer's own top that must be defined for a cell to survive; overrides the plan's `require_top`. `contiguous_from_top` only |
| `--time-window` | | *(all years)* | `YEAR0:YEAR1` — the anomaly baseline and the trend-fit years; separator `:`, `-` or `_` |
| `--no-ensemble` | | off | mean field only: no `_sd`, `<NAME>ENS_` siblings not read |
| `--contract` | | *(none)* | `ME4OH`: the submissions are bare protocol files; infer the quantity table, the layer (from the filename) and mean-only from the protocol, recorded as `inferred_config` |
| `--provenance-link` | | *(none)* | URL/path to the provenance record → `provenance_link` attr |
| `--out` | | `.` | output directory |

## Adding a quantity or a mask

A **quantity**: write a recipe `f(primitives, window) -> DataArray(realization, …)` over the helpers
in [`temporal_transforms.py`](temporal_transforms.py) and register it in
`temporal_transforms.REGISTRY` as `(recipe, primitive)`, naming the step-3 primitive it draws on
(`integral` or `map`) — that sets its `reduction` attr and which kinds it suits. A **mask
prescription**: write `f(level, constituents, reference_bathy, require_top) -> (masked, footprint,
height)`, register it in `masks.REGISTRY`, and declare its kinds in `masks.KINDS` — `footprint` (lat,
lon bool) gives the area, `height` (lat, lon metres, or `None`) gives the volume. In both cases the
runner and the combine do the rest.
