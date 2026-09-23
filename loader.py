"""Step 1 (load) and the final write.

`load_submissions` reads each native-level submission (mean field) and, unless mean-only, its member
sibling — the `<NAME>ENS_` file publish wrote beside the `<NAME>_` submission, NAME being the quantity's
name from the submission's `quantity` attr — and stacks them as `field_value` on a leading
`realization` axis (index 0 = mean field, the rest = members), so downstream steps treat every
realization the same way.

`load_bathy` reads the standard bathymetry as a depth grid.

`write_blob` writes one synthetic level's dataset (each quantity plus its `_sd`) with the provenance
tag and link. `stamp_chain_provenance` rolls every upstream provenance block forward (grouped by
constituent, since derive is a fan-in) and adds this step's own `ohc_derive_*` blocks.
"""
import json
import os

import numpy as np
import xarray as xr

# This step's identity, used to namespace its provenance (`ohc_derive_run_config` / `_run_facts` /
# `_code_version`). Every step rolls all `*_run_config` / `_run_facts` / `_code_version` forward and
# adds its own; at a fan-in they're grouped by constituent, so blocks accrete without collision.
STAGE = "ohc_derive"
_PROV_SUFFIXES = ("_run_config", "_run_facts", "_code_version")


def _maybe_json(v):
    """Parse an upstream block back to JSON so it nests as a real object; leave non-JSON as-is. This is
    structural only — we never read the block's fields, so no coupling to the upstream schema."""
    try:
        return json.loads(v)
    except (TypeError, ValueError):
        return v


def _compact(obj):
    """One-line JSON — reads as a single clean line in `ncdump -h`."""
    return json.dumps(obj, separators=(",", ":"), default=str)


def _to_tlatlon(da):
    return da.transpose("TIME", "LATITUDE", "LONGITUDE").rename(
        {"TIME": "time", "LATITUDE": "lat", "LONGITUDE": "lon"})


def _member_sibling(path, quantity_name):
    """The member file that publish wrote next to a submission: publish names the pair
    `<NAME>_…` / `<NAME>ENS_…` from the quantity's name (upper-cased), so the same name resolves
    the sibling here. Errors if the submission doesn't carry the prefix its quantity attr implies."""
    prefix = quantity_name.upper() + "_"
    base = os.path.basename(path)
    if not base.startswith(prefix):
        raise SystemExit("%s does not carry the %s prefix its quantity attr (name=%r) implies"
                         % (base, prefix, quantity_name))
    return os.path.join(os.path.dirname(path), prefix[:-1] + "ENS_" + base[len(prefix):])


def _quantity(attrs, path):
    """The submission's `quantity` attr (the ingest [quantity] table, compact JSON) -> dict."""
    if "quantity" not in attrs:
        raise SystemExit("%s has no `quantity` attr (expected a submission from the current publish)"
                         % os.path.basename(path))
    return json.loads(attrs["quantity"])


def _load_members(path, quantity_name):
    sib = _member_sibling(path, quantity_name)
    if not os.path.exists(sib):
        raise SystemExit("no member sibling for %s (looked for %s); pass --no-ensemble for mean only"
                         % (os.path.basename(path), os.path.basename(sib)))
    da = xr.open_dataset(sib, decode_times=True)["DATA"]
    return (da.transpose("MEMBER", "TIME", "LATITUDE", "LONGITUDE")
              .rename({"MEMBER": "realization", "TIME": "time", "LATITUDE": "lat", "LONGITUDE": "lon"})
              .astype("float64"))


def _stack(mean_da, member_da):
    """(time, lat, lon) mean + optional (realization, time, lat, lon) members -> one realization stack."""
    mean_r = mean_da.expand_dims(realization=[0])
    if member_da is None:
        return mean_r
    members = member_da.assign_coords(realization=np.arange(1, member_da.sizes["realization"] + 1))
    return xr.concat([mean_r, members], dim="realization")


def load_submissions(paths, with_members=True):
    """paths -> {tag: {"field_value": DataArray(realization, time, lat, lon), "attrs": dict,
    "quantity": dict}} — `quantity` is the ingest [quantity] table the submission carries.

    Submissions are keyed by their native-level tag, and each level a synthetic level needs is selected
    by tag — so passing the whole pool and letting each run pick its constituents is fine. But that only
    works if the pool holds exactly one file per native level: two files with the same tag (a stray
    window / experiment / rerun) is a hard error rather than a silent last-wins that would combine the
    wrong data.
    """
    subs = {}
    seen = {}
    for p in paths:
        ds = xr.open_dataset(p, decode_times=True)
        if "DATA" not in ds.data_vars:
            raise SystemExit("%s has no DATA variable (expected an ME4OH submission)" % p)
        tag = ds.attrs.get("mapped_layer") or ds.attrs.get("layer_m")
        if not tag or "_" not in str(tag):
            raise SystemExit("%s has no usable mapped_layer/layer_m attr (got %r)" % (p, tag))
        tag = str(tag)
        if tag in seen:
            raise SystemExit("two submissions map to native level %s:\n  %s\n  %s\n"
                             "the pool must hold exactly one file per native level." % (tag, seen[tag], p))
        seen[tag] = p
        quantity = _quantity(ds.attrs, p)
        mean_da = _to_tlatlon(ds["DATA"]).astype("float64")
        members = _load_members(p, quantity["name"]) if with_members else None
        subs[tag] = {"field_value": _stack(mean_da, members), "attrs": dict(ds.attrs), "quantity": quantity}
    return subs


def load_bathy(path):
    """The standard bathymetry (etopo60.cdf) -> DataArray(lat, lon), seafloor depth in metres.

    This is the same file ohc_ingest pins the mapping grid to: relief variable ROSE (metres, negative
    below sea level) on dims ETOPO60Y (latitude) and ETOPO60X (longitude), row-major [lat, lon]. We
    read it as-is, so lat/lon come out in the mapping-grid order the submissions also carry, and negate
    the relief to positive-down depth to match ingest's seabed convention.
    """
    da = xr.open_dataset(path)["ROSE"]
    da = da.rename({"ETOPO60Y": "lat", "ETOPO60X": "lon"}).transpose("lat", "lon").astype("float64")
    return -da


def stamp_chain_provenance(blob, level, cfg, submissions):
    """Roll the upstream provenance chain forward and add this step's own blocks.

    derive is a fan-in: N constituents each carry their own `*_run_config` / `_run_facts` /
    `_code_version` (localgp_ingest_*, localgp_publish_*, …). We group each block by the constituent it
    came from — `<block> = {constituent_tag: block}` — so nothing is deduplicated away and downstream
    never has to know who produced what. Blocks are opaque: parsed only to nest cleanly, never read.
    """
    contributors = [c.tag for c in level.contributors]
    forwarded = {}                                            # block_name -> {constituent_tag: block}
    for tag in contributors:
        for k, v in submissions[tag]["attrs"].items():
            if k.endswith(_PROV_SUFFIXES):
                forwarded.setdefault(k, {})[tag] = _maybe_json(v)
    for block_name, per_constituent in forwarded.items():
        blob.attrs[block_name] = _compact(per_constituent)

    require_top = cfg.require_top if cfg.require_top is not None else level.require_top
    blob.attrs["%s_code_version" % STAGE] = cfg.code_version
    blob.attrs["%s_run_config" % STAGE] = _compact(vars(cfg))
    blob.attrs["%s_run_facts" % STAGE] = _compact({
        "level": level.name,
        "quantities": cfg.quantities,
        "mask": cfg.mask,
        "require_top": require_top,
        "time_window": "%d-%d" % cfg.time_window if cfg.time_window else "all",
        "ensemble": not cfg.no_ensemble,
        "area_m2": blob.attrs.get("area_m2"),
        "volume_m3": blob.attrs.get("volume_m3"),
        "constituents": contributors,
        "n_fac": {c.tag: c.n_fac for c in level.contributors},
        "quantity": _maybe_json(blob.attrs.get("quantity")),
    })


def _record_span(blob):
    """(year0, year1) spanned by the blob's own axis — `year` for the annual quantities, `time` for the
    gridded/monthly ones — or None if it carries neither (e.g. a trend-only blob)."""
    if "year" in blob.coords:
        yrs = blob["year"].values.astype(int)
        return int(yrs.min()), int(yrs.max())
    if "time" in blob.coords:
        yrs = blob["time"].values.astype("datetime64[Y]").astype(int) + 1970
        return int(yrs.min()), int(yrs.max())
    return None


def _submissions_span(submissions):
    """(year0, year1) data span from the submissions' shared monthly time axis, or None."""
    for s in submissions.values():
        t = s["field_value"]["time"].values
        if t.size:
            yrs = t.astype("datetime64[Y]").astype(int) + 1970
            return int(yrs.min()), int(yrs.max())
    return None


def _combined_token(data_span, cfg):
    """Filename token carrying both spans: `<data>_tw<baseline>`, each `YYYY_YYYY`. `data` is the years
    actually present; `baseline` is the `--time-window`, defaulting to the whole data span. So a
    windowless run reads as e.g. `2004_2025_tw2004_2025`, and a 2005-2024 baseline as
    `2004_2025_tw2005_2024`."""
    data = "%d_%d" % data_span if data_span else "all"
    baseline = "%d_%d" % cfg.time_window if cfg.time_window else data
    return "%s_tw%s" % (data, baseline)


def _file_token(cfg, blob):
    """The combined `<data>_tw<baseline>` token, data span read from the finished blob's own axis."""
    return _combined_token(_record_span(blob), cfg)


def file_token(cfg, submissions):
    """The combined token, data span read from the submissions' time axis — so the mask/coverage
    auxiliaries (written before the blob exists) share the main file's token."""
    return _combined_token(_submissions_span(submissions), cfg)


def write_blob(blob, level, cfg, token=None):
    """Write one synthetic level's dataset to NetCDF, tagged with cfg.tag and provenance link."""
    os.makedirs(cfg.out, exist_ok=True)
    blob.attrs["level"] = level.name
    # The attr keeps its intent semantics ("all" = whole-record baseline) — the gcos emitter reads it to
    # reject a windowless derive. The filename carries both concrete year ranges (data span and baseline)
    # so runs that differ only in baseline can't collide (same tag+level, different content).
    blob.attrs["time_window"] = "%d-%d" % cfg.time_window if cfg.time_window else "all"
    blob.attrs["provenance_tag"] = cfg.tag
    if cfg.provenance_link is not None:
        blob.attrs["provenance_link"] = cfg.provenance_link
    tok = token if token is not None else _file_token(cfg, blob)
    path = os.path.join(cfg.out, "derive_%s_%s_%s.nc" % (cfg.tag, tok, level.name))
    blob.to_netcdf(path, engine="netcdf4")
    print("wrote", path)
    return path
