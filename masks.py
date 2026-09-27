"""Step 2 — cross-layer masking, pluggable.

A prescription takes the synthetic level, its constituents (each with a `field_value` stack and dbar
bounds), the standard bathy, and the required-top depth, and returns:

    (masked, footprint, height)
    masked    = {tag: {"field_value": DataArray(realization, time, lat, lon), "n_fac": int}}
    footprint = DataArray(lat, lon) bool, True where the level survives (the cells the area counts)
    height    = DataArray(lat, lon) float, the effective column height (metres) per cell (0 outside)

`apply` runs the named prescription, dumps the footprint png and a `coverage_<tag>_<token>_<level>_<mask>_<product_name>_<author>.nc`
(per-cell kept and uncaptured thickness), and returns `(masked, area_m2, volume_m3)` — the summed cell
area of the footprint, and the cell-area-weighted sum of the per-cell height. So land and dropped
columns don't inflate the per-area densities, and the volume tapers with the kept column.

Three prescriptions:

  `fully_wet_nan` — a cell drops where any fully-wet constituent is undefined; a constituent contributes
  where not dry and defined, else 0; the footprint is the surviving wet cells; the height is the level's
  full nominal thickness everywhere in-footprint (a slab). Uses the bathy for wet/dry.

  `contiguous_from_top` — a cell survives only where every constituent covering the layer's top
  `require_top` metres (from `level.low`) is defined; within a surviving cell, keep the constituents
  from the layer top down until the first undefined one, then discard it and everything below (NaN, so
  they vanish from the integral); the height is the n_fac-weighted thickness of the kept run.
  Data-driven; ignores the bathy.

  `as_published` — for the identity level (one constituent) built from our own publish output, asked
  for by name: the mask publish applied is final, so
  the data is untouched; the footprint is the cells finite at every time and member, and the
  prescription verifies that this is the whole story — a cell finite at some times but not others
  (a time-varying footprint) is a hard error, since every stage downstream assumes the footprint is
  time-constant. No height: the level has no vertical extent to report, so no volume and no coverage
  file; the footprint png is still drawn.

Each prescription declares the quantity kinds it is valid for (`KINDS`): `fully_wet_nan` writes 0
into cells that are in the grid but not this level's water, which reads as "contributes nothing"
only for an extensive quantity; `contiguous_from_top` discards with NaN and suits either.

Add a prescription: write `(level, constituents, reference_bathy, require_top) -> (masked, footprint,
height)` and register it, with the kinds it applies to.
"""
import os

import numpy as np
import xarray as xr

import grid


def fully_wet_nan(level, constituents, reference_bathy, require_top=None):
    lat, lon = reference_bathy["lat"], reference_bathy["lon"]
    floor = reference_bathy.values                          # seafloor depth, metres
    nlat, nlon = lat.size, lon.size

    # the one whole-cube reduction: is each cell finite across every member and time step?
    defined = {c["tag"]: c["field_value"].notnull().all(("time", "realization")).values
               for c in constituents}

    exclude = np.zeros((nlat, nlon), dtype=bool)            # cell drops out of the level
    contributes = {c["tag"]: np.zeros((nlat, nlon), dtype=bool) for c in constituents}

    for i in range(nlat):
        for j in range(nlon):
            for c in constituents:
                wet = floor[i, j] >= c["bottom"]            # layer sits entirely above the floor
                dry = np.isnan(floor[i, j]) or floor[i, j] < c["top"]   # entirely below it
                if wet and not defined[c["tag"]][i, j]:
                    exclude[i, j] = True                    # a wet gap kills the column
                elif not dry and defined[c["tag"]][i, j]:
                    contributes[c["tag"]][i, j] = True      # else this cell is a 0 for this level

    masked = {}
    for c in constituents:
        arr = c["field_value"].values.copy()
        arr[:, :, ~contributes[c["tag"]] & ~exclude] = 0.0  # in grid, not this level's water
        arr[:, :, exclude] = np.nan                         # dropped cells
        field_value = xr.DataArray(arr, dims=c["field_value"].dims, coords=c["field_value"].coords)
        masked[c["tag"]] = {"field_value": field_value, "n_fac": c["n_fac"]}

    # footprint: surviving wet cells; height: the level's full nominal thickness there (a slab)
    any_contributes = np.zeros((nlat, nlon), dtype=bool)
    for c in constituents:
        any_contributes |= contributes[c["tag"]]
    footprint = any_contributes & ~exclude
    height = footprint * float(level.nominal_thickness)
    return (masked,
            xr.DataArray(footprint, dims=("lat", "lon"), coords={"lat": lat, "lon": lon}),
            xr.DataArray(height, dims=("lat", "lon"), coords={"lat": lat, "lon": lon}))


def contiguous_from_top(level, constituents, reference_bathy, require_top):
    """Require the layer's top range defined; keep the contiguous defined run down from the top.

    `require_top` is metres of the layer's own top (from `level.low`) that must be present: a cell
    survives only where every constituent covering `level.low .. level.low + require_top` is defined at
    all times and members. The shallowest constituent covers from `level.low` (its n_fac reaches the
    layer top), so any positive `require_top` requires it; a constituent whose `top` sits at or below
    the required depth is not required, and a `require_top` that lands partway into a constituent pulls
    that whole constituent in. In a surviving cell, walk from the layer top down and keep constituents
    until the first undefined one, then discard it and everything below (set to NaN, so they vanish from
    the nan-aware integral). The height is the n_fac-weighted thickness of the kept run. Data-driven —
    the bathy is unused.
    """
    if require_top is None:
        raise SystemExit("contiguous_from_top needs --require-top (metres of the layer top to require)")
    if require_top <= 0:
        raise SystemExit("--require-top must be positive metres, got %g" % require_top)
    if require_top > level.high - level.low:
        raise SystemExit("--require-top %g exceeds the %s layer thickness (%d m)"
                         % (require_top, level.name, level.high - level.low))
    fv0 = constituents[0]["field_value"]
    lat, lon = fv0["lat"], fv0["lon"]
    nlat, nlon = lat.size, lon.size

    defined = {c["tag"]: c["field_value"].notnull().all(("time", "realization")).values
               for c in constituents}
    # required depth is measured from the layer top; the shallowest constituent covers from level.low
    require_depth = level.low + require_top
    required = [c for k, c in enumerate(constituents)
               if (level.low if k == 0 else c["top"]) < require_depth]

    survive = np.ones((nlat, nlon), dtype=bool)             # every required (top) constituent defined
    for c in required:
        survive &= defined[c["tag"]]

    kept = {c["tag"]: np.zeros((nlat, nlon), dtype=bool) for c in constituents}
    for i in range(nlat):
        for j in range(nlon):
            if not survive[i, j]:
                continue
            for c in constituents:                          # layer top first, shallowest-first
                if not defined[c["tag"]][i, j]:
                    break                                   # first gap truncates the column
                kept[c["tag"]][i, j] = True

    masked = {}
    height = np.zeros((nlat, nlon), dtype="float64")
    for c in constituents:
        arr = c["field_value"].values.copy()
        arr[:, :, ~kept[c["tag"]]] = np.nan                 # not kept -> NaN, vanishes in the integral
        field_value = xr.DataArray(arr, dims=c["field_value"].dims, coords=c["field_value"].coords)
        masked[c["tag"]] = {"field_value": field_value, "n_fac": c["n_fac"]}
        height += kept[c["tag"]] * (c["n_fac"] * (c["bottom"] - c["top"]))   # n_fac-weighted kept thickness

    return (masked,
            xr.DataArray(survive, dims=("lat", "lon"), coords={"lat": lat, "lon": lon}),
            xr.DataArray(height, dims=("lat", "lon"), coords={"lat": lat, "lon": lon}))


def as_published(level, constituents, reference_bathy, require_top=None):
    """Take publish's mask as final and verify it is time-constant; see the module docstring."""
    if len(constituents) != 1:
        raise SystemExit("as_published takes exactly one constituent (the identity level); level %s has %d"
                         % (level.name, len(constituents)))
    c = constituents[0]
    finite = c["field_value"].notnull()
    always = finite.all(("time", "realization"))                 # the footprint: defined everywhere in the cube
    ever = finite.any(("time", "realization"))
    partial = int((ever & ~always).sum())
    if partial:
        raise SystemExit("%d cell(s) of %s are defined at some times or members but not all; the footprint "
                         "must be time-constant (publish with incomplete_timeseries and ensemble_incomplete "
                         "honored, or write a mask prescription for a time-varying footprint)"
                         % (partial, c["tag"]))
    masked = {c["tag"]: {"field_value": c["field_value"], "n_fac": c["n_fac"]}}
    return masked, always, None


REGISTRY = {
    "fully_wet_nan": fully_wet_nan,
    "contiguous_from_top": contiguous_from_top,
    "as_published": as_published,
}

# prescription -> the quantity kinds it applies to
KINDS = {
    "fully_wet_nan": {"extensive"},
    "contiguous_from_top": {"extensive", "intensive"},
    "as_published": {"extensive", "intensive"},
}


# The prescription used when --mask is not given, for every level and every input. `as_published` is
# a verification of our own publish step's footprint and is asked for by name.
DEFAULT = "contiguous_from_top"


def apply(name, level, constituents, reference_bathy, out_dir=".", require_top=None, tag=None, token=None,
          product_name="", author="", citation=""):
    """Run the named prescription, dump its footprint png, return (masked, area_m2, volume_m3).
    A prescription that reports no height (`as_published`) yields no volume and no coverage file."""
    if name not in REGISTRY:
        raise SystemExit("unknown mask prescription %r; known: %s" % (name, list(REGISTRY)))
    masked, footprint, height = REGISTRY[name](level, constituents, reference_bathy, require_top)
    _dump_png(footprint, reference_bathy, level.name, name, out_dir, tag, token, product_name, author)
    area = grid.cell_area(footprint["lat"].values, footprint["lon"].values)
    area_m2 = float(area.where(footprint).sum())
    if height is None:
        return masked, area_m2, None
    _dump_coverage(footprint, height, reference_bathy, level, name, out_dir, tag, token,
                   product_name, author, citation)
    return masked, area_m2, float((area * height).sum())


def _dump_coverage(footprint, height, reference_bathy, level, mask_name, out_dir, tag=None, token=None,
                   product_name="", author="", citation=""):
    """Sibling .nc to the mask png, on the product grid, with two per-cell diagnostics (NaN off-footprint):

      kept_thickness       — the n_fac-weighted column thickness kept at the cell (the volume's height).
      uncaptured_thickness — in-layer water the kept column didn't reach: from the bottom of the kept
                             column (level.low + kept_thickness) down to min(bathy, layer bottom), >= 0.
    """
    fp = footprint.values
    column_bottom = level.low + height.values                        # bottom of the deepest kept layer
    water_bottom = np.minimum(reference_bathy.values, float(level.high))
    uncaptured = np.clip(water_bottom - column_bottom, 0.0, None)

    ds = xr.Dataset(
        {"kept_thickness": (("lat", "lon"), np.where(fp, height.values, np.nan)),
         "uncaptured_thickness": (("lat", "lon"), np.where(fp, uncaptured, np.nan))},
        coords={"lat": footprint["lat"], "lon": footprint["lon"]})
    ds["kept_thickness"].attrs = {"units": "m", "long_name": "kept column thickness"}
    ds["uncaptured_thickness"].attrs = {
        "units": "m", "long_name": "in-layer water below the kept column, to min(bathy, layer bottom)"}
    ds.attrs.update({"level": level.name, "mask": mask_name})
    if citation:
        ds.attrs["citation"] = citation                              # full citation sentence, top-level
    os.makedirs(out_dir, exist_ok=True)
    tokens = [t for t in (tag, token) if t]                          # tag/token optional; level+mask always
    tail = [t for t in (product_name, author) if t]                  # publication pair, last before .nc
    stem = "_".join(["coverage", *tokens, level.name, mask_name, *tail])
    path = os.path.join(out_dir, stem + ".nc")
    ds.to_netcdf(path)
    print("wrote", path)


def _dump_png(footprint, reference_bathy, level_name, mask_name, out_dir, tag=None, token=None,
              product_name="", author=""):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap
        from matplotlib.patches import Patch
    except Exception:
        return
    # three classes: in the footprint (white), ocean outside it (black), land (grey, from the bathy)
    category = np.where(footprint.values, 1, 0)            # 1 in footprint, 0 outside
    category[reference_bathy.values < 0] = 2               # negative depth = above sea level -> land
    cmap = ListedColormap(["black", "white", "lightgrey"])

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.pcolormesh(footprint["lon"], footprint["lat"], category, shading="auto",
                  cmap=cmap, vmin=-0.5, vmax=2.5)
    ax.set_title("%s footprint (%s)" % (level_name, mask_name))
    ax.set_xlabel("lon")
    ax.set_ylabel("lat")
    ax.legend(handles=[Patch(facecolor="white", edgecolor="grey", label="in footprint"),
                       Patch(facecolor="black", label="outside footprint"),
                       Patch(facecolor="lightgrey", label="land")],
              loc="lower left", fontsize=8, framealpha=0.9)
    os.makedirs(out_dir, exist_ok=True)
    tokens = [t for t in (tag, token) if t]                          # tag/token optional; level+mask always
    tail = [t for t in (product_name, author) if t]                  # publication pair, last before .png
    stem = "_".join(["mask", *tokens, level_name, mask_name, *tail])
    path = os.path.join(out_dir, stem + ".png")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    print("wrote", path)
