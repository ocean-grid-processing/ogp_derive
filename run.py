#!/usr/bin/env python3
"""ohc_derive factory: build the combined-level analysis quantities for one synthetic level.

Consumes native-level ME4OH submissions (values and NaN on the common grid) plus a standard
bathymetry, and works entirely on the common grid. One invocation handles one synthetic level, so
levels parallelize across jobs. Six steps, in order:

  1. load        the constituents' submissions (mean field + members) and the standard bathy
  2. mask        apply the chosen cross-layer mask prescription over the constituents
  3. primitives  reduce each constituent to its map-level primitives (integral, gridded field)
  4. build       compose the primitives into the requested deliverables (ohca, ohu, trends, map)
  5. collapse    per constituent, collapse the members to a standard deviation; central = mean field
  6. combine     n_fac sum of values, worst-case n_fac sum of standard deviations -> one dataset

Realization axis: the mean field and the members ride together on a leading `realization` axis
(index 0 is the mean field, the rest are members). Steps 3-4 transform every realization the same way;
step 5 reads the central value off index 0 and the spread off the members. So an anomaly demean along
time hits every realization, and each member is referenced to its own window mean.

Members stay per constituent until step 5; step 6 sums values linearly and standard deviations
worst-case. Output is one dataset — each quantity plus its `_sd`, with the footprint area/volume and
the constituents' `quantity` table as attrs; downstream packaging derives the per-area densities and
applies names and layout.
"""
import argparse

import levels
import masks
import map_transforms
import temporal_transforms
import combine
import contracts
import loader


def run(cfg):
    # step 1 — load the constituents' submissions (mean + members) and the standard bathy. Under an input
    # contract the submissions are bare protocol files: no members, quantity and layer inferred.
    if cfg.contract is not None:
        cfg.no_ensemble = True                                   # resolved here so provenance records it
    submissions = loader.load_submissions(cfg.submissions, with_members=not cfg.no_ensemble,
                                          contract=cfg.contract)
    reference_bathy = loader.load_bathy(cfg.bathy)

    # the plan: a table level, or the identity level of a lone submission when --level is omitted.
    level = levels.resolve(cfg.level, submissions)

    token = loader.file_token(cfg, submissions)                    # shared by the .nc and the auxiliaries
    blob = run_level(level, submissions, reference_bathy, cfg, token)
    loader.stamp_chain_provenance(blob, level, cfg, submissions)   # roll upstream chain + stamp ohc_derive_*
    loader.write_blob(blob, level, cfg, token)
    return blob


def run_level(level, submissions, reference_bathy, cfg, token=None):
    """The six steps for one synthetic level -> its dataset."""
    constituents = levels.constituents(level, submissions)          # the native levels this band needs
    quantity = _quantity(submissions, level)
    if cfg.mask is None:
        cfg.mask = masks.default(level, getattr(cfg, "contract", None))   # resolved here so provenance records it
    check_kind(quantity, cfg.quantities, cfg.mask, level)            # the plan must suit the quantity's kind

    # step 2 — apply the cross-layer mask; dumps the mask png and returns the footprint area and volume.
    require_top = cfg.require_top if cfg.require_top is not None else level.require_top
    masked, area_m2, volume_m3 = masks.apply(cfg.mask, level, constituents, reference_bathy,
                                             out_dir=cfg.out, require_top=require_top, tag=cfg.tag,
                                             token=token,
                                             product_name=getattr(cfg, "product_name", ""),
                                             author=getattr(cfg, "author", ""),
                                             citation=getattr(cfg, "citation", ""))

    # step 3 — reduce each constituent to its map-level primitives (integral + gridded field).
    maps = map_transforms.apply(masked, level)

    # step 4 — compose the primitives into the requested deliverables (window sets baseline + trend fit);
    # each is stamped with the field's published units and the primitive it draws on.
    series = temporal_transforms.apply(cfg.quantities, maps, level, window=cfg.time_window,
                                       field_units=quantity["publish_units"])

    # step 5 — collapse each constituent's members to a standard deviation; central from the mean field.
    per_constituent = combine.collapse_sd(series)

    # step 6 — combine constituents: n_fac sum of values, worst-case n_fac sum of standard deviations.
    return combine.combine_synthetic(per_constituent, level, area_m2, volume_m3, quantity)


def check_kind(quantity, quantity_names, mask_name, level):
    """Refuse a plan that would treat the quantity as the wrong kind. `extensive` (a per-area density:
    OHC) sums over area and stacks over layers; `intensive` (a per-cell value: a mixed layer depth) does
    neither. Each primitive and mask prescription declares the kinds it applies to; combining several
    constituents of an intensive quantity would need a thickness-weighted mean, which isn't written.
    """
    kind = quantity["kind"]
    if kind not in ("extensive", "intensive"):
        raise SystemExit("quantity %r has unknown kind %r" % (quantity["name"], kind))
    unknown = [n for n in quantity_names if n not in temporal_transforms.REGISTRY]
    if unknown:
        raise SystemExit("unknown quantity %r; known: %s" % (unknown[0], list(temporal_transforms.REGISTRY)))
    if mask_name not in masks.REGISTRY:
        raise SystemExit("unknown mask prescription %r; known: %s" % (mask_name, list(masks.REGISTRY)))
    problems = []
    for name in quantity_names:
        primitive = temporal_transforms.REGISTRY[name][1]
        if kind not in map_transforms.KINDS[primitive]:
            problems.append("quantity %r draws on the %r primitive, which is not defined for an %s field"
                            % (name, primitive, kind))
    if kind not in masks.KINDS[mask_name]:
        problems.append("mask %r is not defined for an %s field" % (mask_name, kind))
    if kind == "intensive" and len(level.contributors) > 1:
        problems.append("level %s combines %d constituents, but combining an intensive field across "
                        "layers (a thickness-weighted mean) is not implemented"
                        % (level.name, len(level.contributors)))
    if problems:
        raise SystemExit("the %s field (%s) can't be built with this plan:\n  " % (kind, quantity["name"])
                         + "\n  ".join(problems))


def _quantity(submissions, level):
    """The quantity the constituents carry (the ingest [quantity] table). A synthetic level is one
    quantity, so every constituent must carry the same table; a disagreement is a hard error."""
    first = level.contributors[0].tag
    quantity = submissions[first]["quantity"]
    for c in level.contributors[1:]:
        if submissions[c.tag]["quantity"] != quantity:
            raise SystemExit("constituents %s and %s carry different quantity tables; a level is one quantity"
                             % (first, c.tag))
    return quantity


def main():
    ap = argparse.ArgumentParser(description="ohc_derive factory: ME4OH submissions -> one combined level")
    ap.add_argument("submissions", nargs="+", help="the constituent OHC_ submissions (+ OHCENS_ siblings)")
    ap.add_argument("--level", default=None,
                    help="the synthetic level to build (e.g. 0_700). Omit with exactly one submission "
                         "to build its identity level (that native tag, unchanged)")
    ap.add_argument("--bathy", required=True, help="standard bathymetry (NetCDF on the common grid)")
    ap.add_argument("--quantities", required=True,
                    help="comma list of deliverables to build (see temporal_transforms.REGISTRY)")
    ap.add_argument("--mask", default=None,
                    help="mask prescription (see masks.REGISTRY). Default: contiguous_from_top for a "
                         "table level, as_published for an identity level")
    ap.add_argument("--require-top", type=float, default=None,
                    help="metres of the layer's own top that must be defined for a cell to survive; "
                         "overrides the level's own require_top (used by contiguous_from_top)")
    ap.add_argument("--time-window", default=None,
                    help="YEAR0:YEAR1 baseline/trend window (default: all years); separator "
                         "`:`, `-`, or `_` (so the filename token 2004_2025 works too)")
    ap.add_argument("--no-ensemble", action="store_true", help="mean field only; no standard deviations")
    ap.add_argument("--contract", default=None, choices=list(contracts.CONTRACTS),
                    help="treat the submissions as bare protocol files and infer what derive needs from the "
                         "protocol (ME4OH: OHC in TJ/m^2 with the protocol cp0/rho0, layer from the filename, "
                         "no ensemble); recorded as inferred_config. Omit for ingest-published files.")
    ap.add_argument("--tag", required=True, help="provenance tag (filename token + provenance_tag attr)")
    ap.add_argument("--provenance-link", default=None, help="URL/path to the provenance record")
    ap.add_argument("--code-version", required=True,
                    help="URL to the exact ohc_derive code (commit/release); stamped as "
                         "ohc_derive_code_version")
    ap.add_argument("--product-name", required=True,
                    help="product name; trailing filename token on the published mask/coverage auxiliaries "
                         "(placeholder is fine if you don't care)")
    ap.add_argument("--author", required=True,
                    help="author; last filename token on the mask/coverage auxiliaries (e.g. Giglio_etal2026)")
    ap.add_argument("--citation", required=True,
                    help="citation sentence; written to the coverage .nc's top-level `citation` attr")
    ap.add_argument("--out", default=".")
    cfg = ap.parse_args()
    cfg.quantities = [s.strip() for s in cfg.quantities.split(",") if s.strip()]
    if not cfg.quantities:
        raise SystemExit("nothing to build: give --quantities")
    cfg.time_window = _parse_window(cfg.time_window)
    cfg.product_name = "".join(cfg.product_name.split())            # filename tokens: whitespace-stripped,
    cfg.author = "".join(cfg.author.split())                        # case preserved, no other munging
    run(cfg)


def _parse_window(s):
    """YEAR0:YEAR1 -> (int, int); None/empty -> None (all years). Separator may be `:`, `-`, or `_`, so
    the underscore year-range token from the filenames (e.g. `2004_2025`) parses as-is."""
    if not s:
        return None
    y0, y1 = (int(x) for x in s.replace("-", ":").replace("_", ":").split(":"))
    return (y0, y1)


if __name__ == "__main__":
    main()
