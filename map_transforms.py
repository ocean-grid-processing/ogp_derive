"""Step 3 — reduce each masked field to the map-level primitives the deliverables draw on.

Per constituent, across all realizations at once:

    integral   field_value -> (realization, time)            area-weighted global sum
    map        field_value -> (realization, time, lat, lon)   the gridded field, carried through

Step 4 composes these into the deliverables.

Each primitive declares the quantity kinds it is valid for (`KINDS`). `integral` sums a per-area
density over the footprint, which only means something for an extensive quantity; `map` carries
the field through and suits either. An intensive quantity (a mixed layer depth, say) that needs a
footprint-wide value would take an area *mean*, a primitive not written yet — the runner refuses
the plan rather than integrating a depth.
"""
import grid

# primitive -> the quantity kinds it applies to
KINDS = {
    "integral": {"extensive"},
    "map": {"extensive", "intensive"},
}


def integral(field_value):
    """Area-weighted global integral -> (realization, time). NaN cells drop out of the sum."""
    area = grid.cell_area(field_value["lat"].values, field_value["lon"].values)
    return (field_value * area).sum(("lat", "lon"))


def apply(masked, level):
    """Compute the map-level primitives for every constituent.

    -> {tag: {"integral": (realization, time), "map": (realization, time, lat, lon)}}
    """
    return {tag: {"integral": integral(m["field_value"]), "map": m["field_value"]}
            for tag, m in masked.items()}
