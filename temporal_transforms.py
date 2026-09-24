"""Step 4 — build the deliverables from the step-3 primitives.

`primitives` is one constituent's step-3 output: {"integral": (realization, time),
"map": (realization, time, lat, lon)}. `window` is (year0, year1) or None. Each deliverable is a
short recipe over the small helpers below:

    ohca         monthly anomaly (window baseline), then annual mean   -> (realization, year)
    ohu          month-to-month tendency, then annual mean             -> (realization, year)
    ohca_trend   OLS slope of the annual integral over the window      -> (realization,)
    ohu_trend    OLS slope of the annual tendency over the window      -> (realization,)
    map          per-cell monthly anomaly (window baseline)            -> (realization, time, lat, lon)
    field        the masked field itself, per cell and month             -> (realization, time, lat, lon)

`_tendency` drops its leading step (no prior month), so `ohu` takes its annual mean with
`complete=True` — a year missing that step is NaN, not a partial average — and `_slope` skips NaN
years, keeping the dropped step out of the fit. Every rate carries a `per` attr naming the step it
is a rate over — `ohu` "month" (the annual mean of a monthly difference is still per month), the
trends "year" — so packaging divides by the matching seconds-per-step.

Every deliverable is stamped with what it is made of, so packaging converts units from metadata:
`field_units` (the submission's published units, from the `quantity` table) and `reduction` — the
step-3 primitive the recipe draws on: `area_integral` (the field summed over the footprint's cell
areas; units are field units x m2, and dividing by the level's `area_m2` recovers a per-area density)
or `grid` (still per cell; units are the field units). `REGISTRY` maps each name to its recipe and
its primitive.

Add a deliverable: write a recipe over the helpers and register it with the primitive it uses.
"""

# reduction attr value per step-3 primitive
REDUCTION = {"integral": "area_integral", "map": "grid"}


def _in_window(series, dim, window):
    """Clip `series` along `dim` to `window`, an inclusive (lo, hi) pair in that axis's own terms, or
    None for the whole series. Bounds are integer years for the annual `year` axis, year strings (e.g.
    "2005") for the monthly `time` axis, which xarray reads as partial-datetime bounds.
    """
    lo, hi = window or (None, None)
    return series.sel({dim: slice(lo, hi)})


def _anomaly(series, window):
    """Subtract the mean over the window months, per realization."""
    months = None if window is None else (str(window[0]), str(window[1]))
    return series - _in_window(series, "time", months).mean("time")


def _annual(series, complete=False):
    """Calendar-year mean. With complete=True, a year missing any month is NaN instead of a partial
    mean — the treatment a differenced series needs, so its dropped leading step voids that year.
    """
    return series.groupby("time.year").mean("time", skipna=not complete)


def _tendency(series):
    """Backward month-to-month difference; the leading step is NaN (no prior month) by construction."""
    return series - series.shift(time=1)


def _slope(series, dim):
    """OLS slope of `series` against its `dim` coordinate, per realization — per one step of that
    coordinate. NaN entries are skipped: masking x to where the value exists holds the numerator and
    denominator on one valid set, so a hole can't bias the fit. Range selection is the caller's job.
    """
    x = series[dim].astype("float64").where(series.notnull())
    xc = x - x.mean(dim)
    yc = series - series.mean(dim)
    return (xc * yc).sum(dim) / (xc * xc).sum(dim)


def ohca(primitives, window):
    return _annual(_anomaly(primitives["integral"], window))


def ohu(primitives, window):
    return _annual(_tendency(primitives["integral"]), complete=True).assign_attrs(per="month")


def ohca_trend(primitives, window):
    annual = _in_window(_annual(primitives["integral"]), "year", window)
    return _slope(annual, "year").assign_attrs(per="year")


def ohu_trend(primitives, window):
    annual = _in_window(_annual(_tendency(primitives["integral"]), complete=True), "year", window)
    return _slope(annual, "year").assign_attrs(per="year")


def gridded_anomaly(primitives, window):
    return _anomaly(primitives["map"], window)


def field(primitives, window):
    """The masked field as published, per cell and month — no transform, so the collapse yields the
    per-cell mean and member spread. The window is unused."""
    return primitives["map"]


# name -> (recipe, the step-3 primitive it draws on)
REGISTRY = {
    "ohca": (ohca, "integral"),
    "ohu": (ohu, "integral"),
    "ohca_trend": (ohca_trend, "integral"),
    "ohu_trend": (ohu_trend, "integral"),
    "map": (gridded_anomaly, "map"),
    "field": (field, "map"),
}


def _build(name, primitives, window, field_units):
    recipe, primitive = REGISTRY[name]
    out = recipe(primitives, window)
    return out.assign_attrs(field_units=field_units, reduction=REDUCTION[primitive])


def apply(names, maps, level, window=None, field_units=""):
    """Build the named deliverables for every constituent, each stamped with `field_units` (the
    submission's published units) and `reduction` (the primitive it draws on).

    -> {tag: {quantity_name: DataArray(realization, ...)}}
    """
    for name in names:
        if name not in REGISTRY:
            raise SystemExit("unknown quantity %r; known: %s" % (name, list(REGISTRY)))
    return {tag: {name: _build(name, primitives, window, field_units) for name in names}
            for tag, primitives in maps.items()}
