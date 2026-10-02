"""temporal_transforms: the helpers, the deliverable recipes, and the integral<->anomaly invariant."""
import numpy as np
import pytest
import xarray as xr

import temporal_transforms as T
import map_transforms
import conftest


# --- helpers ---------------------------------------------------------------

def test_in_window_selects_years_on_the_monthly_axis():
    s = conftest.series(np.zeros((1, 36)), start_year=2001)          # monthly, 2001-2003
    assert T._in_window(s, "time", ("2002", "2002")).sizes["time"] == 12
    assert T._in_window(s, "time", ("2002", "2003")).sizes["time"] == 24
    assert T._in_window(s, "time", None).sizes["time"] == 36


def test_anomaly_of_constant_is_zero():
    s = conftest.series(np.full((1, 12), 7.0))
    assert np.allclose(T._anomaly(s, None).values, 0.0)


def test_anomaly_removes_the_between_member_level_spread():
    # two members share a shape (a linear ramp) but sit at different levels
    ramp = np.arange(12.0)
    s = conftest.series(np.stack([5.0 + ramp, 20.0 + ramp]))
    raw_spread = s.std("realization")
    anom_spread = T._anomaly(s, None).std("realization")
    assert float(raw_spread.mean()) > 1.0                           # the levels differ
    assert np.allclose(anom_spread.values, 0.0)                     # after demeaning, identical


def test_anomaly_uses_only_the_window_months():
    # year 2001 sits at 0, year 2002 at 10; baseline over 2001 only -> subtract 0
    s = conftest.series(np.concatenate([np.zeros(12), np.full(12, 10.0)])[None, :], start_year=2001)
    out = T._anomaly(s, (2001, 2001))
    assert np.allclose(out.isel(time=slice(0, 12)).values, 0.0)
    assert np.allclose(out.isel(time=slice(12, 24)).values, 10.0)


def test_annual_is_calendar_year_mean():
    vals = np.concatenate([np.full(12, 2.0), np.full(12, 8.0)])     # year1=2, year2=8
    out = T._annual(conftest.series(vals[None, :], start_year=2001))
    assert list(out["year"].values) == [2001, 2002]
    assert np.allclose(out.values, [[2.0, 8.0]])


def test_tendency_is_backward_difference_with_nan_first():
    s = conftest.series(np.array([[1.0, 3.0, 6.0, 10.0]]))
    out = T._tendency(s).values
    assert np.isnan(out[0, 0])
    assert np.allclose(out[0, 1:], [2.0, 3.0, 4.0])


def _annual_da(vals):
    years = np.arange(2001, 2001 + len(vals))
    return xr.DataArray(np.asarray(vals, float)[None, :], dims=("realization", "year"),
                        coords={"realization": [0], "year": years})


def test_in_window_selects_years_on_the_annual_axis():
    annual = _annual_da([1.0, 2.0, 3.0, 4.0])                        # years 2001..2004
    clipped = T._in_window(annual, "year", (2002, 2003))
    assert list(clipped["year"].values) == [2002, 2003]


def test_slope_recovers_a_linear_trend():
    slope = T._slope(_annual_da(3.0 + 2.0 * np.arange(10)), "year")
    assert np.isclose(slope.item(), 2.0)


def test_slope_skips_a_nan_leading_year():
    # a slope-1 line with the leading year voided; the fit must recover 1, not the ~0.25 the old
    # code gave when the denominator counted a year the numerator dropped
    annual = _annual_da([np.nan, 10.0, 11.0])
    assert np.isclose(T._slope(annual, "year").item(), 1.0)


def test_windowing_then_slope_fits_the_selected_years():
    years = np.arange(2001, 2011)
    vals = np.where(years <= 2004, 5.0, 5.0 + (years - 2004) * 2.0)  # flat, then a ramp
    annual = _annual_da(vals)
    flat = T._in_window(annual, "year", (2001, 2004))
    assert np.isclose(T._slope(flat, "year").item(), 0.0)           # flat stretch only
    assert abs(T._slope(annual, "year").item()) > 0.1              # whole record: ramp pulls the fit


# --- recipes ---------------------------------------------------------------

def _primitives():
    integ = conftest.series(np.stack([np.arange(24.0), 100.0 + np.arange(24.0)]), start_year=2001)
    gridded = conftest.const_field(3.0, n_real=2, n_time=24)
    return {"integral": integ, "map": gridded}


def test_ohca_is_annual_of_anomaly():
    p = _primitives()
    got = T.ohca(p, None)
    expect = T._annual(T._anomaly(p["integral"], None))
    assert got.dims == ("realization", "year")
    assert np.allclose(got.values, expect.values)


def test_ohu_voids_its_leading_year():
    p = _primitives()
    got = T.ohu(p, None)
    assert got.dims == ("realization", "year")
    assert np.allclose(got.values, T._annual(T._tendency(p["integral"]), complete=True).values,
                       equal_nan=True)
    assert bool(got.isel(year=0).isnull().all())                   # the dropped leading step voids year 0


def test_trends_reduce_to_one_value_per_realization():
    p = _primitives()
    assert T.ohca_trend(p, None).dims == ("realization",)
    assert T.ohu_trend(p, None).dims == ("realization",)


def test_trends_carry_their_step_as_per():
    p = _primitives()
    assert T.ohca_trend(p, None).attrs["per"] == "year"
    assert T.ohu_trend(p, None).attrs["per"] == "year"


def test_gridded_anomaly_keeps_the_grid_and_demeans_per_cell():
    p = _primitives()
    out = T.gridded_anomaly(p, None)
    assert out.dims == ("realization", "time", "lat", "lon")
    assert np.allclose(out.mean("time").values, 0.0)               # constant field -> zero anomaly


def test_field_is_the_map_primitive_untouched():
    prim = _primitives()
    out = T.field(prim, window=(2001, 2001))
    assert out.dims == ("realization", "time", "lat", "lon")
    assert np.array_equal(out.values, prim["map"].values, equal_nan=True)    # identity, window ignored
    stamped = T.apply(["field"], {"15_20": prim}, level=None, window=None, field_units="m")["15_20"]["field"]
    assert stamped.attrs["reduction"] == "grid" and stamped.attrs["field_units"] == "m"


def test_apply_builds_each_named_quantity_per_constituent():
    maps = {"15_20": _primitives(), "15_300": _primitives()}
    out = T.apply(["ohca", "ohca_trend"], maps, level=None, window=None)
    assert set(out) == {"15_20", "15_300"}
    assert set(out["15_20"]) == {"ohca", "ohca_trend"}


def test_apply_stamps_field_units_and_reduction():
    maps = {"15_20": _primitives()}
    out = T.apply(["ohca", "ohu", "ohca_trend", "map"], maps, level=None, window=None,
                  field_units="TJ/m^2")["15_20"]
    for name in ("ohca", "ohu", "ohca_trend"):                     # integral-based -> area_integral
        assert out[name].attrs["reduction"] == "area_integral"
        assert out[name].attrs["field_units"] == "TJ/m^2"
    assert out["map"].attrs["reduction"] == "grid"                   # stays on the grid
    assert out["map"].attrs["field_units"] == "TJ/m^2"
    assert out["ohca_trend"].attrs["per"] == "year"                  # the recipe's own attrs survive
    assert out["ohu"].attrs["per"] == "month"


def test_apply_rejects_unknown_quantity():
    with pytest.raises(SystemExit):
        T.apply(["bogus"], {"15_20": _primitives()}, level=None, window=None)


# --- invariant: integrate-then-demean == demean-per-cell-then-integrate ----

def test_integral_and_anomaly_commute():
    rng = np.random.default_rng(0)
    fv = conftest.field(rng.normal(size=(2, 12, conftest.NLAT, conftest.NLON)), conftest.months(12))
    left = T._anomaly(map_transforms.integral(fv), None)           # integrate, then demean the series
    right = map_transforms.integral(fv - fv.mean("time"))          # demean each cell, then integrate
    assert np.allclose(left.values, right.values)
