"""combine: the realization split (step 5) and the n_fac folds + stamped geometry (step 6)."""
import json

import numpy as np
import xarray as xr

import combine
import conftest
import levels


def _quantity(mean, members):
    """A (realization, year) quantity: realization 0 = mean field, the rest = members."""
    vals = np.array([mean] + list(members), dtype="float64")[:, None]
    return xr.DataArray(vals, dims=("realization", "year"),
                        coords={"realization": np.arange(len(vals)), "year": [2001]})


def test_split_value_from_mean_and_sd_from_members():
    out = combine._split(_quantity(100.0, [1.0, 2.0, 3.0]))
    assert np.allclose(out["value"].values, 100.0)                  # the mean field, not the member mean
    assert np.allclose(out["sd"].values, np.std([1.0, 2.0, 3.0], ddof=1))
    assert "realization" not in out["value"].coords                # leftover coord dropped


def test_split_mean_only_has_no_sd():
    out = combine._split(_quantity(100.0, []))
    assert np.allclose(out["value"].values, 100.0)
    assert out["sd"] is None


def test_collapse_maps_over_the_nested_dict():
    q = _quantity(100.0, [1.0, 2.0, 3.0])
    out = combine.collapse_sd({"15_20": {"ohca": q}, "15_300": {"ohca": q}})
    assert set(out) == {"15_20", "15_300"}
    assert set(out["15_20"]["ohca"]) == {"value", "sd"}


def _per(v1520, sd1520, v15300, sd15300):
    return {"15_20": {"ohca": {"value": xr.DataArray(v1520), "sd": _maybe(sd1520)}},
            "15_300": {"ohca": {"value": xr.DataArray(v15300), "sd": _maybe(sd15300)}}}


def _maybe(x):
    return None if x is None else xr.DataArray(x)


def test_nfac_sum_value_is_linear_and_sd_is_worst_case():
    per = _per(2.0, 0.5, 7.0, 2.0)
    contribs = levels.get("0_300").contributors                    # 15_20 x3, 15_300 x1
    assert np.isclose(float(combine._nfac_sum(per, contribs, "ohca", "value")), 3 * 2 + 7)      # 13
    assert np.isclose(float(combine._nfac_sum(per, contribs, "ohca", "sd")), 3 * 0.5 + 2)       # 3.5


def test_nfac_sum_none_when_a_part_is_none():
    per = _per(2.0, None, 7.0, None)
    contribs = levels.get("0_300").contributors
    assert combine._nfac_sum(per, contribs, "ohca", "sd") is None


def test_combine_synthetic_builds_dataset_and_stamps_geometry():
    blob = combine.combine_synthetic(_per(2.0, 0.5, 7.0, 2.0), levels.get("0_300"),
                                     area_m2=1000.0, volume_m3=250000.0,
                                     quantity=conftest.OHC_QUANTITY)
    assert np.isclose(float(blob["ohca"]), 13.0)
    assert np.isclose(float(blob["ohca_sd"]), 3.5)
    assert blob.attrs["level"] == "0_300"
    assert np.isclose(blob.attrs["area_m2"], 1000.0)
    assert np.isclose(blob.attrs["volume_m3"], 250000.0)          # passed straight through from the mask
    assert np.isclose(blob.attrs["cp0"], 3989.0)                  # legacy copy of the scale term
    assert json.loads(blob.attrs["quantity"])["name"] == "ohc"


def test_combine_synthetic_mean_only_omits_sd():
    blob = combine.combine_synthetic(_per(2.0, None, 7.0, None), levels.get("0_300"), 1000.0, 250000.0, {})
    assert "ohca" in blob.data_vars and "ohca_sd" not in blob.data_vars


def test_combine_carries_quantity_attrs():
    # a trend's `per` attr must survive the n_fac fold so packaging can read it off the blob
    per = {"15_20": {"ohca_trend": {"value": xr.DataArray(2.0, attrs={"per": "year"}), "sd": None}},
           "15_300": {"ohca_trend": {"value": xr.DataArray(7.0, attrs={"per": "year"}), "sd": None}}}
    blob = combine.combine_synthetic(per, levels.get("0_300"), 1000.0, 250000.0, {})
    assert blob["ohca_trend"].attrs["per"] == "year"


def _grid(vals):
    """A (lat, lon) gridded quantity; np.nan marks a not-kept constituent."""
    return xr.DataArray(np.array(vals, dtype="float64"), dims=("lat", "lon"),
                        coords={"lat": [0.0, 1.0], "lon": [0.0, 1.0]})


def _per_grid(m1520, m15300):
    return {"15_20": {"map": {"value": _grid(m1520), "sd": None}},
            "15_300": {"map": {"value": _grid(m15300), "sd": None}}}


def test_gridded_is_detected_by_lat_lon_dims():
    assert combine._gridded(_grid([[1.0, 1.0], [1.0, 1.0]]))
    assert not combine._gridded(xr.DataArray([1.0], dims=("year",)))


def test_gridded_combine_uses_union_footprint():
    # cell [0,0]: both kept -> 3*1 + 1*2 = 5; cell [0,1]: deep (15_300) not kept -> counts as 0, 3*1 = 3;
    # cell [1,0]: shallow not kept, deep kept -> 1*4 = 4; cell [1,1]: both absent -> NaN
    per = _per_grid([[1.0, 1.0], [np.nan, np.nan]],
                    [[2.0, np.nan], [4.0, np.nan]])
    contribs = levels.get("0_300").contributors                    # 15_20 x3, 15_300 x1
    out = combine._nfac_sum(per, contribs, "map", "value")
    assert np.allclose(out.values[0], [5.0, 3.0])                  # missing constituent contributes 0
    assert np.isclose(out.values[1, 0], 4.0)
    assert np.isnan(out.values[1, 1])                             # NaN only where every constituent absent


def test_series_combine_still_propagates_deliberate_nan():
    # a non-gridded quantity keeps the plain sum, so a deliberate NaN (e.g. OHU's first year) survives
    per = {"15_20": {"ohu": {"value": xr.DataArray([np.nan, 2.0], dims=("year",)), "sd": None}},
           "15_300": {"ohu": {"value": xr.DataArray([np.nan, 7.0], dims=("year",)), "sd": None}}}
    contribs = levels.get("0_300").contributors
    out = combine._nfac_sum(per, contribs, "ohu", "value")
    assert np.isnan(out.values[0])                               # not scrubbed to 0
    assert np.isclose(out.values[1], 3 * 2.0 + 7.0)
