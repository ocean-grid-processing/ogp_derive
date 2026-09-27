"""run: window parsing, constant extraction, and an end-to-end run_level with a known OHCA."""
import json
import types

import numpy as np
import xarray as xr

import run
import levels
import grid
import loader
import conftest


def test_file_token_carries_data_span_and_baseline():
    blob = xr.Dataset({"ohca": ("year", [1.0, 2.0])}, coords={"year": [2004, 2025]})
    # data span from the blob's axis; baseline from --time-window, defaulting to the whole data span
    assert loader._file_token(types.SimpleNamespace(time_window=(2005, 2024)), blob) == "2004_2025_tw2005_2024"
    assert loader._file_token(types.SimpleNamespace(time_window=None), blob) == "2004_2025_tw2004_2025"


def test_load_submissions_rejects_duplicate_native_level(tmp_path):
    import pytest

    def _write(path, tag):
        ds = xr.Dataset({"DATA": (("LONGITUDE", "LATITUDE", "TIME"), np.zeros((2, 2, 1)))})
        ds.attrs["mapped_layer"] = tag
        ds.attrs["quantity"] = json.dumps(conftest.OHC_QUANTITY)
        ds.to_netcdf(path)

    _write(str(tmp_path / "a.nc"), "15_20")
    _write(str(tmp_path / "b.nc"), "15_20")                 # same native level -> collision
    with pytest.raises(SystemExit):
        loader.load_submissions([str(tmp_path / "a.nc"), str(tmp_path / "b.nc")], with_members=False)

    _write(str(tmp_path / "c.nc"), "15_300")               # distinct levels load fine
    subs = loader.load_submissions([str(tmp_path / "a.nc"), str(tmp_path / "c.nc")], with_members=False)
    assert set(subs) == {"15_20", "15_300"}


def test_parse_window():
    assert run._parse_window("2005:2024") == (2005, 2024)
    assert run._parse_window("2005-2024") == (2005, 2024)
    assert run._parse_window("2005_2024") == (2005, 2024)     # filename-token form
    assert run._parse_window(None) is None
    assert run._parse_window("") is None


def test_quantity_is_shared_across_constituents():
    q = conftest.OHC_QUANTITY
    subs = {"15_20": {"quantity": q}, "15_300": {"quantity": dict(q)}}
    assert run._quantity(subs, conftest.level("0_300")) == q


def _kind(kind):
    return dict(conftest.OHC_QUANTITY, name="mld" if kind == "intensive" else "ohc", kind=kind)


def test_check_kind_accepts_the_ohc_plan():
    run.check_kind(_kind("extensive"), ["ohca", "ohu", "ohca_trend", "ohu_trend", "map"],
                   "contiguous_from_top", conftest.level("0_2000"))
    run.check_kind(_kind("extensive"), ["ohca"], "fully_wet_nan", conftest.level("0_300"))


def test_check_kind_refuses_integral_quantities_for_an_intensive_field():
    import pytest
    one = levels.Level("15_20", (levels.Contributor("15_20", 1, 15, 20),), 5)   # a single constituent
    run.check_kind(_kind("intensive"), ["map"], "contiguous_from_top", one)     # per-cell: fine
    with pytest.raises(SystemExit) as e:
        run.check_kind(_kind("intensive"), ["ohca", "map"], "contiguous_from_top", one)
    assert "ohca" in str(e.value) and "integral" in str(e.value)


def test_check_kind_refuses_zero_fill_mask_and_multi_constituent_for_intensive():
    import pytest
    one = levels.Level("15_20", (levels.Contributor("15_20", 1, 15, 20),), 5)
    with pytest.raises(SystemExit) as e:
        run.check_kind(_kind("intensive"), ["map"], "fully_wet_nan", one)
    assert "fully_wet_nan" in str(e.value)
    with pytest.raises(SystemExit) as e:
        run.check_kind(_kind("intensive"), ["map"], "contiguous_from_top", conftest.level("0_300"))
    assert "thickness-weighted" in str(e.value)


def test_check_kind_rejects_unknown_names_and_kinds():
    import pytest
    with pytest.raises(SystemExit):
        run.check_kind(_kind("extensive"), ["bogus"], "contiguous_from_top", conftest.level("0_300"))
    with pytest.raises(SystemExit):
        run.check_kind(_kind("extensive"), ["ohca"], "bogus", conftest.level("0_300"))
    with pytest.raises(SystemExit):
        run.check_kind(dict(conftest.OHC_QUANTITY, kind="sideways"), ["ohca"], "contiguous_from_top",
                       conftest.level("0_300"))


def test_run_level_identity_field_returns_the_input_and_its_spread(tmp_path):
    # one submission, no combine: `field` comes back as the input field; `_sd` is the member spread.
    n_time = 12
    arr = np.empty((4, n_time, conftest.NLAT, conftest.NLON))
    for r, v in enumerate((10.0, 9.0, 10.0, 11.0)):              # realization 0 = mean 10; members 9/10/11 -> sd 1
        arr[r] = v
    field = conftest.field(arr, conftest.months(n_time, start_year=2001))
    lv = levels.identity("15_20")
    subs = {"15_20": {"field_value": field, "attrs": {}, "quantity": conftest.OHC_QUANTITY}}
    cfg = types.SimpleNamespace(quantities=["field"], mask="as_published", require_top=None,
                                time_window=None, out=str(tmp_path), tag="t", no_ensemble=False,
                                product_name="p", author="a", citation="c")
    blob = run.run_level(lv, subs, conftest.bathy(np.full((conftest.NLAT, conftest.NLON), 4000.0)), cfg)
    assert blob.attrs["level"] == "15_20"
    assert "volume_m3" not in blob.attrs and blob.attrs["area_m2"] > 0
    assert np.allclose(blob["field"].values, 10.0)
    assert np.allclose(blob["field_sd"].values, 1.0)
    assert blob["field"].dims == ("time", "lat", "lon")


def test_quantity_disagreement_is_an_error():
    import pytest
    other = dict(conftest.OHC_QUANTITY, name="mld")
    subs = {"15_20": {"quantity": conftest.OHC_QUANTITY}, "15_300": {"quantity": other}}
    with pytest.raises(SystemExit):
        run._quantity(subs, conftest.level("0_300"))


def _ramped(slope, n_time=24):
    """A field uniform in space, equal to slope * month_index in time -> (1, time, lat, lon)."""
    t = slope * np.arange(float(n_time))
    arr = t[None, :, None, None] * np.ones((1, n_time, conftest.NLAT, conftest.NLON))
    return conftest.field(arr, conftest.months(n_time, start_year=2001))


def test_run_level_ohca_matches_hand_computed(tmp_path):
    # 15_20 = t, 15_300 = 2t (uniform in space) over 24 months (2001-2002); deep bathy -> all fully wet.
    # per-cell integral_i(t) = value_i * A; anomaly removes the whole-record mean (month 11.5);
    # annual means over each year give ohca_15_20 = [-6A, 6A], ohca_15_300 = [-12A, 12A];
    # n_fac combine (3, 1) -> [-30A, 30A].
    subs = {
        "15_20": {"field_value": _ramped(1.0), "attrs": {}, "quantity": conftest.OHC_QUANTITY},
        "15_300": {"field_value": _ramped(2.0), "attrs": {}, "quantity": conftest.OHC_QUANTITY},
    }
    reference_bathy = conftest.bathy([[1000.0, 1000.0, 1000.0], [1000.0, 1000.0, 1000.0]])
    cfg = types.SimpleNamespace(mask="fully_wet_nan", quantities=["ohca"], time_window=None,
                                require_top=None, tag="dev", out=str(tmp_path))

    blob = run.run_level(conftest.level("0_300"), subs, reference_bathy, cfg)

    A = float(grid.cell_area(conftest.LAT, conftest.LON).sum())
    assert np.allclose(blob["ohca"].values, [-30.0 * A, 30.0 * A])
    assert "ohca_sd" not in blob.data_vars                          # mean-only run
    assert np.isclose(blob.attrs["area_m2"], A)
    assert np.isclose(blob.attrs["volume_m3"], A * 300)            # nominal thickness of 0_300


def test_run_level_with_members_produces_sd_and_geometry(tmp_path):
    subs = {
        "15_20": {"field_value": conftest.const_field(1.0, n_real=4, n_time=12),
                  "attrs": {}, "quantity": conftest.OHC_QUANTITY},
        "15_300": {"field_value": conftest.const_field(10.0, n_real=4, n_time=12),
                   "attrs": {}, "quantity": conftest.OHC_QUANTITY},
    }
    reference_bathy = conftest.bathy([[1000.0, 1000.0, 1000.0], [1000.0, 1000.0, 1000.0]])
    cfg = types.SimpleNamespace(mask="fully_wet_nan", quantities=["ohca"], time_window=None,
                                require_top=None, tag="dev", out=str(tmp_path))
    blob = run.run_level(conftest.level("0_300"), subs, reference_bathy, cfg)
    A = float(grid.cell_area(conftest.LAT, conftest.LON).sum())
    assert "ohca" in blob.data_vars and "ohca_sd" in blob.data_vars
    for v in ("ohca", "ohca_sd"):                                    # metadata rides through the combine
        assert blob[v].attrs["reduction"] == "area_integral"
        assert blob[v].attrs["field_units"] == conftest.OHC_QUANTITY["publish_units"]
    assert np.allclose(blob["ohca"].values, 0.0, atol=1e-6 * A)     # constant field -> anomaly ~ 0 (machine precision)
    assert np.allclose(blob["ohca_sd"].values, 0.0)                # identical members -> zero spread
