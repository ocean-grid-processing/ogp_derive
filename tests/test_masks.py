"""masks: fully_wet_nan (wet/dry + drops), contiguous_from_top (top-required + contiguous keep), and
the area/volume that apply derives from the footprint and per-cell height.

Bathy laid out so, for the 0_700 constituents (15_20, 15_300, 300_700), the top row exercises the
fully_wet_nan cases: (0,0) all fully wet, (0,1) 300_700 intersecting, (0,2) 300_700 dry.
"""
import numpy as np
import pytest
import xarray as xr

import masks
import levels
import grid
import conftest

LV = levels.get("0_700")
BATHY = conftest.bathy([[1000.0, 500.0, 100.0],
                        [1000.0, 1000.0, 1000.0]])


def cons(f1520, f15300, f300700):
    return [
        conftest.constituent("15_20", 3, 15, 20, f1520),
        conftest.constituent("15_300", 1, 15, 300, f15300),
        conftest.constituent("300_700", 1, 300, 700, f300700),
    ]


def _default():
    return cons(conftest.const_field(1.0), conftest.const_field(10.0), conftest.const_field(100.0))


# --- fully_wet_nan ---------------------------------------------------------

def test_all_defined_no_exclusions_and_dry_is_zeroed():
    masked, footprint, height = masks.fully_wet_nan(LV, _default(), BATHY)
    assert bool(footprint.all())                                    # 15_20 wet everywhere -> all in footprint
    assert np.allclose(height.values, LV.nominal_thickness)         # slab: full nominal thickness in-footprint
    m300 = masked["300_700"]["field_value"]
    assert np.allclose(m300.isel(lat=0, lon=2).values, 0.0)         # dry -> 0
    assert np.allclose(m300.isel(lat=0, lon=0).values, 100.0)       # fully wet -> value


def test_fully_wet_gap_excludes_the_whole_cell():
    f1520 = conftest.const_field(1.0, nan_cells=[(1, 1)])           # 15_20 fully wet at (1,1) but NaN
    masked, footprint, _ = masks.fully_wet_nan(LV, cons(
        f1520, conftest.const_field(10.0), conftest.const_field(100.0)), BATHY)
    assert not bool(footprint.isel(lat=1, lon=1))
    assert int(footprint.sum()) == 5
    for tag in ("15_20", "15_300", "300_700"):
        assert bool(masked[tag]["field_value"].isel(lat=1, lon=1).isnull().all())


def test_intersecting_nan_contributes_zero_and_does_not_exclude():
    f300 = conftest.const_field(100.0, nan_cells=[(0, 1)])          # 300_700 intersecting at (0,1)
    masked, footprint, _ = masks.fully_wet_nan(LV, cons(
        conftest.const_field(1.0), conftest.const_field(10.0), f300), BATHY)
    assert bool(footprint.all())
    assert np.allclose(masked["300_700"]["field_value"].isel(lat=0, lon=1).values, 0.0)


def test_nan_bathy_is_treated_as_dry():
    b = conftest.bathy([[1000.0, 500.0, np.nan], [1000.0, 1000.0, 1000.0]])
    masked, footprint, _ = masks.fully_wet_nan(LV, _default(), b)
    assert not bool(footprint.isel(lat=0, lon=2))                   # all dry there -> outside the footprint
    for tag in ("15_20", "15_300", "300_700"):
        assert np.allclose(masked[tag]["field_value"].isel(lat=0, lon=2).values, 0.0)


# --- contiguous_from_top ---------------------------------------------------

def test_contiguous_all_kept_is_the_full_nominal_height():
    masked, footprint, height = masks.contiguous_from_top(LV, _default(), BATHY, require_top=300)
    assert bool(footprint.all())                                    # every required-top constituent defined
    assert np.allclose(height.values, 700.0)                        # 3*5 + 285 + 400 = full 0_700
    assert np.allclose(masked["300_700"]["field_value"].values, 100.0)


def test_contiguous_required_top_gap_drops_the_cell():
    f1520 = conftest.const_field(1.0, nan_cells=[(1, 1)])           # 15_20 is in the required top range
    masked, footprint, height = masks.contiguous_from_top(LV, cons(
        f1520, conftest.const_field(10.0), conftest.const_field(100.0)), BATHY, require_top=300)
    assert not bool(footprint.isel(lat=1, lon=1))
    assert np.isclose(float(height.isel(lat=1, lon=1)), 0.0)
    for tag in ("15_20", "15_300", "300_700"):
        assert bool(masked[tag]["field_value"].isel(lat=1, lon=1).isnull().all())


def test_contiguous_deep_gap_truncates_the_column():
    f300 = conftest.const_field(100.0, nan_cells=[(0, 1)])          # gap in the deepest constituent
    masked, footprint, height = masks.contiguous_from_top(LV, cons(
        conftest.const_field(1.0), conftest.const_field(10.0), f300), BATHY, require_top=300)
    assert bool(footprint.isel(lat=0, lon=1))                       # top range fine -> cell survives
    assert np.isclose(float(height.isel(lat=0, lon=1)), 300.0)      # kept 15_20 (15) + 15_300 (285)
    assert bool(masked["300_700"]["field_value"].isel(lat=0, lon=1).isnull().all())
    assert np.allclose(masked["15_20"]["field_value"].isel(lat=0, lon=1).values, 1.0)


def test_contiguous_discards_defined_layers_below_a_gap():
    # 0_2000: gap in 300_700 at (0,0); 700_1850 and 1800_1850 are defined but below it -> discarded
    lv = levels.get("0_2000")
    cons5 = [conftest.constituent("15_20", 3, 15, 20, conftest.const_field(1.0)),
             conftest.constituent("15_300", 1, 15, 300, conftest.const_field(10.0)),
             conftest.constituent("300_700", 1, 300, 700, conftest.const_field(100.0, nan_cells=[(0, 0)])),
             conftest.constituent("700_1850", 1, 700, 1850, conftest.const_field(50.0)),
             conftest.constituent("1800_1850", 3, 1800, 1850, conftest.const_field(70.0))]
    masked, footprint, height = masks.contiguous_from_top(lv, cons5, BATHY, require_top=300)
    assert bool(footprint.isel(lat=0, lon=0))                       # top 0-300 present -> survives
    assert np.isclose(float(height.isel(lat=0, lon=0)), 300.0)      # only 15_20 + 15_300 kept
    for tag in ("300_700", "700_1850", "1800_1850"):
        assert bool(masked[tag]["field_value"].isel(lat=0, lon=0).isnull().all())


def test_contiguous_works_for_a_non_surface_layer():
    # 700_2000: require_top is a thickness from the 700 dbar top, so 300 requires 700-1000 -> 700_1850
    lv = levels.get("700_2000")
    cons2 = [conftest.constituent("700_1850", 1, 700, 1850, conftest.const_field(50.0)),
             conftest.constituent("1800_1850", 3, 1800, 1850, conftest.const_field(70.0))]
    _, footprint, height = masks.contiguous_from_top(lv, cons2, BATHY, require_top=300)
    assert bool(footprint.all())                                    # 700_1850 (covers 700-1000) defined
    assert np.allclose(height.values, 1300.0)                       # 1*1150 + 3*50


def test_700_2000_requires_its_shallowest_constituent():
    # require_depth 1000 falls within 700_1850's bounds (700-1850), so it's required. 700_1850 is a
    # single integrated value per cell, so "required" just means that value must be present (not NaN).
    lv = levels.get("700_2000")
    cons2 = [conftest.constituent("700_1850", 1, 700, 1850, conftest.const_field(50.0, nan_cells=[(0, 0)])),
             conftest.constituent("1800_1850", 3, 1800, 1850, conftest.const_field(70.0))]
    _, footprint, _ = masks.contiguous_from_top(lv, cons2, BATHY, require_top=300)
    assert not bool(footprint.isel(lat=0, lon=0))                   # 700_1850 required in full -> drops
    assert bool(footprint.isel(lat=0, lon=1))                      # elsewhere survives


def _cons5(nan_1520=(), nan_15300=(), nan_300700=()):
    return [conftest.constituent("15_20", 3, 15, 20, conftest.const_field(1.0, nan_cells=nan_1520)),
            conftest.constituent("15_300", 1, 15, 300, conftest.const_field(10.0, nan_cells=nan_15300)),
            conftest.constituent("300_700", 1, 300, 700, conftest.const_field(100.0, nan_cells=nan_300700)),
            conftest.constituent("700_1850", 1, 700, 1850, conftest.const_field(50.0)),
            conftest.constituent("1800_1850", 3, 1800, 1850, conftest.const_field(70.0))]


def test_require_top_past_a_constituent_top_pulls_it_into_the_required_set():
    # 0_2000, 300_700 gap at (0,0). require_top 301 crosses the 300 dbar top -> 300_700 becomes required
    # -> the cell drops; require_top 300 stops just above it -> the cell survives (truncated).
    lv, cons5 = levels.get("0_2000"), _cons5(nan_300700=[(0, 0)])
    _, fp_301, _ = masks.contiguous_from_top(lv, cons5, BATHY, require_top=301)
    assert not bool(fp_301.isel(lat=0, lon=0))
    _, fp_300, _ = masks.contiguous_from_top(lv, _cons5(nan_300700=[(0, 0)]), BATHY, require_top=300)
    assert bool(fp_300.isel(lat=0, lon=0))


def test_require_top_inside_the_surface_zone_requires_only_the_top_constituent():
    # require_top 1 m is above the shallowest measured top (15) but the n_fac surface reaches level.low,
    # so only 15_20 is required: a 15_20 gap drops the cell, a 15_300 gap does not.
    lv = levels.get("0_2000")
    _, fp, _ = masks.contiguous_from_top(lv, _cons5(nan_1520=[(0, 0)], nan_15300=[(0, 1)]),
                                         BATHY, require_top=1)
    assert not bool(fp.isel(lat=0, lon=0))                          # 15_20 (the surface) required
    assert bool(fp.isel(lat=0, lon=1))                             # 15_300 not required at 1 m


def test_require_top_in_the_deep_nfac_zone_requires_the_deepest_constituent():
    # 1851 m lands in 1800_1850's n_fac zone (1850-2000); its top 1800 < 1851 so it is required, and a
    # gap in it drops the cell. 1800 (just above its top) stops short, so the cell survives truncated.
    lv = levels.get("0_2000")

    def deep_gap():
        return [conftest.constituent("15_20", 3, 15, 20, conftest.const_field(1.0)),
                conftest.constituent("15_300", 1, 15, 300, conftest.const_field(10.0)),
                conftest.constituent("300_700", 1, 300, 700, conftest.const_field(100.0)),
                conftest.constituent("700_1850", 1, 700, 1850, conftest.const_field(50.0)),
                conftest.constituent("1800_1850", 3, 1800, 1850,
                                     conftest.const_field(70.0, nan_cells=[(0, 0)]))]

    _, fp_1851, _ = masks.contiguous_from_top(lv, deep_gap(), BATHY, require_top=1851)
    assert not bool(fp_1851.isel(lat=0, lon=0))                     # deepest required -> the gap drops it
    _, fp_1800, _ = masks.contiguous_from_top(lv, deep_gap(), BATHY, require_top=1800)
    assert bool(fp_1800.isel(lat=0, lon=0))                        # deepest not required -> survives


def test_contiguous_needs_require_top():
    with pytest.raises(SystemExit):
        masks.contiguous_from_top(LV, _default(), BATHY, require_top=None)


def test_contiguous_rejects_a_nonpositive_require_top():
    with pytest.raises(SystemExit):
        masks.contiguous_from_top(LV, _default(), BATHY, require_top=0)


def test_contiguous_rejects_require_top_past_the_layer_thickness():
    with pytest.raises(SystemExit):
        masks.contiguous_from_top(LV, _default(), BATHY, require_top=701)   # 0_700 is only 700 m thick


# --- apply: area + volume --------------------------------------------------

def _identity_cons(field_value, tag="15_20"):
    lv = levels.identity(tag)
    c = lv.contributors[0]
    return lv, [{"tag": c.tag, "n_fac": 1, "top": c.top, "bottom": c.bottom, "field_value": field_value}]


def test_as_published_leaves_data_alone_and_footprint_is_the_defined_cells():
    fv = conftest.const_field(7.0, n_real=3, nan_cells=[(0, 1)])    # one cell NaN throughout
    lv, cons_ = _identity_cons(fv)
    masked, footprint, height = masks.as_published(lv, cons_, BATHY)
    assert masked["15_20"]["field_value"] is fv                      # untouched, not even copied
    assert masked["15_20"]["n_fac"] == 1
    assert not bool(footprint.isel(lat=0, lon=1)) and int(footprint.sum()) == footprint.size - 1
    assert height is None


def test_as_published_rejects_a_time_varying_footprint():
    fv = conftest.const_field(7.0, n_real=2)
    arr = fv.values.copy()
    arr[0, 1, 0, 0] = np.nan                                         # one realization, one month
    fv = conftest.field(arr, fv["time"].values)
    lv, cons_ = _identity_cons(fv)
    with pytest.raises(SystemExit) as e:
        masks.as_published(lv, cons_, BATHY)
    assert "1 cell" in str(e.value)


def test_as_published_takes_one_constituent_only():
    with pytest.raises(SystemExit):
        masks.as_published(LV, _default(), BATHY)


def test_default_mask_is_contiguous_from_top():
    assert masks.DEFAULT == "contiguous_from_top"
    assert "as_published" in masks.REGISTRY                          # asked for by name, never defaulted


def test_apply_as_published_area_no_volume_no_coverage(tmp_path):
    lv, cons_ = _identity_cons(conftest.const_field(7.0, nan_cells=[(0, 1)]))
    _, area, volume = masks.apply("as_published", lv, cons_, BATHY, out_dir=str(tmp_path), tag="t")
    a = grid.cell_area(conftest.LAT, conftest.LON)
    assert np.isclose(area, float(a.sum()) - float(a.isel(lat=0, lon=1)))
    assert volume is None
    names = [p.name for p in tmp_path.iterdir()]
    assert not any(n.startswith("coverage_") for n in names)        # no vertical extent, no coverage file
    assert any(n.startswith("mask_") for n in names)                # the footprint png is still drawn


def test_apply_fully_wet_area_and_slab_volume(tmp_path):
    _, area, volume = masks.apply("fully_wet_nan", LV, _default(), BATHY, out_dir=str(tmp_path))
    a = grid.cell_area(conftest.LAT, conftest.LON)
    assert np.isclose(area, float(a.sum()))                         # every cell wet
    assert np.isclose(volume, float(a.sum()) * LV.nominal_thickness)   # slab


def test_apply_contiguous_area_and_tapered_volume(tmp_path):
    f300 = conftest.const_field(100.0, nan_cells=[(0, 1)])          # (0,1) truncates to 300 m, rest 700 m
    _, area, volume = masks.apply("contiguous_from_top", LV, cons(
        conftest.const_field(1.0), conftest.const_field(10.0), f300),
        BATHY, out_dir=str(tmp_path), require_top=300)
    a = grid.cell_area(conftest.LAT, conftest.LON)
    cell = float(a.isel(lat=0, lon=1))
    assert np.isclose(area, float(a.sum()))                         # all six survive
    assert np.isclose(volume, (float(a.sum()) - cell) * 700.0 + cell * 300.0)   # one cell tapered


def test_apply_writes_coverage_diagnostics(tmp_path):
    # 0_2000: 300_700 gap at (0,0) truncates the column to 300 m; (0,1) keeps the full 2000 m.
    lv = levels.get("0_2000")
    cons5 = _cons5(nan_300700=[(0, 0)])
    bathy = conftest.bathy([[2500.0, 3000.0, 3000.0], [3000.0, 3000.0, 3000.0]])
    masks.apply("contiguous_from_top", lv, cons5, bathy, out_dir=str(tmp_path), require_top=300, tag="dev")

    cov = xr.open_dataset(str(tmp_path / "coverage_dev_0_2000_contiguous_from_top.nc"))
    # (0,0): kept 15_20+15_300 = 300 m; bathy 2500 capped at layer bottom 2000 -> uncaptured 2000-300
    assert np.isclose(float(cov["kept_thickness"].isel(lat=0, lon=0)), 300.0)
    assert np.isclose(float(cov["uncaptured_thickness"].isel(lat=0, lon=0)), 1700.0)
    # (0,1): full column kept = 2000 m; bathy 3000 capped at 2000 -> nothing uncaptured
    assert np.isclose(float(cov["kept_thickness"].isel(lat=0, lon=1)), 2000.0)
    assert np.isclose(float(cov["uncaptured_thickness"].isel(lat=0, lon=1)), 0.0)


def test_apply_coverage_carries_publication_meta(tmp_path):
    # product_name/author append to the aux filenames (last before the extension); citation lands as a
    # top-level attr on the coverage .nc.
    lv = levels.get("0_2000")
    cons5 = _cons5(nan_300700=[(0, 0)])
    bathy = conftest.bathy([[2500.0, 3000.0, 3000.0], [3000.0, 3000.0, 3000.0]])
    masks.apply("contiguous_from_top", lv, cons5, bathy, out_dir=str(tmp_path), require_top=300,
                tag="dev", product_name="LocalGP", author="Giglio_etal2026",
                citation="Giglio et al. (2026)")

    cov_path = tmp_path / "coverage_dev_0_2000_contiguous_from_top_LocalGP_Giglio_etal2026.nc"
    assert cov_path.exists()                                          # product_name/author suffixed
    assert xr.open_dataset(str(cov_path)).attrs["citation"] == "Giglio et al. (2026)"


def test_apply_unknown_mask_exits(tmp_path):
    with pytest.raises(SystemExit):
        masks.apply("bogus", LV, _default(), BATHY, out_dir=str(tmp_path))
