"""loader: sibling resolution, the realization stack, submission/bathy reads, and the blob write."""
import json
import types

import numpy as np
import pytest
import xarray as xr

import loader
import levels
import conftest


def _write_ohc(path, tag, name="ohc"):
    time = conftest.months(3)
    data = np.arange(conftest.NLON * conftest.NLAT * 3, dtype="float64").reshape(
        conftest.NLON, conftest.NLAT, 3)
    ds = xr.Dataset({"DATA": (("LONGITUDE", "LATITUDE", "TIME"), data)},
                    coords={"LONGITUDE": conftest.LON, "LATITUDE": conftest.LAT, "TIME": time})
    ds.attrs.update({"mapped_layer": tag, "source": "test",
                     "quantity": json.dumps({"name": name, "scale_terms": {"cp0": 3989.0, "rho0": 1030.0}})})
    ds.to_netcdf(path)


def _write_ohcens(path, n_member=3):
    time = conftest.months(3)
    data = np.zeros((n_member, conftest.NLON, conftest.NLAT, 3), dtype="float64")
    for m in range(n_member):
        data[m] = m + 1
    ds = xr.Dataset({"DATA": (("MEMBER", "LONGITUDE", "LATITUDE", "TIME"), data)},
                    coords={"LONGITUDE": conftest.LON, "LATITUDE": conftest.LAT, "TIME": time})
    ds.to_netcdf(path)


def test_member_sibling():
    # publish names the pair <NAME>_ / <NAME>ENS_ from the quantity name; the loader does the same
    assert loader._member_sibling("/x/OHC_a.nc", "ohc") == "/x/OHCENS_a.nc"
    assert loader._member_sibling("/x/MLD_2004_2025_lev0_0.nc", "mld") == "/x/MLDENS_2004_2025_lev0_0.nc"
    with pytest.raises(SystemExit):                                  # prefix must match the quantity
        loader._member_sibling("/x/OHC_a.nc", "mld")
    with pytest.raises(SystemExit):
        loader._member_sibling("/x/bar.nc", "ohc")


def test_load_submissions_members_need_the_quantity_attr(tmp_path):
    # a submission without the attr can't resolve its sibling
    ds = xr.Dataset({"DATA": (("LONGITUDE", "LATITUDE", "TIME"), np.zeros((2, 2, 1)))})
    ds.attrs["mapped_layer"] = "15_20"
    ds.to_netcdf(str(tmp_path / "OHC_x.nc"))
    with pytest.raises(SystemExit):
        loader.load_submissions([str(tmp_path / "OHC_x.nc")], with_members=True)


def test_load_submissions_with_members_any_prefix(tmp_path):
    _write_ohc(str(tmp_path / "MLD_x.nc"), tag="0_0", name="mld")
    _write_ohcens(str(tmp_path / "MLDENS_x.nc"), n_member=2)
    subs = loader.load_submissions([str(tmp_path / "MLD_x.nc")], with_members=True)
    assert subs["0_0"]["field_value"].sizes["realization"] == 3


def test_stack_mean_only():
    mean = conftest.const_field(5.0).isel(realization=0)             # (time, lat, lon)
    out = loader._stack(mean, None)
    assert out.sizes["realization"] == 1
    assert np.allclose(out.isel(realization=0).values, 5.0)


def test_stack_puts_mean_first_then_members():
    mean = conftest.const_field(5.0).isel(realization=0)
    members = conftest.const_field(9.0, n_real=3)                   # (realization, time, lat, lon)
    out = loader._stack(mean, members)
    assert list(out["realization"].values) == [0, 1, 2, 3]
    assert np.allclose(out.isel(realization=0).values, 5.0)         # mean field
    assert np.allclose(out.isel(realization=slice(1, None)).values, 9.0)


def test_load_submissions_mean_only(tmp_path):
    _write_ohc(str(tmp_path / "OHC_x.nc"), tag="15_20")
    subs = loader.load_submissions([str(tmp_path / "OHC_x.nc")], with_members=False)
    assert "15_20" in subs
    fv = subs["15_20"]["field_value"]
    assert fv.dims == ("realization", "time", "lat", "lon")
    assert fv.sizes["realization"] == 1
    assert subs["15_20"]["attrs"]["source"] == "test"              # the submission's attrs ride along
    assert subs["15_20"]["quantity"]["scale_terms"]["cp0"] == 3989.0


def test_load_submissions_with_members(tmp_path):
    _write_ohc(str(tmp_path / "OHC_x.nc"), tag="15_20")
    _write_ohcens(str(tmp_path / "OHCENS_x.nc"), n_member=3)
    subs = loader.load_submissions([str(tmp_path / "OHC_x.nc")], with_members=True)
    fv = subs["15_20"]["field_value"]
    assert fv.sizes["realization"] == 4                            # mean + 3 members
    assert np.allclose(fv.isel(realization=1).values, 1.0)         # member 1's constant


def test_load_submissions_missing_member_sibling_exits(tmp_path):
    _write_ohc(str(tmp_path / "OHC_x.nc"), tag="15_20")            # no OHCENS_ sibling written
    with pytest.raises(SystemExit):
        loader.load_submissions([str(tmp_path / "OHC_x.nc")], with_members=True)


def test_load_bathy_reads_etopo_rose(tmp_path):
    # etopo60.cdf layout: ROSE relief (metres, negative below sea level) on ETOPO60Y/ETOPO60X [lat, lon].
    p = str(tmp_path / "etopo.nc")
    xr.Dataset({"ROSE": (("ETOPO60Y", "ETOPO60X"),
                         [[-1000.0, -2000.0, -3000.0], [-1000.0, -1000.0, -1000.0]])},
               coords={"ETOPO60Y": conftest.LAT, "ETOPO60X": conftest.LON}).to_netcdf(p)
    b = loader.load_bathy(p)
    assert b.dims == ("lat", "lon")                                # Y -> lat, X -> lon
    assert float(b.isel(lat=0, lon=1)) == 2000.0                   # relief negated to positive depth


def test_write_blob_round_trip(tmp_path):
    cfg = types.SimpleNamespace(out=str(tmp_path), tag="TESTTAG", provenance_link="http://prov",
                                time_window=(2005, 2024))
    blob = xr.Dataset({"ohca": ("year", [1.0, 2.0])}, coords={"year": [2001, 2002]})
    path = loader.write_blob(blob, levels.get("0_300"), cfg)
    back = xr.open_dataset(path)
    assert back.attrs["level"] == "0_300"
    assert back.attrs["time_window"] == "2005-2024"
    assert back.attrs["provenance_tag"] == "TESTTAG"
    assert back.attrs["provenance_link"] == "http://prov"


def test_write_blob_window_all_when_none(tmp_path):
    cfg = types.SimpleNamespace(out=str(tmp_path), tag="T", provenance_link=None, time_window=None)
    blob = xr.Dataset({"ohca": ("year", [1.0])}, coords={"year": [2001]})
    back = xr.open_dataset(loader.write_blob(blob, levels.get("0_300"), cfg))
    assert back.attrs["time_window"] == "all"
