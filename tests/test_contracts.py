"""contracts: a bare ME4OH file loads under --contract ME4OH with the protocol quantity, the layer from
the filename, no members, and an inferred_config record; without the flag it is refused."""
import json
import types

import numpy as np
import pytest
import xarray as xr

import contracts
import levels
import loader
import masks
import run
import conftest


def _write_me4oh(path):
    """A protocol-only submission: DATA(LONGITUDE, LATITUDE, TIME) and nothing derive-specific."""
    time = conftest.months(12)
    data = np.full((conftest.NLON, conftest.NLAT, 12), 2.5)
    ds = xr.Dataset({"DATA": (("LONGITUDE", "LATITUDE", "TIME"), data)},
                    coords={"LONGITUDE": conftest.LON, "LATITUDE": conftest.LAT, "TIME": time})
    ds["DATA"].attrs["units"] = "TJ/m^2"
    ds.attrs["experiment"] = "B"
    ds.to_netcdf(path)


def test_me4oh_layer_from_filename():
    assert contracts.me4oh_layer("/x/OHC_2004_2025_lev0_2000_expB.nc") == "0_2000"
    assert contracts.me4oh_layer("OHC_2004_2025_lev15_300_expC_SomeGroup.nc") == "15_300"
    assert contracts.me4oh_layer("OHC_T_2004_2025_lev15_20_expC.nc") == "15_20"     # ours too
    with pytest.raises(SystemExit):
        contracts.me4oh_layer("OHC_2004_2025_expB.nc")
    with pytest.raises(SystemExit):
        contracts.infer("bogus", "OHC_2004_2025_lev0_2000_expB.nc")


def test_bare_me4oh_file_needs_the_contract(tmp_path):
    p = str(tmp_path / "OHC_2004_2025_lev0_2000_expB_group.nc")
    _write_me4oh(p)
    with pytest.raises(SystemExit):
        loader.load_submissions([p], with_members=False)


def test_bare_me4oh_file_loads_under_the_contract(tmp_path):
    p = str(tmp_path / "OHC_2004_2025_lev0_2000_expB_group.nc")
    _write_me4oh(p)
    subs = loader.load_submissions([p], with_members=True, contract="ME4OH")   # members: none, silently
    s = subs["0_2000"]
    assert s["field_value"].sizes["realization"] == 1
    assert s["quantity"] == contracts.ME4OH_QUANTITY
    assert s["inferred"]["contract"] == "ME4OH" and s["inferred"]["ensemble"] is False
    assert s["inferred"]["source_file"].startswith("OHC_2004_2025_lev0_2000")


def test_run_under_contract_builds_ohca_and_records_inferred_config(tmp_path):
    p = str(tmp_path / "OHC_2004_2025_lev0_2000_expB_group.nc")
    _write_me4oh(p)
    cfg = types.SimpleNamespace(
        submissions=[p], level=None, bathy=None, quantities=["ohca"], mask=None, require_top=None,
        time_window=None, no_ensemble=False, contract="ME4OH", tag="T", provenance_link=None,
        code_version="http://c", product_name="p", author="a", citation="c", out=str(tmp_path))
    if cfg.contract is not None:
        cfg.no_ensemble = True
    subs = loader.load_submissions(cfg.submissions, with_members=not cfg.no_ensemble, contract=cfg.contract)
    level = levels.resolve(cfg.level, subs)
    bathy = conftest.bathy(np.full((conftest.NLAT, conftest.NLON), 4000.0))
    blob = run.run_level(level, subs, bathy, cfg)
    loader.stamp_chain_provenance(blob, level, cfg, subs)

    assert cfg.mask == "contiguous_from_top"
    assert "ohca" in blob and "ohca_sd" not in blob                  # mean-only under the protocol
    assert np.allclose(blob["ohca"].values, 0.0)                    # constant field -> zero anomaly
    assert json.loads(blob.attrs["quantity"])["scale_terms"]["cp0"] == 3989.244
    facts = json.loads(blob.attrs["ohc_derive_run_facts"])
    assert facts["inferred_config"]["0_2000"]["contract"] == "ME4OH"
    assert facts["ensemble"] is False
