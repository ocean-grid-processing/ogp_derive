"""Provenance chain: upstream blocks roll forward grouped by constituent; ohc_derive_* stamped.

stamp_chain_provenance is the fan-in courier — it regroups each upstream `*_run_config` / `_run_facts`
/ `_code_version` under the constituent it came from (never reading the blocks' fields) and adds this
step's own blocks.
"""
import json
import types

import xarray as xr

import levels
import loader
import conftest


def _submission_attrs(tag):
    """A couple of upstream provenance blocks (as the submission carries them) plus non-provenance
    attrs that must NOT be forwarded."""
    top, bottom = tag.split("_")
    return {
        "mapped_layer": tag,                               # non-provenance: must not forward
        "provenance_tag": "260824-OP20260507",             # non-provenance: must not forward
        "cp0": 3989.0,
        "localgp_ingest_run_config": json.dumps({"var_name": "potentialTemperature",
                                                  "dir_mean": "/x/%s/FullField" % tag}),
        "localgp_ingest_run_facts": json.dumps({"layer_top": int(top), "layer_bottom": int(bottom)}),
        "localgp_ingest_code_version": "https://github.com/argovis/ohc_ingest/commit/abc",
        "localgp_publish_run_config": json.dumps({"preset": "wmo_wet"}),
        "localgp_publish_run_facts": json.dumps({"mapped_layer": tag}),
        "localgp_publish_code_version": "https://github.com/argovis/ohc_ingest/commit/abc",
    }


def _cfg(tmp_path):
    return types.SimpleNamespace(
        submissions=["OHC_a.nc", "OHC_b.nc"], level="0_300", bathy="etopo60.cdf",
        quantities=["ohca", "ohu"], mask="contiguous_from_top", require_top=None,
        time_window=None, no_ensemble=False, tag="D-260901", provenance_link="http://derive/docs",
        code_version="https://github.com/argovis/ohc_derive/commit/deadbeef", out=str(tmp_path))


def test_chain_forwards_grouped_and_stamps_own(tmp_path):
    level = conftest.level("0_300")
    contributors = [c.tag for c in level.contributors]     # 15_20, 15_300
    submissions = {t: {"attrs": _submission_attrs(t)} for t in contributors}
    blob = xr.Dataset()
    blob.attrs.update({"area_m2": 1.0e12, "volume_m3": 3.0e14, "cp0": 3989.0, "rho0": 1030.0,
                       "quantity": json.dumps(conftest.OHC_QUANTITY)})

    loader.stamp_chain_provenance(blob, level, _cfg(tmp_path), submissions)

    # every upstream block rolled forward, grouped by the constituent it came from
    for block in ("localgp_ingest_run_config", "localgp_ingest_run_facts", "localgp_ingest_code_version",
                  "localgp_publish_run_config", "localgp_publish_run_facts", "localgp_publish_code_version"):
        grouped = json.loads(blob.attrs[block])
        assert set(grouped) == set(contributors)

    # nested as real objects (parsed, not double-escaped strings), content preserved per constituent
    ig = json.loads(blob.attrs["localgp_ingest_run_facts"])
    assert ig["15_20"]["layer_bottom"] == 20
    assert ig["15_300"]["layer_top"] == 15

    # non-provenance submission attrs did not leak through
    assert "provenance_tag" not in blob.attrs
    assert "mapped_layer" not in blob.attrs

    # this step's own blocks
    assert blob.attrs["ohc_derive_code_version"].endswith("deadbeef")
    facts = json.loads(blob.attrs["ohc_derive_run_facts"])
    assert facts["level"] == "0_300"
    assert facts["constituents"] == contributors
    assert facts["area_m2"] == 1.0e12
    assert facts["require_top"] == level.require_top       # resolved from the level (cfg passed None)
    assert facts["quantity"]["name"] == "ohc"                # the [quantity] table, nested as an object
    cfg_blob = json.loads(blob.attrs["ohc_derive_run_config"])
    assert cfg_blob["tag"] == "D-260901"
    assert cfg_blob["quantities"] == ["ohca", "ohu"]
