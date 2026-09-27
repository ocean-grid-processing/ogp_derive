"""Input contracts: what a submission must carry, and what derive can infer when it doesn't.

Derive's native input is what the ingest `publish` step writes: an ME4OH-shaped `DATA` file that also
carries a `quantity` attr (the ingest [quantity] table), a `mapped_layer` attr keying the native level,
its `<NAME>ENS_` member sibling, and the upstream provenance blocks. Ingest stamps all of that.

A file from another group stops at the ME4OH protocol itself, and the protocol fixes everything derive
needs to know: `DATA(LONGITUDE, LATITUDE, TIME)` is ocean heat content density in TJ/m^2, computed with
the protocol's cp0 and rho0, on the common 1 degree grid, missing as NaN, one file per layer named
`OHC_<Y0>_<Y1>_lev<low>_<high>_exp<X>[...].nc`. So under `--contract ME4OH` derive takes the quantity
table from the protocol, the layer from the filename's `lev<low>_<high>` token, and no ensemble (the
protocol has no member files; a `DATA_SD`, if present, is not used). What was inferred is recorded on
the blob as `inferred_config`, keyed by constituent, so the record says the table came from the
specification rather than from the file.
"""
import os
import re

# The ME4OH protocol's quantity: OHC density in TJ/m^2, cp0/rho0 per the protocol (McDougall 2003).
ME4OH_QUANTITY = {
    "name": "ohc",
    "kind": "extensive",
    "units": "J/m2",
    "long_name": "ocean heat content",
    "scale_terms": {"cp0": 3989.244, "rho0": 1030.0},
    "publish_unit_factor": 1e12,
    "publish_units": "TJ/m^2",
}

_LEV = re.compile(r"_lev(\d+)_(\d+)(?=_|\.nc$)")

CONTRACTS = ("ME4OH",)


def me4oh_layer(path):
    """`<low>_<high>` from the protocol filename's `lev<low>_<high>` token."""
    m = _LEV.search(os.path.basename(path))
    if not m:
        raise SystemExit("%s: no `lev<low>_<high>` token in the filename; the ME4OH contract takes the "
                         "layer from the protocol filename" % os.path.basename(path))
    return "%s_%s" % m.groups()


def infer(contract, path):
    """What derive assumes about a submission under `contract` -> (tag, quantity, record)."""
    if contract != "ME4OH":
        raise SystemExit("unknown contract %r; known: %s" % (contract, list(CONTRACTS)))
    tag = me4oh_layer(path)
    record = {
        "contract": "ME4OH",
        "source_file": os.path.basename(path),
        "layer_from": "filename lev<low>_<high> token",
        "quantity_from": "ME4OH protocol (OHC density, TJ/m^2, cp0=3989.244, rho0=1030)",
        "ensemble": False,
    }
    return tag, dict(ME4OH_QUANTITY), record
