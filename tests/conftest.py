"""Shared builders for the factory tests.

Grid choice: two symmetric latitudes (-45, 45) so every cell has the same area, which makes
area-weighted integrals easy to check by hand (integral = value * n_cells * cell_area).
"""
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

LAT = np.array([-45.0, 45.0])
LON = np.array([0.0, 120.0, 240.0])
NLAT, NLON = LAT.size, LON.size


# The ingest [quantity] table as the OHC submissions carry it (the `quantity` attr, parsed).
OHC_QUANTITY = {"name": "ohc", "kind": "extensive", "units": "J/m2", "long_name": "ocean heat content",
                "scale_terms": {"cp0": 3989.0, "rho0": 1030.0},
                "publish_unit_factor": 1e12, "publish_units": "TJ/m^2"}


LOCALGP_LEVELS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                              "levels", "localgp.toml")
_TABLE = None


def table():
    """The LocalGP level plan, loaded once."""
    global _TABLE
    if _TABLE is None:
        import levels
        _TABLE = levels.load(LOCALGP_LEVELS)
    return _TABLE


def level(name):
    return table()[name]


def months(n, start_year=2001):
    """A monthly datetime axis of length n, starting in January of start_year."""
    return pd.date_range("%d-01-01" % start_year, periods=n, freq="MS")


def field(arr, time):
    """(realization, time, lat, lon) ndarray -> DataArray on the standard coords."""
    arr = np.asarray(arr, dtype="float64")
    return xr.DataArray(arr, dims=("realization", "time", "lat", "lon"),
                        coords={"realization": np.arange(arr.shape[0]), "time": time,
                                "lat": LAT, "lon": LON})


def const_field(value, n_real=1, n_time=3, nan_cells=(), member_nan=(), start_year=2001):
    """Field equal to `value` everywhere, with optional NaN columns.

    nan_cells:  iterable of (i, j)      -> NaN at that cell in every realization, all time
    member_nan: iterable of (r, i, j)   -> NaN at that cell in realization r only, all time
    """
    arr = np.full((n_real, n_time, NLAT, NLON), float(value))
    for (i, j) in nan_cells:
        arr[:, :, i, j] = np.nan
    for (r, i, j) in member_nan:
        arr[r, :, i, j] = np.nan
    return field(arr, months(n_time, start_year))


def series(arr, time=None, start_year=2001):
    """(realization, time) ndarray -> DataArray with realization/time coords."""
    arr = np.asarray(arr, dtype="float64")
    time = time if time is not None else months(arr.shape[1], start_year)
    return xr.DataArray(arr, dims=("realization", "time"),
                        coords={"realization": np.arange(arr.shape[0]), "time": time})


def bathy(depth_grid):
    """(lat, lon) seafloor depth in metres -> DataArray."""
    return xr.DataArray(np.asarray(depth_grid, dtype="float64"),
                        dims=("lat", "lon"), coords={"lat": LAT, "lon": LON})


def constituent(tag, n_fac, top, bottom, field_value):
    """A masks/levels-shaped constituent dict."""
    return {"tag": tag, "n_fac": n_fac, "top": top, "bottom": bottom, "field_value": field_value}
