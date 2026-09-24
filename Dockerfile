# Test environment for ohc_derive — deps only, no code copied in (mount the source at runtime).
FROM python:3.12-slim

RUN pip install --no-cache-dir "numpy<2.5" "xarray>=2024.10" netCDF4 matplotlib pytest

WORKDIR /app

