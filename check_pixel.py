#!/usr/bin/env python

import sys
import os
import numpy as np
import xarray as xr
import pandas as pd
import netCDF4 as nc
import datetime
import optparse

def get_data(fn, var, lat, lon):
    ds = xr.open_dataset(fn)
    lats = ds.lat.values
    lons = ds.lon.values
    ii, jj = find_nearest(lats, lat, lons, lon)
    data = ds[var][:,ii,jj].to_dataframe()
    data = data.drop(['lat', 'lon'], axis=1)
    ds.close()

    return data

slice = "1990-2009"
path = "/srv/ccrc/data30/z3393020/NARCliM/postprocess/1990-2009/CCCMA3.1/R1/d01"
nyears = 20
st = int(slice.split("-")[0])

# 72 Eucalyptus grandis,-21.081769,148.4178
lat = -21.081769
lon = 148.4178
for i in range(nyears):
    tag = "%d-%d" % (st, st)
    year = int(st)

    var = "pracc" # precip
    fn = os.path.join(path, "CCRC_NARCliM_01H_%s_%s.nc" % (tag, var))
    df3 = get_data(fn, var, lat, lon)
    df3[var] /= 3600.

    print(df3[var])

    
