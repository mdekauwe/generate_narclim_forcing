#!/usr/bin/env python

import sys
import os
import numpy as np
import xarray as xr
import pandas as pd
import netCDF4 as nc
import datetime
import optparse

def find_nearest(lats, lat, lons, lon):

    # how close each latitude and longitude is to the grid point we want
    abs_lat = np.abs(lats-lat)
    abs_lon = np.abs(lons-lon)

    # combine results and finds the local maximum
    c = np.maximum(abs_lon, abs_lat)

    # fine matching grid point, nb in a flattened array
    #latlon_idx = np.argmin(c)
    #print(latlon_idx)


    # get row, col in non-flattened array
    xx, yy = np.where(c == np.min(c))

    return xx[0], yy[0]

#def find_nearest(a, b):
#    idx = np.argmin(np.abs(a-b))
#
#    return idx

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

    #df3[var] /= 3600.

    for j in range(len(df3[var].values)):
        print((df3[var].values[j])
    p
