#!/usr/bin/env python

"""
For each time slice, GCM, RCM loop through and collect up all the met data
we need to generate a forcing file for CABLE-hydraulics simulations.

TODO:
- will need to get lat/lon pairs from euc species list.
- add netcdf creation for output
- need to interpolate some vars, e.g. tair is hourly, radiation is 3-hourly

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (08.09.2020)"
__email__ = "mdekauwe@gmail.com"

import sys
import os
import numpy as np
import xarray as xr
import pandas as pd

def find_nearest(a, b):
    idx = np.argmin(np.abs(a-b))

    return idx


def get_data(fn, var):
    ds = xr.open_dataset(fn)
    lats = ds.lat[:,0].values # 2D arrays, squeeze
    lons = ds.lon[0,:].values # 2D arrays, squeeze
    ii = find_nearest(lats, lat)
    jj = find_nearest(lons, lon)
    data = ds[var][:,ii,jj].to_dataframe()
    data = data.drop(['lat', 'lon'], axis=1)

    return data

def main(path, slice, GCM, RCM, domain, odir4, var, lat, lon):

    nyears = 19
    df_tas = pd.DataFrame(columns=['tas'])

    st = int(slice.split("-")[0])
    for i in range(nyears):
        st += 1
        tag = "%d-%d" % (st, st)

        if var == "tas":
            fn = os.path.join(path, "CCRC_NARCliM_01H_%s_%s.nc" % (tag, var))
        else:
            raise("not implemented for ... ")

        df = get_data(fn, var)
        df_tas = df_tas.append(df)

    df_tas.to_csv("tas.csv", index=True)

    sys.exit()

if __name__ == "__main__":

    # test loc -> Eucalyptus accedens
    lat = -29.07
    lon = 114.87

    base_path = "/srv/ccrc/data30/z3393020/NARCliM/postprocess/"

    odir = "data"
    if not os.path.exists(odir):
        os.makedirs(odir)

    time_slices = ["1990-2009", "2020-2039", "2060-2079"]
    GCMs = ["CCCMA3.1", "CSIRO-MK3.0", "ECHAM5", "MIROC3.2"]
    RCMs = ["R1", "R2", "R3"]
    domains = ['d01','d02']

    domain = domains[0] # whole of aus
    var = "tas"

    for slice in time_slices:

        odir2 = os.path.join(odir, slice)
        if not os.path.exists(odir2):
            os.makedirs(odir2)

        for GCM in GCMs:

            odir3 = os.path.join(odir2, GCM)
            if not os.path.exists(odir3):
                os.makedirs(odir3)

            for RCM in RCMs:
                print(slice, GCM, RCM)
                odir4 = os.path.join(odir3, RCM)
                if not os.path.exists(odir4):
                    os.makedirs(odir4)

                path = "%s/%s/%s/%s/%s" % (base_path, slice, GCM, RCM, domain)
                main(path, slice, GCM, RCM, domain, odir4, var, lat, lon)
