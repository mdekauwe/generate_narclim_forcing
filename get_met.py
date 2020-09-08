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

def main(path, slice, GCM, RCM, domain, odir4, lat, lon):

    cols = ['tas','huss','pracc', 'wss', 'ps']
    dfx = pd.DataFrame(columns=cols)

    nyears = 20
    st = int(slice.split("-")[0])
    for i in range(nyears):

        print("hourly: %d:%d" % (i, nyears))
        tag = "%d-%d" % (st, st)

        var = "tas" # air temp
        fn = os.path.join(path, "CCRC_NARCliM_01H_%s_%s.nc" % (tag, var))
        print(fn)
        sys.exit()
        df1 = get_data(fn, var)

        var = "huss" # Qair
        fn = os.path.join(path, "CCRC_NARCliM_01H_%s_%s.nc" % (tag, var))
        df2 = get_data(fn, var)

        var = "pracc" # precip
        fn = os.path.join(path, "CCRC_NARCliM_01H_%s_%s.nc" % (tag, var))
        df3 = get_data(fn, var)
        # someone has written the precip with a 30 min timestep (e.g. 06:30:00
        # instead of 06:00:00), even though it is hourly, use the Qair index
        df3.index = df2.index

        var = "wss" # wind
        fn = os.path.join(path, "CCRC_NARCliM_01H_%s_%s.nc" % (tag, var))
        df4 = get_data(fn, var)

        var = "ps" # pressure
        fn = os.path.join(path, "CCRC_NARCliM_01H_%s_%s.nc" % (tag, var))
        df5 = get_data(fn, var)

        frames = [df1, df2, df3, df4, df5]
        result = pd.concat(frames, axis=1)

        dfx = dfx.append(result)

        st += 1

    print("3-hourly")

    # Radiation data is 3-hrly and concatenated into 5 year chunks...
    cols = ['rlds','rsds']
    dfy = pd.DataFrame(columns=cols)

    nyears = 4 # 5 year file segments
    st = int(slice.split("-")[0])
    for i in range(nyears):

        tag = "%d-%d" % (st, st+4)

        var = "rlds" # LWdown
        fn = os.path.join(path, "CCRC_NARCliM_03H_%s_%s.nc" % (tag, var))
        df6 = get_data(fn, var)

        print("3-hourly interp")
        # We need to turn the 3hly data into hrly, linearly interpolate...
        i = pd.DatetimeIndex(df1['tas'].index[0], df1['tas'].index[-1],
                             freq='H')
        df6 = df6.reindex(i).interpolate()

        var = "rsds" # SWdown
        fn = os.path.join(path, "CCRC_NARCliM_03H_%s_%s.nc" % (tag, var))
        df7 = get_data(fn, var)

        # We need to turn the 3hly data into hrly, linearly interpolate...
        df7 = df7.reindex(i).interpolate()

        frames = [df6, df7]
        result = pd.concat(frames, axis=1)

        dfy = dfy.append(result)
        st += 5

    print(df_out2)
    sys.exit()

    # Join the hourly and the interpolated hourly data.
    frames = [dfx, dfy]
    df_out = pd.concat(frames, axis=1)

    df_out['date'] = pd.to_datetime(df_out.index)
    cols = ['date'] + cols
    df_out = df_out[cols]
    df_out.rename(columns={'tas':'Tair', 'huss':'Qair', 'pracc':'Precip',
                           'wss':'Wind', 'ps':'Psurf', 'rlds':'LWdown',
                           'rsds':'SWdown'},
                  inplace=True)
    df_out.to_csv("test.csv", index=False)

    sys.exit()

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
    ds.close()

    return data


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


    for slice in time_slices:

        odir2 = os.path.join(odir, slice)
        if not os.path.exists(odir2):
            os.makedirs(odir2)

        for GCM in GCMs:

            odir3 = os.path.join(odir2, GCM)
            if not os.path.exists(odir3):
                os.makedirs(odir3)

            for RCM in RCMs:

                odir4 = os.path.join(odir3, RCM)
                if not os.path.exists(odir4):
                    os.makedirs(odir4)

                path = "%s/%s/%s/%s/%s" % (base_path, slice, GCM, RCM, domain)
                main(path, slice, GCM, RCM, domain, odir4, lat, lon)
