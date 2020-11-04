#!/usr/bin/env python

"""
Get some NARCLIM precipitation stats...


That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (04.11.2020)"
__email__ = "mdekauwe@gmail.com"

import xarray as xr
import matplotlib.pyplot as plt
import sys
import datetime as dt
import pandas as pd
import numpy as np
from matplotlib.ticker import FixedLocator
import datetime
import os
import glob
from optparse import OptionParser

def main(time_slice, GCMs, RCMs):

    rows = []

    for GCM in GCMs:
        for RCM in RCMs:
            #print(GCM, RCM)

            path = "data/%s/%s/%s" % (time_slice, GCM, RCM)
            files = glob.glob(os.path.join(path, '*.nc'))
            if len(files) > 0:
                for fn in files:
                    print(fn)
                    try:
                        df = read_cable_file(fn)

                        num = os.path.basename(fn).split("_")[-3]

                        map = None
                        cov = None
                        print(fn)
                        print(num)
                        sys.exit()
                        #plc_mean = np.mean(df.resample("D").agg("mean").values)
                        #plc_max = np.max(df.resample("D").agg("max").values)

                        rows.append([GCM, RCM, time_slice, spp, \
                                     num, map, cov])
                    except:
                        rows.append([GCM, RCM, time_slice, spp, \
                                     -999.9, -999.9, -999.9])
            else:
                rows.append([GCM, RCM, time_slice, spp, \
                             -999.9, -999.9, -999.9])
    df_out = pd.DataFrame(rows,
                          columns=['gcm','rcm','time','species',\
                                   'num', 'map','cov'])


    df_out.to_csv("narclim_rainfall.csv", index=False)

def read_cable_file(fname):
    vars_to_keep = ['Rainf']
    ds = xr.open_dataset(fname, decode_times=False)

    time_jump = int(ds.time[1].values) - int(ds.time[0].values)
    if time_jump == 3600:
        freq = "H"
    elif time_jump == 1800:
        freq = "30M"
    else:
        raise("Time problem")

    units, reference_date = ds.time.attrs['units'].split('since')
    df = ds[vars_to_keep].squeeze(dim=["x","y"], drop=True).to_dataframe()
    start = reference_date.strip().split(" ")[0].replace("-","/")
    df['dates'] = pd.date_range(start=start, periods=len(df), freq=freq)
    df = df.set_index('dates')

    return df


if __name__ == "__main__":

    time_slices = ["1990-2009", "2020-2039", "2060-2079"]
    GCMs = ["CCCMA3.1", "CSIRO-MK3.0", "ECHAM5", "MIROC3.2"]
    RCMs = ["R1", "R2", "R3"]

    time_slice = time_slices[0]

    main(time_slice, GCMs, RCMs)
