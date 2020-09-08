#!/usr/bin/env python

"""
Not decided yet...

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (08.09.2020)"
__email__ = "mdekauwe@gmail.com"

import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import cartopy.crs as ccrs
import matplotlib.patches as mpatches
import xarray as xr

def main(path, slice, GCM, RCM, domain, odir4, var):

    st = int(slice.split("-")[0])
    for i in range(19):
        st += 1
        tag = "%d-%d" % (st, st)

        if var == "tas":
            fn = os.path.join(path, "CCRC_NARCliM_01H_%s_%s.nc" % (tag, var))
            print(fn)
            ds = xr.open_dataset(fn)
            print(ds)

            sys.exit()

            

if __name__ == "__main__":

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
                main(path, slice, GCM, RCM, domain, odir4, var)
