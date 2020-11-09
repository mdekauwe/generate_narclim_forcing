#!/usr/bin/env python

"""
For each time slice, GCM, RCM loop through and collect up all the met data
we need to generate a forcing file for CABLE-hydraulics simulations.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (08.09.2020)"
__email__ = "mdekauwe@gmail.com"

import sys
import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import netCDF4 as nc
import datetime
import optparse

def main(path, bias_path, slice, GCM, RCM, domain, opath, spp, lat, lon, df_co2,
         count):


    cols = ['tas','huss','pracc', 'wss', 'ps', 'CO2air']
    dfx = pd.DataFrame(columns=cols)


    # Get all bias corrected PPT files
    files = glob.glob('%s/*_DAY_*_pracc_bc.nc' % (bias_path))

    individual_files = []
    for fn in files:
        dsx = xr.open_dataset(fn)
        lats = dsx.lat[:,0].values # 2D arrays, squeeze
        lons = dsx.lon[0,:].values # 2D arrays, squeeze
        print(lat, lon)
        ii = find_nearest(lats, lat)
        jj = find_nearest(lons, lon)
        print(ii, jj)
        print(dsx['lat'][ii,jj].values)
        print(dsx['lon'][ii,jj].values)
        data = dsx['pracc_bc'][:,ii,jj].to_dataframe()
        data = data.drop(['lat', 'lon'], axis=1)
        print(np.nanmean(data.pracc_bc))
        sys.exit()
        dsx.close()
        individual_files.append(data)
    modis_ds = pd.concat(individual_files, axis=0)

    print(modis_ds)
    sys.exit()


    nyears = 20
    st = int(slice.split("-")[0])

    for i in range(nyears):
        tag = "%d-%d" % (st, st)
        year = int(st)

        var = "tas" # air temp
        fn = os.path.join(path, "CCRC_NARCliM_01H_%s_%s.nc" % (tag, var))
        df1 = get_data(fn, var, lat, lon)

        # There is a time offset issue as it is in UTC, so need to +10 hours
        df1 = df1.shift(periods=10)
        df1[var][0:10] = df1[var][11] # fill the first NaNs we added

        # Add in CO2
        co2 = df_co2[df_co2.Year == year].co2.values[0]
        df1['CO2air'] = co2

        var = "huss" # Qair
        fn = os.path.join(path, "CCRC_NARCliM_01H_%s_%s.nc" % (tag, var))
        df2 = get_data(fn, var, lat, lon)

        # There is a time offset issue as it is in UTC, so need to +10 hours
        df2= df2.shift(periods=10)
        df2[var][0:10] = df2[var][11] # fill the first NaNs we added

        var = "pracc" # precip
        fn = os.path.join(path, "CCRC_NARCliM_01H_%s_%s.nc" % (tag, var))
        df3 = get_data(fn, var, lat, lon)
        # someone has written the precip with a 30 min timestep (e.g. 06:30:00
        # instead of 06:00:00), even though it is hourly, use the Qair index
        df3.index = df2.index

        # units are kg m-2 accumulated over the hour, need to be kg m-2 s-1
        df3[var] /= 3600.

        # is it possible they've aggregated half hour data based on the time
        # stamp and the fact that the rainfall numbers look huge? Will need to
        # check, test this for now.
        #df3[var] /= 7200.

        # There is a time offset issue as it is in UTC, so need to +10 hours
        df3 = df3.shift(periods=10)
        df3[var][0:10] = 0.0 # fill the first NaNs we added

        var = "wss" # wind
        fn = os.path.join(path, "CCRC_NARCliM_01H_%s_%s.nc" % (tag, var))
        df4 = get_data(fn, var, lat, lon)

        # There is a time offset issue as it is in UTC, so need to +10 hours
        df4 = df4.shift(periods=10)
        df4[var][0:10] = df4[var][11] # fill the first NaNs we added

        var = "ps" # pressure
        fn = os.path.join(path, "CCRC_NARCliM_01H_%s_%s.nc" % (tag, var))
        df5 = get_data(fn, var, lat, lon)

        # There is a time offset issue as it is in UTC, so need to +10 hours
        df5 = df5.shift(periods=10)
        df5[var][0:10] = df5[var][11] # fill the first NaNs we added

        frames = [df1, df2, df3, df4, df5]
        result = pd.concat(frames, axis=1)

        # deal with leap year missing data, we will infill as one day isn't
        # a key focus of our analysis
        result = result.fillna(method='ffill')

        dfx = dfx.append(result)

        st += 1

    # Radiation data is 3-hrly and concatenated into 5 year chunks...
    cols = ['rlds','rsds']
    dfy = pd.DataFrame(columns=cols)

    nyears = 4 # 5 year file segments
    st = int(slice.split("-")[0])
    for i in range(nyears):
        tag = "%d-%d" % (st, st+4)

        var = "rlds" # LWdown
        fn = os.path.join(path, "CCRC_NARCliM_03H_%s_%s.nc" % (tag, var))
        df6 = get_data(fn, var, lat, lon)

        # There is a time offset issue as it is in UTC, so need to +10 hours
        df6 = df6.shift(periods=10)
        df6[var][0:10] = df6[var][11] # fill the NaNs we added

        # We need to turn the 3hly data into hrly, linearly interpolate...
        i = pd.date_range(start=df6['rlds'].index[0],
                          end=df6['rlds'].index[-1], freq='H')
        df6 = df6.reindex(i).interpolate()

        var = "rsds" # SWdown
        fn = os.path.join(path, "CCRC_NARCliM_03H_%s_%s.nc" % (tag, var))
        df7 = get_data(fn, var, lat, lon)

        # There is a time offset issue as it is in UTC, so need to +10 hours
        df7 = df7.shift(periods=10)
        df7[var][0:10] = 0.0 # fill the NaNs we added

        # We need to turn the 3hly data into hrly, linearly interpolate...
        df7 = df7.reindex(i).interpolate(method='linear')

        frames = [df6, df7]
        result = pd.concat(frames, axis=1)

        # deal with leap year missing data, we will infill as one day isn't
        # a key focus of our analysis
        result = result.fillna(method='ffill')

        dfy = dfy.append(result)
        st += 5

    # Join the hourly and the interpolated hourly data.
    frames = [dfx, dfy]
    df_out = pd.concat(frames, axis=1)
    df_out = df_out.fillna(method='ffill')

    df_out['date'] = pd.to_datetime(df_out.index)
    cols = ['date','tas','huss','pracc','wss','ps','CO2air','rlds','rsds']
    df_out = df_out[cols]
    df_out.rename(columns={'tas':'Tair', 'huss':'Qair', 'pracc':'Rainf',
                           'wss':'Wind', 'ps':'PSurf', 'rlds':'LWdown',
                           'rsds':'SWdown'},
                  inplace=True)

    #df_out.to_csv("test.csv", index=False)

    out_fname = "narclim_met_%s_%d.nc" % (spp.replace(" ", "_"), count)
    out_fname = os.path.join(opath, out_fname)
    create_cable_nc_file(df_out, lat, lon, out_fname)
    #sys.exit()

def create_cable_nc_file(df, lat, lon, out_fname):

    ndim = 1
    n_timesteps = len(df)
    times = []
    secs = 0.0
    for i in range(n_timesteps):
        times.append(secs)
        secs += 3600.

    # create file and write global attributes
    f = nc.Dataset(out_fname, 'w', format='NETCDF4')
    f.description = 'NARCLIM met data, created by Martin De Kauwe'
    f.history = "Created by: %s" % (os.path.basename(__file__))
    f.creation_date = "%s" % (datetime.datetime.now())
    f.contact = "mdekauwe@gmail.com"

    # set dimensions
    f.createDimension('time', None)
    f.createDimension('z', ndim)
    f.createDimension('y', ndim)
    f.createDimension('x', ndim)
    #f.Conventions = "CF-1.0"

    # create variables
    time = f.createVariable('time', 'f8', ('time',))
    time.units = "seconds since %s 00:00:00" % (df.index[0])
    time.long_name = "time"
    time.calendar = "standard"

    z = f.createVariable('z', 'f8', ('z',))
    z.long_name = "z"
    z.long_name = "z dimension"

    y = f.createVariable('y', 'f8', ('y',))
    y.long_name = "y"
    y.long_name = "y dimension"

    x = f.createVariable('x', 'f8', ('x',))
    x.long_name = "x"
    x.long_name = "x dimension"

    latitude = f.createVariable('latitude', 'f8', ('y', 'x',))
    latitude.units = "degrees_north"
    latitude.missing_value = -9999.
    latitude.long_name = "Latitude"

    longitude = f.createVariable('longitude', 'f8', ('y', 'x',))
    longitude.units = "degrees_east"
    longitude.missing_value = -9999.
    longitude.long_name = "Longitude"

    SWdown = f.createVariable('SWdown', 'f8', ('time', 'y', 'x',))
    SWdown.units = "W/m^2"
    SWdown.missing_value = -9999.
    SWdown.long_name = "Surface incident shortwave radiation"
    SWdown.CF_name = "surface_downwelling_shortwave_flux_in_air"

    Tair = f.createVariable('Tair', 'f8', ('time', 'z', 'y', 'x',))
    Tair.units = "K"
    Tair.missing_value = -9999.
    Tair.long_name = "Near surface air temperature"
    Tair.CF_name = "surface_temperature"

    Rainf = f.createVariable('Rainf', 'f8', ('time', 'y', 'x',))
    Rainf.units = "mm/s"
    Rainf.missing_value = -9999.
    Rainf.long_name = "Rainfall rate"
    Rainf.CF_name = "precipitation_flux"

    Qair = f.createVariable('Qair', 'f8', ('time', 'z', 'y', 'x',))
    Qair.units = "kg/kg"
    Qair.missing_value = -9999.
    Qair.long_name = "Near surface specific humidity"
    Qair.CF_name = "surface_specific_humidity"

    Wind = f.createVariable('Wind', 'f8', ('time', 'z', 'y', 'x',))
    Wind.units = "m/s"
    Wind.missing_value = -9999.
    Wind.long_name = "Scalar windspeed" ;
    Wind.CF_name = "wind_speed"

    PSurf = f.createVariable('PSurf', 'f8', ('time', 'y', 'x',))
    PSurf.units = "Pa"
    PSurf.missing_value = -9999.
    PSurf.long_name = "Surface air pressure"
    PSurf.CF_name = "surface_air_pressure"

    LWdown = f.createVariable('LWdown', 'f8', ('time', 'y', 'x',))
    LWdown.units = "W/m^2"
    LWdown.missing_value = -9999.
    LWdown.long_name = "Surface incident longwave radiation"
    LWdown.CF_name = "surface_downwelling_longwave_flux_in_air"

    CO2air = f.createVariable('CO2air', 'f8', ('time', 'z', 'y', 'x',))
    CO2air.units = "ppm"
    CO2air.missing_value = -9999.
    CO2air.long_name = ""
    CO2air.CF_name = ""

    elevation = f.createVariable('elevation', 'f8', ('y', 'x',))
    elevation.units = "m" ;
    elevation.missing_value = -9999.
    elevation.long_name = "Site elevation above sea level" ;

    za_tq = f.createVariable('za_tq', 'f8', ('y', 'x',))
    za_tq.units = "m"
    za_tq.missing_value = -9999.
    za_tq.long_name = "level of lowest atmospheric model layer"

    za_uv = f.createVariable('za_uv', 'f8', ('y', 'x',))
    za_uv.units = "m"
    za_uv.missing_value = -9999.
    za_uv.long_name = "level of lowest atmospheric model layer"

    # write data to file
    x[:] = ndim
    y[:] = ndim
    z[:] = ndim
    time[:] = times
    latitude[:] = lat
    longitude[:] = lon

    SWdown[:,0,0] = df.SWdown.values.reshape(n_timesteps, ndim, ndim)
    Tair[:,0,0,0] = df.Tair.values.reshape(n_timesteps, ndim, ndim, ndim)
    Rainf[:,0,0] = df.Rainf.values.reshape(n_timesteps, ndim, ndim)
    Qair[:,0,0,0] = df.Qair.values.reshape(n_timesteps, ndim, ndim, ndim)
    Wind[:,0,0,0] = df.Wind.values.reshape(n_timesteps, ndim, ndim, ndim)
    PSurf[:,0,0] = df.PSurf.values.reshape(n_timesteps, ndim, ndim)
    LWdown[:,0,0] = df.LWdown.values.reshape(n_timesteps, ndim, ndim)
    CO2air[:,0,0] = df.CO2air.values.reshape(n_timesteps, ndim, ndim, ndim)
    za_tq[:] = 2.  # temp
    za_uv[:] = 10. # wind
    #elevation[0,0] = elev

    f.close()



def find_nearest(a, b):
    idx = np.argmin(np.abs(a-b))

    return idx

def get_data(fn, var, lat, lon):
    ds = xr.open_dataset(fn)
    lats = ds.lat[:,0].values # 2D arrays, squeeze
    lons = ds.lon[0,:].values # 2D arrays, squeeze
    ii = find_nearest(lats, lat)
    jj = find_nearest(lons, lon)
    data = ds[var][:,ii,jj].to_dataframe()
    data = data.drop(['lat', 'lon'], axis=1)
    ds.close()

    return data

def cmd_line_parser():

    p = optparse.OptionParser()
    p.add_option("-g", default="CCCMA3.1", help="gcm filename")
    options, args = p.parse_args()

    return (options.g)

if __name__ == "__main__":

    (GCM) = cmd_line_parser()

    base_path = "/srv/ccrc/data30/z3393020/NARCliM/postprocess"
    base_path_bias = "/srv/ccrc/data30/z3393020/NARCliM/Bias_corrected"

    df_co2 = pd.read_csv("AmaFACE_co2npdepforcing_1850_2100_AMB.csv", sep=";")
    df_co2.rename(columns={'CO2 [ppm]':'co2'}, inplace=True)

    #df_spp = pd.read_csv("species_locations_sub_sampled.csv")
    df_spp = pd.read_csv("test.csv")

    odir = "data"
    if not os.path.exists(odir):
        os.makedirs(odir)

    time_slices = ["1990-2009", "2020-2039", "2060-2079"]
    time_slices_bias = ["1990-2010", "2020-2040", "2060-2080"]
    #GCMs = ["CCCMA3.1", "CSIRO-MK3.0", "ECHAM5", "MIROC3.2"]
    RCMs = ["R1", "R2", "R3"]
    domains = ['d01','d02']

    domain = domains[0] # whole of aus


    for i,slice in enumerate(time_slices):

        odir2 = os.path.join(odir, slice)
        if not os.path.exists(odir2):
            os.makedirs(odir2)

        #for GCM in GCMs:

        odir3 = os.path.join(odir2, GCM)
        if not os.path.exists(odir3):
            os.makedirs(odir3)

        for RCM in RCMs:

            odir4 = os.path.join(odir3, RCM)
            if not os.path.exists(odir4):
                os.makedirs(odir4)

            path = "%s/%s/%s/%s/%s" % (base_path, slice, GCM, RCM, domain)

            bias_slice = time_slices_bias[i]
            bias_path = "%s/%s/%s/%s/%s" % (base_path_bias, GCM, RCM, bias_slice, domain)

            for i in range(len(df_spp)):
                spp = df_spp.species[i]
                lat = round(df_spp.lat[i], 2)
                lon = round(df_spp.lon[i], 2)
                print(i, spp)
                main(path, bias_path, slice, GCM, RCM, domain, odir4, spp, lat,
                     lon, df_co2, i)
