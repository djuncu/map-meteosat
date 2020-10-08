#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 17:36:22 2018

@author: vincentch
"""
import xarray as xr
import numpy as np
import h5py

def rebuild_lon_lat_vgt(ds, zone):
    if str(zone).lower() == "euro" :
        x_size=10081
        y_size=5601
        lon_x0=-30
        lon_xf=60
        lat_yf=30
        lat_y0=80
    elif str(zone).lower() == "afri" :
        x_size=10081
        y_size=8961
        lon_x0=-30
        lon_xf=60
        lat_yf=-40
        lat_y0=40
    elif str(zone).lower() == "soam" :
        x_size=8961
        y_size=8961
        lon_x0=-110
        lon_xf=-30
        lat_yf=-60
        lat_y0=20
    elif str(zone).lower() == "noam" :
        x_size=13441
        y_size=7841
        lon_x0=-170
        lon_xf=-50
        lat_yf=10
        lat_y0=80
    elif str(zone).lower() == "asia" :
        x_size=13441
        y_size=8961
        lon_x0=60
        lon_xf=180
        lat_yf=0
        lat_y0=80
    elif str(zone).lower() == "ocea" :
        x_size=8961
        y_size=5601
        lon_x0=100
        lon_xf=180
        lat_yf=-50
        lat_y0=0
    dlon = (lon_xf - lon_x0)/(x_size - 1)
    dlat = (lat_yf - lat_y0)/(y_size - 1)       
    lons = np.arange(lon_x0, lon_xf+dlon,dlon)
    lats = np.arange(lat_y0, lat_yf+dlat,dlat)    
    print('dlon: ', dlon,'\n dlat: ', dlat)
    print('len(lons): ', len(lons),'\n len(lats): ', len(lats))
    print('lons: ', lons,'\n lats: ', lats)    
    # check size and end limits
    if len(lons) != x_size: print('ERROR on lons size')
    if len(lats) != y_size: print('ERROR on lats size')
    if lons[-1] != lon_xf: 
        if abs(lons[-1] - lon_xf) < abs(dlon/2):
            lons[-1] = lon_xf # avoid rounding issues
        else:
            print('Error in lons computation')
    if lats[-1] != lat_yf: 
        if abs(lats[-1] - lat_yf ) < abs(dlat/2):
            lats[-1] = lat_yf # avoid rounding issues
        else:
            print('Error in lats computation')
    
    # replace x & y dims by lons,lats and assign values in corresponding coordinate
    ds.rename({'x':'lon','y':'lat'}, inplace=True)
    ds.coords['lon'] = xr.DataArray(lons, dims='lon')
    ds.coords['lat'] = xr.DataArray(lats, dims='lat')
    return ds

def read_msg_lonlats(fillvalue_to_nan=True, 
                     latfile='/cnrm/vegeo/SAT/DATA/MSG/NRT-Operational/INPUTS/LAT-LON/MSG-Disk/HDF5_LSASAF_MSG_LAT_MSG-Disk_4bytesPrecision',
                     lonfile='/cnrm/vegeo/SAT/DATA/MSG/NRT-Operational/INPUTS/LAT-LON/MSG-Disk/HDF5_LSASAF_MSG_LON_MSG-Disk_4bytesPrecision'):
    latmsg=h5py.File(latfile,'r')
    latval=latmsg['LAT'][:]
    latmsg.close()
    latval=np.array(latval,dtype='float')
    latval=latval*0.0001
    if fillvalue_to_nan: 
        latval[latval==91]=np.nan
    else:
        latval[latval==91]=1e5

    lonmsg=h5py.File(lonfile,'r')
    lonval=lonmsg['LON'][:]
    lonmsg.close()
    lonval=np.array(lonval,dtype='float')
    lonval=lonval*0.0001
    if fillvalue_to_nan: 
        lonval[lonval==91]=np.nan
    else:
        lonval[lonval==91]=1e5
    
    # replace negative longitudes 
    lonval[lonval<0.] = 360. + lonval[lonval<0.]	    

    return lonval, latval


def read_mtg_lonlats(lonlatfile='/cnrm/vegeo/SAT/DATA/MTG/NO_SAVE/MTG_TESTING_2019-10-15/QUASI_STATIC/LSA_MTG_LATLON_MTG-Disk_201910090920.nc'):
    ds=xr.open_dataset(lonlatfile, decode_cf = False)
    latval=ds['LAT'].values
    latval = latval*0.0001
    lonval=ds['LON'].values
    lonval = lonval*0.0001
    ds.close()   
    return lonval, latval


def read_msg_iodc_lonlats(fillvalue_to_nan=True, latfile='/cnrm/vegeo/SAT/DATA/MSG-IODC/OperationalChain/LSASAF_MSG-IODC_Inputs/StaticData/LATLON/HDF5_LSASAF_MSG_LAT_IODC-Disk_201711300000',
                          lonfile='/cnrm/vegeo/SAT/DATA/MSG-IODC/OperationalChain/LSASAF_MSG-IODC_Inputs/StaticData/LATLON/HDF5_LSASAF_MSG_LON_IODC-Disk_201711300000'):
    latmsg=h5py.File(latfile,'r')
    latval=latmsg['LAT'][:]
    latmsg.close()
    latval=np.array(latval,dtype='float')
    if fillvalue_to_nan: latval[latval==9100]=np.nan
    latval=latval*0.01
    lonmsg=h5py.File(lonfile,'r')
    lonval=lonmsg['LON'][:]
    lonmsg.close()
    lonval=np.array(lonval,dtype='float')
    if fillvalue_to_nan: lonval[lonval==9100]=np.nan
    lonval=lonval*0.01    
    return lonval, latval


def read_msg_lonlats_to_xr():
    latmsg=h5py.File('/cnrm/vegeo/SAT/DATA/MSG/NRT-Operational/INPUTS/LAT-LON/MSG-Disk/HDF5_LSASAF_MSG_LAT_MSG-Disk_4bytesPrecision','r')
    latval=latmsg['LAT'][:]
    latmsg.close()
    latval=np.array(latval,dtype='float')
    latval[latval==910000]=np.nan
    latval=latval*0.0001
    lonmsg=h5py.File('/cnrm/vegeo/SAT/DATA/MSG/NRT-Operational/INPUTS/LAT-LON/MSG-Disk/HDF5_LSASAF_MSG_LON_MSG-Disk_4bytesPrecision','r')
    lonval=lonmsg['LON'][:]
    lonmsg.close()
    lonval=np.array(lonval,dtype='float')
    lonval[lonval==910000]=np.nan
    lonval=lonval*0.0001
    da_lon = xr.DataArray(lonval,dims=('y','x'))
    ds_geoloc = da_lon.to_dataset(name='lon')
    ds_geoloc['lat'] = xr.DataArray(latval,dims=('y','x'))
    return ds_geoloc

def rebuild_lon_lat_modis(ds):
    x_size=43200
    y_size=21600
    dlon = 1/120
    dlat = -1/120
    lon_x0=-180+dlon/2
    lon_xf=180-dlon/2
    lat_yf=-90-dlat/2
    lat_y0=90+dlat/2
    lons = np.arange(lon_x0, lon_xf+dlon, dlon)
    lats = np.arange(lat_y0, lat_yf+dlat, dlat)    
    print('dlon: ', dlon,'\n dlat: ', dlat)
    print('len(lons): ', len(lons),'\n len(lats): ', len(lats))
    print('lons: ', lons,'\n lats: ', lats)    
    # check size and end limits
    if len(lons) != x_size: print('ERROR on lons size')
    if len(lats) != y_size: print('ERROR on lats size')
    if lons[-1] != lon_xf: 
        if abs(lons[-1] - lon_xf) < abs(dlon/2):
            lons[-1] = lon_xf # avoid rounding issues
        else:
            print('Error in lons computation')
    if lats[-1] != lat_yf: 
        if abs(lats[-1] - lat_yf ) < abs(dlat/2):
            lats[-1] = lat_yf # avoid rounding issues
        else:
            print('Error in lats computation')
    
    # replace x & y dims by lons,lats and assign values in corresponding coordinate
    ds.rename({'x':'lon','y':'lat'}, inplace=True)
    ds.coords['lon'] = xr.DataArray(lons, dims='lon')
    ds.coords['lat'] = xr.DataArray(lats, dims='lat')
    return ds


def rebuild_lon_lat_metop(ds):
    res = 0.01
    lon=np.arange(-180,180,res)
    lat=np.arange(90,-90,-res)
    lat=np.append(lat,[-90],axis=0)

    lons, lats = np.meshgrid(lon,lat)
    lats_sinu = lats
    lons_sinu = lons/np.cos(np.deg2rad(lats))
    
    #ds.rename({'y':'lat'}, inplace=True)
    ds.coords['iy'] = ds.coords['y']
    ds.coords['ix'] = ds.coords['x']
    
    ds.coords['y'] = xr.DataArray(lat, name='y',dims=('y'))

    
    da = xr.DataArray(lons_sinu, name='lon',dims=('y','x'))
    da = da.chunk(chunks={'x':3000,'y':3001}) # chunk defines tiles size for parallelisation of processing & file writing
    ds.coords['lon2D'] = da.astype(np.float32) # set datatype as float 32 istead of default float 64
    del da
    da = xr.DataArray(lats_sinu, name='lat2D',dims=('y','x'))
    ds.coords['lat2D'] = da.chunk(chunks={'x':3000,'y':3001}).astype(np.float32)
    ds['lon2D'] = ds['lon2D'].where(ds['lon2D'] >= -180)
    ds['lon2D'] = ds['lon2D'].where(ds['lon2D'] <= 180)
    # NB: cant use multidim coordinates to filter the dataset
    # but can use groupby to process stats on specific classes/regions http://xarray.pydata.org/en/stable/examples/multidimensional-coords.html
    return ds
