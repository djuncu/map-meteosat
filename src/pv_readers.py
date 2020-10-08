#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 15:27:25 2018

@author: vincentch
"""
import glob
import xarray as xr
import os
import h5py
import pandas as pd
import time
import numpy as np

# homemade tools
import xr_read_sat_files
import utils_xr_process
import my_utils_plot

###########################
# used for map plotting only
# TODO: merge with version used for timeseries
def read_one_file_msg(f_in_path, varname, scaling, offset=0., valid_range = [], missing=[]):
    """
    Read MSG data with h5py
    Parameters : 
        varname : string. Corresponds to the name of the variable to read in the hdf5 files.
        f_in_path: str. Path to one unique file, compatible with glob.glob(). If several files are found, only the first will be used.
        scaling : numeric value. multiplicative factor to apply to the data in da
        offset: numeric value. additive factor to apply to the data in da
        missing: [] or list of values. the missing values provided are expected to refer to data in da BEFORE applying any scaling/offset.
        valid_range: [] or list of two values. min and max (included) thresholds applied on the scaled data from da (threshold are applied AFTER applying scaling/offset to da). Data outside this valid range are returned as NaN.
    """
    # get the list of data available
    f_list = glob.glob(f_in_path)
    if len(f_list) >1: print(f'WARNING: several files found for path {f_in_path}')

    f = f_list[0]
    file = h5py.File(f,'r')
    var = file[varname][:]
    da = xr.DataArray(var, dims = ('y','x'))
    del var
    da = utils_xr_process.apply_scaling_offset_missing(da, scaling=scaling,offset=offset, valid_range = valid_range, missing=missing)
    da = da.assign_attrs(file.attrs)
    file.close()
    return da

###########################
# Used for timeseries extraction & plotting       
# code using xarray reader (only compatible with h5zipped MSG & EPS files))
def read_msg_with_xr(f_list,decode_cf):
    """ Reader for MSG data. 
    NB: the files NEED TO BE PRE-TREATED FIRST WITH H5ZIP to ensure xarray can open these !
    example path to treated data is "/cnrm/vegeo/vincentch/VALIDATION_TOOLS/DATA/H5ZIPPED/MTAL_R/HDF5_LSASAF_MSG_ALBEDO-D10_MSG-Disk_????????0000_withdims.h5"
    """
    f_list.sort()
    datetime_list = get_datetime_from_filename_msg_alb(f_list)
    d = xr_read_sat_files.DatasourcePreprocessingFunctions()
    ds_in = xr.open_mfdataset(f_list,preprocess = d.open_mfdataset_preproc_for_LSAF_msg_eps, concat_dim = 'datetime', engine='h5netcdf',decode_cf=decode_cf) 
#    # verif dates in ds vs dates in filenames
#    for i,d in enumerate(datetime_list):
#        print(i,d,ds_in.isel(datetime=i)['datetime'].values)
#        print(d == ds_in.isel(datetime=i)['datetime'].values)

    # the preprocessing function for MSG (MTAL in particular) assigns the SENSING_START_DATE as the datetime. Yet this can be different from the center datetime (that is used in the filename)
    # we thus overwrite the datetimes with the list of datetime.datetime dates computed from the file_list  
    # force dates with the datelist computed from filenames:
    ds_in['datetime'] = datetime_list
    return ds_in
#def read_and_stack_several_files_msg_xr(f_in_path, varname, scaling, offset=0., valid_range = [], missing=[]):
#    """
#    Read MSG data with xarray
#    Parameters : 
#        varname : string. Corresponds to the name of the variable to read in the hdf5 files.
#        f_in_path: str. Path to one unique file, compatible with glob.glob(). If several files are found, only the first will be used.
#        scaling : numeric value. multiplicative factor to apply to the data in da
#        offset: numeric value. additive factor to apply to the data in da
#        missing: [] or list of values. the missing values provided are expected to refer to data in da BEFORE applying any scaling/offset.
#        valid_range: [] or list of two values. min and max (included) thresholds applied on the scaled data from da (threshold are applied AFTER applying scaling/offset to da). Data outside this valid range are returned as NaN.
#    """
#    # get the list of data available
#    f_list = glob.glob(f_in_path)
#    f_list.sort()
#
#    ds = read_msg_with_xr(f_list,decode_cf=False)
#    da = ds[varname]
#    da = utils_xr_process.apply_scaling_offset_missing(da, scaling=scaling,offset=offset, valid_range = valid_range, missing=missing)
#    return da
def read_and_stack_several_files_msg_xr_v2(f_in_path, sat_key, varnames, f_out_scaled):
    """
    Read MSG data with xarray
    Parameters : 
        varname : string. Corresponds to the name of the variable to read in the hdf5 files.
        f_in_path: str. Path to one unique file, compatible with glob.glob(). If several files are found, only the first will be used.
        scaling : numeric value. multiplicative factor to apply to the data in da
        offset: numeric value. additive factor to apply to the data in da
        missing: [] or list of values. the missing values provided are expected to refer to data in da BEFORE applying any scaling/offset.
        valid_range: [] or list of two values. min and max (included) thresholds applied on the scaled data from da (threshold are applied AFTER applying scaling/offset to da). Data outside this valid range are returned as NaN.
    """
    # get the list of data available
    f_list = glob.glob(f_in_path)
    f_list.sort()
    print('reading')

    ds = read_msg_with_xr(f_list,decode_cf=False)
    v_out = [x for x in ds.data_vars if x not in varnames]
    ds = ds.drop(v_out)
    d = read_input_dtypes(sat_key)
    dtypes_in = {k:v for k,v in d.items() if k in varnames}
    for varname in varnames:
        ds[varname] = ds[varname].astype(dtypes_in[varname])

    d = read_scaled_dtypes(sat_key)
    dtypes_sc = {k:v for k,v in d.items() if k in varnames}
    print('scaling')
    for varname in varnames:
        d = read_scaling_offset_fillvalue(sat_kw='msg_alb_bb', varname=varname)
        ds[varname] = utils_xr_process.apply_scaling_offset_missing_with_valid_before_scaling(ds[varname], dtypes_sc[varname],scaling=d['scaling'],offset=d['offset'], valid_range = d['valid_range_before_scaling'], missing=[d['fillvalue']])
    print('end scaling')
    if f_out_scaled is not None:
        my_utils_plot.ensure_dir(f_out_scaled)
        ds.to_netcdf(f_out_scaled)
    else:
        return ds


# code to read non-h5zipped files with h5py
def read_one_file_msg_v3(f_in_path, varnames, dtypes_dic, keep_attrs, store_read_to_netcdf=False, f_out=None):
    """
    Read MSG data with h5py - no rescaling/filtering
    Parameters : 
        varnames : list of strings. Corresponds to the names of the variables to read in the hdf5 files.
        f_in_path: str. Path to one unique file, compatible with glob.glob(). If several files are found, only the first will be used.
    """
    # get the list of data available
    f_list = glob.glob(f_in_path)
    if len(f_list) >1: print(f'WARNING: several files found for path {f_in_path}')

    f = f_list[0]
    #with h5py.File(f,'r') as file:
    file = h5py.File(f,'r')
    first = None
    for i,varname in enumerate(varnames):       
        print(varname)
        var = file[varname][:]
        if dtypes_dic[varname] is None: da = xr.DataArray(var, dims = ('y','x'),name=varname)
        if dtypes_dic[varname] is not None: da = xr.DataArray(var, dims = ('y','x'),name=varname).astype(dtypes_dic[varname])
        del var
        if keep_attrs: 
            attrs = {k:v for k,v in file[varname].attrs.items()}
            attrs = decode_numpy_bytes_attrs(attrs)
            da = da.assign_attrs(attrs)
        if first is None:
            ds = da.to_dataset()
            first = 'no more first'
        else:
            ds = xr.merge([ds,da])
        del da
    if keep_attrs: 
        attrs = {k:v for k,v in file.attrs.items()}
        attrs = decode_numpy_bytes_attrs(attrs)
        ds = ds.assign_attrs(attrs)
    chunks={'x':1500,'y':1500}
    ds = ds.chunk(chunks=chunks)
    file.close()
    
    if store_read_to_netcdf:
        my_utils_plot.ensure_dir(f_out)
        dict_encode_all={}
#        dict_encode={'zlib':True, 'complevel':4} 
#        for v in ds.data_vars:
#            dict_encode_all[v]=dict_encode
        ds.to_netcdf(f_out, encoding=dict_encode_all)
    else:
        return ds

def gather_on_sub_f_list(sub_f_list, store_read_to_netcdf,f_out_sub, **kw_reader):
    first_t = None
    for i,f in enumerate(sub_f_list):
        print(f)
        ds_one = read_one_file_msg_v3(f_in_path=f, **kw_reader, store_read_to_netcdf=False, f_out=None)
        
        if first_t is None:
            ds_out = ds_one.copy()
            first_t='no more first'
        else:
            ds_out = xr.concat([ds_out,ds_one],dim='datetime')
    if store_read_to_netcdf:
        my_utils_plot.ensure_dir(f_out_sub)
        dict_encode_all={}
#        dict_encode={'zlib':True, 'complevel':4} 
#        for v in ds_out.data_vars:
#            dict_encode_all[v]=dict_encode
        ds_out.to_netcdf(f_out_sub, encoding=dict_encode_all)
    else:
        return ds_out


def read_and_stack_several_files_msg_v4 (f_in_path, varnames, dtypes_dic, keep_attrs = False, store_read_to_netcdf=False, f_out=None, chunk_f_list = 20):
    """
    Read MSG data with h5py - no rescaling/filtering
    Parameters : 
        varnames : list of strings. Corresponds to the names of the variables to read in the hdf5 files.
        f_in_path: str. Path to one unique file, compatible with glob.glob(). If several files are found, only the first will be used.
    """
    # get the list of data available
    f_list = glob.glob(f_in_path)
    f_list.sort()
    
    if chunk_f_list is not None: sub_f_lists = chunk_list(f_list, sub_size=chunk_f_list)
    if chunk_f_list is None: sub_f_lists = [f_list]
    
    first_sub = None
    for it, sub_f_list in enumerate(sub_f_lists):
        kw_reader = dict(varnames=varnames, dtypes_dic=dtypes_dic, keep_attrs = keep_attrs)
        homedir = os.path.expanduser("~")
        f_out_sub = f'{homedir}/tmp/tmp_{time.time()}.nc'
        my_utils_plot.ensure_dir(f_out_sub)
        gather_on_sub_f_list(sub_f_list, store_read_to_netcdf=True, **kw_reader, f_out_sub =f_out_sub )
        ds_sub = xr.open_dataset(f_out_sub)
        os.remove(f_out_sub)
        del f_out_sub
        
        if first_sub is None:
            ds_out = ds_sub.copy()
            first_sub='no more first'
        else:
            ds_out = xr.concat([ds_out,ds_sub],dim='datetime')    

    if keep_attrs: ds_out = ds_out.assign_attrs(ds_sub.attrs)
    del ds_sub
    chunks={'x':1500,'y':1500}
    ds_out = ds_out.chunk(chunks=chunks)

    # add datetime info to the dataset
    datetime_list = get_datetime_from_filename_msg_alb(f_list)
    ds_out['datetime'] = datetime_list    

    print(ds_out)
    if store_read_to_netcdf:
        my_utils_plot.ensure_dir(f_out)
        dict_encode_all={}
#        dict_encode={'zlib':True, 'complevel':4} 
#        for v in ds_out.data_vars:
#            dict_encode_all[v]=dict_encode
#        print(dict_encode_all)
        ds_out.to_netcdf(f_out, engine = 'h5netcdf', encoding=dict_encode_all)
    else:
        return ds_out

def read_and_stack_several_files_msg_no_xr(f_in_path, sat_key, varnames, f_out_scaled):
    d = read_input_dtypes(sat_key)
    dtypes_in = {k:v for k,v in d.items() if k in varnames}
    d = read_scaled_dtypes(sat_key)
    dtypes_sc = {k:v for k,v in d.items() if k in varnames}
    store_read_to_netcdf=True
    chunk_f_list=20 # each nth file, store all to netcdf cache and reload it to remove data from memory
    if store_read_to_netcdf==False:
        f_out=None
        ds_all = read_and_stack_several_files_msg_v4(f_in_path=f_in_path, varnames=varnames, dtypes_dic=dtypes_in, store_read_to_netcdf=store_read_to_netcdf, 
                                                     f_out=f_out, keep_attrs=True, chunk_f_list=chunk_f_list)
    else:
        homedir = os.path.expanduser("~")
        f_out = f'{homedir}/tmp/tmp_{time.time()}.nc'
        my_utils_plot.ensure_dir(f_out)
        # we remove attrs storing in ds (with `keep_attrs=False` as so far can't figure out what type of attrs blocks the writing to netcdf (Unicode Type Error raised)
        read_and_stack_several_files_msg_v4(f_in_path=f_in_path, varnames=varnames, dtypes_dic=dtypes_in, store_read_to_netcdf=store_read_to_netcdf,
                                            f_out=f_out, keep_attrs=False, chunk_f_list=chunk_f_list)
        ds_all = xr.open_dataset(f_out)
        os.remove(f_out)
        
    ds_scaled = ds_all.copy()
    for varname in ds_all.data_vars:
        d = read_scaling_offset_fillvalue(sat_kw='msg_alb_bb', varname=varname)
        ds_scaled[varname] = utils_xr_process.apply_scaling_offset_missing_with_valid_before_scaling(ds_all[varname], dtypes_sc[varname],scaling=d['scaling'],offset=d['offset'], valid_range = d['valid_range_before_scaling'], missing=[d['fillvalue']])
    if f_out_scaled is not None:
        my_utils_plot.ensure_dir(f_out_scaled)
        ds_scaled.to_netcdf(f_out_scaled)
    else:
        return ds_scaled

## aux functions  
# force datetime dates in dataset to the filename/center dates (dates in MSG-MTAL attributes are not the center dates)
def get_datetime_from_filename_msg_alb(f_list):    
    date_list = []
    for f in f_list:
        # file naming examples: 
        # HDF5_LSASAF_MSG_ALBEDO-D10_MSG-Disk_201506150000
        # HDF5_LSASAF_MSG_ALBEDO-D10_MSG-Disk_201506150000_withdims.h5
        datestr = os.path.basename(f).split('MSG-Disk_')[1].split('_')[0][0:10]
        datetime_date = pd.to_datetime(datestr,format='%Y%m%d%H%M')
        date_list.append(datetime_date)    
    return date_list
# codes to read dict of properties for each sat product
def read_scaling_offset_fillvalue(sat_kw, varname):
    from src.all_sat_naming_scaling_missing import AllSatNamingScalingMissing
    d = AllSatNamingScalingMissing()
    dic=vars(d)
    if 'AL' in varname.upper() and (not 'ERR' in varname.upper() or not 'Q-FLAG' in varname.upper()):
        v_key = 'bsa_shwv'
    elif 'Q-FLAG' in varname.upper():
        v_key = 'brdf_alb_qual_flag'
    elif ('AL' in varname.upper()) and ('ERR' in varname.upper()):
        v_key = 'bsa_shwv_err'
    elif ('NMOD' in varname.upper()):
        v_key = 'nmod'
    elif ('AGE' in varname.upper()):
        v_key = 'age_info'
    fillvalue = dic['in_default_fillvalue'][sat_kw][v_key]
    scaling = dic['in_default_scaling'][sat_kw][v_key]
    offset = dic['in_default_offset'][sat_kw][v_key]
    valid_range_before_scaling = dic['in_default_valid_range_before_scaling_and_offset'][sat_kw][v_key]
    valid_range_after_scaling = [None,None]
    for i,v in enumerate(valid_range_before_scaling):
        if v is None:
            valid_range_after_scaling[i]=None
        else:
            valid_range_after_scaling[i] = (float(valid_range_before_scaling[i]) + float(offset)) *float(scaling)     
    d_out = {
            'fillvalue':fillvalue,
            'scaling':scaling,
            'offset':offset,
            'valid_range_before_scaling':valid_range_before_scaling,
            'valid_range_after_scaling':valid_range_after_scaling
            }
    return d_out
def read_varnames(sat_kw):
    from src.all_sat_naming_scaling_missing import AllSatNamingScalingMissing
    d = AllSatNamingScalingMissing()
    dic=vars(d)
    d_out = dic['input_varnames'][sat_kw]
    return d_out
def read_input_dtypes(sat_kw):
    from src.all_sat_naming_scaling_missing import AllSatNamingScalingMissing
    d = AllSatNamingScalingMissing()
    dic=vars(d)
    d_out = dic['input_variables_dtypes'][sat_kw]
    return d_out
def read_scaled_dtypes(sat_kw):
    from src.all_sat_naming_scaling_missing import AllSatNamingScalingMissing
    d = AllSatNamingScalingMissing()
    dic=vars(d)
    d_out = dic['scaled_variables_dtypes'][sat_kw]
    return d_out

def chunk_list(liste, sub_size):
    """ Cuts a list into sub_lists of size "sub_size" """
    # For item i in a range that is a length of liste,
    for i in range(0, len(liste), sub_size):
        # Create an index range for liste of n='sub_size' items:
        yield liste[i:i+sub_size]

def decode_numpy_bytes_attrs(attrs):
    """ required to convert bytes attrs from hdf5 into string for xarray. 
    Useful if we do not use xr_reader directly (where this is already implemented with engine='h5netcdf')
    """
    for k,v in attrs.items():
        if ( type(v) is np.bytes_ ): attrs[k] = v.decode('UTF-8')
        if ( type(v) is np.ndarray): attrs[k] = np.array([ x.decode('UTF-8') for x in v if type(x) is np.bytes_]+[ x for x in v if type(x) is not np.bytes_])
    return attrs
