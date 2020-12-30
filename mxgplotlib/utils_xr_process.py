#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 16:48:07 2018

@author: vincentch
"""
import numpy as np
import xarray as xr
import glob
import h5py


###########################
# used for map plotting only
# TODO: merge with version used for timeseries
def read_one_file_meteosat(f_in_path, varname, scaling, offset=0., valid_range = [], missing=[], meteosatGen = '2'):
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
    if meteosatGen == '2':
        inFile = h5py.File(f,'r')
        var = inFile[varname][:]
        da = xr.DataArray(var, dims = ('y','x'))
        del var
        da = apply_scaling_offset_missing(da, scaling=scaling,offset=offset, valid_range = valid_range, missing=missing)
        da = da.assign_attrs(inFile.attrs)
        inFile.close()
    elif meteosatGen == '3':
        da = xr.open_dataset(f)
        da[varname] = apply_scaling_offset_missing(da[varname], scaling=scaling,offset=offset, valid_range = valid_range, missing=missing)
        da[varname] = da[varname].squeeze()
        if da[varname].dims[0] != 'y':
            da[varname] = da[varname].rename({da[varname].dims[0]:'y', da[varname].dims[1]:'x'})
        da = da[varname]

    return da

def apply_stats_on_xr_da(da,stat,dims,skipna=False, keep_attrs=False):
    """
    Returns da with the statistic stat applied on dimension(s) dims.
    da: xarray DataArray
    dims : str or sequence of str
        Dimension(s) over which to apply `stat`.
    skipna : bool, optional
        If True, skip missing values (as marked by NaN). By default, only
        skips missing values for float dtypes; other dtypes either do not
        have a sentinel missing value (int) or skipna=True has not been
        implemented (object, datetime64 or timedelta64).
    """
    if stat == 'mean':
        da = da.mean(dim=dims, skipna=skipna, keep_attrs=keep_attrs)
    if stat == 'median':
        da = da.load().median(dim=dims, skipna=skipna, keep_attrs=keep_attrs) # load required here as median is not implmeneted for dask dataarrays - so we use the numpy array
    if stat == 'min':
        da = da.min(dim=dims, skipna=skipna, keep_attrs=keep_attrs)
    if stat == 'max':
        da = da.max(dim=dims, skipna=skipna, keep_attrs=keep_attrs)
    if stat == 'std':
        da = da.std(dim=dims, skipna=skipna, keep_attrs=keep_attrs)
    if 'percentile' in stat:
        val = int(stat.split('percentile_')[1])
        da = da.load().quantile(val/100., dim=dims, skipna=skipna, keep_attrs=keep_attrs)
    return da

def apply_spatial_stats_on_full_ds(ds, stat, dims, varlist= None):
    """ Compute spatial stats on each variables from a dataset.
    Parameters:
        ds: xarray dataset
        stats: None or a string among the following : ['mean','min','max','median','std','percentile_xxx'], where _xxx is the value of the xxth percentile (xx in [00-100])
        dims: tuple of the dimensions on which stats are to be processed
    Result:
        ds_out: xarray dataset with the stat being computed on each of the variables from ds + adding the mean location.
    """
    ds_out = xr.Dataset()
    if varlist is None or varlist is []: varlist = list(ds.data_vars)+list(ds.coords)
    for varname in varlist :
        # apply stat on variable, variable-ERR and variable_snow if in ds_extr
        dims_to_average = [ d for d in list(ds[varname].dims) if d in dims ]
        if dims_to_average != []:
            ds_out[varname] = apply_stats_on_xr_da(ds[varname], stat=stat, dims=dims_to_average, skipna=True, keep_attrs=True)
    return ds_out


def filter_on_qflag(da, da_qf, qf_to_drop, fillvalue=np.NaN, remove_matching=True):
    """
    Removes the da pixels where qflag corresponding value is in qf_to_drop.
    Parameters:
        da, da_qf : xarray DataArrays with same shape. da: values to filter; da_qf : quality flag criteria used as the filtering criterion
        qf_to_drop: list of numeric values
        remove_matching: whether to remove or keep the matching pixels
    """
    if qf_to_drop != []:
        for x in qf_to_drop:
            if remove_matching:
                print(x, '=qf removed')
                da = da.where(da_qf != x, fillvalue) # keep info where da_qf != x, if not (ie where da_qf == x), then put fillvalue
            else:
                da = da.where(da_qf == x, fillvalue)
    return da

def filter_bitwise_on_qflag(da, da_qf, qf_bit_mask, qf_bit_results, fillvalue=np.NaN, remove_matching=True):
    """
    Removes the da pixels where qflag compared to bit_mask equals bit_values.
    Parameters:
        da, da_qf : xarray DataArrays with same shape. da: values to filter; da_qf : quality flag criteria used as the filtering criterion
        qf_bit_mask: [] or list of int or of binary values. Used to select bits from da_qf to use in comparison.
        qf_bit_results: [] or list of int or of binary values. Correspond to the results expected in the bitwise comparison of da_qf & bit_mask.
        remove_matching: whether to remove or keep the matching pixels
    """
    if qf_bit_mask != [] and (len(qf_bit_mask) == len(qf_bit_results)):
        for i,x in enumerate(qf_bit_mask):
            if remove_matching:
                print(x, '=qf removed')
                da = da.where(da_qf.astype(int) & qf_bit_mask[i] != qf_bit_results[i], fillvalue)
            else:
                da = da.where(da_qf.astype(int) & qf_bit_mask[i] == qf_bit_results[i], fillvalue)
    return da

def mask_bitwise_on_qflag(da, da_qf, qf_bit_mask, qf_bit_results):
    """
    Keeps the da pixels where qflag compared to bit_mask equals bit_values. Returns a mask with value 1 where ok, NaN elsewhere.
    Parameters:
        da, da_qf : xarray DataArrays with same shape. da: values to filter; da_qf : quality flag criteria used as the filtering criterion
        qf_bit_mask: [] or list of int or of binary values. Used to select bits from da_qf to use in comparison.
        qf_bit_results: [] or list of int or of binary values. Correspond to the results expected in the bitwise comparison of da_qf & bit_mask.
    """
    if qf_bit_mask != [] and (len(qf_bit_mask) == len(qf_bit_results)):
        for i,x in enumerate(qf_bit_mask):
            da = da.where(da_qf.astype(int) & qf_bit_mask[i] != qf_bit_results[i],-9999.8888)
        da = da.where(da == -9999.8888)*0+1
    return da



def filter_on_err(da, da_err, err_max_percent_lim = None, fillvalue=np.NaN):
    """
    Removes the da pixels where da_err/da*100 >  err_max_percent_lim.
    Parameters:
        da, da_err : xarray DataArrays with same shape. da: values to filter, da_err : uncertainty values used as the filtering criterion
        err_max_percent_lim: None or percent value in [0 - 100]. Upper valid limit for data_err/data*100
    """
    if err_max_percent_lim != None:
        da_out = da.where(da_err/da*100 <= err_max_percent_lim, fillvalue)
        da_out = da_out.where(np.isnan(da) != 1) # fills outside with NaNs as default value. > we report initial NAN values in output dataarray (used in pvlike codes when fillvalue = -1 for plotting missing data in grey and NAN as white )

    return da_out

def filter_on_aux_var(da, da_aux, valid_aux = [None,None], fillvalue=np.NaN):
    """
    Removes the da pixels where valid_aux[0]<= da_aux <= valid_aux[1].
    Parameters:
        da, da_aux : xarray DataArrays with same shape. da: values to filter, da_aux : values used as the filtering criterion
        valid_aux: 2-len list of min and max thresholds to apply on da_aux to filter da. Default [None,None]
    """
    if valid_aux != [None,None] or valid_aux !=[]:
        if valid_aux[0]!= None : da = da.where(da_aux >= valid_aux[0])
        if valid_aux[1]!= None : da = da.where(da_aux <= valid_aux[1])
    return da

def apply_scaling_offset_missing(da, scaling=1/10000,offset=0., valid_range = [0, 0.5], missing=[-1]):
    """
    First removes missing values. Then, applies scaling, offset. Finally, removes data outside valid_range provided.
    Parameters:
        da: xarray DataArray
        scaling : numeric value. multiplicative factor to apply to the data in da
        offset: numeric value. additive factor to apply to the data in da
        missing: [] or list of values. the missing values provided are expected to refer to data in da BEFORE applying any scaling/offset.
        valid_range: [] or list of two values. min and max (included) thresholds applied on the scaled data from da (threshold are applied AFTER applying scaling/offset to da). Data outside this valid range are returned as NaN.
    """
    # read missing values as NaN
    for m in missing:
        da = da.where(da != m)
    # apply scaling & offset
    da = da * scaling + offset
    # apply lower & upper limits
    if valid_range != []:
        if valid_range[0]!= None : da = da.where(da >= valid_range[0])
        if valid_range[1]!= None : da = da.where(da <= valid_range[1])
    return da

def apply_scaling_offset_missing_with_valid_before_scaling(da, dtype, scaling=1/10000,offset=0., valid_range = [0, 5000], missing=[-1]):
    """
    First removes missing values. Then, applies scaling, offset. Finally, removes data outside valid_range provided.
    Parameters:
        da: xarray DataArray
        scaling : numeric value. multiplicative factor to apply to the data in da
        offset: numeric value. additive factor to apply to the data in da
        missing: [] or list of values. the missing values provided are expected to refer to data in da BEFORE applying any scaling/offset.
        valid_range: [] or list of two values. min and max (included) thresholds applied on the raw data from da
            (threshold are applied BEFORE applying scaling/offset to da), as done for CF compliant data.
            Data outside this valid range are returned as NaN.
    """
    # read missing values as NaN
    if missing is not None:
        if not isinstance(missing,list):
            missing = [missing]
        for m in missing:
            da = da.where(da != m)
    # apply lower & upper limits
    if valid_range != []:
        if valid_range[0]!= None : da = da.where(da >= valid_range[0])
        if valid_range[1]!= None : da = da.where(da <= valid_range[1])

    # apply scaling & offset
    if scaling != 1.0 or offset != 0.:
        da = da * scaling + offset
        da = da.astype(dtype) #force scaling/offset conversion to keep da type

    return da

def filter_albedo_data(ds, varname, varname_qf, f_out_cache, qf_to_drop, qf_bit_mask, qf_bit_results,
                                       err_max_percent_lim, age_max_lim, nmod_min_lim_included_in_ok, do_snow=False,
                                       snow_qf_bit_mask=[], snow_qf_bit_results=[], varname_err=None,
                                       varname_nmod='nmod', varname_age='age_info'):
    """ Filter albedo data with QFLAG values/bits and ERR percentages,
        and creates a new variable varname+'_snow' which contains the values of variable flagged as snow. (used for snow flag plotting)
    """
    if varname_err is None: varname_err = varname+'_err'
    print(f'filtering on {varname_qf}')
    # filter on qflag if needed
    if qf_to_drop !=[]:
        ds[varname] = filter_on_qflag(ds[varname], ds[varname_qf], qf_to_drop)
    if qf_bit_mask != [] :
        if (len(qf_bit_mask) == len(qf_bit_results) ):
            ds[varname] = filter_bitwise_on_qflag(ds[varname], ds[varname_qf], qf_bit_mask = qf_bit_mask,
                                                                  qf_bit_results = qf_bit_results, fillvalue=np.nan)
        else:
            print('Error: qf_bit_mask & qf_bit_results lists have different lengths')

    # filter on err if needed
    if err_max_percent_lim != None and varname_err in ds.data_vars:
        print(f'filtering on {varname_err}')
        ds[varname] = filter_on_err(ds[varname],ds[varname_err], err_max_percent_lim, fillvalue=np.nan)

    # filter on z_age if needed
    if age_max_lim != None and varname_age in ds.data_vars:
        print(f'filtering on {varname_age} (keep info where age <= {age_max_lim})')
        ds[varname] = filter_on_aux_var(ds[varname], ds[varname_age], valid_aux = [None,age_max_lim], fillvalue=np.NaN)


    # filter on nmod if needed
    if nmod_min_lim_included_in_ok != None and varname_nmod in ds.data_vars:
        print(f'filtering on {varname_nmod} (keep info where nmod >= {nmod_min_lim_included_in_ok})')
        ds[varname] = filter_on_aux_var(ds[varname], ds[varname_nmod], valid_aux = [nmod_min_lim_included_in_ok, None], fillvalue=np.NaN)

    # extract snow data
    if do_snow:
        print(f'creating snow var: {varname}_snow')
        v_snow = varname +'_snow'
        if v_snow not in ds.data_vars: ds[v_snow] = filter_bitwise_on_qflag(ds[varname], ds[varname_qf],
                                         qf_bit_mask = snow_qf_bit_mask, qf_bit_results = snow_qf_bit_results,
                                         fillvalue=np.nan, remove_matching=False)

#    dict_encode={'zlib':True, 'complevel':4}
#    for v in ds.data_vars:
#        ds[v].encoding = dict_encode
    if f_out_cache is not None:
        print('begin storing cache file of scaled/filtered data')
        ds.to_netcdf(f_out_cache)
        print('end storing cache file of scaled/filtered data')
    else:
        return ds
    ds.close()
