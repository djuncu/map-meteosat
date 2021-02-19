#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 11:50:14 2018

@author: vincentch
"""
from datetime import datetime,timedelta #,date
import os
import numpy as np
import PIL
import matplotlib.pyplot as plt

from .MapProj import wgs2local

yyyymmdd_2_yyyyddd = lambda date_yyyymmdd : datetime.strptime(date_yyyymmdd, '%Y%m%d').strftime('%Y%j')
dec2byte = lambda x:np.binary_repr(x,width=8)
byte2dec = lambda x:int(x,2)

def resampleMesh(lonval, latval, da_to_plot, plotXstep=None, plotYstep=None,
                 ixmin=None, ixmax=None, iymin=None, iymax=None, da_ocean=[], oceanVal=-1):
    ## resample values for pcolormesh
    # juncud 09/2020

    if type(da_to_plot).__module__ == 'numpy':
        v = da_to_plot
        if len(da_ocean) != 0:
            ocean  = da_ocean
    else:
        v = da_to_plot.values
        if len(da_ocean) != 0:
            ocean  = da_ocean.values

    if len(lonval.shape) == 2:
        lonval = lonval[iymin:iymax:plotYstep, ixmin:ixmax:plotXstep] # right order?
        latval = latval[iymin:iymax:plotYstep, ixmin:ixmax:plotXstep]
    elif len(lonval.shape) == 3:
        lonval = lonval[0, iymin:iymax:plotYstep, ixmin:ixmax:plotXstep] # right order?
        latval = latval[0, iymin:iymax:plotYstep, ixmin:ixmax:plotXstep]

    shape1 = lonval.shape[0]
    shape2 = lonval.shape[1]
    lon_resamp = np.zeros([shape1-1,shape2-1])
    lat_resamp = np.zeros([shape1-1,shape2-1])

    v_resamp = np.empty([shape1-2,shape2-2])
    #ocean_resamp  = np.zeros([lonval.shape[0]-2,lonval.shape[1]-2])
    for i in range(shape1-1):
        for j in range(shape2-1):
            lon_resamp[i,j] = np.mean(lonval[i:i+2,j:j+2])
            lat_resamp[i,j] = np.mean(latval[i:i+2,j:j+2])

            if j != shape2-2 and i != shape1-2:
                if lon_resamp[i,j] < 1e3:
                    v_resamp[i,j] = v[i+1,j+1]
                    #ocean_resamp[i,j] = ocean[i+1,j+1]

                    if len(da_ocean) != 0 and ocean[i+1,j+1] == 1.:
                        v_resamp[i,j] = oceanVal
                else:
                    v_resamp[i,j] = np.nan
                    #ocean_resamp[i,j]  = np.nan

    v_masked =  np.ma.masked_invalid(v_resamp)

    return lon_resamp, lat_resamp, v_masked

def lonlat2indexLims(lonvals, latvals, lonmin, lonmax, latmin, latmax):
    lonvals = lonvals.squeeze()
    latvals = latvals.squeeze()
    lonvals_r = np.reshape(lonvals, (np.size(lonvals), 1))
    latvals_r = np.reshape(latvals, (np.size(latvals), 1))

    indexArray = np.indices(lonvals_r.shape)[0]

    lonlims = (lonmin, lonmax)
    latlims = (latmin, latmax)

    # only use 'neighborhood'  coordinates
    k = 0
    iLat = np.zeros((4,),dtype=int)
    iLon = np.zeros((4,),dtype=int)
    for i in range(2):
        for j in range(2):
            # lower left - upper left - lower right - upper right
            nhMask = np.logical_and(lonvals_r > lonlims[i]-1., lonvals_r < lonlims[i]+1.) & \
                     np.logical_and(latvals_r > latlims[j]-1., latvals_r < latlims[j]+1.)
            maskedIndexArray = indexArray[nhMask]

            xy = wgs2local(np.hstack((lonvals_r[nhMask, None], latvals_r[nhMask, None])),
                                 [lonlims[i], latlims[j]])
            dist = np.sqrt(np.sum(xy**2,1))

            iLat[k], iLon[k] = np.unravel_index(maskedIndexArray[np.argmin(dist)],
                                                    lonvals.shape)
            k += 1

    iLonMin = np.min((iLon[0], iLon[1]))
    iLonMax = np.max((iLon[2], iLon[3]))
    # lat indices start from top of the map
    iLatMax = np.max((iLat[0], iLat[2]))
    iLatMin = np.min((iLat[1], iLat[3]))
    # make the area a bit bigger so that there are not cut offs due to resampling. careful if this is used for cases where some of the limits are 'None'
    return iLonMin-2, iLonMax+2, iLatMin-2, iLatMax+2


def instanciate_keywords(value, keywords):
    """ This function take a dict, list or string as input, it uses a keyword dictionnary to instanciate the values of some keywords in the object. Note that if the input is a dict, the dict itself is modified.
    **Examples**

    >>> instanciate_keywords('{name}/static-string/{param}.txt', {'name':'foo', 'param':'bar'})
    'foo/static-string/bar.txt'
    """
    if isinstance(value, dict):
       for k in list(value.keys()):
           v = value[k]
           del value[k]
           newk = instanciate_keywords(k, keywords)
           value[newk] = instanciate_keywords(v, keywords)
    if isinstance(value, list):
       value = [instanciate_keywords(v, keywords) for v in value]
       return value
    if isinstance(value, str):
        value = value.format(**keywords)
        return value
    try: # ensure python2.7 compat ugly.
        if isinstance(value, unicode):
            value = value.format(**keywords)
            return value
    except NameError:
        pass # end-of-ugly
    return value

def instanciate_datetime(value, date):
    """ Similar to instanciate_keywords. This function take a dict, list or string as input, it uses a date to instanciate the value in the object. Note that if the input is a dict, the dict itself is modified.
    Examples:

    >>> from datetime import datetime
    >>> instanciate_datetime('%Y/%d/%m/HDF_%Y_%d_%m', datetime(2017,5,12))
    '2017/12/05/HDF_2017_12_05'
    >>> instanciate_datetime({'brdf': 'BRDF_FILENAME_%Y_%d_%m', 'albedo': 'ALBEDO_FILENAME_%Y_%d_%m'}, datetime(2017,5,12))
    {'brdf': 'BRDF_FILENAME_2017_12_05', 'albedo': 'ALBEDO_FILENAME_2017_12_05'}

    """
    if isinstance(value, dict):
       for k in list(value.keys()):
           v = value[k]
           del value[k]
           newk = instanciate_datetime(k, date)
           value[newk] = instanciate_datetime(v, date)
    if isinstance(value, list):
       value = [instanciate_datetime(v, date) for v in value]
       return value
    if isinstance(value, str):
        value = date.strftime(value)
        return value
    try: # ensure python2.7 compat ugly.
        if isinstance(value, unicode):
            value = date.strftime(value)
            return value
    except NameError:
        pass # end-of-ugly
    return value

def ensure_dir(filename):
    """ Creates the basename directory of a file if it does not exist.
    If the argument is a directory, finish with '/' to ensure this directory is created.
    Otherwise only the parent dir will be created"""
    try:
        dirname = os.path.dirname(filename)
    except:
        dirname = ''

    if dirname != '' and not os.path.exists(dirname):
        try:
            os.makedirs(dirname)
        except Exception as e:
            if not os.path.exists(dirname):
                raise(e)

def replace_keywords_in_template_str(str_template, **kwargs):
    """Replaces the keys words given between '{}' in the input string by
    the values corresponding to this keyword in the dictionary given as argument """
    import re
    filebasename = None
    for kk,vv in kwargs.items():
        if kk in str_template :
            if filebasename == None:
                filebasename = re.sub(r'{'+kk+'}', vv, str_template)
            else:
                filebasename = re.sub(r'{'+kk+'}', vv, filebasename)
    return filebasename

def transfo_str_to_fct_ref(str_fct, list_possible_values):
    """ Transforms the text containing a function name into a reference to the function itself
        secured operation as we check first that user input matches with the allowed values.
        In particular we use this to transform user's config file resampling functions names.
        NB: this code shall be called in the same module as where the fcts are defined. """
    if str_fct in list_possible_values:
        return eval(str_fct)

def replace_dict_strings_by_existing_functions_from_a_class(dic,classinit):
    """The user provides the name of the functions to use in the config file. These are given as strings.
    This code is used to replace the function names by their actual address (within the target class object given as argument)."""
    for key,val in dic.items():
        if val != None:
            f= [ getattr(classinit, attr) for attr in dir(classinit) if callable(getattr(classinit, attr)) and (not attr.startswith("__")) and (val in attr)]
            dic[key]=f[0]

def dates_between_start_end_v2(date_start, date_end):
    """ Returns the list of datetime dates between two dates provided in '%Y-%m-%d' string format"""
    date_i0 = datetime.strptime(str(date_start),'%Y-%m-%d') # datetime object
    date_if = datetime.strptime(str(date_end),'%Y-%m-%d') # datetime object
    delta = date_if - date_i0
    date_list=[]
    for i in range(delta.days + 1):
        date_list.append(date_i0 + timedelta(i))
    return date_list

def add_logo_to_figure(fig,ax,logo_path='/cnrm/vegeo/SAT/CODES/internal/general/visualisation_validation_tools/aux_data_for_plots/logos/logo_MF2.png',xycoords='figure fraction', format='png', zoom=0.12, xy_position_tuple= (0.88,0.95)):
    from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
    arr_img = plt.imread(logo_path, format=format)
    imagebox = OffsetImage(arr_img, zoom=zoom)
    ab = AnnotationBbox(imagebox, xy=xy_position_tuple, xycoords= xycoords, frameon=False)
    # see https://matplotlib.org/users/annotations.html for info on annotation coordinates
    ax.add_artist(ab)
    return fig,ax

#def add_logo_to_imagefile(mainfilename, logofilename, outfilename):
#    mimage = PIL.Image.open(mainfilename)
#    limage = PIL.Image.open(logofilename)
#    # resize logo
#    wsize = int(min(mimage.size[0], mimage.size[1]) * 0.27) #0.12
#    wpercent = (wsize / float(limage.size[0]))
#    hsize = int((float(limage.size[1]) * float(wpercent)))
#    simage = limage.resize((wsize, hsize))
#    mbox = mimage.getbbox() # returns left, upper, right, and lower pixel coordinate.
#    sbox = simage.getbbox()
#    # right top corner
#    box = (mbox[0]+int(0.04*abs(mbox[2]-mbox[0]))- sbox[0], mbox[3]-int(0.12*abs(mbox[3]-mbox[1]))- sbox[3])  # bottom left
#    #box = (mbox[2] - sbox[2], mbox[1] - sbox[1]) right top
#    mimage.paste(simage, box)
#    mimage.save(outfilename)
#
#
#def paste_two_logos_horiz(logo1, logo2, logo_merged):
#    im1=PIL.Image.open(logo1)
#    im2=PIL.Image.open(logo2)
#    im1out=im1
#    im2out=im2
#
#   # resize_two_logos(im1, im2, im1out, im2out)
#    widths, heights = zip(*(i.size for i in [im1out, im2out]))
#
#    total_width = sum(widths)
#    max_height = max(heights)
#
#    new_im = PIL.Image.new('RGB', (total_width, max_height))
#    data=np.array(new_im)
#
#    data[:,:]=[255,255,255] # fills background colour with RGB white (instead of black)
#    new_im2=PIL.Image.fromarray(data)    # creates the image from the array
#
#    # we center vertically the smallest image
#    y_offset_smaller=int(max_height/2)
#
#    new_im2.paste(im1out,(0,int(y_offset_smaller-im1out.size[1]/2)))
#    new_im2.paste(im2out,(im1out.size[0],int(y_offset_smaller-im2out.size[1]/2)))
#
#    #new_im2.show()
#    new_im2.save(logo_merged)


## geographic coords extractions

def get_punctual_value_from_da(da,longit,latit,lonname='lon',latname = 'lat'):
    #value = da.sel(lon=[longit], method='nearest', drop=True).sel(lat=[latit], method='nearest', drop=True).values.ravel()
    value = da.sel({lonname:[longit]}, method='nearest', drop=True).sel({latname:[latit]}, method='nearest', drop=True).values.ravel()
    return value

def get_punctual_value_from_da_with_1D_lon_lat_dims(da,longit,latit,lonname='lon',latname = 'lat'):
    #value = da.sel(lon=[longit], method='nearest', drop=True).sel(lat=[latit], method='nearest', drop=True).values.ravel()
    value = da.sel({lonname:[longit]}, method='nearest', drop=True).sel({latname:[latit]}, method='nearest', drop=True).values.ravel()
    return value

def get_point_indices_from_metop(ds_longlat, lonsite, latsite):
    ds_point = ds_longlat.sel(y=latsite, method='nearest').dropna(dim='x')
    ds_point.coords['x'] = ds_point['lon2D'].values
    ds_point = ds_point.sel(x=lonsite, method='nearest')
    #print(ds_point)
    return ds_point.ix.values, ds_point.iy.values

def get_zone_indices_from_metop_macropixel_radius(ds_longlat, lonsite, latsite, radius_x, radius_y):
    """
    radius_x : macropixel size along x dimension is (2*radius_x + 1)
    radius_y : macropixel size along y dimension is (2*radius_y + 1)
    """
    ix, iy = get_point_indices_from_metop(ds_longlat, lonsite, latsite)
    l_ix = np.arange(ix-radius_x, ix+radius_x +1, 1)
    l_iy = np.arange(iy-radius_y, iy+radius_y +1, 1)
    return l_ix, l_iy

def get_zone_2D_mask_from_zone_limits(lons2D, lats2D, lonsite_w, lonsite_e, latsite_s, latsite_n, inorout='inside'):
    """
    Returns a 2D mask of the lons/lats indices which are NOT within the zone desbribed with lonsite_w, lonsite_e, latsite_s, latsite_n limits
    Longitudes and longitude limits (lonsite_w, lonsite_e) shall be in degrees East;
    Latitudes and latitude limits (latsite_s, latsite_n) in degrees north.
    """
    if inorout == 'inside': # keep indices of points inside the zone
        mask_zone = np.where(np.logical_and(np.logical_and(lons2D <= lonsite_e,lons2D >= lonsite_w),np.logical_and(lats2D <= latsite_n,lats2D >= latsite_s)))
    else: # keep indices of points outside the zone
        mask_zone = np.where(np.logical_not(np.logical_and(np.logical_and(lons2D <= lonsite_e,lons2D >= lonsite_w),np.logical_and(lats2D <= latsite_n,lats2D >= latsite_s))))
    return mask_zone

def get_zone_2x1D_mask_from_zone_limits(lons1D, lats1D, lonsite_w, lonsite_e, latsite_s, latsite_n):
    """
    Returns two 1D masks of the lons/lats indices which are NOT within the zone desbribed with lonsite_w, lonsite_e, latsite_s, latsite_n limits
    Longitudes and longitude limits (lonsite_w, lonsite_e) shall be in degrees East;
    Latitudes and latitude limits (latsite_s, latsite_n) in degrees north.
    """
    mask_outside_zone_lons = np.where(np.logical_not(np.logical_and(lons1D <= lonsite_e,lons1D >= lonsite_w)))
    mask_outside_zone_lats = np.where(np.logical_not(np.logical_and(lats1D <= latsite_n,lats1D >= latsite_s)))
    return mask_outside_zone_lons,mask_outside_zone_lats


def get_max_zone_indices_from_macropixel_limits(ds_longlat, lonsite_w, lonsite_e, latsite_s, latsite_n):
    """
    """
    dlat = ds_longlat['y'][1] - ds_longlat['y'][0]
    if dlat <0: latsite_s,latsite_n = latsite_n,latsite_s
    ds_zone = ds_longlat.sel(y=slice(latsite_s,latsite_n))#.dropna(dim='x') # remove nan longitudes
    ds_zone = ds_zone.where(ds_zone['lon2D']<=lonsite_e, drop=True)
    ds_zone = ds_zone.where(ds_zone['lon2D']>=lonsite_w, drop=True)
    # TODO : could also return the list of valid ix, iy
    return ds_zone['ix'].min().values,ds_zone['ix'].max().values,ds_zone['iy'].min().values,ds_zone['iy'].max().values

def filter_metop_on_zone_from_macropixel_limits(ds_geoloc, ds_longlat, lonsite_w, lonsite_e, latsite_s, latsite_n, step_ix=1, step_iy=1, keep_only_points_in_limits=True):
    """
    """
    ix0,ixf,iy0,iyf = get_max_zone_indices_from_macropixel_limits(ds_longlat, lonsite_w, lonsite_e, latsite_s, latsite_n)
    print(ix0,ixf,iy0,iyf)
    ds_zone = ds_geoloc.isel(y=slice(iy0,iyf,step_iy)).isel(x=slice(ix0,ixf,step_ix))
    print(ds_zone)
    if keep_only_points_in_limits:
        # nanify the remaining pixels which are within the max limits but where longitude is not in the required range
        ds_zone = ds_zone.where(ds_zone['lon2D']<=lonsite_e, drop=True)
        ds_zone = ds_zone.where(ds_zone['lon2D']>=lonsite_w, drop=True)
    return ds_zone
