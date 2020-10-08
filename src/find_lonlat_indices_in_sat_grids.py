#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 11:27:20 2018

@author: vincentch
"""
from scipy.spatial import cKDTree
import numpy as np
import h5py
from itertools import product


def find_lonlat_indices_in_2D_np(lon2D,lat2D,lon_min,lon_max,lat_min,lat_max):
    matching_indices = np.where(np.logical_and(lat2D>=lat_min,lat2D<=lat_max) & np.logical_and(lon2D>=lon_min,lon2D<=lon_max))
    return matching_indices

def write_indices_to_hdf(filename_out, site_info, on_zone, matching_indices, matching_lons, matching_lats, npix_x=None, npix_y=None, write_mode='a'):
    
    # verify content in site_info, add default info if not provided
    if 'site_name' not in site_info.keys(): site_info['site_name'] = 'site_name'
    if 'site_network' not in site_info.keys(): site_info['site_network'] = 'unknown_network'
    if 'site_index' not in site_info.keys(): site_info['site_index'] = 0

    # get site info from arguments
    site_name = site_info['site_name']
    site_lon = site_info['site_lon']
    site_lat = site_info['site_lat']
    site_network = site_info['site_network']
    site_index = site_info['site_index']
    
    # write results to file
    h = h5py.File(filename_out,write_mode)
    grp = h.create_group(site_name)
    grp.create_dataset('ilon_msg', data = matching_indices[1], compression="gzip", compression_opts=9, dtype='i')
    grp.create_dataset('ilat_msg', data = matching_indices[0], compression="gzip", compression_opts=9, dtype='i')
    grp.create_dataset('lon_msg', data = matching_lons, compression="gzip", compression_opts=9, dtype='f')
    grp.create_dataset('lat_msg', data = matching_lats, compression="gzip", compression_opts=9, dtype='f')
    dt = h5py.special_dtype(vlen=str) # workaround for h5py issue https://github.com/h5py/h5py/issues/289
    grp.attrs.create(name = 'site_network', data = site_network, dtype=dt) ## TODO better specify type 
    grp.attrs.create(name = 'attrs_missing_value', data = -999, dtype='f')
    grp.attrs.create(name = 'site_index', data = site_index, dtype='i')

    if on_zone:
        site_lon_min,site_lon_max = site_lon
        site_lat_min,site_lat_max = site_lat
        grp.attrs.create(name = 'lon_site', data = -999, dtype='f')
        grp.attrs.create(name = 'lat_site', data = -999, dtype='f')
        grp.attrs.create(name = 'npix_x_msg', data = -999, dtype='i')        
        grp.attrs.create(name = 'npix_y_msg', data = -999, dtype='i')       
        grp.attrs.create(name = 'lon_min_msg', data = min(matching_lons), dtype='f')
        grp.attrs.create(name = 'lat_min_msg', data = min(matching_lats), dtype='f')
        grp.attrs.create(name = 'lon_max_msg', data = max(matching_lons), dtype='f')
        grp.attrs.create(name = 'lat_max_msg', data = max(matching_lats), dtype='f')
        grp.attrs.create(name = 'lon_min_site', data = site_lon_min, dtype='f')
        grp.attrs.create(name = 'lat_min_site', data = site_lat_min, dtype='f')
        grp.attrs.create(name = 'lon_max_site', data = site_lon_max, dtype='f')
        grp.attrs.create(name = 'lat_max_site', data = site_lat_max, dtype='f')
    else:
        grp.attrs.create(name = 'lon_site', data = site_lon, dtype='f')
        grp.attrs.create(name = 'lat_site', data = site_lat, dtype='f')
        grp.attrs.create(name = 'npix_x_msg', data = npix_x, dtype='i')        
        grp.attrs.create(name = 'npix_y_msg', data = npix_y, dtype='i')
        grp.attrs.create(name = 'lon_min_msg', data = min(matching_lons), dtype='f')
        grp.attrs.create(name = 'lat_min_msg', data = min(matching_lats), dtype='f')
        grp.attrs.create(name = 'lon_max_msg', data = max(matching_lons), dtype='f')
        grp.attrs.create(name = 'lat_max_msg', data = max(matching_lats), dtype='f')
        grp.attrs.create(name = 'lon_min_site', data = -999, dtype='f')
        grp.attrs.create(name = 'lat_min_site', data = -999, dtype='f')
        grp.attrs.create(name = 'lon_max_site', data = -999, dtype='f')
        grp.attrs.create(name = 'lat_max_site', data = -999, dtype='f')

    h.close()     

def read_indices_from_hdf_file(filename, site_name):
    h = h5py.File(filename,'r')
    ilon_msg = h[site_name]['ilon_msg'][:]
    ilat_msg = h[site_name]['ilat_msg'][:]
    lon_msg = h[site_name]['lon_msg'][:]
    lat_msg = h[site_name]['lat_msg'][:]
#    #for verif only:
#    print('verif indices extraction & reading:')
#    print(lats_msg[ilat_msg,ilon_msg],lat_msg)
#    print(lons_msg[ilat_msg,ilon_msg],lon_msg)
    h.close()
    return ilat_msg,ilon_msg,lat_msg,lon_msg


def find_lonlat_indices_in_msg(lonval, latval, filename_out, site_info, npix_x=0, npix_y=0, on_zone=False, write_mode='a',other_sats_to_extract=[]): 
    """ 
    Find indices of geolocated sites in MSG grid.
        [this code is based on suman's code best_pixel_identify_ground_stations_userdefined_seviri, in /cnrm/vegeo/SAT/CODES/GROUND_ALBEDO_VALIDATION_git/ ]
    Parameters:
    lonval,latval: 2D numpy arrays of MSG lons and lats for the whole global map. Indices refer to [ilat,ilon].
    filename_out: string. Path and filename for output hdf5 file.
    write_mode: string chosen among [ 'a', 'w', 'a+', 'w+']. Defines the writing mode for output file `filename_out`.
    npix_x, npix_y: int. Used only if 'on_zone'== False. Number of pixels to extract respectively along lon,lat on each side of the central pixel at site_lon/site_lat position.
        e.g. if npix_x = 2 and npix_y = 1, the total number of pixels provided is (2*npix_x+1)*(2*npix_y+1)= 5*3 = 15
    on_zone: boolean. Defines the extraction mode.
        if True: over a zone defined (site_lon_min/site_lon_max in °E, site_lat_min/site_lat_max in °N) coordinates. In this case, site_lon,site_lat are tuples of min/max coordinates.
        if False: extraction of either a single point defined by its (site_lon in °E,site_lat in °N) coordinates or a number of pixels around a single point.
        or a zone defined by (site_lon in °E,site_lat in °N) coordinates and a number of pixels along lon (npix_x), and a number of pixels along lat (npix_y),
    site_info = dict with at least ['site_name','site_lon','site_lat'] keys. 'site_index' and 'site_network' can also be provided. example could be {'site_name':'default_station_name', 'site_lon':0,'site_lat':0,'site_network':'default_network_name'}
    other_sats_to_extract: list of strings within [ 'metop_avhrr','modis','spot_vgt','probav'] # todo implement this function
        enables the search for corresponding indexes in the listed sat grids for ONE MSG pixel. Only activated when on_zone=False and 'npix_x'='npix_y'=0, or when on_zone=True.
    """
    # get site info from arguments
    site_lon = site_info['site_lon']
    site_lat = site_info['site_lat']

    ## find matching pixels
    if on_zone:
        lowlon, uplon = site_lon
        lowlat, uplat = site_lat
        try:
            matching_indices = find_lonlat_indices_in_2D_np(lonval,latval,lon_min=lowlon,lon_max=uplon,lat_min=lowlat,lat_max=uplat)
#            matching_indices = np.where(np.logical_and(latval>=lowlat,latval<=uplat) & np.logical_and(lonval>=lowlon,lonval<=uplon))
            matching_lats = latval[matching_indices]
            matching_lons = lonval[matching_indices]
        except:
            print ("error "+str(IOError))
            print ('Couldnt find MSG pixel(s) on zone \n{site_info}')
            pass                
    else:        
        # define a small zone (+/-0.5°) around required point to have enough points to get a matching MSG pixel but avoid searching the whole globe
        # and get the matching indices
        lowlat = site_lat - 0.5
        uplat = site_lat + 0.5
        lowlon = site_lon - 0.5
        uplon = site_lon + 0.5        
        indx_req_grid=np.where(np.logical_and(latval>=lowlat,latval<=uplat) & np.logical_and(lonval>=lowlon,lonval<=uplon))
        if len(indx_req_grid[0])>1:
            try:
                lat_first_index=min(indx_req_grid[0])
                lat_last_index=max(indx_req_grid[0])
                lon_first_index=min(indx_req_grid[1])
                lon_last_index=max(indx_req_grid[1])

                # create 2D arrays of lon & lat indices of the zone limits (a rectangular area in MSG grid)
                lon_msg_pos=np.arange(lon_first_index,lon_last_index,1)
                lat_msg_pos=np.arange(lat_first_index,lat_last_index,1)            
                latmsg_posmesh,lonmsg_posmesh=np.meshgrid(lat_msg_pos,lon_msg_pos)
                
                latmsg_posmesh=np.transpose(latmsg_posmesh)
                lonmsg_posmesh=np.transpose(lonmsg_posmesh)
                
                latmsgpos_reshp=np.reshape(latmsg_posmesh,[latmsg_posmesh.shape[0]*latmsg_posmesh.shape[1]])
                lonmsgpos_reshp=np.reshape(lonmsg_posmesh,[lonmsg_posmesh.shape[0]*lonmsg_posmesh.shape[1]])                 
                
                # get the 2D lon & lat arrays on the rectangular MSG subset zone
                lon_msg=lonval[lat_first_index:lat_last_index,lon_first_index:lon_last_index]
                lat_msg=latval[lat_first_index:lat_last_index,lon_first_index:lon_last_index]  
               
                latmsg_reshp=np.reshape(lat_msg,[lat_msg.shape[0]*lat_msg.shape[1]])
                lonmsg_reshp=np.reshape(lon_msg,[lon_msg.shape[0]*lon_msg.shape[1]])     
                
                # merge lon & lat in a single array
                latlon_msgbox=np.array([latmsg_reshp,lonmsg_reshp])
                latlon_msgbox=np.transpose(latlon_msgbox)
                
                latlon_posmsgbox=np.array([latmsgpos_reshp,lonmsgpos_reshp])
                latlon_posmsgbox=np.transpose(latlon_posmsgbox)
                
                # remove nan from search (in lons & lats and corresponding indices )
                indx_nan=np.where(np.isnan(latlon_msgbox[:,0])!=1)                               
                latlon_msgbox_removenan=latlon_msgbox[indx_nan]
                latlon_pos_msgbox_removenan=latlon_posmsgbox[indx_nan]
                
                # get nearest point index searching on msg subseted lons&lats
                tree = cKDTree(latlon_msgbox_removenan)
                dists, indexes = tree.query(np.array([site_lat,site_lon]), k=1)
                # find corresponding indices in the global lons & lats arrays
                required_position = latlon_pos_msgbox_removenan[indexes]
                #print('found_pos:',required_position)
                
                # if required, add npix_x and npix_y points around the position found
                # first make sure the indices remain in the lon/lat array lengths
                min_iy = max(0, required_position[0]-npix_y)
                max_iy = min(len(latval)-1,required_position[0]+npix_y)
                min_ix = max(0, required_position[1]-npix_x)
                max_ix = min(len(latval[0])-1,required_position[1]+npix_x)                
#                print(min_iy,required_position[0]-npix_y, max_iy,required_position[0]+npix_y, latval.shape)
#                print(min_ix,required_position[1]-npix_x, max_ix,required_position[1]+npix_x, lonval.shape)
                
                # then get matching indices & corresponding lats/lons
                grid_latpos=np.arange(min_iy,max_iy+1,1)
                grid_lonpos=np.arange(min_ix,max_ix+1,1)
                matching_indices = tuple(np.transpose(list(product(grid_latpos,grid_lonpos)))) 
                matching_lons = lonval[matching_indices]
                matching_lats = latval[matching_indices]
                
                # verify lon & lat extracted in subset
#                temp_lalon_msg=latlon_msgbox_removenan[indexes]
#                print('found pos:',matching_indices)
#                print('matching lon,lats',matching_lons,matching_lats, temp_lalon_msg)
#                print( np.array([matching_indices, matching_lats, matching_lons]))
#                print('found_pos:',tuple(required_position), 'found_coords:',tuple(temp_lalon_msg))#, tuple( latval[required_position[0],required_position[1]],lonval[required_position[0],required_position[1]] )
            except:
                print ("error "+str(IOError))
                print (f'Couldnt find MSG pixel(s) around point \n {site_info}')
                pass
  
    # write results to file
#    print(site_info)
#    print('matching i: ',matching_indices)
    write_indices_to_hdf(filename_out, site_info, on_zone, matching_indices, matching_lons, matching_lats, npix_x=npix_x, npix_y=npix_y, write_mode=write_mode)
    
        
if __name__ == '__main__' :
    
    # test with dummy MSG grid & extraction sites
    lon1D=np.arange(9.5,11.4,0.2)
    lat1D=np.arange(25,26,0.2)
    latval, lonval = np.meshgrid(lat1D, lon1D, indexing='ij') # indexing example in reverse order to match with usual input MSG/hdf5 ordering
    print('latval:',latval)
    print('lonval:',lonval)
    
    print('--- Test on zone ---')
    site_info = {'site_name':'test_station_name', 'site_index':0, 'site_lon':(10,12),'site_lat':(24,25.3),'site_network':'test_network'}    
    find_lonlat_indices_in_msg(lonval=lonval, latval=latval, site_info = site_info, 
                               filename_out='./test_indices_on_zone.h5', npix_x=0, npix_y=0, on_zone=True, write_mode='w')

    print('--- Test on single point ---')
    site_info = {'site_name':'test_station_name', 'site_index':0, 'site_lon':10.002,'site_lat':25.321,'site_network':'test_network'}
    find_lonlat_indices_in_msg(lonval=lonval, latval=latval, site_info = site_info, 
                               filename_out='./test_indices_point.h5', npix_x=0, npix_y=0, on_zone=False, write_mode='w')
    print('--- Test on single point with npix_x, npix_y ---')
    find_lonlat_indices_in_msg(lonval=lonval, latval=latval, site_info = site_info, 
                               filename_out='./test_indices_point_2x3y.h5', npix_x=10, npix_y=3, on_zone=False, write_mode='w')
    
    # test with real sites on MSG files
    import rebuild_lonlat
    lons_msg,lats_msg = rebuild_lonlat.read_msg_lonlats()
    
    print('--- Test on zone / real MSG lonlats ---')
    site_info = {'site_name':'test_station_name', 'site_lon':(12.2,13.457),'site_lat':(49.001,50.3654)}
    f_extr = f'./msg_indices_zone_{list(site_info["site_lon"])[0]}_{list(site_info["site_lon"])[1]}E_{list(site_info["site_lat"])[0]}_{list(site_info["site_lat"])[1]}N.h5'
    find_lonlat_indices_in_msg(lonval=lons_msg, latval=lats_msg, site_info = site_info, 
                               filename_out=f_extr, npix_x=0, npix_y=0, on_zone=True, write_mode='w')
    print('--- Test on single point / real MSG lonlats ---')
    site_info = {'site_name':'test_station_name', 'site_lon':12.2,'site_lat':50.3654}
    f_extr = f'./msg_indices_nearest_{site_info["site_lon"]}E_{site_info["site_lat"]}N.h5'
    find_lonlat_indices_in_msg(lonval=lons_msg, latval=lats_msg, site_info = site_info, 
                               filename_out=f_extr, npix_x=0, npix_y=0, on_zone=False, write_mode='w')
    print('--- Test on single point with npix_x, npix_y / real MSG lonlats ---')
    site_info = {'site_name':'test_station_name', 'site_lon':12.2,'site_lat':50.3654}
    npix_x = 2
    npix_y = 3
    f_extr = f'./msg_indices_nearest_{site_info["site_lon"]}E_{site_info["site_lat"]}N_{npix_x}x{npix_y}y.h5'
    find_lonlat_indices_in_msg(lonval=lons_msg, latval=lats_msg, site_info = site_info, 
                               filename_out=f_extr, npix_x=npix_x, npix_y=npix_y, on_zone=False, write_mode='w')