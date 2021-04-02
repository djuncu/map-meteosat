import numpy as np
import os
from mxgplotlib.MapProj import wgs2local

def getIndex(cfg, lonlatFileName, lonval, latval, gen_):

    if gen_ == '3':
        import xarray as xr

        if lonlatFileName is not None and os.path.isfile(lonlatFileName):
            fn = lonlatFileName
        else:
            fn = cfg['lonlatFileMtg']

        #print(f'Coordinates at y/lat: {args["iy"]} and x/lon: {args["ix"]}')
        #print(f'LON:  {xr.open_dataset(fn)["LON"].values.squeeze()[args["iy"],args["ix"]]}')
        #print(f'LAT:  {xr.open_dataset(fn)["LAT"].values.squeeze()[args["iy"],args["ix"]]}')
    elif gen_ == '2':
        import h5py

        #print(f'Coordinates at y/lat: {args["iy"]} and x/lon: {args["ix"]}')
        with h5py.File(cfg['lonfile'],'r') as lonFileContent:
            lon = lonFileContent['LON'][:] / lonFileContent['LON'].attrs.get('SCALING_FACTOR')
        
        with h5py.File(cfg['latfile'],'r') as latFileContent:
            lat = latFileContent['LAT'][:] / latFileContent['LAT'].attrs.get('SCALING_FACTOR')


    lon_rsmpl = np.reshape(lon, (np.size(lon), 1))
    lat_rsmpl = np.reshape(lat, (np.size(lat), 1))
    indexArray = np.indices(lon_rsmpl.shape)[0]
    
    neighborMask = np.logical_and(lon_rsmpl > lonval-1., lon_rsmpl < lonval+1.) & \
                   np.logical_and(lat_rsmpl > latval-1., lat_rsmpl < latval+1.)
    xy = wgs2local(np.hstack((lon_rsmpl[neighborMask,None], lat_rsmpl[neighborMask,None])), 
                   [lonval, latval])
    
    neighborIndexArray = indexArray[neighborMask]
    dist = np.sqrt(np.sum(xy**2, 1))

    minIndex1d = neighborIndexArray[dist==np.min(dist)]

    iy, ix = np.unravel_index(minIndex1d, lon.shape)
    
    print(f'Nearest indices for coordinates:')
    print(f'ix: {ix[0]}')
    print(f'iy: {iy[0]}')

    return
