import numpy as np
import os
from mxgplotlib.MapProj import wgs2local

def getCoords(cfg, lonlatFileName, ix, iy, gen_):
    if gen_ == '3':
        import xarray as xr

        if lonlatFileName is not None and os.path.isfile(lonlatFileName):
            fn = lonlatFileName
        else:
            fn = cfg['lonlatFileMtg']

        lon = xr.open_dataset(fn)['LON'].isel(time=0).values[iy, ix]
        lat = xr.open_dataset(fn)['LAT'].isel(time=0).values[iy, ix]

    elif gen_ == '2':
        import h5py

        with h5py.File(cfg['lonfile'], 'r') as lonFileContent:
            lon = lonFileContent['LON'][iy, ix] / lonFileContent['LON'].attrs.get('SCALING_FACTOR')
        with h5py.File(cfg['latfile'], 'r') as latFileContent:
            lat = latFileContent['LAT'][iy, ix] / latFileContent['LAT'].attrs.get('SCALING_FACTOR')

    print(f'Coordinates at y/lat: {iy} and x/lon: {ix}')
    print(f'LON:  {lon}')
    print(f'LAT:  {lat}')

def getIndex(cfg, lonlatFileName, lonval, latval, gen_):
    """ get index of coordinates based on coordinates lonval and latval
    """
    if gen_ == '3':
        import xarray as xr

        if lonlatFileName is not None and os.path.isfile(lonlatFileName):
            fn = lonlatFileName
        else:
            fn = cfg['lonlatFileMtg']

        lon = xr.open_dataset(fn)['LON'].isel(time=0).values
        lat = xr.open_dataset(fn)['LAT'].isel(time=0).values

    elif gen_ == '2':
        import h5py

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
