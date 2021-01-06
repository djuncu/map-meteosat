#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 15:20:48 2020

@author: juncud
"""
import numpy as np
import cartopy.crs as ccrs
import os

# homemade libs
from . import my_utils_plot

def makeMsgPlot(val, lon, lat,
                 plot_xstep=None, plot_ystep=None, ixmin=None, ixmax=None, iymin=None, iymax=None, 
                 lonmin=None, lonmax=None, latmin=None, latmax=None,
                 vmin=None, vmax=None, norm=None, f_out_png='./dataarray.png',
                 title="",is_iodc=False, add_logo=False,
                 logo_path_MF=None, logo_path_SAF = None, 
                 figsize=None,cmap='viridis', tick_labels = None, suppressFig=False):
    """ Plot the content of a datarray on a geolocated grid & return the figure or save it to PNG file.
    Parameters: 
        da: Xarray DataArray containing the information to plot
        ixmin,ixmax,iymin,iymax: None or int. Pixel limits to keep for plot. if None, whole zone is used.
        vmin,vmax: limits of the colour palette. Points with values over vmax are plotted in white, points under vmin are in grey (see cmap.set_over/under).
        f_out_png: "None" or path to the output png file that will be created. Intermediate directories will be created if they don't exist. If None, the code returns the fig, ax and input dataArray.
        title: string. Used to define map title.
        is_iodc: Boolean. If True the geolocation and file naming are considered to be those of the MSG-Indian Ocean products.
        add_logo: Boolean : activate or not the addition of logos on the resulting plots
        logo_path_MF/SAF: Path to the logos to add on the figures if needed. These are used only if add_logo=True.    """  
    if suppressFig:
        import matplotlib
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    # define map projection: geostationary
    fig,ax = plt.subplots(figsize=figsize)
    central_lon = 0.0
    if is_iodc: central_lon = 41.5
    proj=ccrs.Geostationary(central_longitude=central_lon, satellite_height=35785831, false_easting=0, false_northing=0, globe=None)
    ax = plt.axes(projection=proj)    

    # plot dataset on the projection based on values of lon/lat (which are given in ccrs.PlateCarree() projection)
    if norm is None:
        p = ax.pcolormesh(lon, lat, val,transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, cmap = cmap)
    else:        
        p = ax.pcolormesh(lon, lat, val,transform=ccrs.PlateCarree(), norm=norm, cmap = cmap)
    
    if tick_labels is not None:        
        cb = fig.colorbar(p, ticks=np.arange(0.5,len(tick_labels)+0.5+1))
        cb.ax.set_yticklabels(tick_labels)
    else:
        cb = fig.colorbar(p)

    ax.coastlines()
    gl = ax.gridlines()
    gl.xlines = True
    gl.xlocator = mticker.FixedLocator([-80, -60, -40, -20,0, 20, 40, 60, 80])
    gl.ylines = True
    gl.ylocator = mticker.FixedLocator([-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75]) 
    if lonmin is not None and latmin is not None and lonmax is not None and latmax is not None:
        try:
            ax.set_extent([lonmin, lonmax, latmin, latmax], crs=ccrs.PlateCarree())
        except Exception as e:
            print(e)
            print("The boundaries provided shall be within the coverage of the geoprojected zone, please adapt parameters lonmin/max & latmin/max.")

    plt.axis('on')
    plt.title(title)
    plt.subplots_adjust(right=0.85)    
    if add_logo and logo_path_MF is not None and logo_path_SAF is not None :
        my_utils_plot.add_logo_to_figure(fig,ax,xy_position_tuple=(0.15,0.1),zoom=0.12, logo_path=logo_path_MF)
        my_utils_plot.add_logo_to_figure(fig,ax,xy_position_tuple=(0.60,0.1),zoom=0.45, logo_path=logo_path_SAF)

    if f_out_png == None and not suppressFig: # returns the figure and the filtered dataset
        return fig,ax
    else: # saves the map to a png file, no info is returned
        f_out_png = os.path.expanduser(f_out_png)
        my_utils_plot.ensure_dir(f_out_png) 
        plt.savefig(f_out_png)
        plt.close()
        print('Plotting finished. File saved to ' + f_out_png)

def makeMtgPlot(da,
                lonmin=None, lonmax=None, latmin=None, latmax=None,
                vmin=None, vmax=None, norm=None, f_out_png=None,
                title="",is_iodc=False, add_logo=False,
                logo_path_MF=None, logo_path_SAF = None, 
                figsize=None,cmap='viridis', cTickLabels = None, 
                mapLabels = False, suppressFig=False):

    if suppressFig:
        import matplotlib
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    # define map projection: geostationary
    fig,ax = plt.subplots(figsize=figsize)
    central_lon = 0.0
    satHeight = 35785831
    if is_iodc: central_lon = 41.5
    proj=ccrs.Geostationary(central_longitude=central_lon, satellite_height=satHeight, 
                            false_easting=0, false_northing=0, globe=None)
    projPc = ccrs.PlateCarree()
    ax = plt.axes(projection=proj)
    
    if norm is None:
        p = ax.pcolormesh(da.x.values*satHeight, da.y.values*satHeight,
                  da.values, transform = proj, vmin=vmin, vmax=vmax, cmap=cmap)
    else:
        p = ax.pcolormesh(da.x.values*satHeight, da.y.values*satHeight,
                  da.values, transform = proj, norm=norm, cmap=cmap)

    #if lonmin is not None or latmin is not None or lonmax is not None or latmax is not None:
    if mapLabels:
        cbPad = 0.15
    else:
        cbPad = 0.05
    
    if cTickLabels is not None:
        cb = fig.colorbar(p, pad=cbPad, ticks=np.arange(0.5,len(cTickLabels)+0.5+1)) # + 1?
        cb.ax.set_yticklabels(cTickLabels)
    else:
        cb = fig.colorbar(p, pad=cbPad)

    ax.coastlines()
    gl = ax.gridlines(crs=projPc, draw_labels=mapLabels)
    gl.xlines = True
    gl.ylines = True
    errMsg = 'The boundaries provided shall be within the coverage of the geoprojected zone, please adapt parameters lonmin/max & latmin/max.'
    if lonmin is not None and latmin is not None and lonmax is not None and latmax is not None:
        try:
            ax.set_extent([lonmin, lonmax, latmin, latmax], crs=projPc)
        except Exception as e:
            print(e)
            print(errMsg)
    elif lonmin is not None and lonmax is not None:
        try:
            ax.set_extent([lonmin, lonmax, -65, 65], crs=projPc)
        except Exception as e:
            print(e)
            print(errMsg)
    elif latmin is not None and latmax is not None:
        try:
            ax.set_extent([-65,65, latmin, latmax], crs=projPc)
        except Exception as e:
            print(e)
            print(errMsg)
    else:
        gl.xlocator = mticker.FixedLocator([-80, -60, -40, -20,0, 20, 40, 60, 80])
        gl.ylocator = mticker.FixedLocator([-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75])
        
    plt.axis('on')
    plt.title(title)
    plt.subplots_adjust(right=0.85)
    if add_logo and logo_path_MF is not None and logo_path_SAF is not None :
        my_utils_plot.add_logo_to_figure(fig,ax,xy_position_tuple=(0.15,0.1),zoom=0.12, logo_path=logo_path_MF)
        my_utils_plot.add_logo_to_figure(fig,ax,xy_position_tuple=(0.60,0.1),zoom=0.45, logo_path=logo_path_SAF)

    if not suppressFig: plt.show()

    if f_out_png == None and not suppressFig: # returns the figure and the filtered dataset
        return fig, ax
    elif f_out_png is not None: # saves the map to a png file, no info is returned
        f_out_png = os.path.expanduser(f_out_png)
        my_utils_plot.ensure_dir(f_out_png) 
        plt.savefig(f_out_png)
        plt.close()
        print('Plotting finished. File saved to ' + f_out_png)
