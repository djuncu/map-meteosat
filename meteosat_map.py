#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: juncud/vincentch
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import numpy as np
import sys
import copy
import yaml

# add the location of the current directory into the python path 
# so that the libs in ./src can be loaded even if the code is called from another location than its own directory
myLibDir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(myLibDir)

# homemade libs
from src import utils_xr_process
from src import my_utils_plot
from src import pv_readers_light as pv_readers
from src import rebuild_lonlat
from src.makeMsgPlot import makeMsgPlot

def plot_msg_geoloc(f_in_tplt="/cnrm/vegeo/SAT/DATA/MSG/Reprocessed-on-2017/MDAL/2004/01/19/HDF5_LSASAF_MSG_ALBEDO_MSG-Disk_200401190000", 
                        varname = 'AL-BB-BH', cfgFile = None, 
                        plot_xstep=10, plot_ystep=10, 
                        ixmin=None, ixmax=None, iymin=None, iymax=None, 
                        lonmin=None, lonmax=None, latmin=None, latmax=None,
                        qf_to_drop=[], err_max_percent_lim = 10,  
                        vmin='default', vmax='default', 
                        f_out_png='None',
                        msg_sce='default',is_iodc=False, 
                        lonfile='default', latfile='default',
                        read_version=False, add_logo=False,  logo_path_MF=None, logo_path_SAF = None, figsize= None,
                        cmap='default', 
                        color_under = 'default', color_over = 'default', color_bad = 'default',
                        convertUnits = False,
                        f_aux = None):
    """ Calls MSG reader for albedo, qflag & err. Apply scaling/offset, filter and remove missing data. Plot & save to PNG file the map of MSG albedo.
    Parameters: 
        f_in_tplt: glob compatible path pointing to the input file to read. (see read_msg in pv_readers). Its basename is also used to get the date/hour of the dataset.
        varname: variable name of the dataset to read, process and plot 
        varname_qf: variable name of the quality flag dataset associated to varame_alb. Used for filtering with qf_to_drop
        qf_to_drop: [] or list of numeric values. The pixels in varname_alb associated to these varname_qf values are removed.
        qf_bit_mask_ocean,missing: [] or list of int or of binary values. Used as a selection mask on da_qf bits (bitwise comparison, where results are compared with qf_bit_results). 
        qf_bit_results_ocean,missing: [] or list of int or of binary values. Correspond to the results expected in the bitwise comparison of da_qf & qf_bit_mask.
            Ocean pixels are i) removed from data, ii) used to plot a light grey background.
            If pixels are already flagged in the dataset, the qf_bit_mask/result_missing arguments are not needed. The missing data pixels will automatically be plotted as grey. 
            If not removed already from variable dataset, use the qf_bit_mask/result_missing arguments to remove these pixels and have them plotted as grey.
            All removed data from variable that are not flagged as ocean, will be plotted in a darker grey layer.
        err_max_percent_lim: None or percent value in [0 - 100]. Upper valid limit for varname_alb_err/varname_alb*100. When above the limit, the pixels in varname_alb are removed.
        plot_xstep, plot_ystep: None or int. Number of pixels to skip when plotting. If None, every pixel is kept for plot.
        ixmin,ixmax,iymin,iymax: None or int. Pixel limits to keep for plot. if None, whole zone is used.
        vmin,vmax: limits of the colour palette. Points with values over vmax are plotted in white, points under vmin are in grey (see cmap.set_over/under).
        f_out_png: "None" or path to the output png file that will be created. Intermediate directories will be created if they don't exist. If None, the code returns the fig, ax and dataArray containing the dataset being plotted.
        msg_sce: None or string. Used to define map title. If None, reads MSG sat number in file attributes & adds 'SEVIRI' suffix and a '/' separator as a prefix. If string, adds separator as prefix and uses full string in title. (if string is empty, adds nothing in title)
        is_iodc: Boolean. If True the geolocation and file naming are considered to be those of the MSG-Indian Ocean products.
        read_version: Boolean. Activates the reading of the algorithm version numeber from the file attributes and uses it in the plot title.
        lonfile, latfile: Path to the longitude and latitude files to pass to the MSG geolocation reader.
        add_logo: Boolean : activate or not the addition of logos on the resulting plots
        logo_path_MF/SAF: Path to the logos to add on the figures if needed. These are used only if add_logo=True.
        convertUnits: True or False. Converts to alternative units, only available for O3 
        others: see in read_msg. 
    """
    # check if input file exists
    if not os.path.isfile(f_in_tplt):
        raise ValueError('ERROR: Input file does not exist.')


    # load config
    if cfgFile is None:
        cfgFile = os.path.join(myLibDir, 'config.yml')
    cfg = yaml.safe_load(open(cfgFile))
    if varname[0:3] == 'AL-':
        vartype = 'Albedo'
    elif varname == 'AOD-CLIM':
        vartype = varname
        varname = 'AER'
    elif varname == 'AER':
        vartype = 'AOD-CLIM'
    elif varname == 'AOD-FC':
        vartype = varname
        varname = 'AOD550'
    elif varname == 'AOD550':
        vartype = 'AOD-FC'
    elif varname == 'BRF-TOC':
        vartype = 'TOC'
    elif varname == 'WV':
        vartype = varname
        varname = 'TCWV'
    elif varname == 'TCWV':
        vartype = 'WV'
    elif varname == 'O3-CLIM':
        vartype = varname
        varname = 'O3'
    elif varname == 'O3':
        vartype = 'O3-CLIM'
    elif varname == 'O3-FC':
        vartype = varname
        varname = 'TCO3'
    elif varname == 'TCO3':
        vartype = 'O3-FC'
    elif varname == 'TOA':
        vartype = varname
        varname ='RADIANCE'
    elif varname == 'TOC-Q':
        vartype = varname
        varname = 'BRF_Q_Flag'
    elif varname == 'ALB-Q':
        vartype = varname
        varname = 'Q-Flag'
    elif varname == 'CMa-Q':
        vartype = varname
        varname = 'CMa_QUALITY'
    else:
        vartype = varname

    # grab parameters
    scaling     = cfg['type'][vartype]['scaling']
    offset      = cfg['type'][vartype]['offset']
    valid_range = cfg['type'][vartype]['valid_range']

    varname_aux = cfg['type'][vartype]['qf_varname']
    qf_bit_mask_missing    = cfg['type'][vartype]['qf_missing']
    qf_bit_results_missing = cfg['type'][vartype]['qf_results_missing']
 
    qf_bit_mask_outside_disk    = cfg['qf_bit_mask']['outside_disk']
    qf_bit_results_outside_disk = cfg['qf_bit_mask']['results_outside_disk']
    qf_bit_mask_ocean           = cfg['qf_bit_mask']['ocean']
    qf_bit_results_ocean        = cfg['qf_bit_mask']['results_ocean']
    
    if isinstance(cfg['type'][vartype]['missing_data_val'], (int, float)):
        missing_data_val = []
        missing_data_val.append(cfg['type'][vartype]['missing_data_val'])
    elif isinstance(cfg['type'][vartype]['missing_data_val'], list):
        missing_data_val = cfg['type'][vartype]['missing_data_val']
    elif cfg['type'][vartype]['missing_data_val'] is None:
        missing_data_val = None
    else:
        raise ValueError('Wrong value for missing_data_val in yml file')
    
    if isinstance(cfg['type'][vartype]['missing_q_val'], (int, float)):
        missing_q_val = []
        missing_q_val.append(cfg['type'][vartype]['missing_q_val'])
    elif isinstance(cfg['type'][vartype]['missing_q_val'], list):
        missing_q_val = cfg['type'][vartype]['missing_q_val']
    elif cfg['type'][vartype]['missing_q_val'] is None:
        missing_q_val = None
    else:
        raise ValueError('Wrong value for missing_q_val in yml file')

    conversion_factor = cfg['type'][vartype]['conversion']['factor']

    if lonfile == 'default': lonfile = cfg['lonfile']
    if latfile == 'default': latfile = cfg['latfile']
    if vmin    == 'default': vmin = cfg['type'][vartype]['vmin']  
    if vmax    == 'default': vmax = cfg['type'][vartype]['vmax']  
    if cmap    == 'default': cmap = cfg['type'][vartype]['cmap'] 
    if color_under == 'default': color_under = cfg['type'][vartype]['color_under'] 
    if color_over  == 'default': color_over = cfg['type'][vartype]['color_over'] 
    if color_bad   == 'default': color_bad = cfg['type'][vartype]['color_bad'] 
    if msg_sce == 'default': msg_sce = cfg['type'][vartype]['msg_sce']

    # grab date from filename
    split_str = 'MSG-Disk_'
    if is_iodc: split_str = 'IODC-Disk_'
    date = os.path.basename(f_in_tplt).split(split_str)[1].split('_')[0][0:12]
    da = pv_readers.read_one_file_msg(f_in_path=f_in_tplt, varname =
                                                varname, scaling=scaling,
                                                offset=offset, valid_range =
                                                valid_range,
                                                missing=missing_data_val)

    if vartype == 'TOA':
        da_aux = pv_readers.read_one_file_msg(f_in_path=f_aux, varname = varname_aux, scaling=1/100., offset=0., valid_range = [0,90], missing=[-20000])

        # radiance to reflectance conversion
        channel = os.path.basename(f_in_tplt).split('HDF5_')[1].split('_')[0]
        da = da / (cfg['type'][vartype]['conversion']['bandfactor'][channel] * np.cos(da_aux *np.pi/180.) )

    if conversion_factor is not None:
        if convertUnits is True:
            da = da * conversion_factor
            units = '(' + cfg['type'][vartype]['conversion']['convertedUnit'] + ')'
        else:
            units = '(' + cfg['type'][vartype]['conversion']['defaultUnit'] + ')'


    # adding geoloc coordinates to dataset
    if not is_iodc:
        lonval, latval = rebuild_lonlat.read_msg_lonlats(fillvalue_to_nan=False,latfile=latfile,lonfile=lonfile)
    else:
        lonval, latval = rebuild_lonlat.read_msg_iodc_lonlats(fillvalue_to_nan=False,latfile=latfile,lonfile=lonfile)
 
    da_to_plot = da.isel(x=slice(ixmin,ixmax,plot_xstep),y=slice(iymin,iymax,plot_ystep))
    #da_to_plot = da_to_plot.fillna(-1)#.astype(np.float32)

    # reading processing version from the product attributes
    if read_version: version_id = str(ascii(da_to_plot.attrs['PRODUCT_ALGORITHM_VERSION'])[2:-1])
    else: version_id= ''

    ## Q-FLAG FILTERING
    # read & subsample qflag data if needed
    if vartype == 'Albedo' or vartype == 'Z_Age' or vartype == 'TOC':
        if qf_to_drop !=[] or qf_bit_mask_missing != [] or qf_bit_mask_outside_disk != [] : 
            da_qf = pv_readers.read_one_file_msg(f_in_path=f_in_tplt, varname =
                                                 varname_aux, scaling=1,
                                                 offset=0., valid_range = [],
                                                 missing = missing_q_val)

            da_qf = da_qf.isel(x=slice(ixmin,ixmax,plot_xstep),y=slice(iymin,iymax,plot_ystep)).astype(np.int16)    

        # filter on qflag if needed
        if qf_to_drop !=[]:
            da_to_plot = utils_xr_process.filter_on_qflag(da_to_plot, da_qf, qf_to_drop)

        if qf_bit_mask_missing != []:
            if len(qf_bit_mask_missing) == len(qf_bit_results_missing):
                da_to_plot = utils_xr_process.filter_bitwise_on_qflag(da_to_plot, da_qf, qf_bit_mask = qf_bit_mask_missing, 
                                                                      qf_bit_results = qf_bit_results_missing,
                                                                      fillvalue=np.NaN)
            else:
                print('Error: qf_bit_mask & qf_bit_results lists have different lengths')

        if qf_bit_mask_ocean != []:
            if len(qf_bit_mask_ocean) == len(qf_bit_results_ocean):
                da_to_plot = utils_xr_process.filter_bitwise_on_qflag(da_to_plot, da_qf, qf_bit_mask = qf_bit_mask_ocean, 
                                                                      qf_bit_results = qf_bit_results_ocean, fillvalue=-1)
            else:
                print('Error: qf_bit_mask & qf_bit_results lists have different lengths')

        # remove pixels outside mask so that the outside disk pixels remain with nan default colour (white)
        if qf_bit_mask_outside_disk != [] :
            if (len(qf_bit_mask_outside_disk) == len(qf_bit_results_outside_disk)): 
                da_to_plot = utils_xr_process.filter_bitwise_on_qflag(da_to_plot, da_qf, qf_bit_mask = qf_bit_mask_outside_disk, 
                                                                      qf_bit_results = qf_bit_results_outside_disk, fillvalue=np.NaN)
            else:
                print('Error: qf_bit_mask_outside_disk & qf_bit_results_outside_disk lists have different length')

        # select ocean pixels to mask on plot
        da_ocean = utils_xr_process.mask_bitwise_on_qflag(da_to_plot, da_qf, qf_bit_mask=qf_bit_mask_ocean, qf_bit_results=qf_bit_results_ocean)
    else:
        da_ocean = []

    ## read, subsample uncertainty, filter 
    if vartype == 'Albedo' and err_max_percent_lim != None:
        da_err = pv_readers.read_one_file_msg(f_in_path=f_in_tplt, varname = varname+'-ERR',scaling=1/10000, offset=0., valid_range = [0,1], missing=[-1])
        da_err = da_err.isel(x=slice(ixmin,ixmax,plot_xstep),y=slice(iymin,iymax,plot_ystep)).astype(np.float32) 
        da_to_plot = utils_xr_process.filter_on_err(da_to_plot, da_err, err_max_percent_lim, fillvalue=-1 )

    # Resample values. For Q-Flags, need to re-assign values first
    if vartype == 'ALB-Q':
        # set bitwise comparison factors : "vmask" defines which bits are considered in the comparison with da_qf, "value" defines the matching bits in the comparison (1 is same, 0 is different)
        # see python bitwise operations for more info
        vmask  = [0b11, 0b11111, 0b10100111, 0b10100111, 0b10100100, 0b10000000, 0b11] 
        value  = [0b00, 0b00001, 0b00000101, 0b00000111, 0b00100100, 0b10000000, 0b10] 
        qVals = da_to_plot.copy().values # i don't know why but it is important to copy()
        for i,val in enumerate(value):
            # we re-number a new variable with contiguous values in [0, len(value)-1] for each qflag category
            mask = np.where(da_to_plot.astype(int) & vmask[i] == val)
            qVals[mask] = i

        lon_resamp, lat_resamp, da_resamp = my_utils_plot.resampleMesh(lonval, latval,
                                                         qVals, plot_xstep, plot_ystep, da_ocean)
    elif vartype == 'TOC-Q':
        # set bitwise comparison factors : "vmask" defines which bits are considered in the comparison with da_qf, "value" defines the matching bits in the comparison (1 is same, 0 is different)
        # see python bitwise operations for more info
        vmask  = [0b11, 0b11, 0b10011111, 0b10011111, 0b10011101, 0b10011101, 0b10011101, 0b10011101, 0b11101, 0b10011001]#, 0b10000000] # equiv to pvwave [3, 3, 159, 159, 157, 157, 157, 157,  29, 157, 157, 128]
        value  = [0b00, 0b10, 0b10000101, 0b10000111, 0b10010001, 0b00000001, 0b00001101, 0b00001001, 0b10101, 0b10011001]#, 0b10000000] # equiv to pvwave [0, 2, 133, 135, 145,   1,  13,   9,  21, 153, 157, 128]
        qVals = da_to_plot.copy().values # get the values in a numpy array. it needs the copy(), no idea why
        for i,val in enumerate(value):
            # we re-number a new variable with contiguous values in [0, len(value)-1] for each qflag category
            mask = np.where(da_to_plot.astype(int) & vmask[i] == val)
            qVals[mask] = i

        lon_resamp, lat_resamp, da_resamp = my_utils_plot.resampleMesh(lonval, latval,
                                                         qVals, plot_xstep, plot_ystep, da_ocean)
    elif vartype == 'CMa':
        qVals = da_to_plot.copy().values # get the values in a numpy array. it needs the copy(), no idea why
   #     for i,val in enumerate(value):
   #         # we re-number a new variable with contiguous values in [0, len(value)-1] for each qflag category
   #         mask = np.where(da_to_plot.astype(int) & vmask[i] == val)
   #         qVals[mask] = i
        qVals[qVals == 0] = 5
        qVals = qVals - 1
        lon_resamp, lat_resamp, da_resamp = my_utils_plot.resampleMesh(lonval, latval,
                                                         qVals, plot_xstep, plot_ystep, da_ocean)
    elif vartype == 'CMa-Q':
        #value  = [0b000000000,0b010000000,0b100000000,0b110000000]# equiv to pvwave [  0, 128, 256, 384]
        #vmask  = [0b110000000,0b110000000,0b110000000,0b110000000]# equiv to pvwave [384, 384, 384, 384]
        
        # not 100 percent sure if these are the right ones
        print('WARNING: This is still in testing, I do not think these results
              are correct at this time.')
        # need to figure out what bits contain the quality information here...
        value  = [0b000000000000, 0b001000000000, 0b010000000000,
                  0b011000000000,
                  0b100000000000]
        vmask  = [0b111000000000, 0b111000000000, 0b111000000000,
                  0b111000000000,
                  0b111000000000]
        qVals = da_to_plot.copy().values # get the values in a numpy array. it needs the copy(), no idea why
        for i,val in enumerate(value):
            # we re-number a new variable with contiguous values in [0, len(value)-1] for each qflag category
            mask = np.where(da_to_plot.astype(int) & vmask[i] == val)
            qVals[mask] = i
        lon_resamp, lat_resamp, da_resamp = my_utils_plot.resampleMesh(lonval, latval,
                                                         qVals, plot_xstep, plot_ystep, da_ocean)

    else:
        lon_resamp, lat_resamp, da_resamp = my_utils_plot.resampleMesh(lonval, latval,
                                                         da_to_plot, plot_xstep, plot_ystep, da_ocean)
    ## Plotting
    # channel dictionary
    dic_channels_wl = {'006':'0.6µm','008':'0.8µm','016':'1.6µm'}

    if msg_sce == 'default':
        if vartype == 'Albedo' or vartype == 'Z_Age':
            sce = ' / '+da.attrs['SATELLITE'][0].astype(str)+'-SEVIRI'
        else:
            sce = ''
    else:
        if msg_sce == '': sce = ''
        if msg_sce != '': sce = msg_sce + version_id

    if vartype =='DEM':
        titleString = f'DEM-USGS (m)'
        clrmap = copy.copy(plt.get_cmap(cmap))
        colors_undersea = clrmap(np.linspace(0, 0.17, 256))
        colors_land = clrmap(np.linspace(0.25, 1, 256))
        all_colors = np.vstack((colors_undersea, colors_land))
        clrmap = mcolors.LinearSegmentedColormap.from_list('terrain_map', all_colors)
        vcenter = 0
        plotnorm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
        tick_labels = None
    elif vartype == 'Z_Age':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]
        titleString = f' Age of Last Observation [days] {datestr}\n {sce}' 
        # plot & customisation
        # for age variable : custom colormap with emphasis on age=0 (defined based on existing 'Purples_r')
        cmap = copy.copy(plt.get_cmap("BuPu_r",15))
        cmap2 = copy.copy(plt.get_cmap("BuPu_r",20))
        cmaplist= [cmap2(i) for i in range(5,5+cmap.N-1)] # we artificially make a gap in the colours decrease to highlight the 0-day age.
        cmaplist[0] = cmap(0) # dark purple
        clrmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
        if color_under is not None: clrmap.set_under(color_under) # plot data < vmin in grey.
        if color_over is not None: clrmap.set_over(color_over) # plot data > vmax in white.
        if color_bad is not None:  clrmap.set_bad(color_bad)
        plotnorm = None
        tick_labels = None
    elif vartype == 'TOC-Q':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]+' '+date[8:10]+':'+date[10:12]+' UTC'
        channel = os.path.basename(f_in_tplt).split('-BRF')[0].split('_')[-1]
        titleString = f'TOC Reflectance Processing Flags at {dic_channels_wl[channel]}\n {datestr} {sce}'
        if '{channel}' in f_out_png: f_out_png = f_out_png.replace('{channel}',channel)
        fname  = ['Ocean', 'Space', 'Land Clear', 'Water Clear', 'Snow Clear','No CMa/Night', 'Cloud filled', 'Cloud Contam.', 'Shadow', 'Bad QA']#, 'Processed']
        colour = ["b","w","chartreuse","skyblue","silver","k","yellow","orange","grey","r"]#,"m"]
        bname  = ['...00', '...10', '1..00101', '1..00111', '1..100.1', '0..000.1', '0..011.1', '0..010.1', '...101.1', '1..11..1']#, '1.......']
        tick_labels = ['\n'.join([x,bname[i]]) for i,x in enumerate(fname) ] # merge fname & bname to legend each class

        n_colors = len(value)
        clrmap = mcolors.ListedColormap(colour)
        plotnorm = mcolors.BoundaryNorm(np.arange(0,n_colors+1), clrmap.N)
    elif vartype == 'ALB-Q':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]
        titleString = f'Albedo Processing Flags \n {datestr} {sce}'
        fname  = ['Ocean','Missing Obs.','Land OK','Water OK','Snow OK','Algo Failed','Space']
        bname  = [ '......00', '...000.1', '0.00.101', '0.00.111', '0.1..1..', '1.......', '......10']
        colour = ["b", "yellow", "chartreuse", "skyblue", "silver", "r","w"]
        tick_labels = ['\n'.join([x,bname[i]]) for i,x in enumerate(fname) ] # merge fname & bname to legend each class
 
        n_colors = len(value)
        clrmap = mcolors.ListedColormap(colour)
        plotnorm = mcolors.BoundaryNorm(np.arange(0,n_colors+1), clrmap.N)
    elif vartype == 'CMa':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]+' '+date[8:10]+':'+date[10:12]+' UTC'
        titleString = f'Cloud Mask (NWC SAF, v2018)\n {datestr}'
        fname  = ['Cloud-Free', 'Cloudy', 'Cloud Contaminated', 'Snow/Ice',
                  'No data / undefined']
        colour = ["dodgerblue","dimgrey","lightgrey","black","w"]
        bname  = ['', '', '', '', '']
        tick_labels = ['\n'.join([x,bname[i]]) for i,x in enumerate(fname) ] # merge fname & bname to legend each class

        n_colors = 5
        clrmap = mcolors.ListedColormap(colour)
        plotnorm = mcolors.BoundaryNorm(np.arange(0,n_colors+1), clrmap.N)
    elif vartype == 'CMa-Q':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]+' '+date[8:10]+':'+date[10:12]+' UTC'
        titleString = f'CMA Processing Flag (NWC SAF, v2012)\n {datestr}'
        fname  = ['No data', 'Good Quality', 'Questionable', 'Bad',
                  'Interpolated / reclassified']
        colour = ["w","chartreuse","orange","red","magenta"]
        bname  = ['...000', '...001', '...010', '...011', '...100']
        tick_labels = ['\n'.join([x,bname[i]]) for i,x in enumerate(fname) ] # merge fname & bname to legend each class

        n_colors = len(value)
        clrmap = mcolors.ListedColormap(colour)
        plotnorm = mcolors.BoundaryNorm(np.arange(0,n_colors+1), clrmap.N)

    else:
        clrmap = copy.copy(plt.get_cmap(cmap))
        if color_under is not None: clrmap.set_under(color_under) # plot data < vmin in grey.
        if color_over is not None: clrmap.set_over(color_over) # plot data > vmax in white.
        if color_bad is not None: clrmap.set_bad(color_bad)
        plotnorm = None
        tick_labels = None

    if vartype == 'Albedo':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]
        titleString = varname[3:]+f' - Albedo {datestr}\n {sce}'
    elif vartype  == 'AOD-CLIM':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]
        titleString = f'AOD550 {datestr}\n {sce}'
    elif vartype  == 'AOD-FC':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]+' '+date[8:10]+':'+date[10:12]+' UTC'
        titleString = f'AOD550 {datestr}\n {sce}'
    elif vartype == 'TOC':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]+' '+date[8:10]+':'+date[10:12]+' UTC'
        channel = os.path.basename(f_in_tplt).split('-BRF')[0].split('_')[-1]
        titleString = f'TOC Reflectance Factor at {dic_channels_wl[channel]}\n {datestr} {sce}'
        if '{channel}' in f_out_png: f_out_png = f_out_png.replace('{channel}',channel)
    elif vartype == 'WV':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]+' '+date[8:10]+':'+date[10:12]+' UTC'
        titleString = f'Water Vapour Total Column Content (kg/m2)\n{datestr} {sce}'
    elif vartype == 'MSL':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]+' '+date[8:10]+':'+date[10:12]+' UTC'
        titleString = f'Pressure at Mean Sea Level (Pa)\n{datestr} {sce}'
    elif vartype == 'O3-CLIM':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]
        titleString = f'TOMS O3 clim {units}\n{datestr} {sce}'
    elif vartype == 'O3-FC':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]+' '+date[8:10]+':'+date[10:12]+' UTC'
        titleString = f'ECMWF O3 forecast {units}\n{datestr} {sce}'
    elif vartype == 'TOA':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]+' '+date[8:10]+':'+date[10:12]+' UTC'
        titleString = f'TOA Reflectance at {dic_channels_wl[channel]}\n {datestr} {sce}'
        if '{channel}' in f_out_png: f_out_png = f_out_png.replace('{channel}',channel)

    makeMsgPlot(da_resamp, lon_resamp, lat_resamp, 
                plot_xstep=plot_xstep, plot_ystep=plot_ystep, ixmin=ixmin, ixmax=ixmax, iymin=iymin, iymax=iymax, 
                lonmin=lonmin, lonmax=lonmax, latmin=latmin, latmax=latmax,
                vmin=vmin, vmax=vmax, norm=plotnorm, f_out_png=f_out_png,
                title=titleString, is_iodc=is_iodc, add_logo=add_logo,
                logo_path_MF=logo_path_MF, logo_path_SAF = logo_path_SAF, 
                figsize=figsize, cmap=clrmap, tick_labels=tick_labels)

