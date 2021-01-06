#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: juncud/vincentch
"""

import argparse
import sys
import os

def main():
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Plot Meteosat datasets.', add_help=False)
    requiredArgs = parser.add_argument_group('required arguments')
    optionalArgs = parser.add_argument_group('optional arguments')
    plotExtent = parser.add_argument_group('(optional) plot extent')
    plotLayout = parser.add_argument_group('(optional) plot title & layout')
    plotColors = parser.add_argument_group('(optional) color scale')
    dataArgs   = parser.add_argument_group('(optional) data modifiers')
    fileNames  = parser.add_argument_group('(optional) files')
    parser.add_argument('file', type=str, help='Input file name')
    requiredArgs.add_argument('-V','--var', help='variable to plot', required='--list' not in sys.argv and '-l' not in sys.argv)
    optionalArgs.add_argument('-h', '--help', action='help', help='show this help message and exit')
    optionalArgs.add_argument('-l','--list', help='list variables in File and exit', action='store_true')
    optionalArgs.add_argument('-G','--gen', help='Meteosat Generation number (default="3")', default='3', choices=['2','3'])
    optionalArgs.add_argument('-S','--suppressfig', help='suppress opening a figure', action='store_true')
    optionalArgs.add_argument('-s','--stride', help='plot only every nth pixel (default=10)', default = 10, type=int)
    fileNames.add_argument('-c', '--config', help='name of config file')
    plotLayout.add_argument('--figsize', help='figure width, height in inches', nargs = 2, type=float, metavar = ('WIDTH','HEIGHT'))
    plotExtent.add_argument('--ixmin', help='minimum plot x-extent py pixel', required='--ixmax' in sys.argv)
    plotExtent.add_argument('--ixmax', help='maximum plot x-extent py pixel', required='--ixmin' in sys.argv)
    plotExtent.add_argument('--iymin', help='mynimum plot y-extent py pixel', required='--iymax' in sys.argv)
    plotExtent.add_argument('--iymax', help='maximum plot y-extent py pixel', required='--iymin' in sys.argv)
    plotExtent.add_argument('--lonmin', help='minimum plot lon. extent', required='--lonmax' in sys.argv)
    plotExtent.add_argument('--lonmax', help='maximum plot lon. extent', required='--lonmin' in sys.argv)
    plotExtent.add_argument('--latmin', help='minimum plot lat. extent', required='--latmax' in sys.argv)
    plotExtent.add_argument('--latmax', help='maximum plot lat. extent', required='--latmin' in sys.argv)
    dataArgs.add_argument('--maxerr', help='None or percent value in [0 - 100]. Upper valid limit for ' + 
                        'varname_alb_err/varname_alb*100. When above the limit, the pixels in varname_alb are removed.',
                        type=int)
    plotColors.add_argument('--vmin', help='minimum value for color scale', type=float)
    plotColors.add_argument('--vmax', help='maximum value for color scale', type=float)
    fileNames.add_argument('-o','--outfile', help='name of output png file')
    plotLayout.add_argument('--msg_sce', help='None or string. Used to define map title. If None, '+
                        'reads MSG sat number in file attributes & adds "SEVIRI" suffix and a '+
                        '"/" separator as a prefix. If string, adds separator as prefix and'+
                        'uses full string in title. (if string is empty, adds nothing in title)',
                        default='default')
    optionalArgs.add_argument('--is_iodc', help='if raised the geolocation and file naming are '+
                        'considered to be those of the MSG-Indian Ocean products.', 
                        action='store_true')
    plotLayout.add_argument('--read_version', help='If raised activates the reading of the algorithm '+
                        'version number from the file attributes and uses it in the plot title.',
                        action='store_true')
    plotLayout.add_argument('--add_logo', help='activates the addition of logos on the resulting plots',
                        action='store_true')
    plotColors.add_argument('--cmap', help='colormap for plot', default='default', metavar='COLORMAP')
    plotColors.add_argument('--color_under', help='color for values under vmin', default='default', metavar='COLOR')
    plotColors.add_argument('--color_over', help='color for values over vmax', default='default', metavar='COLOR')
    plotColors.add_argument('--color_bad', help='color for bad values', default='default', metavar='COLOR')
    plotLayout.add_argument('--maplabels', help='activate coordinate labels on map', action='store_true')
    dataArgs.add_argument('--convert_units', help='converts to alternative units, only available for O3',
                        action='store_true')
    fileNames.add_argument('--auxfile', help='auxiliary file for TOA, needed for conversion from radiance to reflectance')
    dataArgs.add_argument('--qf_to_drop', help='[] or list of numeric values. The pixels in varname_alb '+
                        'associated to these varname_qf values are removed.', nargs='+', default=[])
    args = vars(parser.parse_args())

    # check if input file exists
    if not os.path.isfile(args['file']):
        raise ValueError('ERROR: Input file does not exist.')

    if args['list'] == True:
        from netCDF4 import Dataset

        fileContent = Dataset(args['file'])
        print('File contains the following variables:')
        for var in fileContent.variables:
            if var != 'x' and var != 'y' and var != 'time' and var != 'geostationary': print(var)
        sys.exit()

    
    if args['vmin'] == None: args['vmin'] = 'min'
    if args['vmax'] == None: args['vmax'] = 'max'

    plot_msg_geoloc(f_in_tplt=args['file'], varname=args['var'], cfgFile = args['config'], 
                        plot_xstep=args['stride'], plot_ystep=args['stride'], 
                        ixmin=args['ixmin'], ixmax=args['ixmax'], iymin=args['iymin'], iymax=args['iymax'], 
                        lonmin=args['lonmin'], lonmax=args['lonmax'], latmin=args['latmin'], latmax=args['latmax'],
                        err_max_percent_lim = args['maxerr'],
                        vmin=args['vmin'], vmax=args['vmax'], 
                        f_out_png=args['outfile'], qf_to_drop = args['qf_to_drop'],
                        msg_sce=args['msg_sce'],is_iodc=args['is_iodc'], 
                        read_version=args['read_version'], add_logo=args['add_logo'], figsize= args['figsize'],
                        cmap=args['cmap'], maplabels=args['maplabels'],
                        color_under=args['color_under'], color_over=args['color_over'], color_bad=args['color_bad'],
                        convertUnits = args['convert_units'],
                        f_aux = args['auxfile'], meteosatGen = args['gen'], suppressFig = args['suppressfig'])


def plot_msg_geoloc(f_in_tplt,
                        varname, cfgFile = None, 
                        plot_xstep=10, plot_ystep=10, 
                        ixmin=None, ixmax=None, iymin=None, iymax=None, 
                        lonmin=None, lonmax=None, latmin=None, latmax=None,
                        err_max_percent_lim = 10,  
                        vmin='min', vmax='max', 
                        f_out_png=None, qf_to_drop = [],
                        msg_sce='default',is_iodc=False, 
                        read_version=False, add_logo=False, figsize= None,
                        cmap='default', maplabels=False, 
                        color_under = 'default', color_over = 'default', color_bad = 'default',
                        convertUnits = False,
                        f_aux = None, meteosatGen = '2', suppressFig = False):
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
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import numpy as np
    import copy
    import yaml

    # add the location of the current directory into the python path 
    # so that the libs in ./src can be loaded even if the code is called from another location than its own directory
    myLibDir = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, myLibDir)

    # homemade libs
    from mxgplotlib import utils_xr_process
    from mxgplotlib import my_utils_plot
    from mxgplotlib import rebuild_lonlat
    from mxgplotlib.makeMxgPlot import makeMsgPlot, makeMtgPlot
    
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
    elif varname == 'BRF-TOC' and meteosatGen == '2':
        vartype = 'TOC-MSG'
    elif varname == 'BRF-TOC' and meteosatGen == '3':
        vartype = 'TOC-MTG'
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
        varname = 'RADIANCE'
    elif varname == 'TOC-Q':
        vartype = varname
        varname = 'BRF_Q_Flag'
    elif varname == 'ALB-Q':
        vartype = varname
        varname = 'Q-Flag'
    elif varname == 'CMa_QUALITY':
        vartype = 'CMa-Q'
    elif varname.lower() == 'cma':
        vartype = 'CMa'
    elif varname == 'cma_quality':
        vartype = 'CMa-Q'
    elif varname[0:3] == 'RAD':
        vartype = 'RAD'
    else:
        vartype = varname

    # grab parameters
    if meteosatGen == '2':
        lonfile = cfg['lonfile']
        latfile = cfg['latfile']
    
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
    
    if isinstance(cfg['type'][vartype]['qlike'], bool):
        qlike = cfg['type'][vartype]['qlike']

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

    if qlike and vmin == 'min': vmin = cfg['type'][vartype]['vmin']
    if qlike and vmin == 'max': vmax = cfg['type'][vartype]['vmax']


    if vmin    == 'default': vmin = cfg['type'][vartype]['vmin']  
    if vmax    == 'default': vmax = cfg['type'][vartype]['vmax']  
    if cmap    == 'default': cmap = cfg['type'][vartype]['cmap'] 
    if color_under == 'default': color_under = cfg['type'][vartype]['color_under'] 
    if color_over  == 'default': color_over = cfg['type'][vartype]['color_over'] 
    if color_bad   == 'default': color_bad = cfg['type'][vartype]['color_bad'] 
    if msg_sce == 'default': msg_sce = cfg['type'][vartype]['msg_sce']

    logo_path_MF  = cfg['logo_path_MF']
    logo_path_SAF = cfg['logo_path_SAF']

    # grab date from filename
    split_str = 'MSG-Disk_'
    if is_iodc: split_str = 'IODC-Disk_'
    #baseFileName = os.path
    date = os.path.basename(f_in_tplt).split('.')[0][-12::] 
    #if meteosatGen == '2':
    #    date = os.path.basename(f_in_tplt).split(split_str)[1].split('_')[0][0:12]
    #elif meteosatGen == '3':
    #    date = os.path.basename
 
    da = utils_xr_process.read_one_file_meteosat(f_in_path=f_in_tplt, varname =
                                                 varname, scaling=scaling,
                                                 offset=offset, valid_range =
                                                 valid_range,
                                                 missing=missing_data_val, meteosatGen=meteosatGen)

    if vartype == 'TOA':
        if f_aux is None:
            raise ValueError('no auxfile given')

        da_aux = utils_xr_process.read_one_file_meteosat(f_in_path=f_aux, varname = varname_aux,
                                                         scaling=1/100., offset=0., valid_range = [0,90], 
                                                         missing=[-20000], meteosatGen=meteosatGen)

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
    if meteosatGen == '2':
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
    if vartype == 'Albedo' or vartype == 'Z_Age' or vartype == 'TOC-MSG' or vartype == 'TOC-MTG':
        if qf_to_drop !=[] or qf_bit_mask_missing != [] or qf_bit_mask_outside_disk != [] : 
            da_qf = utils_xr_process.read_one_file_meteosat(f_in_path=f_in_tplt, varname =
                                                 varname_aux, scaling=1,
                                                 offset=0., valid_range = [],
                                                 missing = missing_q_val, meteosatGen=meteosatGen)

            da_qf = da_qf.isel(x=slice(ixmin,ixmax,plot_xstep),y=slice(iymin,iymax,plot_ystep)).astype(np.int16) 

        # filter on qflag if needed
        if qf_to_drop != []:
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
    if vartype == 'Albedo' and err_max_percent_lim is not None:
        da_err = utils_xr_process.read_one_file_meteosat(f_in_path=f_in_tplt, 
                                                         varname = varname+'-ERR',
                                                         scaling=1/10000, offset=0., 
                                                         valid_range = [0,1], missing=[-1],
                                                         meteosatGen=meteosatGen)
        da_err = da_err.isel(x=slice(ixmin,ixmax,plot_xstep),
                             y=slice(iymin,iymax,plot_ystep)).astype(np.float32) 
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

        if meteosatGen == '2':
            lon_resamp, lat_resamp, da_resamp = my_utils_plot.resampleMesh(lonval, latval,
                                                         qVals, plot_xstep, plot_ystep, da_ocean)
        elif meteosatGen == '3':
            da_to_plot.values = qVals

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

        if meteosatGen == '2':
            lon_resamp, lat_resamp, da_resamp = my_utils_plot.resampleMesh(lonval, latval,
                                                         qVals, plot_xstep, plot_ystep, da_ocean)
        elif meteosatGen == '3':
            da_to_plot.values = qVals
    
    elif vartype == 'CMa' and meteosatGen == '2':
        qVals = da_to_plot.copy().values # get the values in a numpy array. it needs the copy(), no idea why
   #     for i,val in enumerate(value):
   #         # we re-number a new variable with contiguous values in [0, len(value)-1] for each qflag category
   #         mask = np.where(da_to_plot.astype(int) & vmask[i] == val)
   #         qVals[mask] = i
        qVals[qVals == 0] = 5
        qVals = qVals - 1
        
        lon_resamp, lat_resamp, da_resamp = my_utils_plot.resampleMesh(lonval, latval, qVals, 
                                                                       plot_xstep, plot_ystep, da_ocean)
    elif vartype == 'CMa-Q' and meteosatGen == '2':
        #value  = [0b000000000,0b010000000,0b100000000,0b110000000]# equiv to pvwave [  0, 128, 256, 384]
        #vmask  = [0b110000000,0b110000000,0b110000000,0b110000000]# equiv to pvwave [384, 384, 384, 384]
        
        # not 100 percent sure if these are the right ones
        print('WARNING: This is still in testing, I do not think these results are correct at this time.')
        # need to figure out what bits contain the quality information here...
        #value  = [0b000000000000, 0b001000000000, 0b010000000000,
        #          0b011000000000,
        #          0b100000000000]
        #vmask  = [0b111000000000, 0b111000000000, 0b111000000000,
        #          0b111000000000,
        #          0b111000000000]
        value  = [0b000000000001, 0b000000000010, 0b000000000100,
                  0b000000001000, 0b000000010000, 0b000000011000, 0b000000100000] # [1, 2, 4, 8, 16, 24, 32]
        vmask  = [0b000000000001, 0b000000000010, 0b000000000100,
                  0b000000111000, 0b000000111000, 0b000000111000, 0b000000111000] 
        qVals = da_to_plot.copy().values # get the values in a numpy array. it needs the copy(), no idea why
        for i,val in enumerate(value):
            # we re-number a new variable with contiguous values in [0, len(value)-1] for each qflag category
            mask = np.where(da_to_plot.astype(int) & vmask[i] == val)
            qVals[mask] = i
        
        
        lon_resamp, lat_resamp, da_resamp = my_utils_plot.resampleMesh(lonval, latval,
                                                         qVals, plot_xstep, plot_ystep, da_ocean)

    elif vartype == 'CMa-Q' and meteosatGen == '3':
        
        value  = [0b000000000001, 0b000000000010, 0b000000000100,
                  0b000000001000, 0b000000010000, 0b000000011000, 0b000000100000] # [1, 2, 4, 8, 16, 24, 32]
        vmask  = [0b000000000001, 0b000000000010, 0b000000000100,
                  0b000000111000, 0b000000111000, 0b000000111000, 0b000000111000] 
        qVals = da_to_plot.copy().values # get the values in a numpy array. it needs the copy(), no idea why
        for i,val in enumerate(value):
            # we re-number a new variable with contiguous values in [0, len(value)-1] for each qflag category
            mask = np.where(da_to_plot.astype(int) & vmask[i] == val)
            qVals[mask] = i
        
        da_to_plot.values = qVals
    
    elif meteosatGen == '2':
        lon_resamp, lat_resamp, da_resamp = my_utils_plot.resampleMesh(lonval, latval,
                                                         da_to_plot, plot_xstep, plot_ystep, da_ocean)
    ## Plotting
    # channel dictionary
    dic_channels_wl = {#'006':'0.6µm','008':'0.8µm','016':'1.6µm', 
                       **dict.fromkeys(['006','VIS06'], '0.6µm'),
                       **dict.fromkeys(['008','VIS08'], '0.8µm'),
                       **dict.fromkeys(['016','NIR16'], '1.6µm')}


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
        if meteosatGen == '2':
            channel = os.path.basename(f_in_tplt).split('-BRF')[0].split('_')[-1]
        if meteosatGen == '3':
            channel = os.path.basename(f_in_tplt).split('_')[-2]
        titleString = f'TOC Reflectance Processing Flags at {dic_channels_wl[channel]}\n {datestr} {sce}'
        if '{channel}' in f_out_png: f_out_png = f_out_png.replace('{channel}',channel)
        fname  = ['Ocean', 'Space', 'Land Clear', 'Water Clear', 'Snow Clear','No CMa/Night', 'Cloud filled', 'Cloud Contam.', 'Shadow', 'Bad QA']#, 'Processed']
        colour = ["xkcd:faded blue","w","xkcd:green","skyblue","silver","k","yellow","orange","grey","r"]#,"m"]
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
    elif vartype == 'CMa' and meteosatGen == '2':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]+' '+date[8:10]+':'+date[10:12]+' UTC'
        titleString = f'Cloud Mask (NWC SAF, v2018)\n {datestr}'
        fname  = ['Cloud-Free', 'Cloudy', 'Cloud Contaminated', 'Snow/Ice',
                  'No data / undefined']
        colour = ["dimgrey","whitesmoke","lightgrey","dodgerblue","red"]
        bname  = ['', '', '', '', '']
        tick_labels = ['\n'.join([x,bname[i]]) for i,x in enumerate(fname) ] # merge fname & bname to legend each class
        n_colors = 5
        clrmap = mcolors.ListedColormap(colour)
        plotnorm = mcolors.BoundaryNorm(np.arange(0,n_colors+1), clrmap.N)
    elif vartype == 'CMa' and meteosatGen == '3':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]+' '+date[8:10]+':'+date[10:12]+' UTC'
        titleString = f'Cloud Mask (NWC SAF, v2018)\n {datestr}'
        fname  = ['Cloud-Free', 'Cloudy']
        colour = ["dimgrey","whitesmoke"]
        bname  = ['', '']
        tick_labels = ['\n'.join([x,bname[i]]) for i,x in enumerate(fname) ] # merge fname & bname to legend each class
        n_colors = len(colour)
        clrmap = mcolors.ListedColormap(colour)
        plotnorm = mcolors.BoundaryNorm(np.arange(0,n_colors+1), clrmap.N)
    elif vartype == 'CMa-Q':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]+' '+date[8:10]+':'+date[10:12]+' UTC'
        titleString = f'CMA Processing Flag (NWC SAF, v2012)\n {datestr}'
        if meteosatGen == '2':
            fname  = ['No data', 'Good Quality', 'Questionable', 'Bad',
                  'Interpolated / reclassified']
            colour = ["w","chartreuse","orange","red","magenta"]
            bname  = ['...000', '...001', '...010', '...011', '...100']
        elif meteosatGen == '3':
            fname  = ['No data', 'Internal Consistency', 'Temporal Consistency', 'Good', 'Questionable', 'Bad', 'Interpolated']
            colour = ['w', 'springgreen' , 'palegreen', 'chartreuse', 'orange', 'firebrick', 'magenta']
            bname  = ['', '', '', '', '','','']
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
    elif vartype == 'TOC-MSG' or vartype == 'TOC-MTG':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]+' '+date[8:10]+':'+date[10:12]+' UTC'
        if meteosatGen == '2':
            channel = os.path.basename(f_in_tplt).split('-BRF')[0].split('_')[-1]
        elif meteosatGen == '3':
            channel = os.path.basename(f_in_tplt).split('_')[-2]
        titleString = f'TOC Reflectance Factor at {dic_channels_wl[channel]}\n {datestr} {sce}'
        if f_out_png is not None and '{channel}' in f_out_png:
            f_out_png = f_out_png.replace('{channel}',channel)
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
        if f_out_png is not None and '{channel}' in f_out_png:
            f_out_png = f_out_png.replace('{channel}',channel)
    elif vartype == 'RAD':
        datestr = date[0:4]+'-'+date[4:6]+'-'+date[6:8]+' '+date[8:10]+':'+date[10:12]+' UTC'
        if meteosatGen == '2':
            channel = os.path.basename(f_in_tplt).split('HDF5_')[1].split('_')[0]
        elif meteosatGen == '3':
            channel = os.path.basename(f_in_tplt).split('RAD-')[1].split('_')[0]

            titleString = f'Radiance at {dic_channels_wl[channel]} \n {datestr} {sce}'
            units = 'mW m$^{-2}$ sr$^{-1}$ cm$^{-1}$'
        if f_out_png is not None and '{channel}' in f_out_png:
            f_out_png = f_out_png.replace('{channel}',channel)


    if meteosatGen == '2':
        if vmin == 'min': vmin = np.min(da_resamp)
        if vmax == 'max': vmax = np.max(da_resamp)
        makeMsgPlot(da_resamp, lon_resamp, lat_resamp, 
                plot_xstep=plot_xstep, plot_ystep=plot_ystep, 
                ixmin=ixmin, ixmax=ixmax, iymin=iymin, iymax=iymax, 
                lonmin=lonmin, lonmax=lonmax, latmin=latmin, latmax=latmax,
                vmin=vmin, vmax=vmax, norm=plotnorm, f_out_png=f_out_png,
                title=titleString, is_iodc=is_iodc, add_logo=add_logo,
                logo_path_MF=logo_path_MF, logo_path_SAF = logo_path_SAF, 
                figsize=figsize, cmap=clrmap, tick_labels=tick_labels, suppressFig=suppressFig)
    elif meteosatGen == '3':
        if vmin == 'min': vmin = np.min(da_to_plot)
        if vmax == 'max': vmax = np.max(da_to_plot)
        makeMtgPlot(da_to_plot,
                 lonmin=lonmin, lonmax=lonmax, latmin=latmin, latmax=latmax,
                 vmin=vmin, vmax=vmax, norm=plotnorm, f_out_png=f_out_png,
                 title=titleString,is_iodc=is_iodc, add_logo=add_logo,
                 logo_path_MF=logo_path_MF, logo_path_SAF = logo_path_SAF, 
                 figsize=figsize,cmap=clrmap, cTickLabels=tick_labels,
                 mapLabels=maplabels, suppressFig=suppressFig, units=units)
if __name__ == '__main__':
    main()
