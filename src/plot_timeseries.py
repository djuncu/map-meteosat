#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 09:41:41 2018

@author: vincentch
"""
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np

def set_date_ticks(fig, ax, freq_major, freq_minor, rotation=None):
    """ Customises the dates ticking spacing and formats on axis 'x'
    (in particular, the weekly ticking step chosen by default when plotting only one or 2 months of data is not nice - end of month and begin of months are too close)
    """
    import matplotlib.dates as mdates
    import locale
    locale.setlocale(locale.LC_ALL,'en_GB.utf8') # set default language to English - required when showing months full name
    
    loc_list=[]
    fmt_list=[]
    
    for i,freq in enumerate([freq_major,freq_minor]):
        if freq == None:
            loc = None
            fmt = None
        if freq == 'year':
            loc = mdates.YearLocator()   # every year
            fmt = mdates.DateFormatter('%Y')
        if freq == 'month':
            loc = mdates.MonthLocator()  # every month
            fmt = mdates.DateFormatter('%Y-%m')
        if freq == 'justmonthtextfull':
            loc = mdates.MonthLocator()  # every month
            fmt = mdates.DateFormatter('%B') # eg. 'September'
        if freq == 'justmonthtextshort':
            loc = mdates.MonthLocator()  # every month
            fmt = mdates.DateFormatter('%b') # eg. 'Sept.'
        if freq == '10days':
            loc = mdates.DayLocator(bymonthday=[1,10,20])
            fmt = mdates.DateFormatter('%Y-%m-%d')
        if freq == 'SAF10days':
            loc = mdates.DayLocator(bymonthday=[5,15,25])
            fmt = mdates.DateFormatter('%Y-%m-%d')
        if freq == 'weekinmonth':
            loc = mdates.DayLocator(bymonthday=[1,8,15,22])
            fmt = mdates.DateFormatter('%Y-%m-%d')
        if freq == 'day':
            loc = mdates.DayLocator()  # every day
            fmt = mdates.DateFormatter('%Y-%m-%d')
        if freq == 'hour_with_days':
            loc = mdates.HourLocator()  # every hour
            fmt = mdates.DateFormatter('%Y-%m-%d\n%H:%M')
        if freq == 'hour_without_days':
            loc = mdates.HourLocator()  # every hour
            fmt = mdates.DateFormatter('%H:%M')
        loc_list.append(loc)
        fmt_list.append(fmt)
    
    if freq_major != None: ax.xaxis.set_major_locator(loc_list[0])
    if freq_major != None: ax.xaxis.set_major_formatter(fmt_list[0])
    if freq_minor != None: ax.xaxis.set_minor_locator(loc_list[1])
    if freq_minor != None: ax.xaxis.set_minor_formatter(fmt_list[1])
    fig.autofmt_xdate() # automatically sets ticks labels rotation and adjust bottom space

    if rotation != None: plt.xticks(rotation=rotation, ha="left")
    

def plot_time_series_from_df(df, datetime_start = None, datetime_end = None, plot_mode = "bars", 
                             vars_to_plot_without_errors = None, vars_to_plot_without_errors_kwargs = {}, 
                             vars_to_plot_with_errors = (None,None,None), vars_to_plot_with_errors_kwargs = {}, 
                             ref_threshold_lines_values = [], ref_threshold_lines_values_kwargs = {'linestyle':'-'},
                             yscale_in_percent = False, title = None, legend_loc = None, grid=None, major_date_ticks=None, minor_date_ticks = None, rotation = None, figsize=None):
    """ Plotting the chosen data with or without error bars or shaded areas
    
    Parameters
        > TODO
        df: Pandas Dataframe
        plot_mode: "bars" for line plots with error bars 
                            none, symetric or up and down
                    "shaded" : fill_between two lines
        
        vars_to_plot_without_errors : None or tuple(size=3) or list of tuples(size 3)
        vars_to_plot_with_errors : None or tuple(size=4) or list of tuples(size 4). eg. ('varname',data_name_in_df,data_lower_limit_to_show,date_upper_limit_to_show,err_plot_kwargs)
        
        errorbar args : see https://matplotlib.org/api/_as_gen/matplotlib.pyplot.errorbar.html
        plt.plot args : > TODO
        grid: None, 'x', 'y' or 'both'. If None, no grid is shown. otherwise, provides the name of the axes where to plot the major ticks grid.
        figsize: tuple. Size of the output figures in inches. used by `plt.subplots()`
        
    Example
        station = 'Agoufou'
        coordinates = '{lon=xx°E; lat=xx°N}'
        plot_time_series_from_df(df[ df['station_name'] == station ], vars_to_plot_without_errors = 'wsa_shwv_etalv4_mean', vars_to_plot_with_errors = ('wsa_shwv_modis6_mean','wsa_shwv_modis6_percentile10','wsa_shwv_modis6_percentile90'), title = f'WSA_SHWV spatial mean on 0.5deg x 0.5deg zone \n around {station} {coordinates}')
        plot_time_series_from_df(df[ df['station_name'] == station ], vars_to_plot_without_errors = 'wsa_shwv_etalv4_mean', vars_to_plot_without_errors_kwargs = [{'linestyle':'None', 'marker':'.','label':'wsa_shwv_etalv4_mean'}], vars_to_plot_with_errors = ('wsa_shwv_modis6_mean','wsa_shwv_modis6_percentile10','wsa_shwv_modis6_percentile90'), vars_to_plot_with_errors_kwargs = [{'linestyle':'None', 'marker':'.','label':'wsa_shwv_modis6_mean'}], title = f'WSA_SHWV spatial mean on 0.5deg x 0.5deg zone \n around {station} {coordinates}')
        plot_time_series_from_df(df[ df['station_name'] == station ], plot_mode = 'shaded', 
                                 vars_to_plot_without_errors = 'wsa_shwv_etalv4_mean', 
                                 vars_to_plot_without_errors_kwargs = [{'linestyle':'None', 'marker':'.','label':'wsa_shwv_etalv4_mean'}], 
                                 vars_to_plot_with_errors = ('wsa_shwv_modis6_mean','wsa_shwv_modis6_percentile10','wsa_shwv_modis6_percentile90'), 
                                 vars_to_plot_with_errors_kwargs = [{'color':'coral','linestyle':'None', 'marker':'.','label':'wsa_shwv_modis6_mean'}], 
                                 title = f'WSA_SHWV spatial mean on 0.5deg x 0.5deg zone \n around {station} {coordinates}')

        df['mean_bias_etal_modis'] = (df['wsa_shwv_etalv4_mean'] - df['wsa_shwv_modis6_mean'])/df['wsa_shwv_modis6_mean']
        plot_time_series_from_df(df[ df['station_name'] == station ], vars_to_plot_without_errors = 'mean_bias_etal_modis', 
                                 vars_to_plot_without_errors_kwargs = [{'linestyle':'None', 'marker':'.','label':'(etal-modis)/modis'}], 
                                 title = f'WSA_SHWV relative bias (etal-modis)/modis on 0.5deg x 0.5deg zone \n around {station} {coordinates}', 
                                 yscale_in_percent = True, legend_loc = 'upper right',
                                 ref_threshold_lines_values = [-0.05,0.05,0.1,-0.1], 
                                 ref_threshold_lines_values_kwargs = [{'linestyle':'--','color':'green'},{'linestyle':'--','color':'green'},{'linestyle':'--','color':'orange'}, {'linestyle':'--','color':'orange'}])
    """
   
    if datetime_start != None:
        df = df[df['datetime']>=datetime_start ]
    if datetime_end != None:
        df = df[df['datetime']<=datetime_end ]
        
    df = df.sort_values(by='datetime')
    
    if vars_to_plot_without_errors == None:
        vars_to_plot_without_errors =[]
    elif isinstance(vars_to_plot_without_errors,str):
        vars_to_plot_without_errors = [ vars_to_plot_without_errors ]
    if vars_to_plot_with_errors == (None,None,None):
        vars_to_plot_with_errors =[]
    elif isinstance(vars_to_plot_with_errors,tuple):
        vars_to_plot_with_errors = [ vars_to_plot_with_errors ]
        
    if isinstance(vars_to_plot_with_errors_kwargs,dict):
        vars_to_plot_with_errors_kwargs = [ vars_to_plot_with_errors_kwargs ]
    if isinstance(vars_to_plot_without_errors_kwargs,dict):
        vars_to_plot_without_errors_kwargs = [ vars_to_plot_without_errors_kwargs ]
        
    allowed_plot_modes = ['bars', 'shaded']
        
    if plot_mode not in allowed_plot_modes: 
        print(f'ERROR: plot_mode given ({plot_mode}) should be within: {allowed_plot_modes}')
    else:
        fig,ax = plt.subplots(figsize=figsize)
        if ref_threshold_lines_values != []:
            for v,kw in zip(ref_threshold_lines_values,ref_threshold_lines_values_kwargs):
                ax.plot(df["datetime"], np.ones(len(df["datetime"]))*v, **kw)
        
        if plot_mode == 'bars':
            for var_nobars, var_nobars_kwargs in zip(vars_to_plot_without_errors,vars_to_plot_without_errors_kwargs): 
                eb=ax.plot(df["datetime"], df[var_nobars], **var_nobars_kwargs) ## TODO correct resample outputs so that the scale_factor is correctly applied for ETAL data
                
            for (var_bars,line_down,line_up),var_bars_kwargs in zip(vars_to_plot_with_errors, vars_to_plot_with_errors_kwargs): 
                if line_down == None:
                    err = df[line_up]
                elif line_up == None:
                    err = df[line_down]
                else:
                    err = np.ones(shape=(2,len(df[line_down])))*np.NaN
                    err[0] = abs(df[var_bars] - df[line_down])
                    err[1] = abs(df[var_bars] - df[line_up])
                eb=ax.errorbar(df["datetime"].dt.to_pydatetime(), df[var_bars], yerr = err, **var_bars_kwargs) # pd.Series().dt.to_pydatetime() required as there is a bug within errorbars which does not support normal datetime64 objects
                eb[-1][0].set_linestyle('--')  #eb1[-1][0] is the LineCollection objects of the errorbar lines
            
        elif plot_mode == 'shaded':
            for var_nobars, var_nobars_kwargs in zip(vars_to_plot_without_errors,vars_to_plot_without_errors_kwargs): 
                ax.plot(df["datetime"].dt.to_pydatetime(),df[var_nobars], **var_nobars_kwargs)
                
            for (var_bars,line_down,line_up),var_bars_kwargs in zip(vars_to_plot_with_errors, vars_to_plot_with_errors_kwargs):
                ax.plot(df["datetime"].dt.to_pydatetime(),df[var_bars],**var_bars_kwargs)
                kw={}
                if 'color' in var_bars_kwargs.keys():
                    kw = {'facecolor':var_bars_kwargs['color']}
                ax.fill_between(df["datetime"].dt.to_pydatetime(), df[line_down], df[line_up], alpha=0.2, **kw)

        if yscale_in_percent:
            vals = ax.get_yticks()
            ax.set_yticklabels(['{:,.1%}'.format(x) for x in vals])

        if title != None: ax.set_title(title)
        ax.legend(loc=legend_loc)
        #fig.autofmt_xdate() # automatically sets ticks spacing, labels rotation and adjust bottom space
        if grid != None: plt.grid(which='major', axis=grid, color='lightgrey', linestyle='-.')
        set_date_ticks(fig=fig, ax=ax, freq_major=major_date_ticks, freq_minor=minor_date_ticks, rotation=rotation)
    return fig,ax


