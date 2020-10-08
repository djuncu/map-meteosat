#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 13:49:51 2018

@author: vincentch
"""
import pandas as pd
import xarray as xr
import coloredlogs, logging
coloredlogs.install(level='INFO')

class DatasourcePreprocessingFunctions:
    """ Class containing all the preprocessing methods to call when reading each source dataset.
    Correspondance between functions and dataset names are given in the preprocessing_functions dictionary in AllSatConfig"""
    ### Preprocessing functions applied when opening the source datasets
    def open_mfdataset_preproc_for_modis6(self,ds):
        """ Preprocessing customised for MODIS collection 6 albedo data & qflag:
              adds date info (from attributes) as a standalone variable called 'datetime'
              renames geographical coordinates to 'x' and 'y' """
        logging.debug(f"open_mfdataset_preproc_for_modis6")
        filename = ds.attrs['CoreMetadata.0'].split("LOCALGRANULEID")[1].split('VALUE')[1].split('=')[1].split("END_OBJECT")[0].split('\n')[0].strip()
        datestr = filename.split('.A')[1].split('.')[0].strip()
        date = pd.to_datetime(datestr,format='%Y%j') # %j=DayOfYear as in modis dates
        ds.coords['datetime'] = [date]

        ds.rename(name_dict={'YDim:Grid_Parameter':'y','XDim:Grid_Parameter':'x'}, inplace=True)

        return ds

    def open_mfdataset_preproc_for_LSAF_msg_eps(self,ds):
        """ Preprocessing customised for LSASAF albedo data (MSG and EPS in HDF5 format):
              adds date info (from attributes) as a standalone variable called 'datetime'
              renames geographical coordinates to 'x' and 'y' """
        logging.debug(f"open_mfdataset_preproc_for_LSAF_msg_eps")
        datestr = ds.attrs['SENSING_START_TIME'][0:8] ### SENSING_START_TIME=20160105000000
        date = pd.to_datetime(datestr,format='%Y%m%d')
        if 'datetime' not in ds.dims : ds = ds.expand_dims(dim='datetime')
        ds.coords['datetime'] = [date]

        ### h5zip creates a different dimension name for each variable: we reassign them to common 'y' and 'x' dimensions
        first = None      
        for v in ds.data_vars:
            da = ds[v]
            dim_x = [x for x in da.dims if 'dim1' in x ][0]
            dim_y = [x for x in da.dims if 'dim0' in x ][0]
            da = da.rename({dim_x:'x',dim_y:'y'})
            if first == None:              
                new_ds = da.to_dataset()
                first = 'no more first'
            else:
                new_ds[da.name] = da
        return new_ds

    def open_mfdataset_preproc_for_CGLS_VGT(self,ds):
        """ Preprocessing customised for VGT GEOV1 albedo data in NC4 format (source CGLS):
              adds date info (from attributes) as a standalone variable called 'datetime'
              renames geographical coordinates to 'x' and 'y' """
        logging.debug(f"open_mfdataset_preproc_for_CGLS_VGT")
        datestr = ds.attrs['TEMPORAL_NOMINAL']        ###  TEMPORAL_NOMINAL=2003-12-24
        date = pd.to_datetime(datestr,format='%Y-%m-%d')
        ds.coords['datetime'] = [date]
        ds.rename(name_dict={'phony_dim_0':'y','phony_dim_1':'x'}, inplace=True)
        return ds


#class ReaderAndPreprocessing:
#    """ Classe qui lit et applique le preprocessing aux donnees d'entree. 
#    > Renvoie un dataset xarray avec des variables aux noms standardises et les dimensions 'datetime','x' et 'y'. """
#
#    def read_and_rename_all_variables_in_ds_withoutclassconfig(self,is_h5zip_needed, preproc_fct, input_varname, file_list):
#        file_list=list(set(file_list)) # returns unique values
#        logging.debug(f'Opening multiple datasets with appropriate reader')
#        if is_h5zip_needed:
#            kw_args= {'preprocess':preproc_fct,'concat_dim':'datetime','engine':'h5netcdf'} ## attention: the 'autoclose' option does not work with h5netcdf engine
#        else:
#            kw_args= {'preprocess':preproc_fct,'concat_dim':'datetime','autoclose':True}
#        ds = xr.open_mfdataset(file_list,**kw_args)
#        ### RENAME SAT VARIABLES WITH THE STANDARD NAMING CHOSEN - only works if the dict corresponds to the actual list of variables in ds
#        logging.debug(f'Renaming variables names')
#        dict_only_ds_vars={}
#        for key,var_in in input_varname.items():
#            if key in ds.data_vars:
#                dict_only_ds_vars[key]=var_in
#        if dict_only_ds_vars != {} : ds.rename(name_dict=dict_only_ds_vars, inplace=True)
#        logging.debug(f'End renaming variables names')
#        return ds
#
#
#    def read_and_rename_all_variables_in_ds(self,instantaneous_config, file_list):
#        file_list=list(set(file_list)) # returns unique values
#        logging.debug(f'Opening multiple datasets with appropriate reader')
#        if instantaneous_config.is_h5zip_needed:
#            kw_args= {'preprocess':instantaneous_config.preproc_fct,'concat_dim':'datetime','engine':'h5netcdf'} ## attention: the 'autoclose' option does not work with h5netcdf engine
#        else:
#            kw_args= {'preprocess':instantaneous_config.preproc_fct,'concat_dim':'datetime','autoclose':True}
#        ds = xr.open_mfdataset(file_list,**kw_args)
#        ### RENAME SAT VARIABLES WITH THE STANDARD NAMING CHOSEN - only works if the dict corresponds to the actual list of variables in ds
#        logging.debug(f'Renaming variables names')
#        dict_only_ds_vars={}
#        for key,var_in in instantaneous_config.input_varname.items():
#            if key in ds.data_vars:
#                dict_only_ds_vars[key]=var_in
#        if dict_only_ds_vars != {} : ds.rename(name_dict=dict_only_ds_vars, inplace=True)
#        logging.debug(f'End renaming variables names')
#        return ds
#
