#!/usr/bin/env python3

import sys

# add the location of the plotting codes to enable import
sys.path.append('/cnrm/vegeo/juncud/codes/map_meteosat')

import meteosat_map as mm

#import os
#import inspect
#print(inspect.getfile(mm))

# Path definitions for example files a number of variables (might not be up to date)
# Albedo file
f_alb      = '/cnrm/vegeo/juncud/examples/MDAL/HDF5_LSASAF_MSG_ALBEDO_MSG-Disk_202009040000'
# BRF file
f_brf      = '/cnrm/vegeo/juncud/examples/MBRDF/HDF5_LSASAF_MSG_006-BRF_MSG-Disk_202009101115'
# Aerosol forecast file
f_aod_fc   = '/cnrm/vegeo/juncud/examples/TestingChain/NWP/AOD550/HDF5_R202009190000_MACC_AOD550_MSG-Disk_202009201500'
# Aerosol climate file
f_aod_clim = '/cnrm/vegeo/juncud/examples/TestingChain/StaticData/MBRDF/HDF5_MF-CAMS_MSG_AOD550_MSG-Disk_200009200000'
# DEM file
f_dem      = '/cnrm/vegeo/vincentch/DATA_DOC_RESULTS/2a__outputs_from_codes_for_analyses/NO_SAVE/MBRDF/v6.2.1__O3_forecasts/from_ipma/0deg_ope/static_inputs/HDF5_LSASAF_USGS_DEM_MSG-Disk_201807110815'
# Water vapour forecast file
f_wv       = '/cnrm/vegeo/vincentch/DATA_DOC_RESULTS/2a__outputs_from_codes_for_analyses/NO_SAVE/MBRDF/v6.2.1__O3_forecasts/from_ipma/0deg_ope/HDF5_R202008131200_ECMWF_TCWV_MSG-Disk_202008141045'
# Pressure at mean sea level file
f_msl      = '/cnrm/vegeo/juncud/examples/TestingChain/NWP/MSL/HDF5_R202009140000_ECMWF_MSL_MSG-Disk_202009141330'
# Ozone climate
f_oc_clim  = '/cnrm/vegeo/vincentch/DATA_DOC_RESULTS/2a__outputs_from_codes_for_analyses/NO_SAVE/MBRDF/v6.2.1__O3_forecasts/from_ipma/0deg_ope/TCO3/HDF5_CLIMA_MF_O3_MSG-Disk_201201010000'
# Ozone forecast file
f_oc_fc    = '/cnrm/vegeo/vincentch/DATA_DOC_RESULTS/2a__outputs_from_codes_for_analyses/NO_SAVE/MBRDF/v6.2.2__iodc_O3-FC_AOD-CAMS_CMa-2016/from_ipma/msg0deg_ope/TCO3/HDF5_R202005091200_ECMWF_TCO3_MSG-Disk_202005100745'
# Top of atmosphere radiance and auxiliary file
f_toa      = '/cnrm/vegeo/juncud/NO_SAVE/RAD/HDF5_006_MSG_RAD_MSG-Disk_202010061545'
f_toa_aux  = '/cnrm/vegeo/SAT/DATA/MSG/NRT-Operational/INPUTS/ANGLES/ANG-20160911/HDF5_LSASAF_MSG_SZA_MSG-Disk_201609111215.h5'

# call the plotting routine
mm.plot_msg_geoloc(f_in_tplt=f_alb, varname = 'AL-BB-BH',  vmin='default', vmax='default', f_out_png='~/map_test.png')
