# Map Meteosat Products
Scripts to map Meteosat products on Cartopy's [Geostationary](https://scitools.org.uk/cartopy/docs/latest/crs/projections.html#geostationary)  projection.

## Usage
Copy the example script **make_map_example.py** and (optionally) the config file **config.yml** to your working directory of choice. 
In that directory, run:

```
./make_map_example.py
```
This will create a map in png format in the specified location.
Adjust the script to your needs.

If you copied the config file and want to use your local version add the ```cfgFile``` argument to the call of mm.plot_msg_geoloc (within make_map_example.py), like so:
```mm.plot_msg_geoloc(..., cfgFile = 'config.yml', ...)```, or whatever the name/location of your own config file is.

Necessary parameters for all types are ```f_in_tplt``` (the data file) and ```varname``` (defining the variable to plot).


## Supported Products
*Product name: necessary parameters*

* Albedo
    - Broadband Bi-Hemispherical: ```varname = 'AL-BB-BH'```
    - Broadband Directional Hemispherical: ```varname = 'AL-BB-DH'```
    - Near-Infrared Directional Hemispherical: ```varname = 'AL-NI-DH```
    - Visual Directional Hemishperical: ```varname = 'AL-VI-DH'``` 
    - Q-Flags: ```varname = 'ALB-Q'```
    - example file name: ```'HDF5_LSASAF_MSG_ALBEDO_MSG-Disk_200401190000'```
 
* Aerosol Optical Depth (AOD)
    - CAMS-Climate: ```varname = 'AOD-CLIM'```
        - location on IPMA server: ```'/home/safpt/OperationalChain/LSASAF_Inputs/StaticData/MBRDF/'```
        - example file name: ```'HDF5_CLIMA_MF_O3_MSG-Disk_201201010000'```
    - CAMS-Forecast: ```varname = 'AOD-FC'```
        - location on IPMA server: ```'/home/safpt/OperationalChain/LSASAF_Inputs/NWP/AOD550/'```
        - example file name: ```'HDF5_R202005140000_MACC_AOD550_MSG-Disk_202005142100'```

* Age of last determined Albedo value: ```varname = 'Z_Age'```
    - in Albedo file

* Digital Elevation Model (DEM): ```varname = 'DEM'```
    - example file name: ```'HDF5_LSASAF_USGS_DEM_MSG-Disk_201807110815'```

* Ozone
    - TOMS O3 Climatology: ```varname = 'O3-CLIM'```
         - example file name: ```'HDF5_CLIMA_MF_O3_MSG-Disk_201201010000'```
    - ECMWF/CAMS O3 Forecast: ```varname = 'O3-FC'```
         - location on IPMA server: ```'/home/safpt/OperationalChain/LSASAF_Inputs/NWP/TCO3/'```
         - example file name: ```'HDF5_R202005091200_ECMWF_TCO3_MSG-Disk_202005100745'```
    - units can be converted with ```convertUnits = True'```

* Pressure at Mean Sea Level: ```varname = 'P-MSL'``` 
    - location on IPMA server: ```'/home/safpt/OperationalChain/LSASAF_Inputs/NWP/MSL/'```
    - example file name: ```'HDF5_R202008131200_ECMWF_MSL_MSG-Disk_202008141045'```

* Top of Atmosphere Radiance: ```varname='TOA', f_aux = '/path/to/sza-file'```
    - location on IPMA server: ```'/home/safpt/OperationalChain/LSASAF_Inputs/MSG/RAD/'```
    - example file name: ```'HDF5_006_MSG_RAD_MSG-Disk_201609111215'```
    - ex. aux file name: ```'HDF5_LSASAF_MSG_SZA_MSG-Disk_201609111215.h5'```

* Top of Canopy: 
    - TOC reflectance: ```varname = 'TOC'```
    - Q-Flags: ```varname = TOC-Q'```
    - location on IPMA server: ```'/home/safpt/OperationalChain/LSASAF_InternalProducts/MBRDF/'```
    - example file name: ```'HDF5_LSASAF_MSG_006-BRF_MSG-Disk_202010021045'```

* Water Vapour Forecast: ```varname = 'WV'```
    - location on IPMA server: ```'/home/safpt/OperationalChain/LSASAF_Inputs/NWP/TCWV/'```
    - example file name: ```'HDF5_R202010011200_ECMWF_TCWV_MSG-Disk_202010021145'```
    
## Configuration file
The scripts rely on a configuration file in YAML format. For your own purposes you can copy the config.yml file to one of your own directories, adjust it to your needs and use it as input. The file contains values for satellite product specific (and general) parameters. Some parameters can be adjusted through the function call as well, so changing the config file may not always be necessary. If you want to do that but you are not familiar with YAML, take a look [here](https://camel.readthedocs.io/en/latest/yamlref.html) and [here](https://rollout.io/blog/yaml-tutorial-everything-you-need-get-started/). 
