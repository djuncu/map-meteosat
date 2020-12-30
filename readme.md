# Map Meteosat Products
Scripts to map Meteosat products on Cartopy's [Geostationary](https://scitools.org.uk/cartopy/docs/latest/crs/projections.html#geostationary)  projection.

## Usage
to use this tool from any directory, link the plotting script ```plotMxg.py``` in a directory that is on your system path, say ```~/bin```:

```bash
# create directory if it  does not exist
mkdir ~/bin
# link script and give execution permission
ln -s /path/to/map-meteosat/plotMxg.py plotMxg && chmod 775 ~/bin/plotMxg
```

Make sure the directory is on your path, by typing this command / adding it to your ~/.bashrc:
```
export PATH="~/bin:$PATH"
```

Now you can run it, e.g.:

```bash
# function help
plotMxg -h

# list variables in a file
plotMxg -l /path/to/file

# example run, TOC file
plotMxg -V BRF-TOC /path/to/file
```

## Config file
By default, the program will use the config.yml file in the tool's directory. If you want to modify it for your purpose, copy it to a different location and call your own config with the -c option.

