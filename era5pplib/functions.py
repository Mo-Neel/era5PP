# Import the libraries

import xarray as xr
import netCDF4 
import matplotlib.pyplot as plt
import earthpy as et
import earthpy.spatial as es
import earthpy.plot as ep
import rioxarray as rxr
from rasterio.plot import plotting_extent
import pandas as pd
import geopandas as gpd
from shapely.geometry import box
import geopandas
from shapely.geometry import mapping
import numpy as np

# Read the shapefile

def shp(directory, type):
    
    '''
    This function takes the directory of a shapefile and it return a visualized plot of the shapefile for fast analysis.
    The function returns two types of plots, either a filled shapefile or a boundary. type: 'filled' or 'boundary'.

    '''
    sh = gpd.read_file(directory)

    if type == 'filled':
        # plot the watershed
        return sh.plot()

    if type == 'boundary':
        # plot the watershed boundary
        return sh.boundary.plot()
    
# Plot a map with layers

def layerplt(raster, shapefile):
    
    '''
    This function takes an uploaded raster image and a shapefile directory, 
    and it returns a layered plot of the shapefile and the raster.

    '''
    # open the shapefile
    sh = gpd.read_file(shapefile)

    # Open the raster data
    data_plotting_extent = plotting_extent(raster, raster.rio.transform())

    # Plot cropped data
    f, ax = plt.subplots()

    ep.plot_bands(raster, extent = data_plotting_extent,
              ax=ax)

    sh.boundary.plot(ax=ax)

    plt.show()

# How to clip the data

def clip(nc4file, coordinates, clippingtype, input):

    '''
    This function takes an uploaded netCDF file, coordinates, and clips the rasters using a shapefile or a geometry  
    The coordinates can be obtained from the following site: https://epsg.io/
    The function also requires the identification of 
    The clippingtype is either 'geometry' or 'shapefile'
    After selecting the clipping type, the input is either a shapefile directory or a geometry list [x1, y1, x2, y2].
    
    '''
    xds = nc4file

    xds.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude", inplace=True)
    
    xds.rio.write_crs(coordinates, inplace=True)

    # Clipping using the shapefile or a geometry

    if clippingtype == 'shapefile':
        
        geodf = gpd.read_file(input)

        clipped = xds.rio.clip(geodf.geometry.apply(mapping), all_touched=True)

        return clipped

    # clipping using geometry

    if clippingtype == 'geometry':

        geodf = geopandas.GeoDataFrame(geometry=[box(input[0], input[1], input[2], input[3])],)

        clipped = xds.rio.clip(geodf.geometry.apply(mapping), all_touched=True)

        return clipped

# Averaging over the clip

def avgclip(clipped, var_name, output):

    '''
    This function takes a clipped file for one variable and calculates the average over the clipped area and returns timeseries data.
    The timeseries data will be saved in a csv file, the output is the name of the output file.
    '''
    avg = clipped[var_name].mean(dim = ('longitude', 'latitude'))

    df = avg.to_dataframe()
    
    df.to_csv(output + '.csv')


# Aggragate the yearly netCDF files into one netCDF

def nc4aggr(directory, listrange, outdir):

    '''
    This function takes the yearly separated netCDF files and aggregates them into one netCDF file. 
    Put the netCDF files in one directory and list them in order from one: 1, 2, ...
    Add the listrange as follows: listrange = (1, last number in the list)
    The function saves the aggregated data into a netCDF file, the outdir is the directory where the file will be downloaded

    '''

    url = [] 
    for i in range(listrange[0], (listrange[1]+1)):
        f = directory + '/' + "{}".format(i) +".nc".format(i)         
        url.append(f)
    
    x = xr.open_dataset(url[0]) 
    for i in range(1, (listrange[1])):                 
        d = xr.open_dataset(url[i])
        x = xr.concat((x,d), dim = 'time')

    x.to_netcdf(outdir)

# Calculate the daily average

def resdaily(directory, type,  outdir):

    '''
    This function resamples the hourly netcdf data into daily netcdf data.
    The function resamples based on the selected type. there are only two types to select. 'mean' or 'sum'

    '''
    x = xr.open_dataset(directory)

    if type == 'mean':

        daily = x.resample(time = 'd').mean(dim = 'time') ## Here we should make sure the parameters are to be averaged not aggregated

        daily.to_netcdf(outdir)

    if type == 'sum':

        daily = x.resample(time = 'd').sum(dim = 'time') ## Here we should make sure the parameters are to be aggregated

        daily.to_netcdf(outdir)

# Calculate the relative humidity
    
def rh(T_directory, Tvar_name, Td_directory, Tdvar_name, out_dir):
    
    '''
    This function takes the temperature as input (T is the air temperature and Td is dew point temperature) in Kelvin (K) and converts it to Degrees Celsius (°C)  
    This function calculates the relative humidity using the Bolton equation which is used for a range of temperature between -30 °C to +35 °C.
    e = 6.112 * np.exp((17.67 * Td)/(Td + 243.5))
    es = 6.112 * np.exp((17.67 * T)/(T + 243.5))
    T = air_temprature in °C
    Td = dew_poin_temprature °C
    The results will be saved in a separate netCDF file using the provided out_dir (example out_dir = '..\..\RH.nc')
    
    '''
    
    t2m = xr.open_dataset(T_directory)
    d2m = xr.open_dataset(Td_directory)
    t2m['t2m_c'] = t2m[Tvar_name] - 273.15
    d2m['d2m_c'] = d2m[Tdvar_name] - 273.15
    # Calculate the actual vapor pressure
    d2m['e'] = 6.112 * np.exp((17.67 * d2m['d2m_c'])/(d2m['d2m_c'] + 243.5))
    # Calculate the standard water vapor pressure
    t2m['es'] = 6.112 * np.exp((17.67 * t2m['t2m_c'])/(t2m['t2m_c'] + 243.5))
    # Calculate the relative humidity
    t2m['rh'] = (d2m['e']/t2m['es'])*100
    # save the results to an output directory
    t2m['rh'].to_netcdf(out_dir)

# Calculate the wind speed

def W_velocity(u_directory, uvar_name, v_directory, vvar_name, out_dir):
    
    '''
    This function takes the U and V components of wind/velocity and calculates the wind speed using the following equation:
    Wind speed = sqrt(u**2 + v**2)
    The results will be saved in a separate netCDF file using the provided out_dir (example out_dir = '..\..\RH.nc')
  
    '''
    u10 = xr.open_dataset(u_directory)
    v10 = xr.open_dataset(v_directory)

    u10['wv'] = np.sqrt((u10[uvar_name]**2) + (v10[vvar_name]**2))

    u10['wv'].to_netcdf(out_dir)
