# This library provides useful tools to analyze and process netCDF climate data such as ERA5 climate data.

In  order to use the library without problems, make sure that the names of the main dimensions are as follows: 'time', 'longitude', 'latitude'.
For more information, you can also check the  uploaded jupyter notebook named "How to use era5pplib".

## Read the shapefile

def shp(directory, type):
    
    '''
    This function takes the directory of a shapefile and it return a visualized plot of the shapefile for fast analysis.
    The function returns two types of plots, either a filled shapefile or a boundary. type: 'filled' or 'boundary'.

    '''
    
## Plot a map with layers

def layerplt(raster, shapefile):
    
    '''
    This function takes an uploaded raster image and a shapefile directory, 
    and it returns a layered plot of the shapefile and the raster.

    '''

## How to clip the data

def clip(nc4file, coordinates, clippingtype, input):

    '''
    This function takes an uploaded netCDF file, coordinates, and clips the rasters using a shapefile or a geometry  
    The coordinates can be obtained from the following site: https://epsg.io/
    The function also requires the identification of 
    The clippingtype is either 'geometry' or 'shapefile'
    After selecting the clipping type, the input is either a shapefile directory or a geometry list [x1, y1, x2, y2].
    
    '''

## Averaging over the clip

def avgclip(clipped, var_name, output):

    '''
    This function takes a clipped file for one variable and calculates the average over the clipped area and returns timeseries data.
    The timeseries data will be saved in a csv file, the output is the name of the output file.
    
    '''


## Aggragate the yearly netCDF files into one netCDF

def nc4aggr(directory, listrange, outdir):

    '''
    This function takes the yearly separated netCDF files and aggregates them into one netCDF file. 
    Put the netCDF files in one directory and list them in order from one: 1, 2, ...
    Add the listrange as follows: listrange = (1, last number in the list)
    The function saves the aggregated data into a netCDF file, the outdir is the directory where the file will be downloaded

    '''

## Calculate the daily average

def resdaily(directory, type,  outdir):

    '''
    This function resamples the hourly netcdf data into daily netcdf data.
    The function resamples based on the selected type. there are only two types to select. 'mean' or 'sum'

    '''

## Calculate the relative humidity
    
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

## Calculate the wind speed

def W_velocity(u_directory, uvar_name, v_directory, vvar_name, out_dir):
    
    '''
    This function takes the U and V components of wind/velocity and calculates the wind speed using the following equation:
    Wind speed = sqrt(u**2 + v**2)
    The results will be saved in a separate netCDF file using the provided out_dir (example out_dir = '..\..\RH.nc')
  
    '''
