# era5PP

`era5PP` is a Python package for processing and analysing netCDF climate datasets, including ERA5 reanalysis data. It provides functions for spatial visualisation, spatial clipping, area-averaged time series, temporal aggregation, and the calculation of derived climate variables such as relative humidity and wind speed.

## Installation

Install the package directly from GitHub using pip:

```bash
pip install "git+https://github.com/Mo-Neel/era5PP.git"
```

## Requirements and Data Conventions

Before using the package, ensure that the main temporal and spatial dimensions in your netCDF datasets are named:

- `time`
- `longitude`
- `latitude`

The relevant functions expect these dimension names. If your dataset uses different names, rename the dimensions before processing it.

For additional guidance, refer to the Jupyter notebook **How to use era5pplib** included in the repository.

## Functions

### 1. Spatial Visualisation

#### `shp(directory, type)`

Reads a shapefile and generates a plot for quick visual inspection.

**Parameters**

- `directory`: Path to the shapefile.
- `type`: Type of plot to generate. Accepted values are:
  - `'filled'`: Displays the shapefile as filled polygons.
  - `'boundary'`: Displays the boundaries of the shapefile.

**Output:** A visualisation of the shapefile.

#### `layerplt(raster, shapefile)`

Creates a layered plot combining a raster dataset with a shapefile, allowing the spatial relationship between the raster data and geographic boundaries to be inspected.

**Parameters**

- `raster`: Raster data to be displayed.
- `shapefile`: Path to the shapefile to overlay on the raster.

**Output:** A plot showing the raster and shapefile layers together.

### 2. Spatial Processing

#### `clip(nc4file, coordinates, clippingtype, input)`

Clips a netCDF dataset to a specified geographic area using either a geometry or a shapefile. This function can be used to extract climate data for a particular region of interest.

**Parameters**

- `nc4file`: Path to the input netCDF file.
- `coordinates`: Coordinate-system information required by the function. Coordinate reference system information can be obtained from [EPSG.io](https://epsg.io/).
- `clippingtype`: Specifies the clipping method. Accepted values are:
  - `'geometry'`: Clips the dataset using a specified geometry.
  - `'shapefile'`: Clips the dataset using a shapefile.
- `input`: Defines the clipping area according to the selected `clippingtype`. This is either a geometry list in the form `[x1, y1, x2, y2]` or a path to a shapefile.

For geometry-based clipping, the coordinate values must correspond to the expected coordinate reference system and ordering.

**Output:** The netCDF data clipped to the specified area.

#### `avgclip(clipped, var_name, output)`

Calculates the spatial average of a specified variable over a previously clipped area and saves the resulting time series to a CSV file. This function is useful for extracting area-averaged climate information for further analysis.

**Parameters**

- `clipped`: The clipped netCDF dataset.
- `var_name`: Name of the variable for which the spatial average is calculated.
- `output`: Name or path of the output CSV file.

**Output:** A CSV file containing the area-averaged time series of the selected variable.

### 3. Temporal Processing

#### `nc4aggr(directory, listrange, outdir)`

Aggregates yearly netCDF files stored in a directory into a single netCDF file. This function is useful when climate data are stored in separate files for individual years and need to be combined into one dataset.

**Parameters**

- `directory`: Path to the directory containing the yearly netCDF files.
- `listrange`: Tuple specifying the range of file indices to aggregate, in the form `(1, n)`, where `n` is the index of the last file to include. The files must be organised and numbered consistently according to the function's expected ordering.
- `outdir`: Output path for the aggregated netCDF file.

**Output:** A single netCDF file containing the aggregated data from the selected yearly files.

#### `resdaily(directory, type, outdir)`

Resamples hourly netCDF data to daily values using the selected aggregation method.

**Parameters**

- `directory`: Path to the input netCDF file or the directory expected by the function.
- `type`: Daily aggregation method. Accepted values are:
  - `'mean'`: Calculates the daily mean.
  - `'sum'`: Calculates the daily sum.
- `outdir`: Output path for the resulting netCDF file.

**Output:** A netCDF file containing daily-resampled data.

### 4. Derived Climate Variables

#### `rh(T_directory, Tvar_name, Td_directory, Tdvar_name, out_dir)`

Calculates relative humidity from air temperature and dew-point temperature using the Bolton saturation vapour pressure formulation.

Both temperature variables must be provided in Kelvin. The function converts the temperatures to degrees Celsius before calculating relative humidity.

The calculation uses the following equations:

e = 6.112 * np.exp((17.67 * Td)/(Td + 243.5))

e_s = 6.112 * np.exp((17.67 * T)/(T + 243.5))

RH = 100 * e/e_s

where:

- (T) is the air temperature in degrees Celsius.
- (T_d) is the dew-point temperature in degrees Celsius.
- (e) is the actual vapour pressure in hPa.
- (e_s) is the saturation vapour pressure in hPa.
- (RH) is the relative humidity expressed as a percentage.

The Bolton formulation is intended for a temperature range of approximately −30 °C to +35 °C, as specified in the original function documentation.

**Parameters**

- `T_directory`: Path to the netCDF file containing air temperature.
- `Tvar_name`: Name of the air temperature variable in the netCDF file.
- `Td_directory`: Path to the netCDF file containing dew-point temperature.
- `Tdvar_name`: Name of the dew-point temperature variable in the netCDF file.
- `out_dir`: Output path for the resulting netCDF file.

**Output:** A netCDF file containing the calculated relative humidity.

#### `W_velocity(u_directory, uvar_name, v_directory, vvar_name, out_dir)`

Calculates wind speed from the two horizontal wind components: the zonal component (\(u\)) and the meridional component (\(v\)).

Wind speed is calculated using:

Wind speed = sqrt(u**2 + v**2)

where:

- (u) is the zonal wind component.
- (v) is the meridional wind component.
- (W) is the resulting wind speed.

The input components must use compatible units. If they are expressed in metres per second, the resulting wind speed will also be in metres per second.

**Parameters**

- `u_directory`: Path to the netCDF file containing the zonal wind component.
- `uvar_name`: Name of the zonal wind variable in the netCDF file.
- `v_directory`: Path to the netCDF file containing the meridional wind component.
- `vvar_name`: Name of the meridional wind variable in the netCDF file.
- `out_dir`: Output path for the resulting netCDF file.

**Output:** A netCDF file containing the calculated wind speed.

## Additional Information

For further information on the library's functionality and practical demonstrations, consult the accompanying Jupyter notebook **How to use era5pplib**.
