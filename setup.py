from setuptools import find_packages, setup
setup(
    name='era5pplib',
    packages=find_packages(),
    version='0.1.0',
    description='Useful library when dealing with netcdf climate data such as era5 data',
    author='Mohamed Elbasheer',
    license='MIT',
    install_requires=["xarray>=2023.1.0", 
                      "rioxarray>=0.13.4",
                      "shapely>=2.0.1",
                      "matplotlib>=3.7.2",
                      "netCDF4>=1.6.4",
                      "rasterio>=1.3.8",
                      "pandas>=2.0.3",
                      "geopandas>=0.13.2"], 
    setup_requires=['pytest-runner'],
    tests_require=['pytest>=4.4.1'],
    test_suite='tests',
)