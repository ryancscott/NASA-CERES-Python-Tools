# ==============================================================================
#
#                   ----***** ERA5 REANALYSIS PYTHON TOOLS *****-----
#
# ==============================================================================
#
# Module:   era5tools.py
#
# Purpose:  This library contains Python3 functions to read, process, analyze
#           and plot data from the European Centre for Medium-Range Weather
#           Forecast's latest ERA5 reanalysis of the global climate.
#           Functions are included to ...
#           See function descriptions below for additional info.
#
# To Use:   import era5tools as era5
#
# Requires: numpy, matplotlib, netcdf4, pyhdf, cartopy, datetime, palettable
#           Recommend installing the above libraries using conda, pip, or an IDE
#
# Author:   Ryan C. Scott, ryan.c.scott@nasa.gov
#
# Last Updated: April 28, 2020
#
# Every attempt is made to conform to PEP-8 standards
# ==============================================================================
# LEVEL 2
# ---------------------
# get_date                    <- get the year, month, day, hour from input file
# get_platform


# ==============================================================================


def print_nc_file_info(file_path):
    """
    ----------------------------------------------------------------------------
    This function prints basic information about the variables in netCDF file
    ----------------------------------------------------------------------------
    :param file_path: path to file  [string]
    ----------------------------------------------------------------------------
    :return: prints information to screen
    """
    print('====================================')
    print('\t\tVariables in file...')
    print('====================================')
    print('\tname - units - dimensions')
    print('====================================')

    from netCDF4 import Dataset

    nc = Dataset(file_path, 'r')
    # print information about variables in netCDF file
    for i in nc.variables:
        print(i, '-', nc.variables[i].units, '-', nc.variables[i].shape)

    return


# ==============================================================================


def read_era5_var(file_path, variable):
    """
    ----------------------------------------------------------------------------
    Extract the data for a particular variable from an ERA5 netCDF file
    ----------------------------------------------------------------------------
    :param file_path: full path to EBAF file
    :param variable: variable to read from file
    ----------------------------------------------------------------------------
    :return: variable array
    """

    from netCDF4 import Dataset

    nc = Dataset(file_path, 'r')

    var = nc.variables[variable]
    var_units = var.units
    var_name = var.standard_name

    print('variable shape...')
    print(var.shape, '\n')

    return var, var_name, var_units


# ==============================================================================


def read_era5_geolocation(file_path):
    """
    ----------------------------------------------------------------------------
    This function reads geolocation information from ERA5 files and then
    returns a 1-D vector and 2-D meshgrid of latitudes and longitudes
    ----------------------------------------------------------------------------
    :param file_path: path to EBAF file
    ----------------------------------------------------------------------------
    :return: (1) latitude  - vector [nlat]
             (2) longitude - vector [nlon]
             (3) lat - matrix [nlat x nlon]
             (4) lon - matrix [nlat x nlon]
    """

    import numpy as np
    from netCDF4 import Dataset

    nc = Dataset(file_path, 'r')

    latitude = nc.variables['latitude'][:]
    longitude = nc.variables['longitude'][:]

    time = nc.variables['time'][:]
    starting_month = str(time[16:18])

    # lat, lon grid
    lon, lat = np.meshgrid(longitude, latitude)

    print('lat array shape...')
    print(lat.shape, '\n')
    print('lon array shape...')
    print(lon.shape, '\n')

    return latitude, longitude, lat, lon, time


# ========================================================================


def polar_map(lon, lat, field,
               varname, varunits,
               nrows, ncols, cen_lon,
               cmap, cmap_lims, ti_str):
    """
    ----------------------------------------------------------------------------
    This function plots a map of gridded data
    ----------------------------------------------------------------------------
    :param lon: longitude matrix
    :param lat: latitude matrix
    :param field: data matrix
    :param varname: variable name
    :param varunits: variable units
    :param nrows: number of axis rows
    :param ncols: number of axis cols
    :param cen_lon: map central longitude
    :param cmap: colormap
    :param cmap_lims: colormap limits
    :param ti_str: title string
    ----------------------------------------------------------------------------
    :return: plots the map
    """

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature
    from mpl_toolkits.axes_grid1 import AxesGrid
    from cartopy.mpl.geoaxes import GeoAxes

    fig = plt.figure(figsize=(8, 6.5))

    # Map projection
    projection = ccrs.SouthPolarStereo(central_longitude=cen_lon)

    # Axis class
    axes_class = (GeoAxes, dict(map_projection=projection))

    # Plot figure
    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(nrows, ncols),
                    axes_pad=(0.4,0.4),
                    share_all = True,
                    cbar_location='right',
                    cbar_mode='each',
                    cbar_pad=0.1,
                    cbar_size='5%',
                    label_mode=1)

    for i, ax in enumerate(axgr):
        ax.add_feature(cartopy.feature.GSHHSFeature('low', edgecolor='grey'), zorder=1, facecolor='white')
        ax.gridlines(color='grey', linestyle='-')
        ax.set_title(ti_str)
        ax.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())
        ax.text(0.5, -0.15, varname + ' \n' + varunits, va='bottom', ha='center',
                rotation='horizontal', rotation_mode='anchor', transform=ax.transAxes, fontsize=10)

    # Use a different colorbar range every time by specifying tuple of tuples
    limits = cmap_lims
    for i in range(1):
        im = axgr[i].pcolor(lon, lat, field, transform=ccrs.PlateCarree(), cmap=cmap,
                            vmin=limits[0], vmax=limits[1])
        axgr.cbar_axes[i].colorbar(im)

    for i, cax in enumerate(axgr.cbar_axes):
        cax.set_yticks((limits[0], limits[1]))

    plt.show()
    return

