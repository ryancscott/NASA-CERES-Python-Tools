# ==============================================================================
#
#                   ----***** NASA CERES PYTHON TOOLS *****-----
#
# ==============================================================================
#
# Module: cerestools.py
#
# Purpose:  This library contains Python 3 functions to read, manipulate, and
#           analyze data from the National Aeronautics and Space Administration
#           (NASA) Clouds and the Earth's Radiant Energy System (CERES)
#           Earth Radiation Budget (ERB) satellite mission, including both
#           footprint-level swath and gridded time-interpolated and spatially
#           averaged (TISA) fields. Functions are provided for both data
#           development purposes (*_dev) as well as the analysis of official
#           release data products. See function descriptions below for more info.
#
# To use: import cerestools as ceres
#
# Dependencies: numpy, matplotlib, netcdf4, pyhdf, cartopy, datetime, palettable
#               I recommend installing the above libraries using conda, pip, or
#               an IDE (e.g., PyCharm).
#
# Author: Ryan C. Scott, ryan.c.scott@nasa.gov
#
# Last Updated: March 18, 2020
#
# ==============================================================================
# FOOTPRINT-LEVEL SWATH
# ---------------------
# get date_dev                <- get the year, month, day, hour from input file
# get platform                <- get the satellite and instrument name
# read_ssf_geolocation        <- get footprint lat, lon, etc. from SSF file
# read_crs_geolocation        <- get footprint lat, lon, etc. from CRS file
# read_ssf_var                <- get variable from SSF file
# read_crs_var                <- get variable from CRS file
# read_crs_dev_var            <- get variable from CRS development file
# swath_difference            <- compute difference between swaths
# set_colormap                <- set colormap from palettable library
# plot_swath                  <- plot SSF, CRS swath
# swath_histogram_scatterplot <- produce histogram & scatter plot of swath diff
# ---------------------
#  TISA GRIDDED FIELDS
# ---------------------
# print_nc_file_info          <- print info about variables in netCDF file
# read_ebaf_geolocation       <- read lat, lon, etc. from EBAF file
# read_ebaf_var               <- read field variable from EBAF file
# compute_monthly_anomalies   <- compute interannual monthly anomalies
# cos_lat_weight              <- compute matrix of cos(lat) weights
# compute_annual_climatology  <- compute long-term mean and std dev (sigma)
# compute_regional_averages   <- compute area-weighted mean for various regions
# composite_difference        <- compute composite mean difference
# global_mean_time_series     <- compute global area-weighted mean time series
# simple_regression           <- regress field on another field
# multiple_regression         <- regress field on multiple other fields
# global_map                  <- plot map of gridded field
# plot_time_series            <- plot time series of field
#
# -------------------
# Under development :
# -------------------
# compute_linear_trend               <- fit linear trend
# ==============================================================================


def get_date(file):
    """
    ----------------------------------------------------------------------------
    This function reads input file yyyymmddhh information and transforms
    the date string into a readable form for plotting. It also outputs date
    information for ingestion by the cartopy NightShade function. Works with
    official and development data files that end with date info. See line below.
    ----------------------------------------------------------------------------
    :param file: e.g., CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.2019010100
    ----------------------------------------------------------------------------
    :return: (1) date for Nightshade
             (2) date_str for plotting [string]
    """
    import datetime

    time_str = file[-10:]
    print(time_str)
    yyyy = time_str[0:4]
    mm = time_str[4:6]
    dd = time_str[6:8]
    hr = time_str[8:10]

    date = datetime.datetime(int(yyyy), int(mm), int(dd), int(hr))

    date_str = mm + '/' + dd + '/' + yyyy + ':' + hr + 'h'

    return date, date_str


# ========================================================================


def get_platform(file):
    """
    ----------------------------------------------------------------------------
    This function retrieves the satellite and flight model info from officially
    released CERES SSF or CRS files.
    ----------------------------------------------------------------------------
    :param file: e.g., CER_CRS_Terra-FM1-MODIS_Edition2G_023034.2010062023
    ----------------------------------------------------------------------------
    :return: (1) platform = satellite and flight model [string]
    """
    terra_aqua = file[8]

    if terra_aqua == "T":
        satellite = "Terra"
    elif terra_aqua == "A":
        satellite = "Aqua"
    else:
        satellite = "ERROR"

    # Get flight model (FM) info
    if satellite == "Terra":
        flight_model = file[14:17]
    elif satellite == "Aqua":
        flight_model = file[13:16]

    platform = satellite + '-' + flight_model

    print("CERES Instrument:", platform)
    return platform


# =====================================================================


def read_ssf_geolocation(file_path):
    """
    ----------------------------------------------------------------------------
    This function reads footprint-level geolocation information from
    CERES Level 2 SSF file
    ----------------------------------------------------------------------------
    :param file_path: path to file
    ----------------------------------------------------------------------------
    :return: (1) fov_lat = FOV latitude          [float]
             (2) fov_lon = FOV longitude         [float]
             (3) sza = FOV SZA at surface        [float]
             (4) obs_time = FOV observation time [float]
    """

    from pyhdf import SD
    hdf = SD.SD(file_path)
    # print(hdf.datasets())

    colatitude = hdf.select('Colatitude of CERES FOV at surface')
    colat = colatitude.get()
    fov_lat = 90 - colat

    longitude = hdf.select('Longitude of CERES FOV at surface')
    fov_lon = longitude.get()

    time_of_obs = hdf.select("Time of observation")
    obs_time = time_of_obs.get()

    solar_zenith = hdf.select("CERES solar zenith at surface")
    sza = solar_zenith.get()

    return fov_lat, fov_lon, obs_time, sza


# =====================================================================


def read_crs_geolocation(file_path):
    """
    ----------------------------------------------------------------------------
    This function reads footprint-level geolocation information from
    CERES Level 2 CRS file
    ----------------------------------------------------------------------------
    :param file_path: path to  file
    ----------------------------------------------------------------------------
    :return: (1) fov_lat = FOV latitude          [float]
             (2) fov_lon = FOV longitude         [float]
             (3) pres_levs = FOV pressure levels [float]
             (3) sza = FOV SZA at surface        [float]
             (4) obs_time = FOV observation time [float]
    """
    from pyhdf import SD
    hdf = SD.SD(file_path)
    # print(hdf.datasets())

    colatitude = hdf.select('Colatitude of CERES FOV at surface')
    colat = colatitude.get()
    fov_lat = 90 - colat

    longitude = hdf.select('Longitude of CERES FOV at surface')
    fov_lon = longitude.get()

    pressure_levels = hdf.select('Pressure levels')
    pres_levs = pressure_levels.get()

    time_of_obs = hdf.select("Time of observation")
    obs_time = time_of_obs.get()

    solar_zenith = hdf.select("CERES solar zenith at surface")
    sza = solar_zenith.get()

    return fov_lat, fov_lon, pres_levs, obs_time, sza


# ==============================================================================

# this function will be obsolete once I rename variables in crs2hdf.f90
def read_crs_geolocation_dev(file_path):
    """
    ----------------------------------------------------------------------------
    This function reads footprint-level geolocation information from
    CERES Level 2 CRS file
    ----------------------------------------------------------------------------
    :param file_path: path to  file
    ----------------------------------------------------------------------------
    :return: (1) fov_lat = FOV latitude          [float]
             (2) fov_lon = FOV longitude         [float]
             (3) pres_levs = FOV pressure levels [float]
             (3) sza = FOV SZA at surface        [float]
             (4) obs_time = FOV observation time [float]
    """

    from pyhdf import SD
    hdf = SD.SD(file_path)
    # print(hdf.datasets())

    colatitude = hdf.select('Colatitude')
    colat = colatitude.get()
    fov_lat = 90 - colat

    longitude = hdf.select('Longitude')
    fov_lon = longitude.get()

    pressure_levels = hdf.select('Pressure')
    pres_levs = pressure_levels.get()

    time_of_obs = hdf.select("JULTIM")
    obs_time = time_of_obs.get()

    return fov_lat, fov_lon, pres_levs, obs_time


# ==============================================================================


def read_ssf_var(file_path, vararg, fill):
    """
    ----------------------------------------------------------------------------
    This function reads variables from CERES Level 2 Single Scanner Footprint
    officially released data
    ----------------------------------------------------------------------------
    :param file_path: path to file
    :param vararg:    variable argument
    :param fill:      fill NaN value option
    ----------------------------------------------------------------------------
    :return: (1) field = variable data   [float]
             (2) name = name of variable [string]
             (3) units = variable units  [string]
    """

    import numpy as np
    from pyhdf import SD
    hdf = SD.SD(file_path)

    # select variable using integer index: vararg
    switch = {
        0: 'CERES downward SW surface flux - Model A',
        1: 'CERES downward LW surface flux - Model A',
        2: 'CERES downward SW surface flux - Model B',
        3: 'CERES downward LW surface flux - Model B'
    }
    var_name = switch.get(vararg)

    print("Getting", switch.get(vararg, "N/A"))

    # select and get the variable
    data = hdf.select(var_name)
    variable = data.get()
    var_units = data.units

    if fill == 1:
        variable[variable == data._FillValue] = np.nan

    # get field at appropriate vertical level
    var_field = variable

    return var_field, var_name, var_units


# ==============================================================================


def read_crs_var(file_path, vararg, levarg, fill):
    """
    ----------------------------------------------------------------------------
    This function reads variables from CERES Level 2 Cloud Radiative Swath
    official-release data files
    ----------------------------------------------------------------------------
    :param file_path: path to file
    :param vararg: variable argument
    :param levarg: level argument
    :param fill: fill NaN value option
    ----------------------------------------------------------------------------
    :return: (1) var_field = desired data    [float]
             (2) var_name = name of variable [string]
             (3) var_units = variable units  [string]
             (4) lev_name = p level name     [string]
    """

    import numpy as np
    from pyhdf import SD
    hdf = SD.SD(file_path)

    # select variable using integer index: vararg
    switch = {
        0: 'Longwave flux - downward - total',
        1: 'Longwave flux - upward - total',
        2: 'Shortwave flux - downward - total',
        3: 'Shortwave flux - upward - total'
    }
    var_name = switch.get(vararg)

    # select p level using integer index: levarg
    switch2 = {
        0: 'TOA',
        1: '70 mb',
        2: '200 mb',
        3: '500 mb',
        4: 'surface'
    }
    lev_name = switch2.get(levarg)

    print("Getting", switch.get(vararg, "N/A"), "at:",
          switch2.get(levarg, "N/A"))

    # select and get the variable
    data = hdf.select(var_name)
    variable = data.get()
    var_units = data.units
    var_fill = data._FillValue

    if fill == 1:
        variable[variable == var_fill] = np.nan

    # get field at appropriate vertical level
    var_field = variable[:, levarg]

    return var_field, var_name, var_units, lev_name


# ==============================================================================

# This function will become obsolete once I rename variables in crs2hdf.f90
def read_crs_var_dev(file_path, vararg, levarg, fill):
    """
    ----------------------------------------------------------------------------
    This function reads data from CERES Level 2 Cloud Radiative Swath
    DEVELOPMENT files - those I produce by running the CRS4 .f90 code on AMI
    ----------------------------------------------------------------------------
    :param file_path: path to data file [string]
    :param vararg: variable argument
    :param levarg: level argument
    :param fill: fill NaN value option
    ----------------------------------------------------------------------------
    :return: (1) var_field = variable field data  [float]
             (2) var_name = variable name         [string]
             (3) var_units = variable units       [string]
             (4) lev_name = p level name          [string]
    """

    import numpy as np
    from pyhdf import SD
    hdf = SD.SD(file_path)

    # select variable using integer index: vararg
    switch = {
        0: 'UT_CLR_SW_DN',
        1: 'UT_CLR_SW_UP',
        2: 'UT_CLR_LW_DN',
        3: 'UT_CLR_LW_UP',

        4: 'UT_TOT_SW_DN',
        5: 'UT_TOT_SW_UP',
        6: 'UT_TOT_LW_DN',
        7: 'UT_TOT_LW_UP',

        8: 'UT_PRS_SW_DN',
        9: 'UT_PRS_SW_UP',
        10: 'UT_PRS_LW_DN',
        11: 'UT_PRS_LW_UP',

        12: 'UT_TNA_SW_DN',
        13: 'UT_TNA_SW_UP',
        14: 'UT_TNA_LW_DN',
        15: 'UT_TNA_LW_UP',

        16: 'MATCH_AOT',
        17: 'SFC_HGT',
        18: 'SZEN_TOA'
    }
    var_name = switch.get(vararg)

    # select level
    switch2 = {
        -1: '-',
        0: 'TOA',
        1: 'surface'
    }
    lev_name = switch2.get(levarg)

    print("Getting", switch.get(vararg, "N/A"), "at:",
          switch2.get(levarg, "N/A"))

    # select the variable and get the data
    data = hdf.select(var_name)
    variable = data.get()
    var_units = data.units

    if fill == 1:
        variable[variable == data._FillValue] = np.nan

    # get field at appropriate vertical level
    if levarg < 0:
        var_field = variable
    elif levarg >= 0:
        var_field = variable[:, levarg]

    return var_field, var_name, var_units, lev_name


# ==============================================================================


def set_colormap(cmap_name, typarg):
    """
    ----------------------------------------------------------------------------
    Selects colormap from palettable library
    ----------------------------------------------------------------------------
    :param cmap_name: name of colormap (from palettable)
    :param typarg: 0 = continuous, 1 = discrete
    ----------------------------------------------------------------------------
    :return: colormap
    """

    switch = {
        0: "continuous",
        1: "discrete"
    }
    print("\nUsing", switch.get(typarg, "N/A"), "colormap")

    if typarg == 0:
        color_map = cmap_name.mpl_colormap
    elif typarg == 1:
        from matplotlib.colors import ListedColormap
        color_map = ListedColormap(cmap_name.mpl_colors)
    elif typarg != 1 or typarg != 0:
        print('Please select 0 or 1')

    return color_map


# ==============================================================================


def plot_swath(lon, lat, field,
               varname, levname, varunits,
               nrows, ncols, cen_lon,
               cmap, cmap_lims, date, nightshade,
               date_str, title_str):
    """
    ----------------------------------------------------------------------------
    This function plots a swath of footprint-level data
    FLASHFlux, SSF, or CRS
    ----------------------------------------------------------------------------
    :param lon: FOV longitude
    :param lat: FOV latitude
    :param field: variable
    :param varname: variable name
    :param levname: level name
    :param varunits: variable units
    :param nrows: number of rows
    :param ncols: number of columns
    :param cen_lon: central longitude
    :param cmap: colormap
    :param cmap_lims: colormap limits
    :param date: date info for nightshade
    :param nightshade: whether to use nighshade feature
    :param date_str: date string
    :param title_str: title string
    ----------------------------------------------------------------------------
    :return: plot of the data
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature
    from mpl_toolkits.axes_grid1 import AxesGrid
    from cartopy.mpl.geoaxes import GeoAxes
    from cartopy.feature.nightshade import Nightshade

    # Map projection
    projection = ccrs.PlateCarree(central_longitude=cen_lon)

    # Axis class
    axes_class = (GeoAxes, dict(map_projection=projection))

    # Create figure
    fig = plt.figure(figsize=(12, 6))
    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(nrows, ncols),
                    axes_pad=(0.4, 0.4),
                    share_all=True,
                    cbar_location='right',
                    cbar_mode='each',
                    cbar_pad=0.1,
                    cbar_size='5%',
                    label_mode=1)

    # Loop over axes
    for i, ax in enumerate(axgr):
        ax.add_feature(cartopy.feature.LAND, zorder=1, facecolor='none',
                       edgecolor='darkgrey')
        ax.gridlines(color='grey', linestyle='--')
        ax.set_title(title_str + ' - ' + date_str, fontsize=10)
        ax.set_extent([-180, 180, -90, 90], projection)
        ax.text(0.5, -0.1, varname + '  -  ' + levname + '  -  ' + varunits,
                va='bottom', ha='center',
                rotation='horizontal', rotation_mode='anchor',
                transform=ax.transAxes, fontsize=10)

        if nightshade == 1:
            ax.add_feature(Nightshade(date, alpha=0.15))

    # To use a different colorbar range every time, can use a tuple of tuples
    limits = cmap_lims
    for i in range(nrows * ncols):
        im = axgr[i].scatter(lon, lat, c=field, s=1,
                             transform=ccrs.PlateCarree(),
                             vmin=limits[0], vmax=limits[1], cmap=cmap)
        axgr.cbar_axes[i].colorbar(im)

    for i, cax in enumerate(axgr.cbar_axes):
        cax.set_yticks(np.linspace(limits[0], limits[1], 5))
        cax.set_yticklabels(np.linspace(limits[0], limits[1], 5),
                            fontsize=8)

    plt.show()
    return


# ==============================================================================


def swath_difference(field2, field1, day_only, sza):
    """
    ----------------------------------------------------------------------------

    ----------------------------------------------------------------------------
    :param field2:
    :param field1:
    :param day_only:
    :param sza:
    ----------------------------------------------------------------------------
    :return:
    """

    import numpy as np

    diff = field2 - field1

    diff[diff == max(diff)] = np.nan
    diff[diff == min(diff)] = np.nan

    # need to num FOVs from file...
    if day_only is True:
        for i in range(field2.shape[0]):
            if sza[i] > 90:
                diff[i] = np.nan

    return diff


# ==============================================================================


def swath_histogram_scatterplot(field1, field2, varname, levname, time_info, platform, day_only, sza):
    """
    ----------------------------------------------------------------------------

    ----------------------------------------------------------------------------
    :param field1:
    :param field2:
    :param varname:
    :param levname:
    :param time_info:
    :param platform:
    :param day_only:
    :param sza:
    ----------------------------------------------------------------------------
    :return:
    """

    import numpy as np
    import matplotlib.pyplot as plt

    swath_diff = swath_difference(field2=field2, field1=field1, day_only=day_only, sza=sza)

    mean_diff = np.nanmean(swath_diff)
    sigma_diff = np.nanstd(swath_diff)

    # number of FOVs compared
    print("Total number of FOVs in swath: ", len(swath_diff))
    print("Total number of FOVs where at least one is NaN: ", np.sum(np.isnan(swath_diff)))
    N = len(swath_diff) - np.sum(np.isnan(swath_diff))
    print("Total number of FOVs compared: ", N)

    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    fig.suptitle('CERES ' + platform + ' - Cloud Radiative Swath (CRS)' + ' - ' + time_info, fontsize=11)

    # histogram of difference
    axs[0].hist(swath_diff, bins=200, align='mid', rwidth=1)
    axs[0].grid()
    axs[0].set_axisbelow("True")
    axs[0].set_title("CRS Ed4 minus CRS Ed2G", fontsize=11)
    axs[0].set_xlabel(varname + ' - ' + levname + "\n" + r"Flux difference ($\Delta$)")
    axs[0].set_xlim([-125, 125])
    axs[0].set_ylabel("Number of CERES Footprints")

    # show descriptive statistics
    textstr = "N = " + str(N) + "\n" + r"Mean: $\bar{\Delta}$ = " + str(np.around(mean_diff, 2)) + "\n" + \
              r"Stdv: $\sigma_{\Delta}$ = " + str(np.around(sigma_diff, 2))
    props = dict(facecolor='white', alpha=0.85)
    axs[0].text(0.05, 0.95, textstr, transform=axs[0].transAxes, fontsize=8,
                verticalalignment='top', bbox=props)

    if np.nanmax(field1) < np.nanmax(field2):
        max = np.nanmax(field2) + 100
    else:
        max = np.nanmax(field1) + 100

    # scatter plot
    axs[1].plot(range(2000), range(2000), color='black', zorder=0)
    axs[1].scatter(field1, field2, s=0.05)
    axs[1].set(xlim=(0, max), ylim=(0, max))
    axs[1].grid()
    axs[1].set_axisbelow("True")
    axs[1].set_title(varname + ' - ' + levname, fontsize=11)
    axs[1].set_xlabel("CRS Edition 2G", fontsize=11)
    axs[1].set_ylabel("CRS Edition 4", fontsize=11)

    plt.show()
    return


# ========================================================================


def print_nc_file_info(filepath):
    """
    ----------------------------------------------------------------------------

    ----------------------------------------------------------------------------
    :param filepath:
    ----------------------------------------------------------------------------
    :return:
    """
    print('====================================')
    print('\t\tVariables in file...')
    print('====================================')
    print('\tname - units - dimensions')
    print('====================================')

    from netCDF4 import Dataset

    nc = Dataset(filepath, 'r')
    # print information about variables in netCDF file
    for i in nc.variables:
        print(i, '-', nc.variables[i].units, '-', nc.variables[i].shape)

    return


# ========================================================================


def read_ebaf_geolocation(filepath):
    """
    ----------------------------------------------------------------------------
    This function reads geolocation information from CERES EBAF files
    ----------------------------------------------------------------------------
    :param filepath: path to EBAF file
    ----------------------------------------------------------------------------
    :return: (1) latitude  - array [nlat]
             (2) longitude - array [nlon]
             (3) lat - matrix [nlat x nlon]
             (4) lon - matrix [nlat x nlon]
    """

    import numpy as np
    from netCDF4 import Dataset

    nc = Dataset(filepath, 'r')

    latitude = nc.variables['lat'][:]
    longitude = nc.variables['lon'][:]

    time = nc.variables['time'].units
    starting_month = str(time[16:18])

    # lat, lon grid
    lon, lat = np.meshgrid(longitude, latitude)

    print('lat array shape...')
    print(lat.shape, '\n')
    print('lon array shape...')
    print(lon.shape, '\n')

    return latitude, longitude, lat, lon


# ========================================================================


def read_ebaf_start_month_year(filepath):
    """
    ----------------------------------------------------------------------------
    This function extracts the starting month and year from EBAF file
    ----------------------------------------------------------------------------
    :param filepath: path to EBAF file
    ----------------------------------------------------------------------------
    :return: (1) month
             (2) year
    """

    import numpy as np
    from netCDF4 import Dataset

    nc = Dataset(filepath, 'r')

    time = nc.variables['time'][:]
    time_units = nc.variables['time'].units
    start_month = str(int(time_units[16:18]))
    start_year = str(time_units[11:15])

    time_series_length = len(time)

    # print(time_series_length)
    # print(start_month)
    # print(start_year)

    return start_month, start_year


# ========================================================================


def read_ebaf_var(filepath, variable):
    """
    ----------------------------------------------------------------------------
    Extract variable data from EBAF file
    ----------------------------------------------------------------------------
    :param filepath: input EBAF file
    :param var: variable to extract from file
    ----------------------------------------------------------------------------
    :return: array of variable var
    """

    from netCDF4 import Dataset

    nc = Dataset(filepath, 'r')

    var = nc.variables[variable]
    var_units = var.units
    var_name = var.standard_name

    print('variable shape...')
    print(var.shape, '\n')

    return var, var_name, var_units


# ========================================================================


def compute_monthly_anomalies(field, fieldname):
    """
    ----------------------------------------------------------------------------
    (1) Compute long-term average for each calendar month, i.e., the monthly mean seasonal cycle
    (2) Subtract the appropriate long-term means from the corresponding monthly fields to get monthly anomalies
    ----------------------------------------------------------------------------
    :param field: variable for which to compute the climatology/anomaly [lat,lon,time]
    ----------------------------------------------------------------------------
    :return: (1) time series of monthly anomalies at each grid box
             (2) monthly mean seasonal cycle at each grid box
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # total length of record in file
    record_len = field.shape[0]

    # initialize array to hold monthly climatology for each month
    monthly_climatology = np.empty((12, 180, 360))

    print('Indices for monthly climatology:\n')
    # loop over each calendar month and average
    for i in range(0, 12):

        # select appropriate indices associated with each month
        indices = list(range(i, record_len, 12))

        print(indices, "Total number of indices used:", len(indices))

        # average using appropriate indices
        monthly_climatology[i, :, :] = np.nanmean(field[indices, :, :], axis=0)

    print('\nComputing # of times to repeat seasonal cycle:\n')
    # compute number of times to repeat seasonal cycle
    i = 0
    num_yr = 0
    while i < record_len:
        i += 12
        num_yr += 1
        print(i, num_yr)

    # repeat 20 times
    seasonal_cycle = np.tile(monthly_climatology, (num_yr, 1, 1))

    # only take data up to the last index
    seasonal_cycle = seasonal_cycle[0:record_len, :, :]

    # calculate anomalies
    monthly_anomalies = np.subtract(field, seasonal_cycle)

    r = 100
    p = 120
    # add lat lon location info to plots...
    fig, (ax1, ax2) = plt.subplots(2, 1)
    # first set of axes
    ax1.plot(np.linspace(0, 0, record_len))
    ax1.plot(monthly_anomalies[:, r, p], label='Anomalies')
    ax1.legend(fontsize=10)
    ax1.set_title(fieldname)
    # second set of axes
    ax2.plot(seasonal_cycle[:, r, p], label='Seasonal Cycle')
    ax2.plot(field[:, r, p], label='Raw Field')
    ax2.legend(fontsize=10)
    plt.show()

    return monthly_anomalies, seasonal_cycle


# ========================================================================


def cos_lat_weight(lat_vector):
    """
    ----------------------------------------------------------------------------
    This function takes a latitude array and computes a matrix of cos(lat)
    weights for performing area-weighted average values of a field
    ----------------------------------------------------------------------------
    :param lat_vector: vector of latitudes
    ----------------------------------------------------------------------------
    :return: matrix of cos(lat) weights
    """

    import numpy as np
    import numpy.matlib
    import matplotlib.pyplot as plt

    # compute cos(lat in degrees to radians)
    weight = np.cos(np.deg2rad(lat_vector))
    print('w array shape...')
    print(weight.shape, '\n')

    # replicate cos(lat) vector 360 times along long dim
    weight = np.matlib.repmat(weight, 360, 1)

    # transpose to be consistent with other fields
    weight = np.transpose(weight)

    print('cos(lat) weight array shape...')
    print(weight.shape, '\n')

    # plt.imshow(weight)
    # plt.title("Cosine of Latitude Weight Matrix")

    return weight


# ========================================================================


def compute_annual_climatology(field):
    """
    ----------------------------------------------------------------------------
    This function calculates the long-term mean and standard deviation of the
    input field
    ----------------------------------------------------------------------------
    :param field: variable under consideration [ntime x nlat x nlon]
    ----------------------------------------------------------------------------
    :return: long-term mean and standard deviation
    """

    import numpy as np

    # average over entire time series
    mean_field = np.nanmean(field, axis=0)

    # standard deviation over entire time series
    standard_deviation_field = np.nanstd(field, axis=0)

    return mean_field, standard_deviation_field


# ========================================================================


def compute_regional_averages(field, latitude, weights):
    """
    ----------------------------------------------------------------------------
    Computes regional area-weighted averages of input field
    ----------------------------------------------------------------------------
    :param field: field to be averaged
    ----------------------------------------------------------------------------
    :return: prints regional averages
    """
    import numpy as np

    # Global mean climatological total cloud fraction
    global_avg = np.average(field, weights=weights)
    print('\n', 'Global average:', global_avg, '\n')

    # 60Sto90S mean climatological total cloud fraction
    sh60to90_avg = np.average(field[0:30, :], weights=weights[0:30, :])
    print('Antarctic (60S-90S) average:', sh60to90_avg)
    print('Latitude range:', latitude[0:30], '\n')

    # SH mid-latitude climatological total cloud fraction
    sh30to60_avg = np.average(field[30:60, :], weights=weights[30:60, :])
    print('S mid-latitude (30S-60S) average:', sh30to60_avg)
    print('Latitude range:', latitude[30:60], '\n')

    # SH tropical mean climatology
    sh0to30_avg = np.average(field[60:90, :], weights=weights[60:90, :])
    print('S tropical (30S-0) average:', sh0to30_avg)
    print('Latitude range:', latitude[60:90], '\n')

    # NH tropical mean climatological total cloud fraction
    nh0to30_avg = np.average(field[90:120, :], weights=weights[90:120, :])
    print('N tropical (0-30N) average:', nh0to30_avg)
    print('Latitude range:', latitude[90:120], '\n')

    # NH tropical mean climatological total cloud fraction
    nh30to60_avg = np.average(field[120:150, :], weights=weights[120:150, :])
    print('N mid-latitude (30N-60N) average:', nh30to60_avg)
    print('Latitude range:', latitude[120:150], '\n')

    # Arctic mean climatological total cloud fraction
    nh60to90_avg = np.average(field[150:180, :], weights=weights[150:180, :])
    print('Arctic 60N-90N average:', nh60to90_avg)
    print('Latitude range:', latitude[150:180], '\n')

    return sh0to30_avg, sh30to60_avg, sh60to90_avg, nh0to30_avg, nh30to60_avg, nh60to90_avg


# ========================================================================


def composite_difference(field, ind1, ind2):
    """
    ----------------------------------------------------------------------------
    This function calculates the composite mean of 'field' for two different
    time period and then takes their difference to illustrate the field change.
    ----------------------------------------------------------------------------
    :param field: field to average
    :param ind1: indices for first time period
    :param ind2: indices for second time period
    ----------------------------------------------------------------------------
    :return: difference between composites
    """
    import numpy as np
    mean_field1 = np.nanmean(field[ind1, :, :], axis=0)
    print('Size of first section to average: ', field[ind1, :, :].shape)

    mean_field2 = np.nanmean(field[ind2, :, :], axis=0)
    print('Size of second section to average: ', field[ind2, :, :].shape)

    composite_diff = np.subtract(mean_field2, mean_field1)

    return composite_diff


# ========================================================================


def global_mean_time_series(field, weight):
    """
    ----------------------------------------------------------------------------
    This function calculates the global area-cos(lat)-weighted mean of the
    input field
    ----------------------------------------------------------------------------
    :param field: input field                [ntime x nlat x nlon]
    :param weight: cos(lat) weight matrix    [nlat x nlon]
    ----------------------------------------------------------------------------
    :return: mean_time_series
    """
    import matplotlib.pyplot as plt
    import numpy as np

    mean_time_series = np.empty(shape=field.shape[0])
    for i in range(field.shape[0]):
        field1 = field[i, :, :]
        mean_time_series[i] = np.average(field1, weights=weight)
        #print(mean_time_series[i])

    return mean_time_series


# ========================================================================


def plot_time_series(var, name, units, start_mo, start_yr):
    """
    ----------------------------------------------------------------------------
    This function plots a time series of arbitrary length, such as one
    produced by the global_mean_time_series function above
    ----------------------------------------------------------------------------
    :param var: time series to plot             [float; ntime]
    :param name: name of the variable           [string]
    :param units: units of the variable         [string]
    :param start_mo: 1st month of time series   [int]
    :param start_yr: 1st year of time series    [int]
    ----------------------------------------------------------------------------
    :return: plot the data
    """
    import matplotlib.pyplot as plt

    record_len = var.shape[0]

    # plot time series
    plt.figure(figsize=(12, 4))
    plt.grid()
    plt.plot(var, label=name, marker='.')
    plt.title(name)
    xticklabels = [str(start_mo) + '/' + str((start_yr+i))[2:] for i in range(record_len)]
    plt.xticks(ticks=range(0, record_len, 12), labels=xticklabels, fontsize=10)
    plt.yticks(fontsize=8)
    plt.ylabel(units)
    plt.legend(fontsize=8)
    plt.show()

    return


# ========================================================================


def simple_regression(y_anomalies, x_anomalies):
    """
    ----------------------------------------------------------------------------
    This function regresses a field of y anomalies onto a field of x anomalies
    ----------------------------------------------------------------------------
    :param x_anomalies: independent variable [ntime x nlat x nlon]
    :param y_anomalies: dependent variable   [ntime x nlat x nlon]
    ----------------------------------------------------------------------------
    :return: matrix of regression coefficients
    """

    import numpy as np

    coefficients = np.zeros([np.shape(x_anomalies)[1], np.shape(x_anomalies)[2]])
    intercept = np.zeros([np.shape(x_anomalies)[1], np.shape(x_anomalies)[2]])

    for i in range(coefficients.shape[0]):
        for j in range(coefficients.shape[1]):
            x = x_anomalies[:, i, j]
            y = y_anomalies[:, i, j]
            a = np.vstack([x, np.ones(len(x))]).T
            coefficients[i, j], intercept[i, j] = np.linalg.lstsq(a, y, rcond=None)[0]

    return coefficients


# ========================================================================


def multiple_regression_2xi(y_anomalies, x1_anomalies, x2_anomalies):
    """
    ----------------------------------------------------------------------------
    This function performs multiple-linear regression analysis of the
    y_anomaly field (dependent variable) onto x1_anomaly and x2_anomaly fields
    (independent variables) at each grid box. It returns regression coefficient
    (beta_i, i = 1,2) maps for x1 and x2. The coefficients represent the partial
    derivative (dy/dx_i) response of y to independent changes in x_i, with the
    other x_i held fixed.
    ----------------------------------------------------------------------------
    :param y_anomalies:   y anomaly field (  dependent variable)
    :param x1_anomalies: x1 anomaly field (independent variable)
    :param x2_anomalies: x2 anomaly field (independent variable)
    :return: matrix of coefficients
    """

    import numpy as np

    coefficients = np.zeros([3, np.shape(x1_anomalies)[1], np.shape(x1_anomalies)[2]])

    for i in range(coefficients.shape[1]):
        for j in range(coefficients.shape[2]):
            x1 = x1_anomalies[:, i, j]
            x2 = x2_anomalies[:, i, j]
            y = y_anomalies[:, i, j]
            a = np.vstack([x1, x2, np.ones(len(x1))]).T
            coefficients[:, i, j] = np.linalg.lstsq(a, y, rcond=None)[0]

    return coefficients

# ========================================================================


def multiple_regression(y_anomalies, *x_anomalies):
    """
    ----------------------------------------------------------------------------
    This function performs multiple-linear regression analysis of y_anomalies
    (the dependent variable) on one or more x_anomalies (independent variables)
    at each grid box. It returns regression coefficient (beta_i) maps for each
    variable x_i. The beta_i represent the partial derivative (dy/dx_i) response
    of y to independent changes in x_i, with all of the other x_i held fixed.
    The star before x_anomalies means that any number of predictors may be input
    to the function.
    ----------------------------------------------------------------------------
    :param y_anomalies:  dependent variable anomaly field [nlat x nlon]
    :param x_anomalies: independent variable anomaly field [*arg x nlat x nlon]
    ----------------------------------------------------------------------------
    :return: coefficient matrix [number of x fields input + 1 x nlat x nlon]
    """

    import numpy as np

    # first dimension of coefficient matrix based on
    # the number of "x_anomalies" variables included
    c = len(x_anomalies) + 1

    # initialize coefficient matrix, size: c x 180 x 360
    coefficients = np.zeros([c, np.shape(x_anomalies[1])[1], np.shape(x_anomalies[1])[2]])

    # loop over each grid box
    # coefficients.shape[0] = c
    for i in range(coefficients.shape[1]):
        for j in range(coefficients.shape[2]):

            # empty array for predictor variable time series
            # size: number of independent variables x length of time series
            x = np.empty([len(x_anomalies), len(x_anomalies[0])])
            # print(len(x_anomalies[0]))

            # add row for each predictor variable
            for k, x_anom in enumerate(x_anomalies):
                x[k, :] = x_anom[:, i, j]
                # if standardized == True:
                #    need to add option

            y = y_anomalies[:, i, j]
            a = np.vstack([x, np.ones(len(x[0, :]))]).T

            # print the input anomalies at the first grid box
            if i == 1 and j == 1:
                print(a)
                print('X matrix shape:', a.shape[0])

            # compute coefficients
            coefficients[:, i, j] = np.linalg.lstsq(a, y, rcond=None)[0]

    return coefficients


# ========================================================================


def global_map(lon, lat, field,
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
    projection = ccrs.PlateCarree(central_longitude=cen_lon)

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
        ax.add_feature(cartopy.feature.LAND, zorder=1, facecolor='none', edgecolor='grey')
        # ax.add_feature(cartopy.feature.GSHHSFeature('low',edgecolor='none'), zorder=1, facecolor='black')
        ax.gridlines(color='black', linestyle=':')
        ax.set_title(ti_str)
        ax.set_extent([-180, 180, -90, 90], projection)
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


# ========================================================================

