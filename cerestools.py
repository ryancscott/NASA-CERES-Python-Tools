# ==============================================================================
#
#                  -----***** NASA CERES PYTHON TOOLS *****-----
#
# ==============================================================================
#
# Module:   cerestools.py
#
# Purpose:  This library contains Python 3 functions to read, process, analyze
#           and visualize data from NASA's Clouds and the Earth's Radiant Energy
#           System (CERES) satellite mission. Functions are included to work
#           with Level 2 swath data & Level 3 time-interpolated and spatially-
#           averaged (TISA) gridded fields. This library is under development
#           for data product development purposes (*_dev) and analysis of
#           official release CERES data products. See function descriptions
#           below for more information.
#
# To Use:   import cerestools as ceres
#
# Requires: numpy, scipy, matplotlib, netcdf4, pyhdf, cartopy, datetime, &
#           palettable. I recommend installing these libraries using conda,
#           pip, or an IDE
#
# Author:   Ryan C. Scott, SSAI
#           ryan.c.scott@nasa.gov
#
# Last Updated: June 1, 2020
#
# This code conforms to Python PEP-8 coding standards.
# ==============================================================================
#       LEVEL 2
# =====================
# get_date                    <- get the year, month, day, hour from input file
# get_platform                <- get the satellite and instrument name
# get_platform_dev            <- same as above but for dev files
# read_ssf_geolocation        <- get footprint lat, lon, etc. from SSF file
# read_crs_geolocation        <- get footprint lat, lon, etc. from CRS file
# read_crs_geolocation_dev    <- same as above but for dev files
# read_ssf_var                <- get variable from SSF file
# read_crs_var                <- get variable from CRS file
# read_crs_var_dev            <- get variable from CRS development file
# read_day_of_crs_files       <- get 24 hr of data from CRS file
# read_day_of_ssf_files       <- get 24 hr of data from SSF file
# read_month_of_ssf_files     <- get month of data from SSF file
# read_month_of_crs_files     <- get month of data from CRS file
# swath_lon_360_to_180        <- convert FOV lon range 0 to 360 to -180 to 180
# swath_daytime_only          <- extract only daytime FOVs
# swath_difference            <- compute difference between swaths
# set_colormap                <- set colormap from palettable library
# plot_swath                  <- plot SSF, CRS swath
# swath_histogram_scatterplot <- produce histogram of swath diff & scatter plot
# grid_to_1x1_degree          <- bin FOVs in 1x1 grid boxes & compute statistics
# plot_gridded_fields         <- plot gridded fields, will replace global_map
# =====================
#       LEVEL 3
# =====================
# read_syn1deg_hdf            <- read variable from SYN1deg HDF file
#
# print_nc_file_info          <- print info about data in netCDF file
# read_ebaf_geolocation       <- read lat, lon, etc. from EBAF file
# read_ebaf_var               <- read field variable from EBAF file
# compute_monthly_anomalies   <- compute interannual monthly anomalies
# cos_lat_weight              <- compute matrix of cos(lat) weights
# compute_annual_climatology  <- compute long-term mean and std dev (sigma)
# compute_regional_averages   <- compute area-weighted mean for various regions
# composite_difference        <- compute composite difference of field
# global_mean_time_series     <- compute global area-weighted mean time series
# simple_regression           <- regress a field onto another field
# multiple_regression         <- regress field on multiple other fields
# global_map                  <- plot map of gridded field
# plot_time_series            <- plot time series of field
# =====================
# Under development :
# =====================
# compute_linear_trend        <- calculate linear trend
# ==============================================================================


def get_date(file):
    """
    ----------------------------------------------------------------------------
    This function reads input file "yyyymmddhh" information and transforms
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


# ==============================================================================


def get_platform(file):
    """
    ----------------------------------------------------------------------------
    This function retrieves the satellite and flight model information from
    official release CERES SSF or CRS HDF files.
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
        satellite = "Error"

    # Get flight model (FM) info
    if satellite == "Terra":
        flight_model = file[14:17]
    elif satellite == "Aqua":
        flight_model = file[13:16]

    platform = satellite + '-' + flight_model

    print("CERES Instrument:", platform)
    return platform


# ==============================================================================


def get_platform_dev(file):
    """
    ----------------------------------------------------------------------------
    This function retrieves the satellite and flight model info from CERES
    CRS *DEVELOPMENT* HDF files...
    ----------------------------------------------------------------------------
    :param file: e.g., CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.2019011212
    ----------------------------------------------------------------------------
    :return: (1) platform = satellite and flight model [string]
    """
    terra_aqua = file[9]

    if terra_aqua == "T":
        satellite = "Terra"
    elif terra_aqua == "A":
        satellite = "Aqua"
    else:
        satellite = "ERROR"

    # Get flight model (FM) info
    if satellite == "Terra":
        flight_model = file[15:18]
        print(flight_model)
    elif satellite == "Aqua":
        flight_model = file[14:17]

    platform = satellite + '-' + flight_model

    print("CERES Instrument:", platform)
    return platform


# ==============================================================================


def read_ssf_geolocation(file_path):
    """
    ----------------------------------------------------------------------------
    This function reads footprint-level geolocation information from
    CERES Level 2 SSF HDF files
    ----------------------------------------------------------------------------
    :param file_path: path to SSF file
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


# ==============================================================================


def read_crs_geolocation(file_path):
    """
    ----------------------------------------------------------------------------
    This function reads footprint-level geolocation information from
    CERES Level 2 CRS HDF files
    ----------------------------------------------------------------------------
    :param file_path: path to CRS file
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


def read_crs_geolocation_dev(file_path):
    """
    ----------------------------------------------------------------------------
    This function reads footprint-level geolocation information from CERES
    Level 2 CRS HDF *DEVELOPMENT* files
    ----------------------------------------------------------------------------
    :param file_path: path to file
    ----------------------------------------------------------------------------
    :return: (1) fov_lat = FOV latitude          [float]
             (2) fov_lon = FOV longitude         [float]
             (3) pres_levs = FOV pressure levels [float]
             (4) obs_time = FOV observation time [float]
             (5) sfc_ind = surface level index   [float]
             (6) sza = FOV SZA at surface        [float]
    """

    from pyhdf import SD
    hdf = SD.SD(file_path)
    # print(hdf.datasets())

    colatitude = hdf.select('Colatitude of CERES FOV at surface')
    colat = colatitude.get()
    fov_lat = 90 - colat

    longitude = hdf.select('Longitude of CERES FOV at surface')
    fov_lon = longitude.get()

    pressure_levels = hdf.select('Pressure profile')  # or 'Pressure profile'
    pres_levs = pressure_levels.get()

    time_of_obs = hdf.select("Julian Time")
    obs_time = time_of_obs.get()

    surface_indices = hdf.select('Surface level index')
    sfc_ind = surface_indices.get()

    solar_zenith = hdf.select("CERES solar zenith at surface")
    sza = solar_zenith.get()

    return fov_lat, fov_lon, pres_levs, obs_time, sfc_ind, sza


# ==============================================================================


def read_ssf_var(file_path, var_name, index, fill):
    """
    ----------------------------------------------------------------------------
    This function reads data from CERES Level 2 Single Scanner Footprint
    (SSF) official release HDF files
    ----------------------------------------------------------------------------
    :param file_path: path to file               [string]
    :param var_name:  variable name              [string]
    :param index:     index of desired column    [integer]
    :param fill:      option to fill NaN values  [boolean]
    ----------------------------------------------------------------------------
    :return: (1) field = variable data           [float]
             (2) name = name of variable         [string]
             (3) units = variable units          [string]
    """

    import numpy as np
    from pyhdf import SD
    hdf = SD.SD(file_path)

    # select and get the variable
    data = hdf.select(var_name)
    variable = data.get()
    var_units = data.units

    # replace fill values with NaN
    if fill is True:
        variable[variable == data._FillValue] = np.nan

    var_field = np.empty([variable.shape[0]])
    # if the data being read is of size (numFOVs x 1), set index = -1
    # get data from appropriate column of HDF file
    if index < 0:
        var_field = variable
    elif index >= 0:
        var_field = variable[:, index]

    return var_field, var_name, var_units


# ==============================================================================


def read_crs_var(file_path, var_name, lev_arg, fill):
    """
    ----------------------------------------------------------------------------
    This function reads data from CERES Level 2 Cloud Radiative Swath
    (CRS) official release HDF files
    ----------------------------------------------------------------------------
    :param file_path: path to file           [string]
    :param var_name: variable name           [string]
    :param lev_arg: level argument           [integer]
    :param fill: option to fill NaN values   [boolean]
    ----------------------------------------------------------------------------
    :return: (1) var_field = desired data    [float]
             (2) var_name = name of variable [string]
             (3) var_units = variable units  [string]
             (4) lev_name = p level name     [string]
    """

    import numpy as np
    from pyhdf import SD
    hdf = SD.SD(file_path)

    # Select p level using integer index: lev_arg
    switch2 = {
        -1: ' ',    # not associated with a particular level
        0: 'TOA',
        1: '70 mb',
        2: '200 mb',
        3: '500 mb',
        4: 'surface'
    }
    # retrieve name of the level
    lev_name = switch2.get(lev_arg)

    print("Getting", var_name, "at",
          switch2.get(lev_arg, "N/A"))

    # select and get the variable
    data = hdf.select(var_name)
    variable = data.get()
    var_units = data.units
    var_fill = data._FillValue

    # option to replace fill values with NaN
    if fill is True:
        variable[variable == var_fill] = np.nan

    var_field = np.empty([variable.shape[0]])
    # get field at appropriate vertical level
    if lev_arg < 0:
        var_field = variable
    elif lev_arg >= 0:
        var_field = variable[:, lev_arg]

    return var_field, var_name, var_units, lev_name


# ==============================================================================


def read_crs_var_dev(file_path, var_name, lev_arg, fill):
    """
    ----------------------------------------------------------------------------
    This function reads data from CERES Level 2 Cloud Radiative Swath
    *DEVELOPMENT* HDF files - produce by running the CRS4 .f90 code on AMI
    ----------------------------------------------------------------------------
    :param file_path: path to data file  [string]
    :param var_name: variable name       [string]
    :param lev_arg: level argument       [integer]
    :param fill: fill NaN value option   [boolean]
    ----------------------------------------------------------------------------
    :return: (1) var_field = variable field data  [float]
             (2) var_name = variable name         [string]
             (3) var_units = variable units       [string]
             (4) lev_name = p level name          [string]
    """

    import numpy as np
    from pyhdf import SD
    hdf = SD.SD(file_path)

    # Select p level using integer index: lev_arg
    switch2 = {
        -1: ' ',
        0: 'TOA',
        1: '70',
        2: '200',
        3: '500',
        4: '850',
        5: 'surface'
    }
    lev_name = switch2.get(lev_arg)

    print("Getting", var_name, "at:",
          switch2.get(lev_arg, "N/A"))

    # select the variable and get the data
    data = hdf.select(var_name)
    variable = data.get()
    var_units = data.units

    # option to replace fill values with NaN
    if fill is True:
        variable[variable == data._FillValue] = np.nan

    var_field = np.empty([variable.shape[0]])
    # get field at appropriate vertical level
    if lev_arg < 0:
        var_field = variable
    elif lev_arg >= 0:
        var_field = variable[:, lev_arg]

    # elif lev_arg == 5:
    #     surface_indices = hdf.select('Surface level index')
    #     sfc_ind = surface_indices.get()
    #     sfc_ind[sfc_ind > 10] = 5  # if footprint is bad
    #     sfc_ind = list(sfc_ind)
    #     sfc_ind = [int(el) for el in sfc_ind]
    #     for i in range(variable.shape[0]):
    #         ind = np.int(sfc_ind[i] - 1.0)
    #         var_field[i] = variable[i, ind]

    return var_field, var_name, var_units, lev_name


# ==============================================================================


def read_day_of_ssf_files(path, file_struc, variable, index, fill):
    """
    ----------------------------------------------------------------------------
    This function loops over and reads an entire day (24 hr) of SSF data.
    ----------------------------------------------------------------------------
    :param path: path to directory containing data files
    :param file_struc: file structure (without the hour portion at the end)
    :param variable: variable to be read from file
    :param index:
    :param fill: option to replace fill values as NaN
    :return: variable, lat, lon, sza
    ----------------------------------------------------------------------------
    """
    import numpy as np

    print('============================================')
    print('\tReading SSF Files...\t\t\t')
    print('============================================')

    len_tot = []
    sza_all = np.empty([])
    tim_all = np.empty([])
    lat_all = np.empty([])
    lon_all = np.empty([])
    var_all = np.empty([])

    for k in range(24):        # Should be 24 for CERES SSF, 23 for FLASHFlux since
        if k < 10:             # the final swath length differs between FF and CERES
            k = '0' + str(k)   # Still need to look into why... and fix?

        file = file_struc + str(k)

        file_path = path + file
        print(file_path)

        # read geolocation info from file
        lat, lon, tim, sza = read_ssf_geolocation(file_path)
        # read variable from file
        var, var_name, var_units = \
            read_ssf_var(file_path=file_path,
                         var_name=variable,
                         index=index,
                         fill=fill)

        # len_tot contains the length (num FOVs) of each individual swath
        len_tot.append(lat.shape[0])

        # all of the swaths combined
        sza_all = np.concatenate((sza_all, sza), axis=None)
        tim_all = np.concatenate((tim_all, tim), axis=None)
        lat_all = np.concatenate((lat_all, lat), axis=None)
        lon_all = np.concatenate((lon_all, lon), axis=None)
        var_all = np.concatenate((var_all, var), axis=None)

    print(len_tot)
    print(var_all.shape)

    return var_all, lon_all, lat_all, sza_all


# ==============================================================================


def read_day_of_crs_files(path, file_struc, variable, lev_arg, fill):
    """
    ----------------------------------------------------------------------------
    This function loops over and reads an entire day (24 hr) of CRS data
    ----------------------------------------------------------------------------
    :param path: path to files
    :param file_struc: file structure (without the hour portion at the end)
    :param variable: variable to be read from file
    :param lev_arg: level argument (0 = TOA, 5 = sfc)
    :param fill: option to replace fill value as NaN
    :return: variable, lat, lon, sza
    ----------------------------------------------------------------------------
    """
    import numpy as np

    print('============================================')
    print('\tReading CRS Files...\t\t\t')
    print('============================================')

    len_tot = []
    sza_all = np.empty([])
    lat_all = np.empty([])
    lon_all = np.empty([])
    var_all = np.empty([])

    for k in range(24):
        if k < 10:
            k = '0' + str(k)

        file = file_struc + str(k)

        file_path = path + file
        print(file_path)

        lat, lon, pres_levs, obs_tim, sfc_ind, sza = read_crs_geolocation_dev(file_path)

        var, var_name, var_units, var_lev = \
            read_crs_var_dev(file_path=file_path,
                             var_name=variable,
                             lev_arg=lev_arg,
                             fill=fill)

        len_tot.append(lat.shape[0])
        sza_all = np.concatenate((sza_all, sza), axis=None)
        lat_all = np.concatenate((lat_all, lat), axis=None)
        lon_all = np.concatenate((lon_all, lon), axis=None)
        var_all = np.concatenate((var_all, var), axis=None)

    print(len_tot)
    print(var_all.shape)

    return var_all, lon_all, lat_all, sza_all


# =====================================================================


def read_month_of_ssf_files(path, file_struc, variable):
    """
    ----------------------------------------------------------------------------
    This function loops over and reads an entire month of SSF data
    ----------------------------------------------------------------------------
    :param path: path to files
    :param file_struc: file structure (without the day & hour portion at the end)
    :param variable: variable to be read from file
    :return: variable, lat, lon, sza
    ----------------------------------------------------------------------------
    """
    import numpy as np

    print('============================================')
    print('\tReading SSF Files...\t\t\t')
    print('============================================')

    len_tot = []
    sza_all = np.empty([])
    tim_all = np.empty([])
    lat_all = np.empty([])
    lon_all = np.empty([])
    var_all = np.empty([])

    for d in range(1, 16, 1):
        if d < 10:
            d = '0' + str(d)

        for k in range(24):
            if k < 10:
                k = '0' + str(k)

            file = file_struc + str(d) + str(k)

            file_path = path + file
            print(file_path)

            # read geolocation info from file
            lat, lon, tim, sza = read_ssf_geolocation(file_path)
            # read variable from file
            var, var_name, var_units = \
            read_ssf_var(file_path=file_path, var_name=variable, fill=True)

            # len_tot contains the length (num FOVs) of each individual swath
            len_tot.append(lat.shape[0])

            # all of the swaths combined
            sza_all = np.concatenate((sza_all, sza), axis=None)
            tim_all = np.concatenate((tim_all, tim), axis=None)
            lat_all = np.concatenate((lat_all, lat), axis=None)
            lon_all = np.concatenate((lon_all, lon), axis=None)
            var_all = np.concatenate((var_all, var), axis=None)

    print(len_tot)
    print(var_all.shape)

    return var_all, lon_all, lat_all, sza_all


# =============================================================================


def read_month_of_crs_files(path, file_struc, variable, lev_arg):
    """
    ----------------------------------------------------------------------------
    This function loops over and reads an entire month of CRS data
    ----------------------------------------------------------------------------
    :param path: path to files
    :param file_struc: file structure (without day & hour portion at the end)
    :param variable: variable to be read from file
    :param lev_arg: level argument (0 = TOA, 5 = sfc)
    :return: variable, lat, lon, sza
    ----------------------------------------------------------------------------
    """
    import numpy as np

    print('============================================')
    print('\tReading CRS Files...\t\t\t')
    print('============================================')

    len_tot = []
    sza_all = np.empty([])
    lat_all = np.empty([])
    lon_all = np.empty([])
    var_all = np.empty([])

    for d in range(1, 16, 1):
        if d < 10:
            d = '0' + str(d)

        for k in range(24):
            if k < 10:
                k = '0' + str(k)

            file = file_struc + str(d) + str(k)

            file_path = path + file
            print(file_path)

            lat, lon, pres_levs, obs_tim, sfc_ind, sza = read_crs_geolocation_dev(file_path)

            var, var_name, var_units, var_lev = \
            read_crs_var_dev(file_path=file_path,
                                   var_name=variable,
                                   lev_arg=lev_arg, fill=False)

            len_tot.append(lat.shape[0])
            sza_all = np.concatenate((sza_all, sza), axis=None)
            lat_all = np.concatenate((lat_all, lat), axis=None)
            lon_all = np.concatenate((lon_all, lon), axis=None)
            var_all = np.concatenate((var_all, var), axis=None)

    print(len_tot)
    print(var_all.shape)

    return var_all, lon_all, lat_all, sza_all


# ==============================================================================


def swath_lon_360_to_180(lon):
    """
    ----------------------------------------------------------------------------
    This function converts the range of footprint longitude values from
    0 to 360 deg to -180 to 180 deg
    ----------------------------------------------------------------------------
    :param lon: swath longitudes     ranging from 0 to 360      [float]
    :return: (1) lon_ new longitudes ranging from -180 to 180   [float]
    ----------------------------------------------------------------------------
    """
    import numpy as np

    lon_ = np.empty([len(lon)])
    print(lon_.shape)
    for i in range(len(lon)):
        if lon[i] > 180:
            lon_[i] = lon[i] - 360
        elif lon[i] <= 180:
            lon_[i] = lon[i]

    return lon_


# ==============================================================================


def set_colormap(cmap_name, typ_arg):
    """
    ----------------------------------------------------------------------------
    Selects colormap from palettable library
    ----------------------------------------------------------------------------
    :param cmap_name: name of colormap from palettable [string]
    :param typ_arg:   0 = continuous, 1 = discrete     [integer]
    ----------------------------------------------------------------------------
    :return: colormap object
    """

    switch = {
        0: "continuous",
        1: "discrete"
    }
    print("\nUsing a", switch.get(typ_arg, "N/A"), "colormap...")

    if typ_arg == 0:
        colormap = cmap_name.mpl_colormap
    elif typ_arg == 1:
        from matplotlib.colors import ListedColormap
        colormap = ListedColormap(cmap_name.mpl_colors)
    elif typ_arg != 1 or typ_arg != 0:
        print('Please input 0 for continuous or 1 for discrete ...')

    return colormap


# ==============================================================================


def plot_swath(lon, lat, field,
               varname, levname, varunits,
               nrows, ncols, cen_lon,
               cmap, cmap_lims, date, nightshade,
               date_str, title_str):
    """
    ----------------------------------------------------------------------------
    This function plots a swath of footprint-level CERES data
    e.g., FLASHFlux, SSF, or CRS
    ----------------------------------------------------------------------------
    :param lon: FOV longitude               [float]
    :param lat: FOV latitude                [float]
    :param field: variable                  [float]
    :param varname: variable name           [string]
    :param levname: level name              [string]
    :param varunits: variable units         [string]
    :param nrows: number of rows            [integer]
    :param ncols: number of columns         [integer]
    :param cen_lon: central longitude       [float]
    :param cmap: colormap                   []
    :param cmap_lims: colormap limits       [tuple]
    :param date: date info for nightshade   []
    :param nightshade: use nightshade?      [boolean]
    :param date_str: date string            [string]
    :param title_str: title string          [string]
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
    projection = ccrs.Robinson(central_longitude=cen_lon)

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
        #ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
        ax.set_global()
        ax.text(0.5, -0.1, varname + ' - ' + levname + '\n' + varunits,
                va='bottom', ha='center',
                rotation='horizontal', rotation_mode='anchor',
                transform=ax.transAxes, fontsize=10)

        if nightshade is True:
            ax.add_feature(Nightshade(date, alpha=0.15))

    # To use a different colorbar range each time, use a tuple of tuples
    limits = cmap_lims
    for i in range(nrows * ncols):
        im = axgr[i].scatter(lon, lat, c=field, s=1, transform=ccrs.PlateCarree(),
                             vmin=limits[0], vmax=limits[1], cmap=cmap)
        axgr.cbar_axes[i].colorbar(im)

    for i, cax in enumerate(axgr.cbar_axes):
        cax.set_yticks(np.linspace(limits[0], limits[1], 5))
        cax.set_yticklabels(np.linspace(limits[0], limits[1], 5),
                            fontsize=8)

    plt.show()
    return


# ==============================================================================


def swath_daytime_only(lat, lon, var, sza, sza_cutoff):
    """
    ----------------------------------------------------------------------------
    This function extracts daytime footprints from a time series of CERES FOVs
    using a solar zenith angle cut-off value. Nighttime values more than the
    threshold (e.g., 90 degrees) are ignored.
    ----------------------------------------------------------------------------
    :param lat: FOV latitudes                   [float]
    :param lon: FOV longitudes                  [float]
    :param var: variable under consideration    [float]
    :param sza: FOV solar zenith angle          [float]
    :param sza_cutoff: SZA cut-off value        [float]
    :return:
    ----------------------------------------------------------------------------
    """
    import numpy as np

    # if footprint SZA > cutoff value, ignore
    for i in range(len(sza)):
        if sza[i] >= sza_cutoff:
            sza[i] = np.nan
            var[i] = np.nan
            lat[i] = np.nan
            lon[i] = np.nan

    # ignore and remove NaNs
    bad_indices = np.isnan(var)
    good_indices = ~bad_indices
    lat = lat[good_indices]
    lon = lon[good_indices]
    sza = sza[good_indices]
    var = var[good_indices]

    return lat, lon, var, sza


# ==============================================================================


def swath_nighttime_only(lat, lon, var, sza, sza_cutoff):
    """
    ----------------------------------------------------------------------------
    This function extracts nighttime footprints from a time series of CERES FOVs
    using a solar zenith angle cut-off value. Daytime values less than the
    threshold (e.g., 90 degrees) are ignored.
    ----------------------------------------------------------------------------
    :param lat: FOV latitudes                   [float]
    :param lon: FOV longitudes                  [float]
    :param var: variable under consideration    [float]
    :param sza: FOV solar zenith angle          [float]
    :param sza_cutoff: SZA cut-off value        [float]
    :return:
    ----------------------------------------------------------------------------
    """
    import numpy as np

    # if footprint SZA < cutoff value, ignore
    for i in range(len(sza)):
        if sza[i] <= sza_cutoff:
            sza[i] = np.nan
            var[i] = np.nan
            lat[i] = np.nan
            lon[i] = np.nan

    # ignore and remove NaNs
    bad_indices = np.isnan(var)
    good_indices = ~bad_indices
    lat = lat[good_indices]
    lon = lon[good_indices]
    sza = sza[good_indices]
    var = var[good_indices]

    return lat, lon, var, sza


# ==============================================================================


def swath_difference(field2, field1, day_only, sza):
    """
    ----------------------------------------------------------------------------
    This function calculates the difference between two swaths, which requires
    them to have the same footprints co-located in time and space. Its purpose
    is to compare different versions of CERES products, such as SSF-
    parameterized vs CRS-computed surface radiative fluxes.
    ----------------------------------------------------------------------------
    :param field2: swath 2 variable                      [float]
    :param field1: swath 1 variable                      [float]
    :param day_only: option to only use daytime FOVs     [boolean]
    :param sza: solar zenith angle                       [float]
    ----------------------------------------------------------------------------
    :return: difference = difference between swaths      [float]
    """

    import numpy as np

    difference = field2 - field1

    # a footprint might be good in one but bad in another...
    difference[difference == max(difference)] = np.nan
    difference[difference == min(difference)] = np.nan

    # field2.shape[0] = field1.shape[0]
    if day_only is True:
        for i in range(field2.shape[0]):
            if sza[i] > 90:
                difference[i] = np.nan

    # night time only
    elif day_only is False:
        for i in range(field2.shape[0]):
            if sza[i] < 90:
                difference[i] = np.nan

    return difference


# ==============================================================================


def swath_histogram_scatterplot(field2, field1, var_name, lev_name,
                                suptitle, ti_str2, ti_str1,
                                date_str, platform, day_only, sza):
    """
    ----------------------------------------------------------------------------
    This function creates a histogram of the difference between two CERES L2
    footprint-level swath fields and also plots a scatterplot of the data. It
    also yields basic descriptive statistics about the data.
    ----------------------------------------------------------------------------
    :param field2:                                   [float]
    :param field1:                                   [float]
    :param var_name: name of variable being compared [string]
    :param lev_name: name of vertical level          [string]
    :param date_str: time/data information           [string]
    :param platform: satellite and flight model      [string]
    :param day_only: day only?                       [boolean]
    :param sza: FOV solar zenith angle               [float]
    ----------------------------------------------------------------------------
    :return: plot of histogram and scatterplot
    """

    import numpy as np
    import matplotlib.pyplot as plt

    swath_diff = swath_difference(field2=field2, field1=field1, day_only=day_only, sza=sza)

    mean_diff = np.nanmean(swath_diff)             # mean of the difference, i.e., bias
    sigma_diff = np.nanstd(swath_diff)             # standard deviation of the difference
    rms_diff = np.sqrt(np.nanmean(swath_diff**2))  # RMS difference = RMSD

    # number of FOVs compared
    print("Number of FOVs in swath: ", len(swath_diff))
    print("Number of FOVs where at least one is NaN: ", np.sum(np.isnan(swath_diff)))
    N = len(swath_diff) - np.sum(np.isnan(swath_diff))
    print("Number of FOVs compared: ", N)

    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    fig.suptitle('CERES ' + platform + ' - ' + suptitle + ' - ' + date_str, fontsize=11)

    # histogram of difference
    axs[0].hist(swath_diff, bins=200, align='mid', rwidth=1)
    axs[0].grid()
    axs[0].set_axisbelow("True")
    axs[0].set_title(ti_str2 + ' minus ' + ti_str1 + ' Flux', fontsize=11)
    axs[0].set_xlabel(var_name + ' - ' + lev_name + "\n" + r"Flux difference ($\Delta$)")
    axs[0].set_xlim([-125, 125])
    axs[0].set_ylabel("Number of CERES Footprints")

    # show basic descriptive statistics
    textstr = "N = " + str(N) + "\n" + \
              r"Mean: $\bar{\Delta}$ = " + str(np.around(mean_diff, 2)) + "\n" + \
              r"Stdv: $\sigma_{\Delta}$ = " + str(np.around(sigma_diff, 2)) + "\n" + \
              "RMSD = " + str(np.around(rms_diff, 2))
    props = dict(facecolor='white', alpha=0.85)
    axs[0].text(0.05, 0.95, textstr, transform=axs[0].transAxes, fontsize=8,
                verticalalignment='top', bbox=props)

    if day_only is True:
        day_only_str = "Daytime only"
        axs[0].text(0.05, 0.8, day_only_str, transform=axs[0].transAxes, fontsize=8,
                    verticalalignment='top', bbox=props)
    elif day_only is False:
        night_only_str = "Nighttime only"
        axs[0].text(0.05, 0.8, night_only_str, transform=axs[0].transAxes, fontsize=8,
                    verticalalignment='top', bbox=props)

    field1[field1 > 1.e6] = np.nan
    field2[field2 > 1.e6] = np.nan

    # remove NaN values
    bad_indices = np.isnan(field1) | np.isnan(field2)
    good_indices = ~bad_indices
    field1 = field1[good_indices]
    field2 = field2[good_indices]

    if np.nanmax(field1) < np.nanmax(field2):
        max = np.nanmax(field2) + 100
    else:
        max = np.nanmax(field1) + 100

    # scatter plot
    axs[1].plot(range(2000), range(2000), color='black', zorder=0)
    axs[1].scatter(field1, field2, s=0.05)
    # axs[1].hist2d(field1, field2, bins=50, cmap='ocean_r')
    axs[1].set(xlim=(0, max), ylim=(0, max))
    axs[1].grid()
    axs[1].set_axisbelow("True")
    axs[1].set_title(var_name + ' - ' + lev_name, fontsize=11)
    axs[1].set_xlabel(ti_str1, fontsize=11)
    axs[1].set_ylabel(ti_str2, fontsize=11)

    # show descriptive statistics
    str1 = "r = " + str(np.around(np.corrcoef(field1, field2)[0][1], decimals=3))
    props = dict(facecolor='white', alpha=0.85)
    axs[1].text(0.05, 0.95, str1, transform=axs[1].transAxes, fontsize=8,
                verticalalignment='top', bbox=props)

    plt.show()
    return

# =============================================================================


def grid_to_1x1_deg_equal_angle(lat_data, lon_data, variable, lon_360=True):
    """
    ----------------------------------------------------------------------------
    This functions bins & grids CERES footprints to a 1 deg x 1 deg equal angle
    latitude-longitude grid using the SciPy stats routine binned_statistic_2d.
    After FOVs are aggregated into 1 deg x 1 deg regions it computes the mean
    of the input "variable" - alternatively, it can compute the # of footprints,
    the median, or other statistics (including user-defined functions).
    ----------------------------------------------------------------------------
    :param lat_data: FOV latitude array                            [float]
    :param lon_data: FOV longitude array                           [float]
    :param variable: for which gridded statistic will be computed  [float]
    :param lon_360:  use 0 to 360 or -180 to 180 longitude bins    [boolean]
    :return: the field of gridded FOVs                             [float]
    ----------------------------------------------------------------------------
    """
    import time
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats

    # consider generalizing and adding bins as parameters of the function?
    lat_bins = np.arange(-90, 91)    # lat: -90 to 90 deg

    if lon_360 is True:
        lon_bins = np.arange(0, 361)     # lon: 0 to 360
    elif lon_360 is False:
        lon_bins = np.arange(-180, 181)  # lon: -180 to 180

    print('Griding and averaging footprints to 1 x 1',
          'degree equal angle grid boxes...\n')

    print('Lat bins:\n')
    print(lat_bins)
    print('Lon bins:\n')
    print(lon_bins)

    # # each FOV has lat and lon indices that map to each grid box
    # lat_ind = np.digitize(lat_data, lat_bins)
    # lon_ind = np.digitize(lon_data, lon_bins)
    #
    # # loop over FOVs and show their index
    # for n in range(lat_data.size):
    #     print(lat_bins[lat_ind[n]-1], "<=", lat_data[n], "<", lat_bins[lat_ind[n]], '...', lat_ind[n])
    #     print(lon_bins[lon_ind[n]-1], "<=", lon_data[n], "<", lon_bins[lon_ind[n]], '...', lon_ind[n])

    # start timing it
    tic = time.time()

    # compute mean in each grid box - statistics: 'count', 'mean', 'median'
    gridded = stats.binned_statistic_2d(lon_data, lat_data, variable,
                                        statistic=np.nanmean,
                                        bins=[lon_bins, lat_bins])

    gridded_stat = np.rot90(gridded.statistic)

    # finish timing it, print relevant info
    toc = time.time()
    print(toc - tic, 'seconds elapsed during grid_to_1x1_deg_equal_angle\n')

    print("Shape of 1 x 1 gridded field:")
    print(gridded_stat.shape)

    # quick & dirty plot of the result
    # plt.pcolor(gridded_stat)
    # plt.colorbar()
    # plt.show()

    # might be nice to write the result to a netCDF or HDF file
    # in cases where this takes a long time to run...

    return gridded_stat


# =============================================================================


def grid_to_1x1_deg_ceres_nested(lat_data, lon_data, variable, lon_360=True):
    """
    ----------------------------------------------------------------------------
    This functions bins & grids CERES footprints to a 1 deg x 1 deg equal angle
    latitude-longitude grid using the SciPy stats routine binned_statistic_2d.
    After FOVs are aggregated into 1 deg x 1 deg regions it computes the mean
    of the input "variable" - alternatively, it can compute the # of footprints,
    the median, or other statistics (including user-defined functions).
    ----------------------------------------------------------------------------
    :param lat_data: FOV latitude array                            [float]
    :param lon_data: FOV longitude array                           [float]
    :param variable: for which gridded statistic will be computed  [float]
    :param lon_360:  use 0 to 360 or -180 to 180 longitude bins    [boolean]
    :return: the field of gridded FOVs                             [float]
    ----------------------------------------------------------------------------
    """
    import time
    import copy
    import numpy as np
    from scipy import stats

    # 1 degree latitude bins, -90 to 90 deg
    lat_bins = np.arange(-90, 91)

    # different longitude bin extents for different latitude zones
    lon_ext = [1, 2, 4, 8]

    # loop over nested grid zones
    for i, el in enumerate(lon_ext):

        lon_bins = np.empty([])
        # set appropriate longitude extent
        if lon_360 is True:                      # lon ranges 0 to 360
            lon_bins = np.arange(0, 361, el)
        elif lon_360 is False:                   # lon ranges -180 to 180
            lon_bins = np.arange(-180, 181, el)  # different bin sizes

        print('Zone {}, lon width {}'.format(i+1, el))
        print('Lon bins:\n', lon_bins)

        # time the griding procedure
        tic = time.time()

        # Zone 1
        if i == 0:
            var1 = copy.copy(variable)
            var1[abs(lat_data) > 45] = np.nan
            zone_1 = stats.binned_statistic_2d(lon_data, lat_data, var1,
                                               statistic=np.nanmean,
                                               bins=[lon_bins, lat_bins])
            zone1 = np.rot90(zone_1.statistic)

        # Zone 2
        if i == 1:
            var2 = copy.copy(variable)
            var2[abs(lat_data) < 45] = np.nan
            var2[abs(lat_data) > 70] = np.nan
            zone_2 = stats.binned_statistic_2d(lon_data, lat_data, var2,
                                               statistic=np.nanmean,
                                               bins=[lon_bins, lat_bins])
            zone2 = np.rot90(zone_2.statistic)

        # Zone 3
        if i == 2:
            var3 = copy.copy(variable)
            var3[abs(lat_data) < 70] = np.nan
            var3[abs(lat_data) > 80] = np.nan
            zone_3 = stats.binned_statistic_2d(lon_data, lat_data, var3,
                                               statistic=np.nanmean,
                                               bins=[lon_bins, lat_bins])
            zone3 = np.rot90(zone_3.statistic)

        # Zone 4
        if i == 3:
            var4 = copy.copy(variable)
            var4[abs(lat_data) < 80] = np.nan
            var4[abs(lat_data) > 89] = np.nan
            zone_4 = stats.binned_statistic_2d(lon_data, lat_data, var4,
                                               statistic=np.nanmean,
                                               bins=[lon_bins, lat_bins])
            zone4 = np.rot90(zone_4.statistic)

        # finish timing it, print relevant info
        toc = time.time()
        print(toc - tic, 'seconds elapsed during grid_to_1x1_ceres_nested\n')

    print('Populating CERES 1 deg x 1 deg Nested Grid\n')
    nested_grid = np.empty([180, 360])

    for k in range(180):
        if k == 0:  # zone 5 : pole, 1 box
            nested_grid[k, :] = 0
        elif 1 <= k <= 9:  # zone 4 : 80 to 89 deg, 9 boxes
            nested_grid[k, :] = np.repeat(zone4[k, :], 8)
        elif 10 <= k <= 19:  # zone 3 : 70 to 80 deg, 10 boxes
            nested_grid[k, :] = np.repeat(zone3[k, :], 4)
        elif 20 <= k <= 44:  # zone 2 : 45 to 70 deg, 25 boxes
            nested_grid[k, :] = np.repeat(zone2[k, :], 2)
        elif 45 <= k <= 134:  # zone 1 : -45 to 45 deg, 90 boxes
            nested_grid[k, :] = zone1[k, :]
        elif 135 <= k <= 159:  # zone 2
            nested_grid[k, :] = np.repeat(zone2[k, :], 2)
        elif 160 <= k <= 169:  # zone 3
            nested_grid[k, :] = np.repeat(zone3[k, :], 4)
        elif 170 <= k <= 178:  # zone 4
            nested_grid[k, :] = np.repeat(zone4[k, :], 8)
        elif k == 179:  # zone 5 - pole
            nested_grid[k, :] = 0

    return nested_grid


# =============================================================================


def plot_gridded_fields(nrows, ncols, cen_lon,
                        date_str, title_str,
                        cmap, cmap_lims,
                        varname, levname, varunits,
                        lon, lat, field):
    """
    ----------------------------------------------------------------------------
    This function provides a framework for plotting multiple different gridded
    fields.
    ----------------------------------------------------------------------------
    :param nrows: number of rows
    :param ncols: number of columns
    :param cen_lon: central longitude
    :param date_str: date string
    :param title_str: title string
    :param cmap: colormap
    :param cmap_lims: colormap limits
    :param varname: variable name
    :param levname: level name
    :param varunits: variable units
    :param lon: longitude
    :param lat: latitude
    :param field: variable to plot
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
    projection = ccrs.Robinson(central_longitude=cen_lon)

    # Axis class
    axes_class = (GeoAxes, dict(map_projection=projection))

    # Create figure
    fig = plt.figure(figsize=(10, 8))
    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(nrows, ncols),
                    axes_pad=(0.4, 0.6),
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
        #ax.set_extent([-180, 180, -90, 90], projection)
        ax.text(0.5, -0.1, varname + ' - ' + levname + '\n' + varunits,
                va='bottom', ha='center',
                rotation='horizontal', rotation_mode='anchor',
                transform=ax.transAxes, fontsize=10)

    # To use a different color bar range each time use a tuple of tuples
    #        = ((0, 120), (0, 120), (0, 120), (0, 120), (0, 120), (0, 120), (0, 120))
    # limits = (-30, 30)
    limits = cmap_lims
    for i in range((nrows * ncols)):
        im = axgr[i].pcolor(lon, lat, field[:, :, i], transform=ccrs.PlateCarree(),
                            vmin=limits[i][0], vmax=limits[i][1], cmap=cmap)
        axgr.cbar_axes[i].colorbar(im)
        axgr[i].set_title(title_str[i], fontsize=10)

    plt.show()
    return


# ========================================================================


def read_syn1deg_hdf(file_path, var_name, fill):
    """
    ----------------------------------------------------------------------------
    This function reads data from CERES Level 3 Synoptic 1-degree (SYN1deg)
    HDF data files. Presumably, it should also work for other L3 files...
    ----------------------------------------------------------------------------
    :param file_path:
    :param var_name: variable name [string]
    :return:
    ----------------------------------------------------------------------------
    """
    import numpy as np
    from pyhdf import SD
    hdf = SD.SD(file_path)

    # select and get the variable
    data = hdf.select(var_name)
    variable = data.get()
    var_name2 = data.long_name
    var_units = data.units

    print('Reading... ' + var_name2)

    # replace fill values with NaN
    if fill is True:
        print('Replacing fill values with NaN...')
        variable[variable == data._FillValue] = np.nan

    var_field = variable

    return var_field, var_name2, var_units


# ========================================================================


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


# ========================================================================


def read_ebaf_geolocation(file_path):
    """
    ----------------------------------------------------------------------------
    This function reads geolocation information from CERES EBAF files and then
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

    latitude = nc.variables['lat'][:]
    longitude = nc.variables['lon'][:]

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


def read_ebaf_start_month_year(file_path):
    """
    ----------------------------------------------------------------------------
    This function extracts the starting month and year from EBAF file
    ----------------------------------------------------------------------------
    :param file_path: path to EBAF file
    ----------------------------------------------------------------------------
    :return: (1) month
             (2) year
    """

    import numpy as np
    from netCDF4 import Dataset

    nc = Dataset(file_path, 'r')

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


def read_ebaf_var(file_path, variable):
    """
    ----------------------------------------------------------------------------
    Extract the data for a particular variable from an EBAF netCDF file
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


# ========================================================================


def compute_monthly_anomalies(field, field_name):
    """
    ----------------------------------------------------------------------------
    (1) Compute long-term average for each calendar month, i.e., the monthly mean seasonal cycle
    (2) Subtract the appropriate long-term means from the corresponding monthly fields to get monthly anomalies
    ----------------------------------------------------------------------------
    :param field: variable for which to compute the climatology/anomaly [lat,lon,time]
    :param field_name: name of variable
    ----------------------------------------------------------------------------
    :return: (1) monthly_anomalies : monthly anomaly time series at each grid box
             (2) seasonal_cycle    : monthly mean seasonal cycle at each grid box
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
        # print(i, num_yr)

    # repeat 20 times
    seasonal_cycle = np.tile(monthly_climatology, (num_yr, 1, 1))

    # only take data up to the last index
    seasonal_cycle = seasonal_cycle[0:record_len, :, :]

    # calculate anomalies
    monthly_anomalies = np.subtract(field, seasonal_cycle)

    # add which lat/lon this corresponds to on the plots...
    r = 100
    p = 120
    # add lat lon location info to plots...
    fig, (ax1, ax2) = plt.subplots(2, 1)
    # first set of axes
    ax1.plot(np.linspace(0, 0, record_len))
    ax1.plot(monthly_anomalies[:, r, p], label='Anomalies')
    ax1.legend(fontsize=10)
    ax1.set_xticks(ticks=range(0, record_len, 12))
    ax1.set_xticklabels(range(0, int(record_len/12+1)))
    ax1.set_title(field_name)
    ax1.grid()
    # second set of axes
    ax2.plot(seasonal_cycle[:, r, p], label='Seasonal Cycle')
    ax2.plot(field[:, r, p], label='Raw Field')
    ax2.legend(fontsize=10)
    ax2.set_xticks(ticks=range(0, record_len, 12))
    ax2.set_xticklabels(range(0, int(record_len/12+1)))
    ax2.grid()
    plt.show()

    return monthly_anomalies, seasonal_cycle


# ========================================================================


def cos_lat_weight(lat_vector):
    """
    ----------------------------------------------------------------------------
    This function takes a latitude array and creates a matrix of cos(lat)
    weights in order to perform area-weighted averaging of a field.
    ----------------------------------------------------------------------------
    :param lat_vector: vector of latitudes
    ----------------------------------------------------------------------------
    :return: (1) weight: matrix of cos(lat) weights
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
    This function calculates the long-term temporal mean and standard deviation
    of the input field.
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

    # Global mean
    global_avg = np.ma.average(np.ma.masked_array(field, np.isnan(field)), weights=weights)
    print('\n', 'Global average:', global_avg, '\n')

    # 60Sto90S mean
    sh60to90_avg = np.ma.average(np.ma.masked_array(field[0:30, :], np.isnan(field[0:30, :])),
                                 weights=weights[0:30, :])
    print('Antarctic (60S-90S) average:', sh60to90_avg)
    print('Latitude range:', latitude[0:30], '\n')

    # SH mid-latitude
    sh30to60_avg = np.ma.average(np.ma.masked_array(field[30:60, :], np.isnan(field[30:60, :])),
                                 weights=weights[30:60, :])
    print('S mid-latitude (30S-60S) average:', sh30to60_avg)
    print('Latitude range:', latitude[30:60], '\n')

    # SH tropical
    sh0to30_avg = np.ma.average(np.ma.masked_array(field[60:90, :], np.isnan(field[60:90, :])),
                                weights=weights[60:90, :])
    print('S tropical (30S-0) average:', sh0to30_avg)
    print('Latitude range:', latitude[60:90], '\n')

    # NH tropical mean
    nh0to30_avg = np.ma.average(np.ma.masked_array(field[90:120, :], np.isnan(field[90:120, :])),
                                weights=weights[90:120, :])
    print('N tropical (0-30N) average:', nh0to30_avg)
    print('Latitude range:', latitude[90:120], '\n')

    # NH mid-latitude
    nh30to60_avg = np.ma.average(np.ma.masked_array(field[120:150, :], np.isnan(field[120:150, :])),
                                 weights=weights[120:150, :])
    print('N mid-latitude (30N-60N) average:', nh30to60_avg)
    print('Latitude range:', latitude[120:150], '\n')

    # Arctic mean climatological total cloud fraction
    nh60to90_avg = np.ma.average(np.ma.masked_array(field[150:180, :], np.isnan(field[150:180, :])),
                                 weights=weights[150:180, :])
    print('Arctic 60N-90N average:', nh60to90_avg)
    print('Latitude range:', latitude[150:180], '\n')

    return global_avg, sh0to30_avg, sh30to60_avg, sh60to90_avg, nh0to30_avg, nh30to60_avg, nh60to90_avg


# ========================================================================


def composite_difference(field, ind1, ind2):
    """
    ----------------------------------------------------------------------------
    This function calculates the composite mean of 'field' for two different
    time periods and then takes the difference to illustrate the field change.
    ----------------------------------------------------------------------------
    :param field: field to average                [float; ntime x nlat x nlon]
    :param ind1: indices for first time period    [integer range]
    :param ind2: indices for second time period   [integer range]
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
    import numpy as np

    mean_time_series = np.empty(shape=field.shape[0])
    for i in range(field.shape[0]):
        field1 = field[i, :, :]
        mean_time_series[i] = np.average(field1, weights=weight)
        # print(mean_time_series[i])

    return mean_time_series


# ========================================================================


def plot_monthly_time_series(var, name, units, start_mo, start_yr):
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
    plt.yticks(fontsize=10)
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
    :return: matrix of linear regression coefficients
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
    projection = ccrs.Robinson(central_longitude=cen_lon)

    # Axis class
    axes_class = (GeoAxes, dict(map_projection=projection))

    # Plot figure
    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(nrows, ncols),
                    axes_pad=(0.4, 0.4),
                    share_all=True,
                    cbar_location='right',
                    cbar_mode='each',
                    cbar_pad=0.1,
                    cbar_size='5%',
                    label_mode=1)

    for i, ax in enumerate(axgr):
        ax.add_feature(cartopy.feature.LAND, zorder=1, facecolor='none', edgecolor='grey')
        ax.gridlines(color='black', linestyle=':')
        ax.set_title(ti_str)
        #ax.set_extent([-180, 180, -60, 60], ccrs.PlateCarree())
        ax.set_global()
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

