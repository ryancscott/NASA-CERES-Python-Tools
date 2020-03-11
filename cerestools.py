# =====================================================================
#                           CERES PYTHON TOOLS
# =====================================================================
#
# Author: Ryan C. Scott, ryan.c.scott@nasa.gov
#
# Purpose: This library contains functions to process and manipulate
#          satellite data from the NASA Clouds and the Earth's Radiant
#          Energy System (CERES) satellite mission, including
#          footprint-level and gridded fields.
#
# =====================================================================


def read_crs_geolocation(file_path):
    """
    Reads geolocation information from file and converts
    colatitude to latitude for each CERES footprint/FOV
    :param file_path: path to input file
    :return: lat/lon of each FOV, pressure levels, and
    the time of each observation
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
    pres_lev = pressure_levels.get()

    time_of_obs = hdf.select("Time of observation")
    time_obs = time_of_obs.get()

    return fov_lat, fov_lon, pres_lev, time_obs


# ==============================================================================


def read_crs_var(file_path, vararg, levarg, fill):
    """Select field at appropriate vertical level
    from the official CERES CRS HDF file.

    Function arguments:

    (1) file_path = path to file
    (2) vararg = variable argument and
    (3) levarg = level argument
    (4) fill = option to fill missing values
    """

    from pyhdf import SD
    hdf = SD.SD(file_path)

    # select variable
    switch = {
        0: 'Longwave flux - downward - total',
        1: 'Longwave flux - upward - total',
        2: 'Shortwave flux - downward - total',
        3: 'Shortwave flux - upward - total'
    }
    var_name = switch.get(vararg)

    # select level
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


def get_date(file):

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


def read_crs_var2(file_path, vararg, levarg, fill):
    """
    Select field at desired vertical level from
    my run of the CERES CRS4 code.

    Function arguments:

    (1) file_path = path to file   [string]
    (2) vararg = variable argument [int]
    (3) levarg = level argument    [int]
    (4) fill = fillvalues -> NaN   [logical]
    """

    from pyhdf import SD
    hdf = SD.SD(file_path)

    # select variable
    switch = {
        0: 'UT_TOT_LW_DN',
        1: 'UT_TOT_LW_UP',
        2: 'UT_TOT_SW_DN',
        3: 'UT_TOT_SW_UP'
    }
    var_name = switch.get(vararg)

    # select level
    switch2 = {
        0: 'TOA',
        1: 'surface'
    }
    lev_name = switch.get(levarg)

    print("Reading", switch.get(vararg, "N/A"), "at:",
          switch2.get(levarg, "N/A"))

    # select the variable and get the data
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


def compute_diff(field2, field1):

    import numpy as np
    diff = field2 - field1

    diff[diff == max(diff)] = np.nan
    diff[diff == min(diff)] = np.nan

    return diff


# ==============================================================================


def set_colormap(cmap_name, typarg):
    """
    Selects colormap for plotting using palettable library

    Function arguments:

    (1) colormap name from palettable
    (2) colormap type: 0 for continuous, 1 for discrete
    """

    switch = {
        0: "continuous",
        1: "discrete"
    }
    print("\nUsing", switch.get(typarg, "N/A"), "colormap")

    if typarg == 0:
        c_map = cmap_name.mpl_colormap
    elif typarg == 1:
        from matplotlib.colors import ListedColormap
        c_map = ListedColormap(cmap_name.mpl_colors)
    elif typarg != 1 or typarg != 0:
        print('Please select 0 or 1')

    return c_map


# ==============================================================================


def plot_swath(nrows, ncols, cen_lon,
                field, varname, levname, varunits,
                cmap, cmap_lims,
                date, nightshade, title_str):
    """
    Plots map of CERES footprint swath data.

    Function arguments:

    (1) number of rows     [int]
    (2) number of columns  [int]
    (3) central longitude  [float]
    (4) field to plot      [string]
    (5) variable name      [string]
    (6) level name         [string]
    (7) variable units     [string]
    (8) colormap           [object]
    (9) colormap limits    [tuple]
    (10) date info         [object]
    (11) nightshade option [object]
    (12) title string      [string]
    """

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
        ax.text(0.5, -0.1, varname + ' - ' + levname + ' \n' + varunits,
                va='bottom', ha='center',
                rotation='horizontal', rotation_mode='anchor',
                transform=ax.transAxes, fontsize=10)

        if nightshade == 1:
            ax.add_feature(Nightshade(date, alpha=0.15))

    # To use a different colorbar range every time, use a tuple of tuples
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


def histogram_scatterplot(field1, field2, varname, levname, time_info):

    import numpy as np
    import matplotlib.pyplot as plt

    swath_diff = compute_diff(field2=field2, field1=field1)

    mean_diff = np.nanmean(swath_diff)
    stdv_diff = np.nanstd(swath_diff)

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
              r"Stdv: $\sigma_{\Delta}$ = " + str(np.around(stdv_diff, 2))
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
    Function prints out information about variables
    in the input netCDF file
    :param filepath: path to netCDF file
    :return: print to screen
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
    Reads lat, lon info from file
    :param filepath: path to input file
    :return:
    (1) latitude vector
    (2) longitude vector
    (3) lat grid
    (4) lon grid
    """

    import numpy as np
    from netCDF4 import Dataset

    nc = Dataset(filepath, 'r')

    latitude = nc.variables['lat'][:]
    longitude = nc.variables['lon'][:]

    # lat, lon grid
    lon, lat = np.meshgrid(longitude, latitude)

    print('lat array shape...')
    print(lat.shape, '\n')
    print('lon array shape...')
    print(lon.shape, '\n')

    return latitude, longitude, lat, lon


# ========================================================================


def read_ebaf_var(filepath, variable):
    """
    Extract variable data from EBAF file
    :param filepath: input EBAF file
    :param var: variable to extract from file
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


def compute_monthly_anomalies(field):
    """
    Compute long-term average for each calendar month, i.e., the monthly seasonal cycle,
    and then subtract the long-term means from the raw fields to get monthly anomalies
    :param field: variable for which to compute the climatology/anomaly [lat,lon,time]
    :return:
    (1) time series of monthly anomalies
    (2) monthly mean seasonal cycle
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

    r = 100
    p = 120
    fig, (ax1, ax2) = plt.subplots(2, 1)
    # first set of axes
    ax1.plot(np.linspace(0, 0, record_len))
    ax1.plot(monthly_anomalies[:, r, p], label='Anomalies')
    ax1.legend(fontsize=10)
    # second set of axes
    ax2.plot(seasonal_cycle[:, r, p], label='Seasonal Cycle')
    ax2.plot(field[:, r, p], label='Raw Field')
    ax2.legend(fontsize=10)

    return monthly_anomalies, seasonal_cycle


# ========================================================================


def cos_lat_weight(lat_vector):
    """
    Returns cos(lat) weight array [size nlat x nlon]
    :param lat_vector: vector of latitudes
    :return: matrix of cos(lat) weights
    """

    import numpy as np
    import numpy.matlib

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

    import matplotlib.pyplot as plt
    plt.imshow(weight)

    return weight


# ========================================================================


def compute_climatology(field):
    """
    Compute mean and standard deviation of field
    :param field: field under consideration
    :return: the mean and standard deviation
    """

    import numpy as np

    # average over entire time series
    mean_field = np.nanmean(field, axis=0)

    # standard deviation over entire time series
    stdv_field = np.nanstd(field, axis=0)

    return mean_field, stdv_field


# ========================================================================


def compute_regional_averages(field, latitude, w):
    """
    Computes weighted regional averages of input field
    :param field: field to be averaged
    :return: prints regional averages
    """
    import numpy as np

    # Global mean climatological total cloud fraction
    global_avg = np.average(field, weights=w)
    print('\n', 'Global average:', global_avg, '\n')

    # 60Sto90S mean climatological total cloud fraction
    sh60to90_avg = np.average(field[0:30, :], weights=w[0:30, :])
    print('Antarctic (60S-90S) average:', sh60to90_avg)
    print('Latitude range:', latitude[0:30], '\n')

    # SH mid-latitude climatological total cloud fraction
    sh30to60_avg = np.average(field[30:60, :], weights=w[30:60, :])
    print('S mid-latitude (30S-60S) average:', sh30to60_avg)
    print('Latitude range:', latitude[30:60], '\n')

    # SH tropical mean climatology
    sh0to30_avg = np.average(field[60:90, :], weights=w[60:90, :])
    print('S tropical (30S-0) average:', sh0to30_avg)
    print('Latitude range:', latitude[60:90], '\n')

    # NH tropical mean climatological total cloud fraction
    nh0to30_avg = np.average(field[90:120, :], weights=w[90:120, :])
    print('N tropical (0-30N) average:', nh0to30_avg)
    print('Latitude range:', latitude[90:120], '\n')

    # NH tropical mean climatological total cloud fraction
    nh30to60_avg = np.average(field[120:150, :], weights=w[120:150, :])
    print('N mid-latitude (30N-60N) average:', nh30to60_avg)
    print('Latitude range:', latitude[120:150], '\n')

    # Arctic mean climatological total cloud fraction
    nh60to90_avg = np.average(field[150:180, :], weights=w[150:180, :])
    print('Arctic 60N-90N average:', nh60to90_avg)
    print('Latitude range:', latitude[150:180], '\n')

    return sh0to30_avg, sh30to60_avg, sh60to90_avg, nh0to30_avg, nh30to60_avg, nh60to90_avg


# ========================================================================


def composite_diff(field, ind1, ind2):
    """
    Calculates composite means and then takes
    their difference to illustrate change between
    two time periods
    :param field: field to average
    :param ind1: indices for first time period
    :param ind2: indices for second time period
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


def global_map(lon, lat, field,
               varname, varunits,
               nrows, ncols, cen_lon,
               cmap, cmap_lims, ti_str):

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature
    from mpl_toolkits.axes_grid1 import AxesGrid
    from cartopy.mpl.geoaxes import GeoAxes

    # Map projection
    projection = ccrs.PlateCarree(central_longitude=cen_lon)

    # Axis class
    axes_class = (GeoAxes, dict(map_projection=projection))

    # Plot figure
    fig = plt.figure(figsize=(8, 6.5))
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
        ax.add_feature(cartopy.feature.LAND, zorder=1, facecolor='none',edgecolor='grey')
        # ax.add_feature(cartopy.feature.GSHHSFeature('low',edgecolor='none'), zorder=1, facecolor='black')
        ax.gridlines(color='black', linestyle=':')
        ax.set_title(ti_str)
        ax.set_extent([-180, 180, -90, 90], projection)
        ax.text(0.5, -0.15, varname + ' \n' + varunits, va='bottom', ha='center',
                rotation='horizontal', rotation_mode='anchor', transform=ax.transAxes, fontsize=10)

    # Use a different colorbar range every time by specifying tuple of tuples
    limits=cmap_lims
    for i in range(1):
        im = axgr[i].pcolor(lon, lat, field,transform=ccrs.PlateCarree(),cmap=cmap,
                              vmin=limits[0], vmax=limits[1])
        axgr.cbar_axes[i].colorbar(im)

    for i, cax in enumerate(axgr.cbar_axes):
        cax.set_yticks((limits[0], limits[1]))

    plt.show()
    return

