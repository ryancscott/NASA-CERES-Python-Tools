import math
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
import datetime
from netCDF4 import Dataset
import cartopy
import cartopy.crs as ccrs
import cartopy.feature
from cartopy.mpl.geoaxes import GeoAxes
#from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from mpl_toolkits.axes_grid1 import AxesGrid
from palettable.cmocean.diverging import Balance_19
from palettable.cmocean.sequential import Thermal_20
from palettable.cmocean.sequential import Ice_20
from palettable.mycarta import Cube1_20
from palettable.cubehelix import cubehelix2_16


print('====================================')
print('\t\t\tEBAF file...\t\t\t')
print('====================================')


path = '/Users/rcscott2/Desktop/CERES/EBAF/'
file = 'CERES_EBAF-TOA_Ed4.1_Subset_200003-201910.nc'

file_path = path + file


# ========================================================================


def print_nc_file_info(filepath):
    """
    Function prints out information about variables
    in the input netCDF file
    :param filepath: path to netcdf file
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

    print('\n')
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

    print('field shape...')
    print(var.shape, '\n')

    return var, var_name, var_units


# ========================================================================


def compute_monthly_anomalies(field):
    """
    Compute average for each calendar month, i.e., the monthly mean seasonal cycle,
    and then subtract the long-term means from the raw fields to get monthly anomalies
    :param field: variable to compute the climatology/anomaly for
    :return:
    (1) time series of monthly anomalies
    (2) monthly mean seasonal cycle
    """

    # total length of record in file
    record_len = field.shape[0]

    # initialize array to hold monthly climatology for each month
    monthly_climatology = np.empty((12, 180, 360))

    print('Indices for monthly climatology:\n')
    # loop over each calendar month and average
    for i in range(0, 12):

        # select appropriate indices associated with each month
        indices = list(range(i,record_len, 12))

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


# def get_date(file):
#     """
#     Get date/time info from input data file for
#     Nightshade and displaying it on plots
#     :param file: input file
#     :return: date object and time string
#     """
#
#     time_str = file[-10:]
#     print(time_str)
#     yyyy = time_str[0:4]
#     mm = time_str[4:6]
#     dd = time_str[6:8]
#     hr = time_str[8:10]
#
#     date = datetime.datetime(int(yyyy), int(mm), int(dd), int(hr))
#
#     date_str = mm + '/' + dd + '/' + yyyy + ':' + hr + 'h'
#
#     return date, date_str


# ========================================================================


def cos_lat_weight(lat_vector):
    """
    Returns cos(lat) weight array [size nlat x nlon]
    :param lat_vector: vector of latitudes
    :return: matrix of cos(lat) weights
    """

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

    plt.imshow(weight)

    return weight


# ========================================================================


def compute_climatology(field):
    """
    Compute mean and standard deviation of field
    :param field: field under consideration
    :return: the mean and standard deviation
    """

    # average over entire time series
    mean_field = np.nanmean(field, axis=0)

    # standard deviation over entire time series
    stdv_field = np.nanstd(field, axis=0)

    return mean_field, stdv_field


# ========================================================================


def compute_regional_averages(field):
    """
    Computes weighted regional averages of input field
    :param field: field to be averaged
    :return: prints regional averages
    """

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
    mean_field1 = np.nanmean(field[ind1, :, :], axis=0)
    print('Size of first section to average: ', field[ind1, :, :].shape)

    mean_field2 = np.nanmean(field[ind2, :, :], axis=0)
    print('Size of second section to average: ', field[ind2, :, :].shape)

    composite_diff = np.subtract(mean_field2, mean_field1)

    return composite_diff


# ========================================================================


def set_colormap(cmap_name, argument):
    """Sets colormap to be used for mapping geospatial data. Inputs include the
     colormap name from palettable, 0/1 for continuous/discrete colormap, and
     an optional number of colors for the discrete case"""
    switch = {
        0: "continuous",
        1: "discrete"
    }
    print("\nUsing", switch.get(argument, "NO"), "colormap")

    if argument == 0:
        color_map = cmap_name.mpl_colormap
    elif argument == 1:
        from matplotlib.colors import ListedColormap
        color_map = ListedColormap(cmap_name.mpl_colors, N = num_colors)

    return color_map


# ========================================================================


def global_map(nrows, ncols, cen_lon,
                field, varname, varunits,
                cmap, cmap_lims, ti_str):

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


#====================================

# print information about the files
print_nc_file_info(file_path)

# read latitude and longitude information from file
latitude, _, lat, lon = read_ebaf_geolocation(filepath=file_path)

# compute weights for area averaging
w = cos_lat_weight(latitude)

# read variable including its name and units
field, name, units = read_ebaf_var(filepath=file_path, variable='cldarea_total_daynight_mon')

# compute the long-term mean and standard deviation
mean_field, stdv_field = compute_climatology(field)

# compute area weighted averages of the long-term mean
compute_regional_averages(mean_field)

# compute anomalies by removing long-term monthly means
monthly_anomalies, seasonal_cycle = compute_monthly_anomalies(field)

# set the colormap
cmap = set_colormap(Balance_19, 0)

# plot global map
global_map(nrows=1, ncols=1, cen_lon=180, field=monthly_anomalies[100, :, :],
           varname=name, varunits=units, cmap=cmap,
           cmap_lims=(-30, 30), ti_str='CERES-EBAF Ed4.1')

