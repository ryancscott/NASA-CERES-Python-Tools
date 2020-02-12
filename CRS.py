import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import datetime
import cartopy
import cartopy.crs as ccrs
import cartopy.feature
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.feature.nightshade import Nightshade
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from palettable.cmocean.sequential import Thermal_20
from palettable.cubehelix import cubehelix2_16
from palettable.cmocean.diverging import Balance_20
from palettable.mycarta import Cube1_20


print('============================================')
print('\t\t\tCRS official file...\t\t\t')
print('============================================')


path1 = '/Users/rcscott2/Desktop/CERES/CRS/'
file1 = 'CER_CRS_Terra-FM1-MODIS_Edition2G_023034.2010062023'
file_path1 = path1 + file1
print(file_path1)


print('============================================')
print('\t\t\tCRS development file...\t\t\t')
print('============================================')


file2 = 'CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.2010062023'
path2 = '/Users/rcscott2/Desktop/CRS/output/'
file_path2 = path2 + file2
print(file_path2)


# =======================================================


def read_crs_geolocation(file_path):
    """
    Reads geolocation information from file
    :param file_path: file to be read
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


# =======================================================


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
        1: 'Shortwave flux - downward - total',
        2: 'Longwave flux - upward - total',
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

    print("Getting", switch.get(vararg, "N/A"), "at:", switch2.get(levarg, "N/A"))

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


# =======================================================


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


# =======================================================


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
        1: 'UT_TOT_SW_DN',
        2: 'UT_TOT_LW_UP',
        3: 'UT_TOT_SW_UP'
    }
    var_name = switch.get(vararg)

    # select level
    switch2 = {
        0: 'TOA',
        1: 'surface'
    }
    lev_name = switch.get(levarg)

    print("Getting", switch.get(vararg, "N/A"), "at:", switch2.get(levarg, "N/A"))

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


# =======================================================


def compute_diff(field2, field1):

    diff = field2 - field1

    diff[diff == max(diff)] = np.nan
    diff[diff == min(diff)] = np.nan

    # num_good = np.sum(np.isnan(diff))
    # print(num_good)

    return diff


# =======================================================


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


# =======================================================


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
    fig = plt.figure(figsize=(10, 6.5))
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
        ax.add_feature(cartopy.feature.LAND, zorder=1, facecolor='none', edgecolor='darkgrey')
        ax.gridlines(color='grey', linestyle='--')
        ax.set_title(title_str + ' - ' + date_str, fontsize=10)
        ax.set_extent([-180, 180, -90, 90], projection)
        ax.text(0.5, -0.1, varname + ' - ' + levname + ' \n' + varunits, va='bottom', ha='center',
                rotation='horizontal', rotation_mode='anchor', transform=ax.transAxes, fontsize=10)

        if nightshade == 1:
            ax.add_feature(Nightshade(date, alpha=0.15))

    # To use a different colorbar range every time, use a tuple of tuples
    limits = cmap_lims
    for i in range(nrows * ncols):
        im = axgr[i].scatter(lon, lat, c=field, s=1, transform=ccrs.PlateCarree(),
                             vmin=limits[0], vmax=limits[1], cmap=cmap)
        axgr.cbar_axes[i].colorbar(im)

    for i, cax in enumerate(axgr.cbar_axes):
        cax.set_yticks(np.linspace(limits[0], limits[1], 5))
        cax.set_yticklabels(np.linspace(limits[0], limits[1], 5), fontsize=8)

    plt.show()
    return


# ========================================================================


def histogram_scatterplot(difference, varname, levname, time_info):

    # Ed 2G
    field1, var1, units1, lev1 = read_crs_var(file_path=file_path1, vararg=0, levarg=4, fill=1)
    # Ed 4
    field2, var2, units2, lev2 = read_crs_var2(file_path=file_path2, vararg=0, levarg=1, fill=1)

    difference = compute_diff(field2=field2, field1=field1)

    print(np.nanmean(difference))

    # number of FOVs compared
    print("Total number of FOVs in swath: ", len(difference))
    print("Total number of FOVs where at least one is NaN: ", np.sum(np.isnan(difference)))
    N = len(difference) - np.sum(np.isnan(difference))
    print("Total number of FOVs compared: ", N)

    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    fig.suptitle('CERES Terra Cloud Radiative Swath' + ' - ' + time_info, fontsize=11)

    # histogram of difference
    axs[0].hist(difference, bins=200, align='mid', rwidth=1)
    axs[0].grid()
    axs[0].set_axisbelow("True")
    axs[0].set_title("CRS Ed4 minus CRS Ed2G", fontsize=11)
    axs[0].set_xlabel(varname + ' - ' + levname + "\n" + r"Flux difference ($\Delta$LW)")
    axs[0].set_xlim([-125, 125])
    axs[0].set_ylabel("Number of CERES FOVs, Total: " + str(N))

    # scatter plot
    axs[1].plot(range(500), range(500), color='black', zorder=0)
    axs[1].scatter(field1, field2, s=0.5)
    axs[1].set(xlim=(0, 500), ylim=(0, 500))
    axs[1].grid()
    axs[1].set_axisbelow("True")
    axs[1].set_title(varname + ' - ' + levname, fontsize=11)
    axs[1].set_xlabel("CRS Edition 2G")
    axs[1].set_ylabel("CRS Edition 4")

    plt.show()
    return


# =========================================================================


# def scatterplot(field1, field2):
#     plt.scatter(field1, field2)
#     plt.
#     plt.grid(b=1, zorder=0)
#     plt.title("LW Flux", fontsize=11)
#     plt.xlabel("CRS Ed2G")
#     plt.xlim([0, 500])
#     plt.ylabel("CRS Ed4")
#     plt.show()


print('============================================')
print('\t\t\tReading geolocation...\t\t\t')
print('============================================')

lat, lon, p_levels, time_obs = read_crs_geolocation(file_path=file_path1)

print('============================================')
print('\t\t\tReading data...\t\t\t')
print('============================================')

field1, var1, units1, lev1 = read_crs_var(file_path=file_path1, vararg=0, levarg=4, fill=1)
field2, var2, units2, lev2 = read_crs_var2(file_path=file_path2, vararg=0, levarg=1, fill=1)

print('============================================')
print('\t\t\tReading time/date info...\t\t\t')
print('============================================')

date, date_str = get_date(file=file1)

print('============================================')
print('\t\t\tComputing difference...\t\t\t')
print('============================================')

difference = compute_diff(field2=field2, field1=field1)

print('============================================')
print('\t\t\tSetting colormap...\t\t\t')
print('============================================')

colormap = set_colormap(cmap_name=Balance_20, typarg=0)

print('============================================')
print('\t\t\tPlotting data...\t\t\t')
print('============================================')

print("Plotting: ", field1)
a, b = input("Enter colormap limits: ").split(',')
print("Specified colormap range [{}, {}]".format(a, b))
cmap_lim = (float(a), float(b))

#title_str = 'CERES Terra CRS Ed4'
title_str = r'CERES Terra CRS Ed4 - Ed2G difference ($\Delta$)'

plot_swath(nrows=1, ncols=1, cen_lon=0, field=difference,
           varname=var1, levname=lev1, varunits=units1,
           cmap=colormap, cmap_lims=(-50,50), date=date,
           nightshade=1, title_str=title_str)

histogram_scatterplot(difference, var1, lev1, date_str)
