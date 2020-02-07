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
print('\t\t\tFLASHFLux official file...\t\t\t')
print('============================================')


path1 = '/Users/rcscott2/Desktop/CERES/FLASHFlux/'
file1 = 'FLASH_SSF_Aqua-FM3-MODIS_Version3C_232103.2020020523'
file_path1 = path1 + file1
print(file_path1)


print('============================================')
print('\t\t\tSSF official file...\t\t\t')
print('============================================')


file2 = 'CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.2010062023'
path2 = '/Users/rcscott2/Desktop/CRS/output/'
file_path2 = path2 + file2
print(file_path2)


print('============================================')
print('\t\t\tDefining function...\t\t\t')
print('============================================')


def read_ssf_geolocation(file_path):
    """Reads geolocation info from file"""

    from pyhdf import SD
    hdf = SD.SD(file_path)
    # print(hdf.datasets())

    colatitude = hdf.select('Colatitude of CERES FOV at surface')
    colat = colatitude.get()
    fov_lat = 90 - colat

    longitude = hdf.select('Longitude of CERES FOV at surface')
    fov_lon = longitude.get()

    time_of_obs = hdf.select("Time of observation")
    time_obs = time_of_obs.get()

    return fov_lat, fov_lon, time_obs


print('============================================')
print('\t\t\tDefining function...\t\t\t')
print('============================================')


def read_ssf_var(file_path, vararg, fill):
    """
    Read CERES SSF footprint variable

    :param file_path: path to file
    :param vararg: variable argument
    :param fill: fill option
    :return: field, name, units
    """

    from pyhdf import SD
    hdf = SD.SD(file_path)

    # select variable
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
    var_fill = data._FillValue

    if fill == 1:
        variable[variable == var_fill] = np.nan

    # get field at appropriate vertical level
    var_field = variable

    return var_field, var_name, var_units


print('============================================')
print('\t\t\tDefining function...\t\t\t')
print('============================================')


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


print('============================================')
print('\t\t\tDefining function...\t\t\t')
print('============================================')


def select_var2(file_path, vararg, levarg, fill):
    """
    Select field variable at desired vertical level from
    official release CERES CRS HDF file.

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


print('============================================')
print('\t\t\tDefining function...\t\t\t')
print('============================================')


def compute_diff(field2, field1):

    diff = field2 - field1

    diff[diff == max(diff)] = np.nan
    diff[diff == min(diff)] = np.nan

    return diff


print('============================================')
print('\t\t\tDefining function...\t\t\t')
print('============================================')


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
    print("\tConstructing", switch.get(typarg, "N/A"), "colormap")

    if typarg == 0:
        c_map = cmap_name.mpl_colormap
    elif typarg == 1:
        from matplotlib.colors import ListedColormap
        c_map = ListedColormap(cmap_name.mpl_colors)
    elif typarg != 1 or typarg != 0:
        print('Please select 0 or 1')

    return c_map


print('============================================')
print('\t\t\tDefining function...\t\t\t')
print('============================================')


def plot_ssf_swath(nrows, ncols, cen_lon,
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

    # Creat figure
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
        ax.set_title(title_str + ' ' + date_str, fontsize=10)
        ax.set_extent([-180, 180, -90, 90], projection)
        ax.text(0.5, -0.1, varname + ' - ' + ' \n' + varunits, va='bottom', ha='center',
                rotation='horizontal', rotation_mode='anchor', transform=ax.transAxes, fontsize=10)

        if nightshade == 1:
            ax.add_feature(Nightshade(date, alpha=0.15))

    # To use a different colorbar range every time, use a tuple
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


print('============================================')
print('\t\t\tReading geolocation...\t\t\t')
print('============================================')

lat, lon, time_obs = read_ssf_geolocation(file_path=file_path1)

print('============================================')
print('\t\t\tReading data...\t\t\t')
print('============================================')

field1, var1, units1 = read_ssf_var(file_path=file_path1, vararg=3, fill=1)

print('============================================')
print('\t\t\tReading time/date info...\t\t\t')
print('============================================')

date, date_str = get_date(file=file1)

print('============================================')
print('\t\t\tComputing difference...\t\t\t')
print('============================================')

#difference = compute_diff(field2=field2, field1=field1)

print('============================================')
print('\t\t\tSetting colormap...\t\t\t')
print('============================================')

colormap = set_colormap(cmap_name=Balance_20, typarg=0)

print('============================================')
print('\t\t\tPlotting data...\t\t\t')
print('============================================')

a, b = input("Enter colormap limits: ").split(',')
print("First number is {} and second number is {}".format(a, b))
cmap_lim = (float(a), float(b))

title_str = 'CERES FLASHFlux Edition 3C' + ' -'

# title_str = r'Difference ($\Delta$) between CERES Terra CRS Ed4 - Ed2G' + ' for ' + date_str

plot_ssf_swath(nrows=1, ncols=1, cen_lon=0, field=field1,
               varname=var1, levname="", varunits=units1,
               cmap=colormap, cmap_lims=cmap_lim, date=date,
               nightshade=1, title_str=title_str)

