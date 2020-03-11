import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import datetime
import cartopy
import cartopy.crs as ccrs
import shapely.geometry as sgeom
import cartopy.feature
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.feature.nightshade import Nightshade
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from palettable.cmocean.sequential import Thermal_20
from palettable.cubehelix import cubehelix2_16
from palettable.cmocean.diverging import Balance_19
from palettable.colorbrewer.sequential import Greys_8
from palettable.mycarta import Cube1_20


print('====================================')
print('\t\t\tMODIS file...\t\t\t')
print('====================================')


path = '/Users/rcscott2/Desktop/CERES/MODIS/'
file = 'CLDPROP_L2_MODIS_Aqua.A2019359.2140.011.2019360180828.nc'

file_path = path + file


# ========================================================================


def print_nc_file_info_groups(filepath):
    """
    This function prints out info about the variables
    within the input netCDF file
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
    geolocation_group = nc.groups['geolocation_data']
    geophysical_group = nc.groups['geophysical_data']

    # print information about variables in netCDF file
    print("\nGeolocation:\n")
    for i in geolocation_group.variables:
        print(i, '-', geolocation_group.variables[i].units, '-', geolocation_group.variables[i].shape)

    # print information about variables in netCDF file
    print("\nGeophysical:\n")
    for i in geophysical_group.variables:
        print(i, '-', geophysical_group.variables[i].units, '-', geophysical_group.variables[i].shape)

    return


# =============================================================================


def read_modis_geoloc(filepath):
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
    geolocation_group = nc.groups['geolocation_data']

    lat = geolocation_group.variables['latitude'][:, :]
    lon = geolocation_group.variables['longitude'][:, :]

    #lon[lon<0] = lon[lon<0]+360

    return lat, lon


# =============================================================================


def read_modis_var(filepath, variable):
    """"""

    from netCDF4 import Dataset

    nc = Dataset(filepath, 'r')
    geophysical_group = nc.groups['geophysical_data']

    var = geophysical_group.variables[variable]

    print("===========================")
    print("\t\t Variable:")
    print("===========================")

    var_name = str(variable).split("_")
    var_name = " ".join(var_name)

    print("Name: ", var_name)

    var_sf = var.scale_factor
    print("Scale factor: ", var_sf)

    var_units = var.units
    print("Units: ", var_units)

    var_vmin = var.valid_min
    print("Valid min: ", var_vmin)

    var_vmax = var.valid_max
    print("Valid max: ", var_vmax)

    var_fill = var._FillValue
    print("Fill value: ", var_fill)

    # get the actual data and fill missing values
    var = var[:, :]
    var[var == var_fill] = np.nan

    return var, var_name, var_units, var_vmin, var_vmax

# =============================================================================


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


def modis_swath_map(nrows, ncols, field, varname, varunits, cmap, cmap_lims, dlatlon, title, inset):

    print('\nPlot limits:')
    min_lon = np.min(lon)-dlatlon
    max_lon = np.max(lon)+dlatlon
    min_lat = np.min(lat)-dlatlon
    max_lat = np.max(lat)+dlatlon
    print("Minimum longitude: ", min_lon)
    print("Maximum longitude: ", max_lon)
    print("Minimum latitude: ", min_lat)
    print("Maximum latitude: ", max_lat)

    cen_lon = (min_lon+max_lon)/2
    print("Central longitude: ", cen_lon)

    # map projection
    projection = ccrs.PlateCarree()

    # axis class
    axes_class = (GeoAxes, dict(map_projection=projection))

    # plot figure
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

    stamen_terrain = cimgt.Stamen('terrain-background')

    for i, ax in enumerate(axgr):
        ax.add_feature(cartopy.feature.GSHHSFeature('intermediate', edgecolor='black'), zorder=0, facecolor='none')
        #ax.stock_img()
        ax.add_image(stamen_terrain, 8)
        ax.gridlines(color='black', linestyle=':')
        ax.set_title(title, fontsize=10)
        ax.set_extent([min_lon, max_lon, min_lat, max_lat], projection)
        ax.text(0.5, -0.1, varname + ' \n' + "Units: " + varunits, va='bottom', ha='center',
                rotation='horizontal', rotation_mode='anchor', transform=ax.transAxes, fontsize=10)

    limits = cmap_lims
    for i in range(1):
        im = axgr[i].pcolor(lon, lat, field, transform=projection, cmap=cmap, vmin=limits[0], vmax=limits[1])
        axgr.cbar_axes[i].colorbar(im)

    # colorbar axis ticks and labels
    for i, cax in enumerate(axgr.cbar_axes):
        cax.set_yticks(np.linspace(np.round(limits[0]), np.round(limits[1]), 5))
        cax.set_yticklabels(np.linspace(np.round(limits[0]), np.round(limits[1]), 5), fontsize=8)

    if inset == 1:
        # create an inset GeoAxes showing the location of the MODIS swath
        geodetic = ccrs.Geodetic(globe=ccrs.Globe(datum='WGS84'))
        sub_ax = fig.add_axes([0.65, 0.65, 0.18, 0.18], projection=projection)
        sub_ax.set_extent([min_lon-30, max_lon+30, min_lat-20, max_lat+20], geodetic)

        # add land and extent box to inset
        sub_ax.add_feature(cartopy.feature.GSHHSFeature('intermediate', edgecolor='none'), zorder=0, facecolor='black')
        extent = [min_lon, max_lon, min_lat, max_lat]
        extent_box = sgeom.box(extent[0], extent[2], extent[1], extent[3])
        sub_ax.add_geometries([extent_box], ccrs.PlateCarree(), facecolor='none', edgecolor='blue', linewidth=1)

    print('\nPlotting map of the data...')

    plt.show()

    return


# ====================================


# print information about file contents
print_nc_file_info_groups(file_path)

# read geo-location data
lat, lon = read_modis_geoloc(file_path)

# read desired field and attributes
data, name, units, valid_min, valid_max = read_modis_var(file_path, "Cloud_Optical_Thickness")

# set the colormap
cmap = set_colormap(Balance_19, 0)

# plot regional map with inset
modis_swath_map(nrows=1, ncols=1, field=data,
                varname=name, varunits=units, cmap=cmap,
                cmap_lims=(np.min(data), np.max(data)),
                dlatlon=5, title=file, inset=1)
