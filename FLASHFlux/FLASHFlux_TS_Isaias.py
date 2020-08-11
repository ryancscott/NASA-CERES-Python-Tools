# ==========================================================================
# This script grids surface & TOA radiative flux retrievals from FLASHFlux
# as tropical storm Isaias approached and made landfall over North Carolina
# during early August 2020.
#
# Author: Ryan Scott, SSAI
#         ryan.c.scott@nasa.gov
# ==========================================================================

import sys
import numpy as np
import cerestools as ceres
import matplotlib.pyplot as plt

from palettable.cmocean.diverging import Balance_20
from palettable.scientific.sequential import LaJolla_20_r
from palettable.cmocean.sequential import Deep_20_r
from palettable.cmocean.sequential import Thermal_20

from palettable.scientific.sequential import Nuuk_3

# ==================================================================
# Read entire day of swaths
# ==================================================================
# terra_var_all, terra_lon_all, terra_lat_all, terra_sza_all = \
#     ceres.read_day_of_ssf_files(
#         path='/Users/rcscott2/Desktop/FLASHFlux/ASDC_archive/SSF/TS_Isaias/Terra/',
#         file_struc='FLASH_SSF_Terra-FM1-MODIS_Version3C_232103.20200802',
#         variable='CERES LW TOA flux - upwards',
#         index=-1,
#         fill=True
#     )
#
# terra_var_all, terra_lon_all, terra_lat_all, terra_sza_all = \
#     ceres.swath_nighttime_only(
#         var=terra_var_all,
#         lon=terra_lon_all,
#         lat=terra_lat_all,
#         sza=terra_sza_all,
#         sza_cutoff=90
#     )
# ==================================================================
# Plot a single swath
# ==================================================================
# mo = '08'
# day = '02'
# hr = '16'
#
# path = '/Users/rcscott2/Desktop/FLASHFlux/ASDC_archive/SSF/Terra/'
# file = 'FLASH_SSF_Terra-FM1-MODIS_Version3C_232103.2020' + mo + day + hr
# filepath = path + file
#
# terra_lat_all, terra_lon_all, _, terra_sza_all = \
#     ceres.read_ssf_geolocation(file_path=filepath)
#
# terra_var_all, _, _ = \
#     ceres.read_ssf_var(
#         file_path=filepath,
#         var_name='CERES LW TOA flux - upwards',
#         index=-1,
#         fill=True)
#
# gridded_field, grid_lon, grid_lat = ceres.grid_to_equal_angle_grid(
#     grid_res_lat=0.5,
#     grid_res_lon=0.5,
#     variable=terra_var_all,
#     lon_data=terra_lon_all,
#     lat_data=terra_lat_all,
#     lon_360=True)
#
# date_str = mo + '/'+day+'/2020:'+hr+'h'
#
# # set colormap
# cmap = ceres.set_colormap(Deep_20_r, 0)
# # cmap = 'viridis'
#
# # gridded_field[gridded_field > 180] = np.nan
# #
# ceres.plot_gridded_fields(nrows=1, ncols=1, cen_lon=0,
#                           date_str='',
#                           title_str=['FLASHFlux v3C Terra FM1 - ' + date_str],
#                           cmap=cmap, cmap_lims=((130, 330), (130, 330)),
#                           varname='Outgoing Longwave Radiation', levname='TOA', varunits=r'W m$^{-2}$',
#                           lon=grid_lon, lat=grid_lat,
#                           field=gridded_field[:, :],
#                           itr = 1)
#
# sys.exit()

# ==================================================================
# Terra FM1 OLR
# ==================================================================
# make all necessary plots
# Terra
# m = ['07', '07', '08', '08', '08', '08', '08', '08', '08', '08']
# d = ['31', '31', '01', '01', '02', '02', '03', '03', '04', '04']
# h = ['02', '15', '03', '15', '02', '16', '03', '15', '02', '16']
#
#
# for j, (mo, day, hr) in enumerate(list(zip(m, d, h))):
#
#     path = '/Users/rcscott2/Desktop/FLASHFlux/ASDC_archive/SSF/Terra/'
#     file = 'FLASH_SSF_Terra-FM1-MODIS_Version3C_232103.2020' + mo + day + hr
#     filepath = path + file
#
#     terra_lat_all, terra_lon_all, _, terra_sza_all = \
#         ceres.read_ssf_geolocation(file_path=filepath)
#
#     terra_var_all, _, _ = \
#         ceres.read_ssf_var(
#             file_path=filepath,
#             var_name='CERES LW TOA flux - upwards',
#             index=-1,
#             fill=True)
#
#     gridded_field, grid_lon, grid_lat = \
#         ceres.grid_to_equal_angle_grid(
#             grid_res_lat=0.5,
#             grid_res_lon=0.5,
#             variable=terra_var_all,
#             lon_data=terra_lon_all,
#             lat_data=terra_lat_all,
#             lon_360=True)
#
#     # Plot fluxes below a certain threshold
#     # gridded_field[gridded_field > 180] = np.nan
#
#     date_str = mo + '/' + day + '/2020:' + hr + 'h'
#
#     # set colormap
#     cmap = ceres.set_colormap(Deep_20_r, 0)
#     # cmap = ceres.set_colormap(Nuuk_3, 0)
#     # cmap = 'viridis'
#     #
#     ceres.plot_gridded_fields(
#         nrows=1, ncols=1, cen_lon=0,
#         date_str='',
#         title_str=['FLASHFlux v3C Terra FM1 - ' + date_str],
#         cmap=cmap, cmap_lims=((130, 320), (130, 320)),
#         varname='Outgoing Longwave Radiation', levname='TOA', varunits=r'W m$^{-2}$',
#         lon=grid_lon, lat=grid_lat,
#         field=gridded_field[:, :],
#         itr=j)

#
# Aqua
# m = ['07', '08', '08', '08']
# d = ['31', '01', '02', '03']
# h = ['18', '18', '18', '18']

m = ['08']
d = ['03']
h = ['18']

for j, (mo, day, hr) in enumerate(list(zip(m, d, h))):

    path = '/Users/rcscott2/Desktop/FLASHFlux/ASDC_archive/SSF/Aqua/'
    file = 'FLASH_SSF_Aqua-FM3-MODIS_Version3C_232103.2020' + mo + day + hr
    filepath = path + file

    terra_lat_all, terra_lon_all, _, terra_sza_all = \
        ceres.read_ssf_geolocation(file_path=filepath)

    terra_var_all, _, _ = \
        ceres.read_ssf_var(
            file_path=filepath,
            var_name='CERES SW TOA flux - upwards',
            index=-1,
            fill=True)

    gridded_field, grid_lon, grid_lat = \
        ceres.grid_to_equal_angle_grid(
            grid_res_lat=0.5,
            grid_res_lon=0.5,
            variable=terra_var_all,
            lon_data=terra_lon_all,
            lat_data=terra_lat_all,
            lon_360=True)

    # Plot fluxes below a certain threshold
    # gridded_field[gridded_field > 180] = np.nan

    date_str = mo + '/' + day + '/2020:' + hr + 'h'

    # set colormap
    cmap = ceres.set_colormap(Deep_20_r, 0)
    # cmap = 'viridis'
    #
    ceres.plot_gridded_fields(
        nrows=1, ncols=1, cen_lon=0,
        date_str='',
        title_str=['FLASHFlux v3C Aqua FM3 - ' + date_str],
        cmap=cmap, cmap_lims=((0, 700), (0, 700)),
        varname='Reflected Shortwave Radiation', levname='TOA', varunits=r'W m$^{-2}$',
        lon=grid_lon, lat=grid_lat,
        field=gridded_field[:, :],
        itr=j)

