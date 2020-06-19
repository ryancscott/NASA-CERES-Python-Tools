# ==========================================================================
# Author: Ryan C. Scott
#         ryan.c.scott@nasa.gov
#
# This script compares CRS1deg_beta to SYN1deg-Hour at the Terra/Aqua
# overpass time.
# ==========================================================================

import numpy as np
import cerestools as ceres
import matplotlib.pyplot as plt

from palettable.cmocean.diverging import Balance_20
# from palettable.scientific.sequential import LaJolla_20_r
# from palettable.cmocean.sequential import Thermal_20


path2 = '/Users/rcscott2/Desktop/CERES/SYN1deg/'
file2 = 'CER_SYN1deg-1Hour_Terra-Aqua-MODIS_Edition4A_407406.20190101'
file_path2 = path2 + file2


var_syn1deg, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
                                           var_name='obs_clr_toa_lw',
                                           fill=False)

lat_syn1deg, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
                                           var_name='latitude',
                                           fill=False)

lon_syn1deg, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
                                           var_name='longitude',
                                           fill=False)

num_sw_obs, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
                                          var_name='num_sw_obs',
                                          fill=False)

num_lw_obs, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
                                          var_name='num_lw_obs',
                                          fill=False)

num_geo_sw_obs, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
                                              var_name='num_geo_sw_obs',
                                              fill=False)

num_geo_lw_obs, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
                                              var_name='num_geo_lw_obs',
                                              fill=False)


print("SYN1deg Terra-Aqua SW obs shape: ", num_sw_obs.shape)
print("SYN1deg GEO SW obs shape: ", num_geo_sw_obs.shape)
print("SYN1deg Terra-AquaLW obs shape: ", num_lw_obs.shape)
print("SYN1deg GEO LW obs shape: ", num_geo_lw_obs.shape)


terra_var_gridded = np.zeros([24, 180, 360])
aqua_var_gridded = np.zeros([24, 180, 360])

terra_mask = np.zeros([24, 180, 360])
aqua_mask = np.zeros([24, 180, 360])
both_mask = np.zeros([24, 180, 360])


for k in range(24):

    if k < 10:
        i = '0' + str(k)
    elif k >= 10:
        i = str(k)

    terra_crs_file = '/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/' \
                     'CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20190101' + str(i)
    aqua_crs_file = '/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/' \
                    'CER_CRS4_Aqua-FM3-MODIS_GH4_1111TH.20190101' + str(i)

    date, date_str = ceres.get_date(terra_crs_file)

    print(terra_crs_file)
    print(aqua_crs_file)

    terra_lat_all, terra_lon_all, _, _, _, terra_sza_all = \
        ceres.read_crs_geolocation_dev(file_path=terra_crs_file)

    terra_var_all, _, _, _ = \
        ceres.read_crs_var_dev(file_path=terra_crs_file,
                               var_name='Longwave flux - upward - total sky',
                               lev_arg=0,
                               fill=True)

    aqua_lat_all, aqua_lon_all, _, _, _, aqua_sza_all = \
        ceres.read_crs_geolocation_dev(file_path=aqua_crs_file)

    aqua_var_all, _, _, _ = \
        ceres.read_crs_var_dev(file_path=aqua_crs_file,
                               var_name='Longwave flux - upward - total sky',
                               lev_arg=0,
                               fill=True)

    terra_lon_all = ceres.swath_lon_360_to_180(terra_lon_all)
    aqua_lon_all = ceres.swath_lon_360_to_180(aqua_lon_all)

    terra_lat_all, terra_lon_all, terra_var_all, terra_sza_all = \
        ceres.swath_daytime_only(lat=terra_lat_all,
                                 lon=terra_lon_all,
                                 var=terra_var_all,
                                 sza=terra_sza_all,
                                 sza_cutoff=90)

    aqua_lat_all, aqua_lon_all, aqua_var_all, aqua_sza_all = \
        ceres.swath_daytime_only(lat=aqua_lat_all,
                                 lon=aqua_lon_all,
                                 var=aqua_var_all,
                                 sza=aqua_sza_all,
                                 sza_cutoff=90)

    terra_var_gridded[k, :, :] = \
        ceres.grid_to_1x1_deg_equal_angle(lat_data=terra_lat_all,
                                          lon_data=terra_lon_all,
                                          variable=terra_var_all,
                                          lon_360=False)

    aqua_var_gridded[k, :, :] = \
        ceres.grid_to_1x1_deg_equal_angle(lat_data=aqua_lat_all,
                                          lon_data=aqua_lon_all,
                                          variable=aqua_var_all,
                                          lon_360=False)
    print('...', k, ' ...')

    for i in range(180):
        for j in range(360):
            if terra_var_gridded[k, i, j] > 0:
                terra_mask[k, i, j] = 1
            if aqua_var_gridded[k, i, j] > 0:
                aqua_mask[k, i, j] = 1
            if terra_var_gridded[k, i, j] > 0 and aqua_var_gridded[k, i, j] > 0:
                both_mask[k, i, j] = 1


# isolate grid boxes observed by one satellite
terra_only_mask = terra_mask - both_mask
aqua_only_mask = aqua_mask - both_mask
terra_aqua_mask = terra_mask + aqua_mask


for k in range(0):

    if k < 10:
        j = '0' + str(k)
    elif k >= 10:
        j = str(k)

    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(11, 7))
    for i, ax in enumerate(axes.flat):
        if i == 0:
            im = ax.imshow(terra_mask[k, :, :], vmin=0, vmax=1)
            ax.set_title(r'CRS1deg$_{\beta}$ Terra FM1, 1-1-2019:' + str(j) + 'h')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
        elif i == 1:
            im = ax.imshow(aqua_only_mask[k, :, :], vmin=0, vmax=1)
            ax.set_title(r'CRS1deg$_{\beta}$ Aqua - (Terra$\bigcap$Aqua), 1-1-2019:' + str(j) + 'h')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
        elif i == 2:
            im = ax.imshow(aqua_mask[k, :, :], vmin=0, vmax=1)
            ax.set_title(r'CRS1deg$_{\beta}$ Aqua FM3, 1-1-2019:' + str(j) + 'h')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
        elif i == 3:
            im = ax.imshow(num_sw_obs[k, :, :], vmin=0, vmax=1)
            ax.set_title('SYN1deg Terra & Aqua SW mask, 1-1-2019:' + str(j) + 'h')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
        elif i == 4:
            im = ax.imshow(terra_aqua_mask[k, :, :], vmin=0, vmax=2)
            ax.set_title(r'CRS1deg$_{\beta}$ Terra & Aqua, 1-1-2019:' + str(j) + 'h')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
        elif i == 5:
            im = ax.imshow(num_lw_obs[k, :, :], vmin=0, vmax=1)
            ax.set_title('SYN1deg Terra & Aqua LW mask, 1-1-2019:' + str(j) + 'h')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])

    fig.subplots_adjust(right=0.8)
    cb_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cb_ax)

    plt.show()


# take the difference CRS1deg_beta minus SYN1deg hour-by-hour
aqua_diff = aqua_only_mask*aqua_var_gridded - aqua_only_mask*var_syn1deg
terra_diff = terra_only_mask*terra_var_gridded - terra_only_mask*var_syn1deg

# ignore grid boxes observed by Terra and Aqua
aqua_diff[terra_aqua_mask == 2] = np.nan
terra_diff[terra_aqua_mask == 2] = np.nan

#
fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.hist(np.reshape(terra_diff, 24*180*360), bins=200, align='mid', rwidth=1)
ax2.hist(np.reshape(aqua_diff, 24*180*360), bins=200, align='mid', rwidth=1)
plt.show()


for k in range(0):

    if k < 10:
        j = '0' + str(k)
    elif k >= 10:
        j = str(k)
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(6, 7))
    for i, ax in enumerate(axes.flat):
        if i == 0:
            im = ax.imshow(terra_diff[k, :, :], vmin=-30, vmax=30)
            ax.set_title(r'Terra Only CRS1deg$_{\beta}$ minus SYN1deg' + '\n' +
                         r'Outgoing LW Radiation [W m$^{-2}$], 1-1-2019:' + str(j) + 'h')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
        elif i == 1:
            im = ax.imshow(aqua_diff[k, :, :], vmin=-30, vmax=30)
            ax.set_title(r'Aqua Only CRS1deg$_{\beta}$ minus SYN1deg' + '\n' +
                         r'Outgoing LW Radiation [W m$^{-2}$], 1-1-2019:' + str(j) + 'h')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])

    fig.subplots_adjust(right=0.8)
    cb_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cb_ax)
    plt.show()


# take the mean difference at each grid box
terra_mean_diff = np.nanmean(terra_diff, axis=0)
aqua_mean_diff = np.nanmean(aqua_diff, axis=0)


# combine fields into a single matrix to loop over figure axes
field = np.stack((terra_mean_diff, aqua_mean_diff), axis=2)
title_str = [r'Mean Daytime $Terra$ $Only$ CRS1deg$_{\beta}$ minus SYN1deg' + '\n' +
             r'Clear-Sky TOA Outgoing LW Radiation [W m$^{-2}$], 1-1-2019:00-23h',
             r'Mean Daytime $Aqua$ $Only$ CRS1deg$_{\beta}$ minus SYN1deg' + '\n' +
             r'Clear-Sky TOA Outgoing LW Radiation [W m$^{-2}$], 1-1-2019:00-23h'
             ]


# select colormap
cmap = ceres.set_colormap(Balance_20, 0)

ceres.plot_gridded_fields(nrows=2, ncols=1, cen_lon=0,
                          date_str='', title_str=title_str,
                          cmap=cmap, cmap_lims=((-30, 30), (-30, 30)),
                          varname='', levname='', varunits='',
                          lon=lon_syn1deg, lat=lat_syn1deg,
                          field=field)


# compute regional averages
weights = ceres.cos_lat_weight(np.flipud(lat_syn1deg))
global_avg, sh0to30_avg, sh30to60_avg, sh60to90_avg, nh0to30_avg, nh30to60_avg, nh60to90_avg = \
ceres.compute_regional_averages(field=terra_mean_diff,
                                latitude=np.flipud(lat_syn1deg),
                                weights=weights)


# -----------------
#
# field2 = np.stack((terra_mask[0, :, :], terra_only_mask[0, :, :], aqua_mask[0, :, :],
#                    num_sw_obs[0, :, :], both_mask[0, :, :], num_lw_obs[0, :, :]), axis=2)
# title_str2 = [r'CRS1deg$_{\beta}$ Terra FM1, 1-1-2019:00h',
#               r'CRS1deg$_{\beta}$ Terra - (Terra$\bigcap$Aqua) 1-1-2019:00h',
#               r'CRS1deg$_{\beta}$ Aqua FM1, 1-1-2019:00h',
#               r'SYN1deg Terra & Aqua SW mask, 1-1-2019:00h',
#               r'CRS1deg$_{\beta}$ Terra & Aqua, 1-1-2019:00h',
#               r'SYN1deg Terra & Aqua LW mask, 1-1-2019:00h']
#
# print(field2.shape)
#
# cmap = ceres.set_colormap(LaJolla_20_r, 0)
#
# ceres.plot_gridded_fields(nrows=3, ncols=2, cen_lon=0,
#                           date_str='', title_str=title_str2,
#                           cmap=cmap,
#                           cmap_lims=((0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1)),
#                           varname='', levname='', varunits='',
#                           lon=lon_syn1deg, lat=lat_syn1deg,
#                           field=field2)
#

# lw_syn1deg_obs, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
#                                               var_name='obs_clr_toa_lw',
#                                               fill=True)
#
# lw_syn1deg_comp, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
#                                                var_name='init_clr_toa_lw_up',
#                                                fill=True)
#
# field2 = np.stack((lw_syn1deg_obs[0, :, :], lw_syn1deg_comp[0, :, :]), axis=2)
# title_str2 = [r'Observed OLR SYN1deg Terra FM1, 1-1-2019:00h',
#               r'Computed OLR SYN1deg Terra FM1, 1-1-2019:00h']
#
# print(field2.shape)
#
# cmap = ceres.set_colormap(Thermal_20, 0)
#
# ceres.plot_gridded_fields(nrows=2, ncols=1, cen_lon=0,
#                           date_str='', title_str=title_str2,
#                           cmap=cmap,
#                           cmap_lims=((150, 350), (-30, 30)),
#                           varname='', levname='', varunits='',
#                           lon=lon_syn1deg, lat=lat_syn1deg,
#                           field=field2)

