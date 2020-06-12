
import numpy as np
import cerestools as ceres
import matplotlib.pyplot as plt
from palettable.cmocean.sequential import Thermal_20


path2 = '/Users/rcscott2/Desktop/CERES/SYN1deg/'
file2 = 'CER_SYN1deg-1Hour_Terra-Aqua-MODIS_Edition4A_407406.20190101'
file_path2 = path2 + file2


var_syn1deg, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
                                           var_name='init_all_toa_sw_up',
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

# terra_mask[:] = np.nan
# aqua_mask[:] = np.nan
# both_mask[:] = np.nan

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

    terra_lat_all, terra_lon_all, _, _, _, terra_sza_all = ceres.read_crs_geolocation_dev(file_path=terra_crs_file)
    terra_var_all, _, _, _ = ceres.read_crs_var_dev(file_path=terra_crs_file,
                                                    var_name='Shortwave flux - upward - total sky',
                                                    lev_arg=0,
                                                    fill=True)

    aqua_lat_all, aqua_lon_all, _, _, _, aqua_sza_all = ceres.read_crs_geolocation_dev(file_path=aqua_crs_file)
    aqua_var_all, _, _, _ = ceres.read_crs_var_dev(file_path=aqua_crs_file,
                                                   var_name='Shortwave flux - upward - total sky',
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

    terra_var_gridded[k, :, :] = ceres.grid_to_1x1_deg_equal_angle(terra_lat_all, terra_lon_all, terra_var_all, False)
    aqua_var_gridded[k, :, :] = ceres.grid_to_1x1_deg_equal_angle(aqua_lat_all, aqua_lon_all, aqua_var_all, False)

    print(k)

    for i in range(180):
        for j in range(360):
            if terra_var_gridded[k, i, j] > 0:
                terra_mask[k, i, j] = 1
            if aqua_var_gridded[k, i, j] > 0:
                aqua_mask[k, i, j] = 1
            if terra_var_gridded[k, i, j] > 0 and aqua_var_gridded[k, i, j] > 0:
                both_mask[k, i, j] = 1


terra_only_mask = terra_mask - both_mask
aqua_only_mask = aqua_mask - both_mask
terra_aqua_mask = terra_mask + aqua_mask


# for k in range(24):
#
#     if k < 10:
#         j = '0' + str(k)
#     elif k >= 10:
#         j = str(k)
#
#     fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(11, 7))
#     for i, ax in enumerate(axes.flat):
#         if i == 0:
#             im = ax.imshow(terra_mask[k, :, :], vmin=0, vmax=1)
#             ax.set_title(r'CRS1deg$_{\beta}$ Terra FM1, 1-1-2019:' + str(j) + 'h')
#             ax.set_xticklabels([])
#             ax.set_yticklabels([])
#             ax.set_xticks([])
#             ax.set_yticks([])
#         elif i == 3:
#             im = ax.imshow(num_sw_obs[k, :, :], vmin=0, vmax=1)
#             ax.set_title('SYN1deg Terra & Aqua SW mask, 1-1-2019:' + str(j) + 'h')
#             ax.set_xticklabels([])
#             ax.set_yticklabels([])
#             ax.set_xticks([])
#             ax.set_yticks([])
#         elif i == 2:
#             im = ax.imshow(aqua_mask[k, :, :], vmin=0, vmax=1)
#             ax.set_title(r'CRS1deg$_{\beta}$ Aqua FM3, 1-1-2019:' + str(j) + 'h')
#             ax.set_xticklabels([])
#             ax.set_yticklabels([])
#             ax.set_xticks([])
#             ax.set_yticks([])
#         elif i == 5:
#             im = ax.imshow(num_lw_obs[k, :, :], vmin=0, vmax=1)
#             ax.set_title('SYN1deg Terra & Aqua LW mask, 1-1-2019:' + str(j) + 'h')
#             ax.set_xticklabels([])
#             ax.set_yticklabels([])
#             ax.set_xticks([])
#             ax.set_yticks([])
#         elif i == 4:
#             im = ax.imshow(terra_aqua_mask[k, :, :], vmin=0, vmax=2)
#             ax.set_title(r'CRS1deg$_{\beta}$ Terra & Aqua, 1-1-2019:' + str(j) + 'h')
#             ax.set_xticklabels([])
#             ax.set_yticklabels([])
#             ax.set_xticks([])
#             ax.set_yticks([])
#         elif i == 1:
#             im = ax.imshow(aqua_only_mask[k, :, :], vmin=0, vmax=1)
#             ax.set_title(r'CRS1deg$_{\beta}$ Aqua - (Terra$\bigcap$Aqua), 1-1-2019:' + str(j) + 'h')
#             ax.set_xticklabels([])
#             ax.set_yticklabels([])
#             ax.set_xticks([])
#             ax.set_yticks([])
#
#     fig.subplots_adjust(right=0.8)
#     cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
#     fig.colorbar(im, cax=cbar_ax)
#
#     plt.show()


# take the difference CRS1deg_beta minus SYN1deg hour-by-hour
aqua_diff = aqua_only_mask*aqua_var_gridded - aqua_only_mask*var_syn1deg
terra_diff = terra_only_mask*terra_var_gridded - terra_only_mask*var_syn1deg
aqua_diff[terra_aqua_mask == 2] = np.nan
terra_diff[terra_aqua_mask == 2] = np.nan

# for k in range(24):
#
#     if k < 10:
#         j = '0' + str(k)
#     elif k >= 10:
#         j = str(k)
#
#     fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(6, 7))
#     for i, ax in enumerate(axes.flat):
#         if i == 0:
#             im = ax.imshow(terra_diff[k, :, :], vmin=-30, vmax=30)
#             ax.set_title(r'Terra Only CRS1deg$_{\beta}$ minus SYN1deg' + '\n' +
#                          r'Outgoing LW Radiation [W m$^{-2}$], 1-1-2019:' + str(j) + 'h')
#             ax.set_xticklabels([])
#             ax.set_yticklabels([])
#             ax.set_xticks([])
#             ax.set_yticks([])
#         elif i == 1:
#             im = ax.imshow(aqua_diff[k, :, :], vmin=-30, vmax=30)
#             ax.set_title(r'Aqua Only CRS1deg$_{\beta}$ minus SYN1deg' + '\n' +
#                          r'Outgoing LW Radiation [W m$^{-2}$], 1-1-2019:' + str(j) + 'h')
#             ax.set_xticklabels([])
#             ax.set_yticklabels([])
#             ax.set_xticks([])
#             ax.set_yticks([])
#
#     fig.subplots_adjust(right=0.8)
#     cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
#     fig.colorbar(im, cax=cbar_ax)
#
#     plt.show()

terra_mean_diff = np.nanmean(terra_diff, axis=0)
aqua_mean_diff = np.nanmean(aqua_diff, axis=0)


print(terra_diff[0, 0:100, 0:100])


fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(6, 7))
for i, ax in enumerate(axes.flat):
    if i == 0:
        im = ax.imshow(terra_mean_diff, vmin=-100, vmax=100)
        ax.set_title(r'Mean Daytime $Terra$ $Only$ CRS1deg$_{\beta}$ minus SYN1deg' + '\n' +
                     r'Reflected SW Radiation [W m$^{-2}$], 1-1-2019:00-23h')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
    elif i == 1:
        im = ax.imshow(aqua_mean_diff, vmin=-100, vmax=100)
        ax.set_title(r'Mean Daytime $Aqua$ $Only$ CRS1deg$_{\beta}$ minus SYN1deg' + '\n' +
                     r'Reflected SW Radiation [W m$^{-2}$], 1-1-2019:00-23h')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)

plt.show()


field = np.stack((terra_mean_diff, aqua_mean_diff), axis=2)
title_str = [r'Mean Daytime $Terra$ $Only$ CRS1deg$_{\beta}$ minus SYN1deg' + '\n' +
             r'Reflected SW Radiation [W m$^{-2}$], 1-1-2019:00-23h',
             r'Mean Daytime $Aqua$ $Only$ CRS1deg$_{\beta}$ minus SYN1deg' + '\n' +
             r'Reflected SW Radiation [W m$^{-2}$], 1-1-2019:00-23h'
             ]


def plot_gridded_fields(nrows, ncols, cen_lon,
                        date_str, title_str,
                        varname, levname, varunits,
                        lon, lat, input):
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
    :param nightshade: whether to use nighshade feature    [boolean]
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
    fig = plt.figure(figsize=(10, 8))
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
        ax.set_extent([-180, 180, -90, 90], projection)
        ax.text(0.5, -0.1, varname + ' - ' + levname + '\n' + varunits,
                va='bottom', ha='center',
                rotation='horizontal', rotation_mode='anchor',
                transform=ax.transAxes, fontsize=10)

    # To use a different color bar range each time use a tuple of tuples
    #        ((0, 120), (0, 120), (0, 120), (0, 120), (0, 120), (0, 120), (0, 120))
    limits = (-100, 100)
    for i in range((nrows * ncols)):
        im = axgr[i].pcolor(lon, lat, field[:, :, i],
                            vmin=limits[0], vmax=limits[1])
        axgr.cbar_axes[i].colorbar(im)
        axgr[i].set_title(title_str[i])

    #plt.tight_layout(pad=3.0)
    plt.show()
    return


plot_gridded_fields(nrows=2, ncols=1, cen_lon=0,
                    date_str='', title_str=title_str,
                    varname='', levname='', varunits='',
                    lon=lon_syn1deg, lat=lat_syn1deg,
                    input=field)





# ==============================================================================
#
# # # CERES Ed4 CRS data
# # var_all, lon_all, lat_all, sza_all = ceres.read_day_of_crs_files(
# #                                path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
# #                                file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20190101',
# #                                variable='Longwave flux - upward - total sky',
# #                                lev_arg=0,
# #                                fill=True)
#
#
# terra_crs_file = '/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.2019010101'
# aqua_crs_file = '/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/CER_CRS4_Aqua-FM3-MODIS_GH4_1111TH.2019010101'
#
#
# date, date_str = ceres.get_date(terra_crs_file)
#
#
# terra_lat_all, terra_lon_all, _, _, _, terra_sza_all = ceres.read_crs_geolocation_dev(file_path=terra_crs_file)
# terra_var_all, _, _, _ = ceres.read_crs_var_dev(file_path=terra_crs_file,
#                                                 var_name='Longwave flux - upward - total sky',
#                                                 lev_arg=0,
#                                                 fill=True)
#
# aqua_lat_all, aqua_lon_all, _, _, _, aqua_sza_all = ceres.read_crs_geolocation_dev(file_path=aqua_crs_file)
# aqua_var_all, _, _, _ = ceres.read_crs_var_dev(file_path=aqua_crs_file,
#                                                var_name='Longwave flux - upward - total sky',
#                                                lev_arg=0,
#                                                fill=True)
#
#
# # convert footprint longitude range from 0 to 360 deg to -180 to 180 deg
# terra_lon_all = ceres.swath_lon_360_to_180(terra_lon_all)
# aqua_lon_all = ceres.swath_lon_360_to_180(aqua_lon_all)
#
# # terra_lat_all, terra_lon_all, terra_var_all, terra_sza_all \
# #     = ceres.swath_daytime_only(lat=terra_lat_all,
# #                                lon=terra_lon_all,
# #                                var=terra_var_all,
# #                                sza=terra_sza_all,
# #                                sza_cutoff=90)
# #
# # aqua_lat_all, aqua_lon_all, aqua_var_all, aqua_sza_all \
# #     = ceres.swath_daytime_only(lat=aqua_lat_all,
# #                                lon=aqua_lon_all,
# #                                var=aqua_var_all,
# #                                sza=aqua_sza_all,
# #                                sza_cutoff=90)
#
#
# # grid and average footprints in 1 deg x 1 deg grid boxes
# terra_var_gridded = ceres.grid_to_1x1_deg_equal_angle(terra_lat_all, terra_lon_all, terra_var_all, False)
# aqua_var_gridded = ceres.grid_to_1x1_deg_equal_angle(aqua_lat_all, aqua_lon_all, aqua_var_all, False)
#
#
# terra_mask = np.zeros([180, 360])
# aqua_mask = np.zeros([180, 360])
# both_mask = np.zeros([180, 360])
#
#
# for i in range(180):
#     for j in range(360):
#         if terra_var_gridded[i, j] > 0:
#             terra_mask[i, j] = 1
#         if aqua_var_gridded[i, j] > 0:
#             aqua_mask[i, j] = 1
#         if terra_var_gridded[i, j] > 0 and aqua_var_gridded[i, j] > 0:
#             both_mask[i, j] = 1
#
#
# plt.imshow(terra_mask)
# plt.title("Terra CERES FM1 mask")
# plt.colorbar()
# plt.show()
#
# plt.imshow(aqua_mask)
# plt.title("Aqua CERES FM3 mask")
# plt.colorbar()
# plt.show()
#
# plt.imshow(terra_mask + aqua_mask)
# plt.title("Terra CERES FM1 + Aqua CERES FM3 mask")
# plt.colorbar()
# plt.show()
#
# plt.imshow(both_mask)
# plt.title(r"Terra $\bigcap$ Aqua")
# plt.colorbar()
# plt.show()
#
# terra_only_mask = terra_mask-both_mask
# plt.imshow(terra_only_mask)
# plt.title(r"Terra minus (Terra $\bigcap$ Aqua)")
# plt.colorbar()
# plt.show()
#
#
# # ==============================================================================
#
#
# path2 = '/Users/rcscott2/Desktop/CERES/SYN1deg/'
# file2 = 'CER_SYN1deg-1Hour_Terra-Aqua-MODIS_Edition4A_407406.20190101'
# file_path2 = path2 + file2
#
#
# toa_lw_up_all, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
#                                              var_name='init_all_toa_lw_up',
#                                              fill=False)
#
# lat_syn1deg, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
#                                            var_name='latitude',
#                                            fill=False)
#
# lon_syn1deg, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
#                                            var_name='longitude',
#                                            fill=False)
#
# num_sw_obs, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
#                                           var_name='num_sw_obs',
#                                           fill=False)
#
# num_lw_obs, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
#                                           var_name='num_lw_obs',
#                                           fill=False)
#
# num_geo_sw_obs, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
#                                               var_name='num_geo_sw_obs',
#                                               fill=False)
#
# num_geo_lw_obs, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
#                                               var_name='num_geo_lw_obs',
#                                               fill=False)
#
# print("SYN1deg Terra-Aqua SW obs shape: ", num_sw_obs.shape)
# print("SYN1deg GEO SW obs shape: ", num_geo_sw_obs.shape)
# print("SYN1deg Terra-AquaLW obs shape: ", num_lw_obs.shape)
# print("SYN1deg GEO LW obs shape: ", num_geo_lw_obs.shape)
#
#
# for k in range(7):
#     plt.imshow(num_sw_obs[k, :, :])
#     plt.title('SYN1deg Terra-Aqua SW obs ' + date_str)
#     plt.show()
#
# for k in range(7):
#     plt.imshow(num_lw_obs[k, :, :])
#     plt.title('SYN1deg Terra-Aqua LW obs ' + date_str)
#     plt.show()
#
# # for k in range(24):
# #     plt.imshow(num_geo_sw_obs[k, :, :]+num_sw_obs[k, :, :])
# #     plt.title('SYN1deg Terra-Aqua + GEO SW obs ' + str(k) + 'hr')
# #     plt.show()
# #
# # for k in range(24):
# #     plt.imshow(num_geo_lw_obs[k, :, :])
# #     plt.title('SYN1deg Terra-Aqua + GEO LW obs ' + str(k) + 'hr')
# #     plt.show()
#
#
# # select colormap
# cmap = ceres.set_colormap(Thermal_20, 0)
#
# # plot the newly gridded field
# ceres.global_map(lon=lon_syn1deg, lat=lat_syn1deg, field=terra_var_gridded,
#                  varname='CERES LW TOA flux - upwards', varunits=r'W m$^{-2}$',
#                  nrows=1, ncols=1, cen_lon=0,
#                  cmap=cmap, cmap_lims=(150, 350), ti_str=r'Terra CERES FM1 CRS1deg$_{\beta}$ ' + date_str)
#
# # plot the newly gridded field
# ceres.global_map(lon=lon_syn1deg, lat=lat_syn1deg, field=aqua_var_gridded,
#                  varname='CERES LW TOA flux - upwards', varunits=r'W m$^{-2}$',
#                  nrows=1, ncols=1, cen_lon=0,
#                  cmap=cmap, cmap_lims=(150, 350), ti_str=r'Aqua CERES FM3 CRS1deg$_{\beta}$ ' + date_str)
#
#
# # ========================================================
# # # CRS1deg-beta
# # plt.imshow(terra_var_gridded)
# # plt.colorbar()
# # plt.title(r'Terra CRS1deg$_{\beta}$ OLR [W m$^{-2}$] JAN-1-2019')
# # plt.clim(150, 350)
# # plt.show()
# #
# #
# # # SYN1deg
# # plt.imshow(np.nanmean(toa_lw_up_all, axis=0))
# # plt.colorbar()
# # plt.title('SYN1deg OLR [W m$^{-2}$] JAN-1-2019')
# # plt.clim(150, 350)
# # plt.show()
# #
# #
# # # plot of the difference between daily means
# # plt.imshow(gridded_var_all-np.nanmean(toa_lw_up_all, axis=0))
# # plt.colorbar()
# # plt.title(r'Terra CRS1deg$_{\beta}$ - SYN1deg OLR [W m$^{-2}$] JAN-1-2019')
# # plt.clim(-50, 50)
# # plt.show()
# #
# # for k in range(1):
# #     plt.imshow(num_lw_obs[k, :, :]*toa_lw_up_all[k, :, :])
# #     plt.title('TOA LW ' + str(k) + 'hr')
# #     plt.colorbar()
# #     plt.clim(150, 350)
# #     plt.show()
# #
#
#
#
#
