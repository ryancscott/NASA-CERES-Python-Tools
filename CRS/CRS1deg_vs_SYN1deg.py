
import numpy as np
import cerestools as ceres
import matplotlib.pyplot as plt
from palettable.cmocean.sequential import Thermal_20


# ==============================================================================

# # CERES Ed4 CRS data
# var_all, lon_all, lat_all, sza_all = ceres.read_day_of_crs_files(
#                                path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
#                                file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20190101',
#                                variable='Longwave flux - upward - total sky',
#                                lev_arg=0,
#                                fill=True)


terra_crs_file = '/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.2019010100'
aqua_crs_file = '/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/CER_CRS4_Aqua-FM3-MODIS_GH4_1111TH.2019010100'


date, date_str = ceres.get_date(terra_crs_file)


terra_lat_all, terra_lon_all, _, _, _, terra_sza_all = ceres.read_crs_geolocation_dev(file_path=terra_crs_file)
terra_var_all, _, _, _ = ceres.read_crs_var_dev(file_path=terra_crs_file,
                                                var_name='Longwave flux - upward - total sky',
                                                lev_arg=0,
                                                fill=True)

aqua_lat_all, aqua_lon_all, _, _, _, aqua_sza_all = ceres.read_crs_geolocation_dev(file_path=aqua_crs_file)
aqua_var_all, _, _, _ = ceres.read_crs_var_dev(file_path=aqua_crs_file,
                                               var_name='Longwave flux - upward - total sky',
                                               lev_arg=0,
                                               fill=True)


# convert footprint longitude range from 0 to 360 deg to -180 to 180 deg
terra_lon_all = ceres.swath_lon_360_to_180(terra_lon_all)
aqua_lon_all = ceres.swath_lon_360_to_180(aqua_lon_all)

# terra_lat_all, terra_lon_all, terra_var_all, terra_sza_all \
#     = ceres.swath_daytime_only(lat=terra_lat_all,
#                                lon=terra_lon_all,
#                                var=terra_var_all,
#                                sza=terra_sza_all,
#                                sza_cutoff=90)
#
# aqua_lat_all, aqua_lon_all, aqua_var_all, aqua_sza_all \
#     = ceres.swath_daytime_only(lat=aqua_lat_all,
#                                lon=aqua_lon_all,
#                                var=aqua_var_all,
#                                sza=aqua_sza_all,
#                                sza_cutoff=90)


# grid and average footprints in 1 deg x 1 deg grid boxes
terra_var_gridded = ceres.grid_to_1x1_deg_equal_angle(terra_lat_all, terra_lon_all, terra_var_all, False)
aqua_var_gridded = ceres.grid_to_1x1_deg_equal_angle(aqua_lat_all, aqua_lon_all, aqua_var_all, False)


terra_mask = np.zeros([180, 360])
aqua_mask = np.zeros([180, 360])
both_mask = np.zeros([180, 360])

for i in range(180):
    for j in range(360):
        if terra_var_gridded[i, j] > 0:
            terra_mask[i, j] = 1
        if aqua_var_gridded[i, j] > 0:
            aqua_mask[i, j] = 1
        if terra_var_gridded[i, j] > 0 and aqua_var_gridded[i, j] > 0:
            both_mask[i, j] = 1


plt.imshow(terra_mask)
plt.title("Terra CERES FM1 mask")
plt.colorbar()
plt.show()

plt.imshow(aqua_mask)
plt.title("Aqua CERES FM3 mask")
plt.colorbar()
plt.show()

plt.imshow(aqua_mask)
plt.title("Terra CERES FM1 + Aqua CERES FM3 mask")
plt.colorbar()
plt.show()

plt.imshow(both_mask)
plt.title(r"Terra $\bigcap$ Aqua")
plt.colorbar()
plt.show()

plt.imshow(terra_mask-both_mask)
plt.title(r"Terra minus (Terra $\bigcap$ Aqua)")
plt.colorbar()
plt.show()


# ==============================================================================


path2 = '/Users/rcscott2/Desktop/CERES/SYN1deg/'
file2 = 'CER_SYN1deg-1Hour_Terra-Aqua-MODIS_Edition4A_407406.20190101'
file_path2 = path2 + file2


toa_lw_up_all, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
                                             var_name='init_all_toa_lw_up',
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

print("Num SW obs shape: ", num_sw_obs.shape)
print("Num GEO SW obs shape: ", num_geo_sw_obs.shape)
print("Num LW obs shape: ", num_lw_obs.shape)
print("Num GEO LW obs shape: ", num_geo_lw_obs.shape)

# for k in range(24):
#     plt.imshow(num_lw_obs[k, :, :]*toa_lw_up_all[k, :, :])
#     plt.title('TOA LW flux ' + str(k) + 'hr' + ' times Terra-Aqua mask')
#     plt.show()

for k in range(2):
    plt.imshow(num_sw_obs[k, :, :])
    plt.title('Number of SW obs ' + date_str)
    plt.show()

# for k in range(24):
#     plt.imshow(num_geo_sw_obs[k, :, :]+num_sw_obs[k, :, :])
#     plt.title('Number of SW + GEO SW obs ' + str(k) + 'hr')
#     plt.show()

for k in range(2):
    plt.imshow(num_lw_obs[k, :, :])
    plt.title('Number of LW obs ' + date_str)
    plt.show()
#
# for k in range(24):
#     plt.imshow(num_geo_lw_obs[k, :, :])
#     plt.title('Number of GEO LW obs ' + str(k) + 'hr')
#     plt.show()


# select colormap
cmap = ceres.set_colormap(Thermal_20, 0)
# plot the newly gridded field
ceres.global_map(lon=lon_syn1deg, lat=lat_syn1deg, field=terra_var_gridded,
                 varname='CERES LW TOA flux - upwards', varunits=r'W m$^{-2}$',
                 nrows=1, ncols=1, cen_lon=0,
                 cmap=cmap, cmap_lims=(150, 350), ti_str=r'Terra CERES FM1 CRS1deg$_{\beta}$ ' + date_str)

# select colormap
cmap = ceres.set_colormap(Thermal_20, 0)
# plot the newly gridded field
ceres.global_map(lon=lon_syn1deg, lat=lat_syn1deg, field=aqua_var_gridded,
                 varname='CERES LW TOA flux - upwards', varunits=r'W m$^{-2}$',
                 nrows=1, ncols=1, cen_lon=0,
                 cmap=cmap, cmap_lims=(150, 350), ti_str=r'Aqua CERES FM3 CRS1deg$_{\beta}$ ' + date_str)


# # CRS1deg-beta
# plt.imshow(terra_var_gridded)
# plt.colorbar()
# plt.title(r'Terra CRS1deg$_{\beta}$ OLR [W m$^{-2}$] JAN-1-2019')
# plt.clim(150, 350)
# plt.show()
#
#
# # SYN1deg
# plt.imshow(np.nanmean(toa_lw_up_all, axis=0))
# plt.colorbar()
# plt.title('SYN1deg OLR [W m$^{-2}$] JAN-1-2019')
# plt.clim(150, 350)
# plt.show()
#
#
# # plot of the difference between daily means
# plt.imshow(gridded_var_all-np.nanmean(toa_lw_up_all, axis=0))
# plt.colorbar()
# plt.title(r'Terra CRS1deg$_{\beta}$ - SYN1deg OLR [W m$^{-2}$] JAN-1-2019')
# plt.clim(-50, 50)
# plt.show()
#
# for k in range(1):
#     plt.imshow(num_lw_obs[k, :, :]*toa_lw_up_all[k, :, :])
#     plt.title('TOA LW ' + str(k) + 'hr')
#     plt.colorbar()
#     plt.clim(150, 350)
#     plt.show()
#




