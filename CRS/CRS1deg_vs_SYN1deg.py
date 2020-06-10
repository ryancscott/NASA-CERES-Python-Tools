
import time
import numpy as np
import cerestools as ceres
import matplotlib.pyplot as plt
from palettable.cmocean.sequential import Thermal_20
from palettable.colorbrewer.sequential import BuPu_9_r

# ==============================================================================

# path = '/Users/rcscott2/Desktop/CERES/EBAF/'
# file = 'CERES_EBAF-TOA_Ed4.1_Subset_200003-201910.nc'
# file_path = path + file
#
# # ceres.print_nc_file_info(file_path)
#
# lat_ebaf, lat_name, lat_units = ceres.read_ebaf_var(file_path=file_path, variable='lat')
# lon_ebaf, lon_name, lon_units = ceres.read_ebaf_var(file_path=file_path, variable='lon')

# --------------------------------------------------

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

# ==============================================================
# # CERES Ed4 CRS data
# var_all, lon_all, lat_all, sza_all = ceres.read_day_of_crs_files(
#                                path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
#                                file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20190101',
#                                variable='Longwave flux - upward - total sky',
#                                lev_arg=0,
#                                fill=True)


terra_crs_file = '/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.2019010100'
aqua_crs_file = '/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/CER_CRS4_Aqua-FM3-MODIS_GH4_1111TH.2019010100'


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


# convert FOV longitudinal range from 0 to 360 deg to -180 to 180 deg
terra_lon_all = ceres.swath_lon_360_to_180(terra_lon_all)
aqua_lon_all = ceres.swath_lon_360_to_180(aqua_lon_all)


# grid and average footprints in 1 deg x 1 deg grid boxes
terra_var_gridded = ceres.grid_to_1x1_deg_equal_angle(terra_lat_all, terra_lon_all, terra_var_all, False)
aqua_var_gridded = ceres.grid_to_1x1_deg_equal_angle(aqua_lat_all, aqua_lon_all, aqua_var_all, False)


# select colormap
cmap = ceres.set_colormap(Thermal_20, 0)

# plot the newly gridded field
ceres.global_map(lon=lon_syn1deg, lat=lat_syn1deg, field=gridded_var_all,
                 varname='CERES LW TOA flux - upwards', varunits=r'W m$^{-2}$',
                 nrows=1, ncols=1, cen_lon=0,
                 cmap=cmap, cmap_lims=(150, 350), ti_str='CERES Terra FM1 "CRS1deg" JAN-01-2019')


# CRS1deg-beta
plt.imshow(gridded_var_all)
plt.colorbar()
plt.title(r'Terra CRS1deg$_{\beta}$ OLR [W m$^{-2}$] JAN-1-2019')
plt.clim(150, 350)
plt.show()


# SYN1deg
plt.imshow(np.nanmean(toa_lw_up_all, axis=0))
plt.colorbar()
plt.title('SYN1deg OLR [W m$^{-2}$] JAN-1-2019')
plt.clim(150, 350)
plt.show()


# plot of the difference between daily means
plt.imshow(gridded_var_all-np.nanmean(toa_lw_up_all, axis=0))
plt.colorbar()
plt.title(r'Terra CRS1deg$_{\beta}$ - SYN1deg OLR [W m$^{-2}$] JAN-1-2019')
plt.clim(-50, 50)
plt.show()


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


for k in range(1):
    plt.imshow(num_lw_obs[k, :, :]*toa_lw_up_all[k, :, :])
    plt.title('TOA LW ' + str(k) + 'hr')
    plt.colorbar()
    plt.clim(150, 350)
    plt.show()





