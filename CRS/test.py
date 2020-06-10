# import numpy as np
#
# a = np.array([[0.5, 1, 6, 7],
#               [0.6, 4, 5, 8],
#               [1, 4, 6, 9],
#               [0.5, 0.3, 6, 7],
#               ])
#
# b = np.array([3., 2., 2., 1.])
# b = list(b)
# b = [int(el) for el in b]
#
# print(b)
#
# c = np.empty([4])
# for i in range(4):
#     print(a[i, :])
#     print(a[i, b[i]])
#     c[i] = a[i, b[i]]
#
# print(c)
#

# ==============================================================================

# import matplotlib.pyplot as plt
# import cerestools as ceres
#
# path = '/Users/rcscott2/Desktop/CERES/SYN1deg/'
# file = 'CER_SYN1deg-1Hour_Terra-Aqua-MODIS_Edition4A_407406.20190101'
#
# file_path = path + file
#
# def read_syn1deg_hdf(file_path, var_name, fill):
#     """
#     ----------------------------------------------------------------------------
#     This function reads data from CERES Level 3 Synoptic 1-degree (SYN1deg)
#     HDF data files.
#     ----------------------------------------------------------------------------
#     :param file_path:
#     :param var_name: variable name [string]
#     :return:
#     ----------------------------------------------------------------------------
#     """
#     import numpy as np
#     from pyhdf import SD
#     hdf = SD.SD(file_path)
#
#     # select and get the variable
#     data = hdf.select(var_name)
#     variable = data.get()
#     var_name_ = data.long_name
#     var_units = data.units
#
#     print('Reading... ' + var_name_)
#
#     # replace fill values with NaN
#     if fill is True:
#         print('Replacing fill values with NaN...')
#         variable[variable == data._FillValue] = np.nan
#
#     var_field = variable
#
#     return var_field, var_name_, var_units
#
#
# cld_lwp, _, _ = read_syn1deg_hdf(file_path=file_path, var_name='obs_cld_lwp', fill=False)
#
# print(cld_lwp.shape)
#
# toa_lw_up_all, _, _ = read_syn1deg_hdf(file_path=file_path, var_name='init_all_toa_lw_up', fill=False)
#
# print(toa_lw_up_all[0, :, :].shape)
#
# for k in range(24):
#     plt.imshow(toa_lw_up_all[k, :, :])
#     plt.show()
#
# print(toa_lw_up_all[0, 0:10, 0:10])

#
# path2 = '/Users/rcscott2/Desktop/CERES/SYN1deg/'
# file2 = 'CER_SYN1deg-1Hour_Terra-Aqua-MODIS_Edition4A_407406.20190101'
# file_path2 = path2 + file2
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
# print(num_sw_obs.shape)
#
# for k in range(24):
#     plt.imshow(num_lw_obs[k, :, :]*toa_lw_up_all[k, :, :])
#     plt.title('Number of SW obs ' + str(k) + 'hr')
#     plt.show()
#
# for k in range(24):
#     plt.imshow(num_geo_sw_obs[k, :, :]+num_sw_obs[k, :, :])
#     plt.title('Number of SW + GEO SW obs ' + str(k) + 'hr')
#     plt.show()
#
# for k in range(24):
#     plt.imshow(num_lw_obs[k, :, :])
#     plt.title('Number of LW obs ' + str(k) + 'hr')
#     plt.show()
#
# for k in range(24):
#     plt.imshow(num_geo_lw_obs[k, :, :])
#     plt.title('Number of GEO LW obs ' + str(k) + 'hr')
#     plt.show()

# ====================================================================================

import cerestools as ceres
import matplotlib.pyplot as plt
import numpy as np


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


terra_lon_all = ceres.swath_lon_360_to_180(terra_lon_all)
aqua_lon_all = ceres.swath_lon_360_to_180(aqua_lon_all)


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
plt.title("Terra mask")
plt.colorbar()
plt.show()

plt.imshow(aqua_mask)
plt.title("Aqua mask")
plt.colorbar()
plt.show()

plt.imshow(both_mask)
plt.title(r"Terra $\bigcap$ Aqua")
plt.colorbar()
plt.show()

# CERES NESTED GRID....
# import numpy as np
#
# # longitude bins
# # -45 to 45 deg latitude
# a = np.arange(0, 361)
# # 45 to 70 deg latitude
# b = np.arange(0, 361, 2)
# # 70 to 80 deg latitude
# c = np.arange(0, 361, 4)
# # 80 to 89 deg latitude
# d = np.arange(0, 361, 8)
#
# print(a)
# print(b)
# print(c)
# print(d)