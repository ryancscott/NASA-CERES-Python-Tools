
import cerestools as ceres
from palettable.cmocean.sequential import Thermal_20
from palettable.colorbrewer.sequential import BuPu_9_r

path = '/Users/rcscott2/Desktop/CERES/EBAF/'
file = 'CERES_EBAF-TOA_Ed4.1_Subset_200003-201910.nc'
file_path = path + file

# ceres.print_nc_file_info(file_path)

lat, lat_name, lat_units = ceres.read_ebaf_var(file_path=file_path, variable='lat')
lon, lon_name, lon_units = ceres.read_ebaf_var(file_path=file_path, variable='lon')


# CERES Ed4A SSF data
var_all4, lon_all4, lat_all4, sza_all4 = \
     ceres.read_day_of_ssf_files(path='/Users/rcscott2/Desktop/CRS/ASDC_archive/SSF_Ed4A/',
                                 file_struc='CER_SSF_Terra-FM1-MODIS_Edition4A_404405.20190101',
                                 variable='CERES LW TOA flux - upwards',
                                 fill=True)

# CERES Ed4 CRS data
# var_all4, lon_all4, lat_all4, sza_all4 = ceres.read_day_of_crs_files(
#                                path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
#                                file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20190101',
#                                variable='Shortwave flux - upward - total sky',
#                                lev_arg=0,
#                                fill_nan=True)


gridded_var_all4 = ceres.grid_to_1x1_degree(lat_all4, lon_all4, var_all4)

print("Shape of 1 deg x 1 deg gridded field:\n")
print(gridded_var_all4.shape)


cmap = ceres.set_colormap(Thermal_20, 0)

ceres.global_map(lon=lon, lat=lat, field=gridded_var_all4,
                 varname='CERES LW TOA flux - upwards', varunits=r'W m$^{-2}$',
                 nrows=1, ncols=1, cen_lon=0,
                 cmap=cmap, cmap_lims=(0, 350), ti_str='CERES Terra FM1 "CRS1deg" JAN-01-2019')


# import numpy as np
# import cerestools as ceres
# from scipy import stats
# import matplotlib.pyplot as plt
#
# path = '/Users/rcscott2/Desktop/CERES/EBAF/'
# file = 'CERES_EBAF-TOA_Ed4.1_Subset_200003-201910.nc'
# file_path = path + file
#
# ceres.print_nc_file_info(file_path)
#
# #lat, lat_name, lat_units = ceres.read_ebaf_var(file_path=file_path, variable='lat')
# #lon, lon_name, lon_units = ceres.read_ebaf_var(file_path=file_path, variable='lon')
#
# # lat = np.array(lat)
# # lat_bins = lat - 0.5
# # print(lat_bins)
# lat_bins = np.arange(-90, 91)
# lon_bins = np.arange(0, 361)
#
# print(lat_bins)
# print(lon_bins)
#
# # fake FOV lat data
# lat_data = np.random.random_sample(2000000)
# lat_data = 180*lat_data-90
# print(lat_data)
#
# # fake FOV lon data
# lon_data = np.random.random_sample(2000000)
# lon_data = 360*lon_data
# print(lon_data)
#
# lat_ind = np.digitize(lat_data, lat_bins)
# lon_ind = np.digitize(lon_data, lon_bins)
#
# # print(lat_ind)
# # print(lon_ind)
#
# for n in range(lat_data.size):
#     print(lat_bins[lat_ind[n]-1], "<=", lat_data[n], "<", lat_bins[lat_ind[n]], '...', lat_ind[n])
#     print(lon_bins[lon_ind[n]-1], "<=", lon_data[n], "<", lon_bins[lon_ind[n]], '...', lon_ind[n])
#
#
# print(len(lat_bins))
# print(len(lon_bins))
#
# ret = stats.binned_statistic_2d(lon_data, lat_data, lon_data, 'mean', bins=[lon_bins, lat_bins])
# ret.statistic
#
# plt.imshow(ret.statistic)
# plt.show()
#
# gridded = np.empty([180, 360, 200])
# for n in range(200):
#     # i = lat_ind[n]-1
#     # j = lon_ind[n]-1
#     # print(i, j, lat_data[n])
#     gridded[lat_ind[n]-1, lon_ind[n]-1, n] = lat_data[n]
#
# print(gridded[9, 9, :].reshape((1, 200)))
#
#
#
#
# # a = []
# # b = [[a for n in range(360)] for m in range(180)]
# #
# # print(b)
# #
# # for n in range(2000):
# #     i = lat_ind[n]-1
# #     j = lon_ind[n]-1
# #     print(i, j, lat_data[n])
# #     b[i][j].append(lat_data[n])
# #
# # print(b[0][0])
#

# # print(np.sort(lat_ind))
# # print(np.sort(lat_data))
# #
# # print(np.sort(lon_ind))
# # print(np.sort(lon_data))
#
# a1 = [[[1, 2, 3], [], []],
#      [[], [], []]]
#
# a1[0][0].append(1)
#
# print(a1[0][0])
# #
#
# # x = np.array([0.2, 6.4, 3.0, 1.6])
# # bins = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])
# # inds = np.digitize(x, bins)
# #
# # for n in range(x.size):
# #     print(bins[inds[n]-1], "<=", x[n], "<", bins[inds[n]])
# #



