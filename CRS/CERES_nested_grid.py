# ==========================================================================
# Author: Ryan C. Scott
#         ryan.c.scott@nasa.gov
#
# This script grids CERES FOVs to the CERES 1 deg x 1 deg nested grid.
# ==========================================================================

import cerestools as ceres
import matplotlib.pyplot as plt

terra_var_all, terra_lon_all, terra_lat_all, terra_sza_all = ceres.read_day_of_crs_files(
                               path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
                               file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20190101',
                               variable='Shortwave flux - upward - total sky',
                               lev_arg=0,
                               fill=True)

terra_lon_all = ceres.swath_lon_360_to_180(terra_lon_all)

# terra_lat_all, terra_lon_all, terra_var_all, terra_sza_all = \
#     ceres.swath_daytime_only(lat=terra_lat_all,
#                              lon=terra_lon_all,
#                              var=terra_var_all,
#                              sza=terra_sza_all,
#                              sza_cutoff=90)


def grid_to_1x1_deg_ceres_nested(lat_data, lon_data, variable, lon_360=True):
    """
    ----------------------------------------------------------------------------
    This functions bins & grids CERES footprints to a 1 deg x 1 deg equal angle
    latitude-longitude grid using the SciPy stats routine binned_statistic_2d.
    After FOVs are aggregated into 1 deg x 1 deg regions it computes the mean
    of the input "variable" - alternatively, it can compute the # of footprints,
    the median, or other statistics (including user-defined functions).
    ----------------------------------------------------------------------------
    :param lat_data: FOV latitude array                            [float]
    :param lon_data: FOV longitude array                           [float]
    :param variable: for which gridded statistic will be computed  [float]
    :param lon_360:  use 0 to 360 or -180 to 180 longitude bins    [boolean]
    :return: the field of gridded FOVs                             [float]
    ----------------------------------------------------------------------------
    """
    import time
    import copy
    import numpy as np
    from scipy import stats

    # 1 degree latitude bins, -90 to 90 deg
    lat_bins = np.arange(-90, 91)

    # different lon bin extents for different lat zones
    lon_ext = [1, 2, 4, 8]

    # loop over nested grid zones
    for i, el in enumerate(lon_ext):

        # set appropriate longitude extent
        if lon_360 is True:
            lon_bins = np.arange(0, 361, el)     # lon: 0 to 360
        elif lon_360 is False:
            lon_bins = np.arange(-180, 181, el)

        print(i, el)
        print('Lon bins:\n', lon_bins)

        # time the griding procedure
        tic = time.time()

        # Zone 1
        if i == 0:
            var1 = copy.copy(variable)
            var1[abs(lat_data) > 45] = np.nan
            zone_1 = stats.binned_statistic_2d(lon_data, lat_data, var1,
                                               statistic=np.nanmean,
                                               bins=[lon_bins, lat_bins])
            zone1 = np.rot90(zone_1.statistic)

        # Zone 2
        if i == 1:
            var2 = copy.copy(variable)
            var2[abs(lat_data) < 45] = np.nan
            var2[abs(lat_data) > 70] = np.nan
            zone_2 = stats.binned_statistic_2d(lon_data, lat_data, var2,
                                               statistic=np.nanmean,
                                               bins=[lon_bins, lat_bins])
            zone2 = np.rot90(zone_2.statistic)

        # Zone 3
        if i == 2:
            var3 = copy.copy(variable)
            var3[abs(lat_data) < 70] = np.nan
            var3[abs(lat_data) > 80] = np.nan
            zone_3 = stats.binned_statistic_2d(lon_data, lat_data, var3,
                                               statistic=np.nanmean,
                                               bins=[lon_bins, lat_bins])
            zone3 = np.rot90(zone_3.statistic)

        # Zone 4
        if i == 3:
            var4 = copy.copy(variable)
            var4[abs(lat_data) < 80] = np.nan
            var4[abs(lat_data) > 89] = np.nan
            zone_4 = stats.binned_statistic_2d(lon_data, lat_data, var4,
                                               statistic=np.nanmean,
                                               bins=[lon_bins, lat_bins])
            zone4 = np.rot90(zone_4.statistic)

        # finish timing it, print relevant info
        toc = time.time()
        print(toc - tic, 'seconds elapsed during grid_to_1x1_ceres_nested\n')

    print('==============================================================\n')
    print(zone1.shape)
    print(zone2.shape)
    print(zone3.shape)
    print(zone4.shape)

    # for i in range(180):
    #     for j in range(360):
    #         if i




    return zone1, zone2, zone3, zone4


z1, z2, z3, z4 = grid_to_1x1_deg_ceres_nested(lat_data=terra_lat_all,
                                              lon_data=terra_lon_all,
                                              variable=terra_var_all,
                                              lon_360=False)

plt.imshow(z1)
plt.show()

plt.imshow(z2)
plt.show()

plt.imshow(z3)
plt.show()

plt.imshow(z4)
plt.show()


# variable = terra_var_all
# lat_data = terra_lat_all
# lon_data = terra_lon_all
# lon_360 = False
#
# lat_bins = np.arange(-90, 91)  # lat: -90 to 90 deg
#

# # each FOV has lat and lon indices that map to each grid box
# lat_ind = np.digitize(lat_data, lat_bins)
# lon_ind = np.digitize(lon_data, lon_bins)
#
# # loop over FOVs and show their index
# for n in range(lat_data.size):
#     print(lat_bins[lat_ind[n]-1], "<=", lat_data[n], "<", lat_bins[lat_ind[n]], '...', lat_ind[n])
#     print(lon_bins[lon_ind[n]-1], "<=", lon_data[n], "<", lon_bins[lon_ind[n]], '...', lon_ind[n])

# start timing it
# tic = time.time()
#
# # compute mean in each grid box - statistics: 'count', 'mean', 'median'
# gridded = stats.binned_statistic_2d(lon_data, lat_data, variable,
#                                     statistic=np.nanmean,
#                                     bins=[lon_bins, lat_bins])
#
# gridded_stat = np.rot90(gridded.statistic)

# finish timing it, print relevant info
# toc = time.time()
# print(toc - tic, 'seconds elapsed during grid_to_1x1_deg_equal_angle\n')
#
# print("Shape of 1 x 1 gridded field:")
# print(gridded_stat.shape)
#
# # quick & dirty plot of the result
# plt.pcolor(gridded_stat)
# plt.colorbar()
# plt.show()

# # grid each zone
# for i, el in enumerate(lon_ext):
#
#     if lon_360 is True:
#         lon_bins = np.arange(0, 361, el)     # lon: 0 to 360
#     elif lon_360 is False:
#         lon_bins = np.arange(-180, 181, el)
#
#     if el == 1:
#         print("Griding ZONE 1: abs(lat) < 45")
#         var1 = variable
#         var1[abs(lat_data) > 45] = np.nan
#
#             gridded_stat1 = stats.binned_statistic_2d(lon_data, lat_data, var1,
#                                                   statistic=np.nanmean,
#                                                   bins=[lon_bins, lat_bins])
#
#             gridded1[:, :] = np.rot90(gridded_stat1.statistic)
#
#         elif el == 2:
#             print("Griding ZONE 2: 45 < abs(lat) < 70")
#             var2 = variable
#             var2[abs(lat_data) < 45] = np.nan
#             var2[abs(lat_data) > 70] = np.nan
#
#             gridded_stat2 = stats.binned_statistic_2d(lon_data, lat_data, var2,
#                                                   statistic=np.nanmean,
#                                                   bins=[lon_bins, lat_bins])
#
#             gridded2[:, :] = np.rot90(gridded_stat2.statistic)
#
#         elif el == 4:
#             print("Griding ZONE 3: 70 < abs(lat) < 80")
#             var3 = variable
#             var3[abs(lat_data) < 70] = np.nan
#             var3[abs(lat_data) > 80] = np.nan
#
#             gridded_stat3 = stats.binned_statistic_2d(lon_data, lat_data, var3,
#                                                   statistic=np.nanmean,
#                                                   bins=[lon_bins, lat_bins])
#
#             gridded3[:, :] = np.rot90(gridded_stat3.statistic)
#
#         elif el == 8:
#             print("Griding ZONE 4: 80 < abs(lat) < 89")
#             var4 = variable
#             var4[abs(lat_data) < 80] = np.nan
#             var4[abs(lat_data) > 89] = np.nan
#
#             gridded_stat4 = stats.binned_statistic_2d(lon_data, lat_data, var4,
#                                                   statistic=np.nanmean,
#                                                   bins=[lon_bins, lat_bins])
#
#             gridded4[:, :] = np.rot90(gridded_stat4.statistic)

