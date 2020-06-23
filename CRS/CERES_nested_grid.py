# ==========================================================================
# Author: Ryan C. Scott
#         ryan.c.scott@nasa.gov
#
# This script grids CERES FOVs to the CERES 1 deg x 1 deg nested grid.
# ==========================================================================

import cerestools as ceres
import matplotlib.pyplot as plt
import numpy as np

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

    # different longitude bin extents for different latitude zones
    lon_ext = [1, 2, 4, 8]

    # loop over nested grid zones
    for i, el in enumerate(lon_ext):

        # set appropriate longitude extent
        if lon_360 is True:                      # lon ranges 0 to 360
            lon_bins = np.arange(0, 361, el)
        elif lon_360 is False:                   # lon ranges -180 to 180
            lon_bins = np.arange(-180, 181, el)  # different bin sizes

        print('Zone {}, lon width {}'.format(i+1, el))
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

        print('Populating CERES 1 deg x 1 deg Nested Grid\n')
        nested_grid = np.empty([180, 360])

        for k in range(180):
            if k == 0:  # zone 5 : pole, 1 box
                nested_grid[k, :] = 0
            elif 1 <= k <= 9:  # zone 4 : 80 to 89 deg, 9 boxes
                nested_grid[k, :] = np.repeat(zone4[k, :], 8)
            elif 10 <= k <= 19:  # zone 3 : 70 to 80 deg, 10 boxes
                nested_grid[k, :] = np.repeat(zone3[k, :], 4)
            elif 20 <= k <= 44:  # zone 2 : 45 to 70 deg, 25 boxes
                nested_grid[k, :] = np.repeat(zone2[k, :], 2)
            elif 45 <= k <= 134:  # zone 1 : -45 to 45 deg, 90 boxes
                nested_grid[k, :] = zone1[k, :]
            elif 135 <= k <= 159:  # zone 2
                nested_grid[k, :] = np.repeat(zone2[k, :], 2)
            elif 160 <= k <= 169:  # zone 3
                nested_grid[k, :] = np.repeat(zone3[k, :], 4)
            elif 170 <= k <= 178:  # zone 4
                nested_grid[k, :] = np.repeat(zone4[k, :], 8)
            elif k == 179:  # zone 5 - pole
                nested_grid[k, :] = 0

    return nested_grid


gridded_field = grid_to_1x1_deg_ceres_nested(lat_data=terra_lat_all,
                                             lon_data=terra_lon_all,
                                             variable=terra_var_all,
                                             lon_360=False)


plt.imshow(gridded_field)
plt.show()

