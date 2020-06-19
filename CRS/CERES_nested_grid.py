# ==========================================================================
# Author: Ryan C. Scott
#         ryan.c.scott@nasa.gov
#
# This script is for griding CRS FOVs to CERES NESTED grid
# ==========================================================================

import cerestools as ceres

terra_var_all, terra_lon_all, terra_lat_all, _ = ceres.read_day_of_crs_files(
                               path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
                               file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20190101',
                               variable='Shortwave flux - upward - total sky',
                               lev_arg=0,
                               fill=True)

terra_lon_all = ceres.swath_lon_360_to_180(terra_lon_all)


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
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats

    lat_bins = np.arange(-90, 91)  # lat: -90 to 90 deg

    # consider generalizing and adding bins as parameters of the function?
    lat_bins = np.arange(-90, 91)    # lat: -90 to 90 deg

    lon_ext = [1, 2, 4, 8]

    # gridded1 = np.empty([180, 360])
    # gridded2 = np.empty([180, 180])
    # gridded3 = np.empty([180, 90])
    # gridded4 = np.empty([180, 45])
    for el in lon_ext:

        if lon_360 is True:
         lon_bins = np.arange(0, 361, el)     # lon: 0 to 360
        elif lon_360 is False:
         lon_bins = np.arange(-180, 181, el)

        if el == 1:
            variable[abs(lat_data) > 45] = np.nan

            gridded1 = stats.binned_statistic_2d(lon_data, lat_data, variable,
                                                 statistic=np.nanmean,
                                                 bins=[lon_bins, lat_bins])
            gridded_stat1 = gridded1.stat

        elif el == 2:
            variable[abs(lat_data) < 45] = np.nan
            variable[abs(lat_data) > 70] = np.nan

            gridded2 = stats.binned_statistic_2d(lon_data, lat_data, variable,
                                                 statistic=np.nanmean,
                                                 bins=[lon_bins, lat_bins])
            gridded_stat2 = gridded2.stat

        elif el == 4:
            variable[abs(lat_data) < 70] = np.nan
            variable[abs(lat_data) > 80] = np.nan

            gridded3 = stats.binned_statistic_2d(lon_data, lat_data, variable,
                                                 statistic=np.nanmean,
                                                 bins=[lon_bins, lat_bins])
            gridded_stat3 = gridded3.stat

        elif el == 8:
            variable[abs(lat_data) < 80] = np.nan
            variable[abs(lat_data) > 89] = np.nan

            gridded4 = stats.binned_statistic_2d(lon_data, lat_data, variable,
                                                 statistic=np.nanmean,
                                                 bins=[lon_bins, lat_bins])
            gridded_stat4 = gridded4.stat

    return gridded_stat1, gridded_stat2, gridded_stat3, gridded_stat4


grid_to_1x1_deg_ceres_nested(lat_data=terra_lat_all,
                             lon_data=terra_lon_all,
                             variable=terra_var_all,
                             lon_360=False)

# -------------------------------
# ZONE 1
# if lon_360 is True:
#     lon_bins = np.arange(0, 361)     # lon: 0 to 360
# elif lon_360 is False:
#     lon_bins = np.arange(-180, 181)  # lon: -180 to 180

# zone 1 : equatorward of 45 - set to np.nan
# variable[abs(lat_data) > 45] = np.nan

# -------------------------------
# # ZONE 2
# if lon_360 is True:
#     lon_bins = np.arange(0, 361, 2)  # lon: 0 to 360
# elif lon_360 is False:
#     lon_bins = np.arange(-180, 181, 2)  # lon: -180 to 180
#
# # zone 2 : use 2 deg lon bins between 45 and 70 - set to np.nan
# variable[abs(lat_data) < 45] = np.nan
# variable[abs(lat_data) > 70] = np.nan

# -------------------------------
# ZONE 3
# if lon_360 is True:
#     lon_bins = np.arange(0, 361, 4)  # lon: 0 to 360
# elif lon_360 is False:
#     lon_bins = np.arange(-180, 181, 4)  # lon: -180 to 180
#
# # zone 3 : use 2 deg lon bins between 70 and 80 - set to np.nan
# variable[abs(lat_data) < 70] = np.nan
# variable[abs(lat_data) > 80] = np.nan

# -------------------------------
# ZONE 4
# if lon_360 is True:
#     lon_bins = np.arange(0, 361, 8)  # lon: 0 to 360
# elif lon_360 is False:
#     lon_bins = np.arange(-180, 181, 8)  # lon: -180 to 180

# zone 4 : use 2 deg lon bins between 70 and 80 - set to np.nan
# variable[abs(lat_data) < 80] = np.nan
# variable[abs(lat_data) > 89] = np.nan

# print('Griding and averaging footprints to 1 x 1',
#       'degree CERES nested grid...\n')
#
# print('Lat bins:\n')
# print(lat_bins)
# print('Lon bins:\n')
# print(lon_bins)

# gridded1 = stats.binned_statistic_2d(lon_data, lat_data, variable,
#                                      statistic=np.nanmean,
#                                      bins=[lon_bins, lat_bins])
# gridded_stat4 = gridded4.stat


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