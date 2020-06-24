# ==========================================================================
# Author: Ryan C. Scott
#         ryan.c.scott@nasa.gov
#
# This script grids CERES FOVs to the CERES 1 deg x 1 deg nested grid.
# ==========================================================================

import cerestools as ceres
import matplotlib.pyplot as plt

terra_var_all, terra_lon_all, terra_lat_all, terra_sza_all =\
    ceres.read_day_of_crs_files(
        path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
        file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20190101',
        variable='Longwave flux - upward - total sky',
        lev_arg=0,
        fill=True
    )

terra_lon_all = ceres.swath_lon_360_to_180(lon=terra_lon_all)

# terra_lat_all, terra_lon_all, terra_var_all, terra_sza_all = \
#     ceres.swath_daytime_only(
#         lat=terra_lat_all,
#         lon=terra_lon_all,
#         var=terra_var_all,
#         sza=terra_sza_all,
#         sza_cutoff=90
#     )

gridded_field = ceres.grid_to_1x1_deg_ceres_nested(
    lat_data=terra_lat_all,
    lon_data=terra_lon_all,
    variable=terra_var_all,
    lon_360=False
)


plt.imshow(gridded_field)
plt.show()

