

import time
import cerestools as ceres
from palettable.cmocean.sequential import Thermal_20
from palettable.colorbrewer.sequential import BuPu_9_r
import matplotlib.pyplot as plt

# ==============================================================================

path = '/Users/rcscott2/Desktop/CERES/EBAF/'
file = 'CERES_EBAF-TOA_Ed4.1_Subset_200003-201910.nc'
file_path = path + file


# ceres.print_nc_file_info(file_path)


# Read lat, lon from EBAF file
lat, lat_name, lat_units = ceres.read_ebaf_var(file_path=file_path, variable='lat')
lon, lon_name, lon_units = ceres.read_ebaf_var(file_path=file_path, variable='lon')


# CERES Ed4A SSF data
var_all4, lon_all4, lat_all4, sza_all4 = \
     ceres.read_day_of_ssf_files(path='/Users/rcscott2/Desktop/CRS/ASDC_archive/SSF_Ed4A/',
                                 file_struc='CER_SSF_Terra-FM1-MODIS_Edition4A_404405.20190101',
                                 variable='CERES LW TOA flux - upwards',
                                 index=-1,
                                 fill=True)

# CERES Ed4 CRS data
# var_all4, lon_all4, lat_all4, sza_all4 = ceres.read_day_of_crs_files(
#                                path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
#                                file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20190101',
#                                variable='Shortwave flux - upward - total sky',
#                                lev_arg=0,
#                                fill=True)


# Grid and average footprints into 1 deg x 1 deg grid boxes
gridded_var_all4 = ceres.grid_to_1x1_deg_equal_angle(lat_all4, lon_all4, var_all4)


# select colormap
cmap = ceres.set_colormap(Thermal_20, 0)

# plot the gridded field
ceres.global_map(lon=lon, lat=lat, field=gridded_var_all4,
                 varname='CERES LW TOA flux - upwards', varunits=r'W m$^{-2}$',
                 nrows=1, ncols=1, cen_lon=0,
                 cmap=cmap, cmap_lims=(150, 350), ti_str='CERES Terra FM1 "CRS1deg" JAN-01-2019')

plt.imshow(gridded_var_all4)
plt.colorbar()
plt.clim(vmin=150, vmax=350)
plt.show()


