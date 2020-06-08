
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
# for k in range(24):
#     plt.imshow(toa_lw_up_all[k, :, :])
#     plt.colorbar()
#     plt.show()

lat_syn1deg, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
                                           var_name='latitude',
                                           fill=False)
lon_syn1deg, _, _ = ceres.read_syn1deg_hdf(file_path=file_path2,
                                           var_name='longitude',
                                           fill=False)

# CERES Ed4 CRS data
var_all, lon_all, lat_all, sza_all = ceres.read_day_of_crs_files(
                               path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
                               file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20190101',
                               variable='Longwave flux - upward - total sky',
                               lev_arg=0,
                               fill=True)


# convert longitude range from 0 to 360 deg to -180 to 180 deg
for i in range(len(lon_all)):
    if lon_all[i] > 180:
        lon_all[i] = lon_all[i] - 360


# grid and average footprints in 1 deg x 1 deg grid boxes
print('Gridding and averaging footprints to 1 x 1 degree grid boxes')
tic = time.time()
gridded_var_all = ceres.grid_to_1x1_degree(lat_all, lon_all, var_all)
toc = time.time()
print(toc-tic, 'seconds elapsed during ceres.grid_to_1x1_degree...')

# select colormap
cmap = ceres.set_colormap(Thermal_20, 0)

# plot the gridded field
ceres.global_map(lon=lon_syn1deg, lat=lat_syn1deg, field=gridded_var_all,
                 varname='CERES LW TOA flux - upwards', varunits=r'W m$^{-2}$',
                 nrows=1, ncols=1, cen_lon=0,
                 cmap=cmap, cmap_lims=(150, 350), ti_str='CERES Terra FM1 "CRS1deg" JAN-01-2019')


print("Shape of 1 deg x 1 deg gridded field:\n")
print(gridded_var_all.shape)


plt.imshow(gridded_var_all)
plt.colorbar()
plt.title(r'Terra CRS1deg$_{\beta}$ OLR [W m$^{-2}$] JAN-1-2019')
plt.clim(150, 350)
plt.show()


plt.imshow(np.nanmean(toa_lw_up_all, axis=0))
plt.colorbar()
plt.title('SYN1deg OLR [W m$^{-2}$] JAN-1-2019')
plt.clim(150, 350)
plt.show()


# plot of the difference
plt.imshow(gridded_var_all-np.nanmean(toa_lw_up_all, axis=0))
plt.colorbar()
plt.title(r'Terra CRS1deg$_{\beta}$ - SYN1deg OLR [W m$^{-2}$] JAN-1-2019')
plt.clim(-50, 50)
plt.show()


















