
import cerestools as ceres
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from palettable.cmocean.diverging import Balance_20


path1 = '/Users/rcscott2/Desktop/CRS/ASDC_archive/FLASHFlux/'
path2 = '/Users/rcscott2/Desktop/CRS/ASDC_archive/SSF_Ed4A/'

file1 = 'FLASH_SSF_Terra-FM1-MODIS_Version4A_400400.2019010112'
file2 = 'CER_SSF_Terra-FM1-MODIS_Edition4A_404405.2019010112'

file_path1 = path1 + file1
file_path2 = path2 + file2

# # ===============================================================================================
# # Section 1 - single swath
# # ===============================================================================================
#
# print('============================================')
# print('\t\t\tReading time/date info...\t\t\t')
# print('============================================')
#
# date, date_str = ceres.get_date(file=file_path1)
#
# lat1, lon1, time_obs1, sza1 = ceres.read_ssf_geolocation(file_path=file_path1)
#
# lat2, lon2, time_obs2, sza2 = ceres.read_ssf_geolocation(file_path=file_path2)
#
# print('============================================')
# print('\t\t\tReading data...\t\t\t')
# print('============================================')
#
# field1, var1, units1 = ceres.read_ssf_var(file_path=file_path1,
#                                           var_name='CERES SW TOA flux - upwards',
#                                           fill=True)
#
# field2, var2, units2 = ceres.read_ssf_var(file_path=file_path2,
#                                           var_name='CERES SW TOA flux - upwards',
#                                           fill=True)
#
# vza, vza_name, vza_units = ceres.read_ssf_var(file_path=file_path1,
#                                               var_name='CERES viewing zenith at surface',
#                                               fill=True)
#
# sza, sza_name, sza_units = ceres.read_ssf_var(file_path=file_path1,
#                                               var_name='CERES solar zenith at surface',
#                                               fill=True)
#
# print('============================================')
# print('\t\t\tComputing swath difference...\t\t\t')
# print('============================================')
#
# difference = ceres.swath_difference(field2=field1, field1=field2, day_only=True, sza=sza)
#
# print('============================================')
# print('\t\t\tPlotting swath...\t\t\t')
# print('============================================')
#
# colormap = ceres.set_colormap(cmap_name=Balance_20, typ_arg=0)
#
# title_str = 'FLASHFlux minus CERES SSF'
#
# ceres.plot_swath(lon=lon1, lat=lat1, field=difference, nrows=1, ncols=1, cen_lon=0,
#            varname=var1, levname='', varunits='',
#            cmap=colormap, cmap_lims=(-30, 30), date=date, date_str=date_str,
#            nightshade=True, title_str=title_str)
#
#
# print('============================================')
# print('\t\t\tGenerating 2D histogram...\t\t\t')
# print('============================================')
#
# vza_bins = np.arange(0, 71, 1)
# sza_bins = np.arange(0, 91, 1)
#
# print('VZA bins:')
# print(vza_bins)
# print('SZA bins:')
# print(sza_bins)
#
# print('Length of VZA, SZA:')
# print(len(vza))
# print(len(sza))
#
# gridded = stats.binned_statistic_2d(vza, sza, difference, statistic=np.nanmean, bins=[vza_bins, sza_bins])
# gridded_stat = np.flipud(np.rot90(gridded.statistic))
#
# # Plot of the result
# plt.pcolor(gridded_stat, vmin=-10, vmax=10)
# plt.ylabel('SZA')
# plt.xlabel('VZA')
# plt.title('Average ' + 'TOA SW flux difference \n (FLASHFlux - CERES SSF) binned by VZA, SZA \n' + date_str)
# # plt.title('Number of CERES FOVs binned by VZA, SZA')
# plt.colorbar()
# plt.show()

# ===============================================================================================
# Section 2 - full day (well, 23 hr) of swaths
# ===============================================================================================


print('============================================')
print('\t\t\tReading day of FF & CERES data...\t\t ')
print('============================================')


var_all1, lon_all1, lat_all1, sza_all1 = ceres.read_day_of_ssf_files(
                                        path=path1,
                                        file_struc='FLASH_SSF_Terra-FM1-MODIS_Version4A_400400.20190101',
                                        variable='CERES SW TOA flux - upwards',
                                        index=-1,
                                        fill=True)

var_all2, lon_all2, lat_all2, sza_all2 = ceres.read_day_of_ssf_files(
                                        path=path2,
                                        file_struc='CER_SSF_Terra-FM1-MODIS_Edition4A_404405.20190101',
                                        variable='CERES SW TOA flux - upwards',
                                        index=-1,
                                        fill=True)

vza_all3, lon_all3, lat_all3, sza_all3 = ceres.read_day_of_ssf_files(
                                        path=path2,
                                        file_struc='CER_SSF_Terra-FM1-MODIS_Edition4A_404405.20190101',
                                        variable='CERES viewing zenith at surface',
                                        index=-1,
                                        fill=True)

diff = ceres.swath_difference(field2=var_all1, field1=var_all2, day_only=True, sza=sza_all2)

# # Quick & dirty plot of the difference
# plt.plot(diff)
# plt.title('Difference, daytime only')
# plt.show()

print('============================================')
print('\t\t\tGetting daytime data only...\t\t\t')
print('============================================')

# # Day time only
# # (SZA 86.5 deg cutoff used since ADMs can't be
# # applied beyond this SZA)
for i in range(len(sza_all3)):
    if sza_all3[i] >= 86.5:
        sza_all3[i] = np.nan
        vza_all3[i] = np.nan
        diff[i] = np.nan
        lat_all3[i] = np.nan
        lon_all3[i] = np.nan

# ignore/remove NaNs
bad_indices = np.isnan(diff)
good_indices = ~bad_indices
lat_all3 = lat_all3[good_indices]
lon_all3 = lon_all3[good_indices]
sza_all3 = sza_all3[good_indices]
vza_all3 = vza_all3[good_indices]
diff = diff[good_indices]

print('Length of SZA, VZA, diff:')
print(len(sza_all3))
print(len(vza_all3))
print(len(diff))

print('Length of lon, lat:')
print(len(lon_all3))
print(len(lat_all3))

print('============================================')
print('\t\t\tPlotting swaths...\t\t\t')
print('============================================')

diff_colormap = ceres.set_colormap(cmap_name=Balance_20, typ_arg=0)

ceres.plot_swath(lon=lon_all3, lat=lat_all3, field=sza_all3,
                 nrows=1, ncols=1, cen_lon=0,
                 varname='CERES solar zenith angle', levname='surface',
                 varunits='degrees',
                 cmap=diff_colormap, cmap_lims=(0, 90), date='20190101', date_str='01/01/2019:00-22h',
                 nightshade=False, title_str='CERES Terra FM1 Ed4A')

ceres.plot_swath(lon=lon_all3, lat=lat_all3, field=vza_all3,
                 nrows=1, ncols=1, cen_lon=0,
                 varname='CERES viewing zenith angle', levname='surface',
                 varunits='degrees',
                 cmap=diff_colormap, cmap_lims=(0, 70), date='20190101', date_str='01/01/2019:00-22h',
                 nightshade=False, title_str='CERES Terra FM1 Ed4A')

ceres.plot_swath(lon=lon_all3, lat=lat_all3, field=diff,
                 nrows=1, ncols=1, cen_lon=0,
                 varname='CERES SW TOA flux - upwards', levname='',
                 varunits='Watts per square meter',
                 cmap=diff_colormap, cmap_lims=(-30, 30), date='20190101', date_str='01/01/2019:00-22h',
                 nightshade=False, title_str='Terra FM1 FLASHFlux V4A minus CERES SSF Ed4A')

# # Quick & dirty plot of each time series
# plt.plot(diff)
# plt.title('TOA flux difference')
# plt.show()
#
# plt.plot(sza_all3)
# plt.title('SZA')
# plt.show()
#
# plt.plot(vza_all3)
# plt.title('VZA')
# plt.show()

ceres.swath_histogram_scatterplot(field2=var_all1, field1=var_all2,
                                  var_name='CERES SW TOA flux - upwards',
                                  lev_name='',
                                  suptitle=' ',
                                  ti_str2='FLASHFlux',
                                  ti_str1='CERES SSF',
                                  date_str='01/01/2019:00-22h',
                                  platform='Terra FM1',
                                  day_only=True, sza=sza_all2)

print('============================================')
print('\t\t\tGenerating 2D histogram...\t\t\t')
print('============================================')

vza_bins = np.arange(0, 71, 1)
sza_bins = np.arange(0, 91, 1)

print('VZA bins')
print(vza_bins)
print('SZA bins')
print(sza_bins)

print('Length of VZA, SZA, Difference')
print(len(vza_all3))
print(len(sza_all3))
print(len(diff))

# Compute mean diff binned by VZA, SZA
gridded = stats.binned_statistic_2d(vza_all3, sza_all3, diff, statistic=np.nanmean, bins=[vza_bins, sza_bins])
gridded_stat = np.flipud(np.rot90(gridded.statistic))

# Plot of 2D histogram
plt.pcolor(gridded_stat, vmin=-5, vmax=5)
plt.ylabel('SZA')
plt.xlabel('VZA')
plt.title('Average ' + 'TOA SW flux difference \n (FLASHFlux - CERES SSF) binned by VZA, SZA \n JAN-01-2019:00-22h')
# plt.title('Number of CERES FOVs binned by VZA, SZA')
plt.colorbar()
plt.show()

