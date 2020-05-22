import cerestools as ceres
from palettable.cmocean.diverging import Balance_20

path1 = '/Users/rcscott2/Desktop/CRS/ASDC_archive/FLASHFlux/'
path2 = '/Users/rcscott2/Desktop/CRS/ASDC_archive/SSF_Ed4A/'

file1 = 'FLASH_SSF_Terra-FM1-MODIS_Version4A_400400.2019010112'
file2 = 'CER_SSF_Terra-FM1-MODIS_Edition4A_404405.2019010112'

file_path1 = path1 + file1
file_path2 = path2 + file2

print('============================================')
print('\t\t\tReading time/date info...\t\t\t')
print('============================================')

date, date_str = ceres.get_date(file=file_path1)


lat1, lon1, time_obs1, sza1 = ceres.read_ssf_geolocation(file_path=file_path1)

lat2, lon2, time_obs2, sza2 = ceres.read_ssf_geolocation(file_path=file_path2)

print('============================================')
print('\t\t\tReading data...\t\t\t')
print('============================================')

field1, var1, units1 = ceres.read_ssf_var(file_path=file_path1,
                                          var_name='CERES downward SW surface flux - Model B',
                                          fill=True)

field2, var2, units2 = ceres.read_ssf_var(file_path=file_path2,
                                          var_name='CERES downward SW surface flux - Model B',
                                          fill=True)

print('============================================')
print('\t\t\tComputing difference...\t\t\t')
print('============================================')

difference = ceres.swath_difference(field2=field1, field1=field2, day_only=False, sza=sza1)


print('============================================')
print('\t\t\tSetting colormap...\t\t\t')
print('============================================')

colormap = ceres.set_colormap(cmap_name=Balance_20, typ_arg=0)

title_str = 'FLASHFlux minus CERES SSF'

# ceres.plot_swath(lon=lon1, lat=lat1, field=difference, nrows=1, ncols=1, cen_lon=0,
#            varname=var1, levname='', varunits='',
#            cmap=colormap, cmap_lims=(-100, 100), date=date, date_str=date_str,
#            nightshade=True, title_str=title_str)


print('============================================')
print('\t\t\tREADING DAY OF FF + SSF data...\t\t\t')
print('============================================')


var_all1, lon_all1, lat_all1, sza_all1 = ceres.read_day_of_ssf_files(
                                        path=path1,
                                        file_struc='FLASH_SSF_Terra-FM1-MODIS_Version4A_400400.20190101',
                                        variable='CERES LW TOA flux - upwards',
                                        fill_nan=True)

var_all2, lon_all2, lat_all2, sza_all2 = ceres.read_day_of_ssf_files(
                                        path=path2,
                                        file_struc='CER_SSF_Terra-FM1-MODIS_Edition4A_404405.20190101',
                                        variable='CERES LW TOA flux - upwards',
                                        fill_nan=True)

diff = ceres.swath_difference(field2=var_all1, field1=var_all2, day_only=True, sza=sza_all2)

print(diff)

diff_colormap = ceres.set_colormap(cmap_name=Balance_20, typ_arg=0)

ceres.plot_swath(lon=lon_all2, lat=lat_all2, field=diff,
                 nrows=1, ncols=1, cen_lon=0,
                 varname='CERES LW TOA flux - upwards', levname='difference',
                 varunits=' ',
                 cmap=diff_colormap, cmap_lims=(-30, 30), date='20190101', date_str='01/01/2019:00-22h',
                 nightshade=False, title_str='Terra FM1 FLASHFlux V4A minus CERES SSF Ed4A')


ceres.swath_histogram_scatterplot(field2=var_all1, field1=var_all2,
                                  var_name='CERES LW TOA flux - upwards',
                                  lev_name='',
                                  suptitle='FF vs SSF',
                                  ti_str2='FLASHFlux',
                                  ti_str1='CERES SSF',
                                  time_info='01/01/2019:00-22h',
                                  platform='Terra FM1',
                                  day_only=True, sza=sza_all2)




