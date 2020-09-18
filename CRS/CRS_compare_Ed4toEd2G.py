

import cerestools as ceres
from palettable.cmocean.diverging import Balance_20
from palettable.cmocean.sequential import Thermal_20
from palettable.scientific.sequential import LaPaz_20


print('============================================')
print('\t\t\tCRS official file...\t\t\t')
print('============================================')


path1 = '/Users/rcscott2/Desktop/CERES/CRS/'
file1 = 'CER_CRS_Terra-FM1-MODIS_Edition2G_023034.2010062023'
file_path1 = path1 + file1
print(file_path1)

ceres.get_platform(file1)


print('============================================')
print('\t\t\tCRS development file...\t\t\t')
print('============================================')


file2 = 'CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.2010062023'
path2 = '/Users/rcscott2/Desktop/CRS/my_output/First_runs/'
file_path2 = path2 + file2
print(file_path2)


print('============================================')
print('\t\t\tReading geolocation...\t\t\t')
print('============================================')


lat, lon, p_levels, time_obs, sza = ceres.read_crs_geolocation(file_path=file_path1)


print('============================================')
print('\t\t\tReading data...\t\t\t')
print('============================================')

field0, var0, units0, lev0 = ceres.read_crs_var(file_path=file_path1,
                                                var_name='Shortwave flux adjustment at TOA - upward - total',
                                                lev_arg=-1, fill=True)

field1, var1, units1, lev1 = ceres.read_crs_var(file_path=file_path1,
                                                var_name='Longwave flux - downward - total',
                                                lev_arg=4, fill=True)

field1 = field1# -field0

field2, var2, units2, lev2 = ceres.read_crs_var_dev(file_path=file_path2,
                                                    var_name='Longwave flux - downward - total sky',
                                                    lev_arg=5, fill=True)

field3, var3, units3, lev3 = ceres.read_crs_var_dev(file_path=file_path2,
                                                    var_name='Cloud temperature',
                                                    lev_arg=1, fill=True)


print('============================================')
print('\t\t\tReading time/date info...\t\t\t')
print('============================================')

date, date_str = ceres.get_date(file=file1)

print('============================================')
print('\t\t\tComputing difference...\t\t\t')
print('============================================')

difference, _ = ceres.swath_difference(field2=field2, field1=field1, day=-1, sza=sza, lat=lat)

print('============================================')
print('\t\t\tSetting colormap...\t\t\t')
print('============================================')

colormap = ceres.set_colormap(cmap_name=Balance_20, typ_arg=0)
colormap1 = ceres.set_colormap(cmap_name=LaPaz_20, typ_arg=0)

print('============================================')
print('\t\t\tPlotting data...\t\t\t')
print('============================================')

title_str = 'CERES ' + ceres.get_platform(file1) + ' CRS Ed2G'
cen_lon = 180

ceres.plot_swath(lon=lon, lat=lat, field=field1, nrows=1, ncols=1, cen_lon=cen_lon,
                 varname=var1, levname=lev1, varunits=units1,
                 cmap=colormap, cmap_lims=(150, 800), date=date, date_str=date_str,
                 nightshade=True, title_str=title_str,
                 marker='o', msize=1)

title_str = 'CERES ' + ceres.get_platform(file1) + ' CRS EdX (new)'

ceres.plot_swath(lon=lon, lat=lat, field=field2, nrows=1, ncols=1, cen_lon=cen_lon,
                 varname=var2, levname=lev2, varunits='Watts per square meter',
                 cmap=colormap, cmap_lims=(150, 800), date=date, date_str=date_str,
                 nightshade=True, title_str=title_str,
                 marker='o', msize=1)


title_str = 'CERES ' + ceres.get_platform(file1) + ' latest CRS'

ceres.plot_swath(lon=lon, lat=lat, field=field3, nrows=1, ncols=1, cen_lon=cen_lon,
                 varname='Cloud temperature', levname='', varunits='K',
                 cmap=colormap, cmap_lims=(200, 300), date=date, date_str=date_str,
                 nightshade=True, title_str=title_str,
                 marker='o', msize=1)


# print("Plotting ", var1, "difference")
# a, b = input("Enter colormap limits: ").split(',')
# print("Specified colormap range [{}, {}]".format(a, b))
# cmap_lim = (float(a), float(b))

colormap = ceres.set_colormap(cmap_name=Balance_20, typ_arg=0)

#title_str = 'CERES ' + ceres.get_platform(file1) + r' calculated minus observed flux ($\Delta$)'
title_str = 'CERES ' + ceres.get_platform(file1) + r' - New CRS minus CRS Ed2G ($\Delta$, untuned fluxes)'

ceres.plot_swath(lon=lon, lat=lat, field=difference, nrows=1, ncols=1, cen_lon=cen_lon,
                 varname='Downward Longwave Radiation', levname=lev1, varunits=units1,
                 cmap=colormap, cmap_lims=(-50, 50), date=date, date_str=date_str,
                 nightshade=True, title_str=title_str,
                 marker='o', msize=1)

# ceres.swath_histogram_scatterplot(field2=field1, field1=field2,
#                                   var_name=var1,
#                                   lev_name=lev1,
#                                   f1_str='',
#                                   f2_str='',
#                                   date_str=date_str,
#                                   platform=ceres.get_platform(file1),
#                                   day=False, sza=sza)

