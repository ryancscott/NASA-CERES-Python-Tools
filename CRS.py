

import cerestools as ceres
from palettable.cmocean.diverging import Balance_20


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


file2 = 'CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.2010062023'
path2 = '/Users/rcscott2/Desktop/CRS/my_output/'
file_path2 = path2 + file2
print(file_path2)


print('============================================')
print('\t\t\tReading geolocation...\t\t\t')
print('============================================')


lat, lon, p_levels, time_obs, sza = ceres.read_crs_geolocation(file_path=file_path1)


print('============================================')
print('\t\t\tReading data...\t\t\t')
print('============================================')

field1, var1, units1, lev1 = ceres.read_crs_var(file_path=file_path1, vararg=0, levarg=4, fill=1)
field2, var2, units2, lev2 = ceres.read_crs_var_dev(file_path=file_path2, vararg=0, levarg=1, fill=1)

print('============================================')
print('\t\t\tReading time/date info...\t\t\t')
print('============================================')

date, date_str = ceres.get_date(file=file1)

print('============================================')
print('\t\t\tComputing difference...\t\t\t')
print('============================================')

difference = ceres.swath_difference(field2=field2, field1=field1, day_only=0, sza=sza)

print('============================================')
print('\t\t\tSetting colormap...\t\t\t')
print('============================================')

colormap = ceres.set_colormap(cmap_name=Balance_20, typarg=0)

print('============================================')
print('\t\t\tPlotting data...\t\t\t')
print('============================================')

title_str = 'CERES ' + ceres.get_platform(file1) + ' CRS Ed2G'

ceres.plot_swath(lon=lon, lat=lat, field=field1, nrows=1, ncols=1, cen_lon=0,
           varname=var1, levname=lev1, varunits=units1,
           cmap=colormap, cmap_lims=(0, 500), date=date, date_str=date_str,
           nightshade=1, title_str=title_str)

title_str = 'CERES ' + ceres.get_platform(file1) + ' CRS Ed4'

ceres.plot_swath(lon=lon, lat=lat, field=field2, nrows=1, ncols=1, cen_lon=0,
           varname=var1, levname=lev1, varunits=units1,
           cmap=colormap, cmap_lims=(0, 500), date=date, date_str=date_str,
           nightshade=1, title_str=title_str)


# print("Plotting ", var1, "difference")
# a, b = input("Enter colormap limits: ").split(',')
# print("Specified colormap range [{}, {}]".format(a, b))
# cmap_lim = (float(a), float(b))


title_str = 'CERES ' + ceres.get_platform(file1) + r' CRS Ed4 - Ed2G difference ($\Delta$)'


ceres.plot_swath(lon=lon, lat=lat, field=difference, nrows=1, ncols=1, cen_lon=0,
           varname=var1, levname=lev1, varunits=units1,
           cmap=colormap, cmap_lims=(-100, 100), date=date, date_str=date_str,
           nightshade=1, title_str=title_str)

ceres.swath_histogram_scatterplot(field1, field2, var1, lev1, date_str, ceres.get_platform(file1), day_only=0, sza=sza)
