

import cerestools as ceres
from palettable.cmocean.diverging import Balance_20


print('============================================')
print('\t\t\tFLASHFlux official file...\t\t\t')
print('============================================')


path1 = '/Users/rcscott2/Desktop/CERES/FLASHFlux/'
file1 = 'FLASH_SSF_Aqua-FM3-MODIS_Version3C_232103.2020020523'
file_path1 = path1 + file1
print(file_path1)


# ========================================================================


print('============================================')
print('\t\t\tReading geolocation...\t\t\t')
print('============================================')

lat, lon, time_obs, sza = ceres.read_ssf_geolocation(file_path=file_path1)

print('============================================')
print('\t\t\tReading data...\t\t\t')
print('============================================')

field1, var1, units1 = ceres.read_ssf_var(file_path=file_path1,
                                          var_name='CERES downward SW surface flux - Model B',
                                          fill=True)

print('============================================')
print('\t\t\tReading time/date info...\t\t\t')
print('============================================')

date, date_str = ceres.get_date(file=file1)

print('============================================')
print('\t\t\tComputing difference...\t\t\t')
print('============================================')

#difference = compute_diff(field2=field2, field1=field1)

print('============================================')
print('\t\t\tSetting colormap...\t\t\t')
print('============================================')

colormap = ceres.set_colormap(cmap_name=Balance_20, typ_arg=0)

print('============================================')
print('\t\t\tPlotting data...\t\t\t')
print('============================================')

a, b = input("Enter colormap limits: ").split(',')
print("First number is {} and second number is {}".format(a, b))
cmap_lim = (float(a), float(b))

title_str = 'CERES FLASHFlux Edition 3C'

ceres.plot_swath(lon=lon, lat=lat, field=field1,
                 varname=var1, levname="", varunits=units1,
                 nrows=1, ncols=1, cen_lon=0,
                 cmap=colormap, cmap_lims=cmap_lim, date=date,
                 nightshade=1, date_str=date_str, title_str=title_str)

