from palettable.cmocean.diverging import Balance_20
import cerestools as ceres

print('============================================')
print('\tCRS File...\t\t\t')
print('============================================')


# Aqua-FM3
#file = 'CER_CRS4_Aqua-FM3-MODIS_GH4_2222TH.2019010100'

# Terra-FM1
#file = 'CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.2010062023'
#file = 'CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.2019010100'
#file = 'CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.2019011200'
file = 'CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.2019011212'

path = '/Users/rcscott2/Desktop/CRS/my_output/'
file_path = path + file

print(file_path)

print('============================================')
print('\tReading Data...\t\t\t')
print('============================================')


# Geolocation
lat, lon, pres_levs, obs_tim = ceres.read_crs_geolocation_dev(file_path)

# Top-of-atmosphere fluxes
lwut, lwut_name, lwut_units, lwut_lev = ceres.read_crs_var_dev(file_path=file_path, vararg=1, levarg=0, fill=1)
swdt, swdt_name, swdt_units, swdt_lev = ceres.read_crs_var_dev(file_path=file_path, vararg=2, levarg=0, fill=1)
swut, swut_name, swut_units, swut_lev = ceres.read_crs_var_dev(file_path=file_path, vararg=3, levarg=0, fill=1)

# Surface fluxes
lwds, lwds_name, lwds_units, lwds_lev = ceres.read_crs_var_dev(file_path=file_path, vararg=0, levarg=1, fill=1)
lwus, lwus_name, lwus_units, lwus_lev = ceres.read_crs_var_dev(file_path=file_path, vararg=1, levarg=1, fill=1)
swds, swds_name, swds_units, swds_lev = ceres.read_crs_var_dev(file_path=file_path, vararg=2, levarg=1, fill=1)
swus, swus_name, swus_units, swus_lev = ceres.read_crs_var_dev(file_path=file_path, vararg=3, levarg=1, fill=1)

# No particular level OR column integrated parameters
aot, aot_name, aot_units, aot_lev = ceres.read_crs_var_dev(file_path=file_path, vararg=4, levarg=-1, fill=1)
hgts, hgts_name, hgts_units, hgts_lev = ceres.read_crs_var_dev(file_path=file_path, vararg=5, levarg=-1, fill=1)

# Date info
date, date_str = ceres.get_date(file=file)

print('============================================')
print('\tPlotting Data...\t\t\t')
print('============================================')

cmap = ceres.set_colormap(cmap_name=Balance_20, typarg=0)

title_str = 'CERES Cloud Radiative Swath (CRS) Ed4 Development'

ceres.plot_swath(lon=lon, lat=lat, field=hgts, nrows=1, ncols=1, cen_lon=0,
           varname=hgts_name, levname=hgts_lev, varunits=hgts_units,
           cmap=cmap, cmap_lims=(0, 4000), date=date, date_str=date_str,
           nightshade=1, title_str=title_str)