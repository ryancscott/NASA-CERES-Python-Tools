

import cerestools as ceres

# cmocean
from palettable.cmocean.sequential import Deep_20_r
from palettable.cmocean.diverging import Balance_20
from palettable.cmocean.sequential import Thermal_20
# scientific
from palettable.scientific.sequential import LaJolla_20_r
from palettable.scientific.sequential import LaPaz_20
from palettable.scientific.sequential import Nuuk_20
from palettable.scientific.sequential import Acton_20
from palettable.scientific.sequential import Bilbao_20
from palettable.scientific.sequential import Devon_20
# colorbrewer
from palettable.colorbrewer.sequential import BuPu_9_r

# ===============================================================================================
#
# READ A FULL DAY OF DATA
#
# ===============================================================================================
#
# COMPARE CRS4 to CERES SSF
#
# ===============================================================================================


# var_all3, lon_all3, lat_all3, sza_all3 = ceres.read_day_of_crs_files(
#                                path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
#                                file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20190101',
#                                variable='Longwave flux - upward - clear sky',
#                                lev_arg=0,
#                                fill_nan=False)
#
#
# # CERES SW TOA flux - upwards
# # CERES downward SW surface flux - Model B
# # CERES LW TOA flux - upwards
# #  CERES downward LW surface flux - Model B
# var_all4, lon_all4, lat_all4, sza_all4 = \
#     ceres.read_day_of_ssf_files(path='/Users/rcscott2/Desktop/CRS/ASDC_archive/SSF_Ed4A/',
#                                 file_struc='CER_SSF_Terra-FM1-MODIS_Edition4A_404405.20190101',
#                                 variable='CERES LW TOA flux - upwards',
#                                 fill_nan=False)
#
# # computes var_all3 minus var_all4
# diff = ceres.swath_difference(field2=var_all3, field1=var_all4, day_only=True, sza=sza_all3)
#
#
# diff_colormap = ceres.set_colormap(cmap_name=Balance_20, typ_arg=0)
#
#
# ceres.plot_swath(lon=lon_all3, lat=lat_all3, field=diff, nrows=1, ncols=1, cen_lon=0,
#                  varname='Longwave flux - upward - clear sky', levname='TOA', varunits='Watts per square meter',
#                  cmap=diff_colormap, cmap_lims=(-100, 100), date='20190101', date_str='01/01/2019:00-23h',
#                  nightshade=False, title_str='Terra FM1 CERES CRS Ed4 (computed) minus CERES SSF Ed4A (observed)')


# ceres.swath_histogram_scatterplot(field2=var_all3, field1=var_all4,
#                                   var_name='Longwave flux - upward - total sky',
#                                   lev_name='TOA',
#                                   ti_str1='SSF Model B',
#                                   ti_str2='CRS computed',
#                                   time_info='JAN 1 2019',
#                                   platform='Terra FM1',
#                                   day_only=False, sza=sza_all3)


# ===============================================================================================
#
# COMPARE CRS4 FLUXES for different sky conditions
#
# ===============================================================================================

var_all3, lon_all3, lat_all3, sza_all3 = ceres.read_day_of_crs_files(
                               path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
                               file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20190101',
                               variable='Shortwave flux - downward - total sky no aerosol',
                               lev_arg=5,
                               fill_nan=False)

var_all4, lon_all4, lat_all4, sza_all4 = ceres.read_day_of_crs_files(
                               path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
                               file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20190101',
                               variable='Shortwave flux - downward - pristine sky',
                               lev_arg=5,
                               fill_nan=False)

# computes var_all3 minus var_all4
diff = ceres.swath_difference(field2=var_all3, field1=var_all4, day_only=True, sza=sza_all3)

# set the colormap
diff_colormap = ceres.set_colormap(cmap_name=Balance_20, typ_arg=0)

# plot the swath
ceres.plot_swath(lon=lon_all3, lat=lat_all3, field=diff, nrows=1, ncols=1, cen_lon=0,
                 varname='Shortwave flux - downward', levname='surface', varunits='Watts per square meter',
                 cmap=diff_colormap, cmap_lims=(-100, 100), date='20190101', date_str='01/01/2019:00-23h',
                 nightshade=False, title_str='Terra FM1 CERES CRS Ed4 total-sky no aerosol minus pristine-sky')


# ===============================================================================================
#
# READ A FULL MONTH OF DATA
#
# ===============================================================================================
#
# var_all2, lon_all2, lat_all2, sza_all2 = \
#     read_full_day_of_crs_files(path='/Users/rcscott2/Desktop/CRS/my_output/JAN-01-2019/',
#                                file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.20190101',
#                                variable='Longwave flux - downward - total sky',
#                                lev_arg=5)
#
#
# # CERES LW TOA flux - upwards
# # CERES downward LW surface flux - Model B
# var_all1, lon_all1, lat_all1, sza_all1 =\
#     read_full_day_of_ssf_files(path='/Users/rcscott2/Desktop/CRS/ASDC_archive/',
#                                file_struc='CER_SSF_Terra-FM1-MODIS_Edition4A_404405.20190101',
#                                variable='CERES downward SW surface flux - Model B')


# diff = ceres.swath_difference(field2=var_all4, field1=var_all3, day_only=True, sza=sza_all2)
#
# #var_all2[sza_all2 > 90] = np.nan
#
# lw_colormap = ceres.set_colormap(cmap_name=Balance_20, typ_arg=0)
# ceres.plot_swath(lon=lon_all3, lat=lat_all3, field=diff, nrows=1, ncols=1, cen_lon=0,
#                  varname='Longwave flux - downward - total sky', levname='surface', varunits='Watts per square meter',
#                  cmap=lw_colormap, cmap_lims=(-50, 50), date='20190101', date_str='01/01/2019:00-23h',
#                  nightshade=False, title_str='Terra FM1 CERES CRS Ed4 (computed) minus CERES SSF Ed4A (observed)')


# ceres.swath_histogram_scatterplot(field2=var_all2, field1=var_all1,
#                                   var_name='Shortwave flux - downward',
#                                   lev_name='surface',
#                                   ti_str1='CRS Computed',
#                                   ti_str2='SSF Model B',
#                                   time_info='01/01/2019:00-23h',
#                                   platform='Terra FM-1',
#                                   day_only=True, sza=sza_all2)
#
# var_all2[sza_all2 > 90] = np.nan
# var_all1[sza_all1 > 90] = np.nan
#
# LW colormap - nice options = Deep_20_r, Thermal_20, Bilbao_20
# lw_colormap = ceres.set_colormap(cmap_name=Nuuk_20, typ_arg=0)
# ceres.plot_swath(lon=lon_all2, lat=lat_all2, field=var_all2, nrows=1, ncols=1, cen_lon=0,
#                  varname='Shortwave flux - downward - total sky', levname='TOA', varunits='Watts per square meter',
#                  cmap=lw_colormap, cmap_lims=(0, 1000), date='20190101', date_str='01/01/2019:00-23h',
#                  nightshade=False, title_str='CERES Terra-FM1 CRS Edition 4')
#
# # LW colormap - nice options = Deep_20_r, Thermal_20, Bilbao_20
# lw_colormap = ceres.set_colormap(cmap_name=Nuuk_20, typ_arg=0)
# ceres.plot_swath(lon=lon_all1, lat=lat_all1, field=var_all1, nrows=1, ncols=1, cen_lon=0,
#                  varname='Shortwave flux - downward - total sky', levname='TOA', varunits='Watts per square meter',
#                  cmap=lw_colormap, cmap_lims=(0, 1000), date='20190101', date_str='01/01/2019:00-23h',
#                  nightshade=False, title_str='CERES Terra-FM1 SSF Ed4A parameterized')

