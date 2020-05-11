

import cerestools as ceres
import numpy as np


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


def read_full_day_of_crs_files(path, file_struc, variable, lev_arg):
    """

    :param path:
    :param file_struc:
    :param variable:
    :return:
    """
    print('============================================')
    print('\tReading CRS Files...\t\t\t')
    print('============================================')

    len_tot = []
    sza_all = np.empty([])
    lat_all = np.empty([])
    lon_all = np.empty([])
    var_all = np.empty([])

    for k in range(24):
        if k < 10:
            k = '0' + str(k)

        file = file_struc + str(k)

        file_path = path + file
        print(file_path)

        lat, lon, pres_levs, obs_tim, sfc_ind, sza = ceres.read_crs_geolocation_dev(file_path)

        var, var_name, var_units, var_lev = \
            ceres.read_crs_var_dev(file_path=file_path,
                                   var_name=variable,
                                   lev_arg=lev_arg, fill=0)

        len_tot.append(lat.shape[0])
        sza_all = np.concatenate((sza_all, sza), axis=None)
        lat_all = np.concatenate((lat_all, lat), axis=None)
        lon_all = np.concatenate((lon_all, lon), axis=None)
        var_all = np.concatenate((var_all, var), axis=None)

    print(len_tot)
    print(var_all.shape)

    # # LW colormap - nice options = Deep_20_r, Thermal_20, Bilbao_20
    # lw_colormap = ceres.set_colormap(cmap_name=Deep_20_r, typ_arg=0)
    #
    # ceres.plot_swath(lon=lon_all, lat=lat_all, field=var_all, nrows=1, ncols=1, cen_lon=0,
    #              varname=var_name, levname=var_lev, varunits=var_units,
    #              cmap=lw_colormap, cmap_lims=(150, 450), date='20190101', date_str='01/01/2019:00-23h',
    #              nightshade=False, title_str='CERES Terra CRS Edition 4')

    return var_all, lon_all, lat_all, sza_all


# ============================================================================


def read_full_day_of_ssf_files(path, file_struc, variable):
    print('============================================')
    print('\tReading SSF Files...\t\t\t')
    print('============================================')

    len_tot = []
    sza_all = np.empty([])
    tim_all = np.empty([])
    lat_all = np.empty([])
    lon_all = np.empty([])
    var_all = np.empty([])

    for k in range(24):
        if k < 10:
            k = '0' + str(k)

        file = file_struc + str(k)

        file_path = path + file
        print(file_path)

        # read geolocation info from file
        lat, lon, tim, sza = ceres.read_ssf_geolocation(file_path)
        # read variable from file
        var, var_name, var_units = \
            ceres.read_ssf_var(file_path=file_path, var_name=variable, fill=True)

        # len_tot contains the length (num FOVs) of each individual swath
        len_tot.append(lat.shape[0])

        # all of the swaths combined
        sza_all = np.concatenate((sza_all, sza), axis=None)
        tim_all = np.concatenate((tim_all, tim), axis=None)
        lat_all = np.concatenate((lat_all, lat), axis=None)
        lon_all = np.concatenate((lon_all, lon), axis=None)
        var_all = np.concatenate((var_all, var), axis=None)

    print(len_tot)
    print(var_all.shape)

    # # LW colormap - nice options = Deep_20_r, Thermal_20, Bilbao_20
    # lw_colormap = ceres.set_colormap(cmap_name=Deep_20_r, typ_arg=0)
    #
    # ceres.plot_swath(lon=lon_all, lat=lat_all, field=var_all, nrows=1, ncols=1, cen_lon=0,
    #              varname=var_name, levname='', varunits=var_units,
    #              cmap=lw_colormap, cmap_lims=(150, 450), date=date, date_str='01/01/2019:00-23h',
    #              nightshade=False, title_str='CERES Terra SSF Edition 4A')

    return var_all, lon_all, lat_all, sza_all


# =============================================================================


def read_full_month_of_crs_files(path, file_struc, variable, lev_arg):
    """
    :param path:
    :param file_struc:
    :param variable:
    :return:
    """
    print('============================================')
    print('\tReading CRS Files...\t\t\t')
    print('============================================')

    len_tot = []
    sza_all = np.empty([])
    lat_all = np.empty([])
    lon_all = np.empty([])
    var_all = np.empty([])

    for d in range(1, 7, 1):
        if d < 10:
            d = '0' + str(d)

        for k in range(24):
            if k < 10:
                k = '0' + str(k)

            file = file_struc + str(d) + str(k)

            file_path = path + file
            print(file_path)

            lat, lon, pres_levs, obs_tim, sfc_ind, sza = ceres.read_crs_geolocation_dev(file_path)

            var, var_name, var_units, var_lev = \
            ceres.read_crs_var_dev(file_path=file_path,
                                   var_name=variable,
                                   lev_arg=lev_arg, fill=0)

            len_tot.append(lat.shape[0])
            sza_all = np.concatenate((sza_all, sza), axis=None)
            lat_all = np.concatenate((lat_all, lat), axis=None)
            lon_all = np.concatenate((lon_all, lon), axis=None)
            var_all = np.concatenate((var_all, var), axis=None)

    print(len_tot)
    print(var_all.shape)

    return var_all, lon_all, lat_all, sza_all


# =====================================================================


def read_full_month_of_ssf_files(path, file_struc, variable):
    print('============================================')
    print('\tReading SSF Files...\t\t\t')
    print('============================================')

    len_tot = []
    sza_all = np.empty([])
    tim_all = np.empty([])
    lat_all = np.empty([])
    lon_all = np.empty([])
    var_all = np.empty([])

    for d in range(1,7,1):
        if d < 10:
            d = '0' + str(d)

        for k in range(24):
            if k < 10:
                k = '0' + str(k)

            file = file_struc + str(d) + str(k)

            file_path = path + file
            print(file_path)

            # read geolocation info from file
            lat, lon, tim, sza = ceres.read_ssf_geolocation(file_path)
            # read variable from file
            var, var_name, var_units = \
            ceres.read_ssf_var(file_path=file_path, var_name=variable, fill=True)

            # len_tot contains the length (num FOVs) of each individual swath
            len_tot.append(lat.shape[0])

            # all of the swaths combined
            sza_all = np.concatenate((sza_all, sza), axis=None)
            tim_all = np.concatenate((tim_all, tim), axis=None)
            lat_all = np.concatenate((lat_all, lat), axis=None)
            lon_all = np.concatenate((lon_all, lon), axis=None)
            var_all = np.concatenate((var_all, var), axis=None)

    print(len_tot)
    print(var_all.shape)

    return var_all, lon_all, lat_all, sza_all



# =====================================================================



var_all3, lon_all3, lat_all3, sza_all3 = read_full_month_of_crs_files(
                               path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019/',
                               file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.201901',
                               variable='Longwave flux - downward - total sky',
                               lev_arg=5)

var_all4, lon_all4, lat_all4, sza_all4 = \
    read_full_month_of_ssf_files(path='/Users/rcscott2/Desktop/CRS/ASDC_archive/',
                               file_struc='CER_SSF_Terra-FM1-MODIS_Edition4A_404405.201901',
                               variable='CERES downward LW surface flux - Model B')


ceres.swath_histogram_scatterplot(field2=var_all3, field1=var_all4,
                                  var_name='Longwave flux - downward - total sky',
                                  lev_name='surface',
                                  ti_str1='CRS Computed',
                                  ti_str2='SSF Model B',
                                  time_info='January 1-7, 2019',
                                  platform='Terra FM-1',
                                  day_only=False, sza=sza_all3)


# ===============================================================================================


# var_all2, lon_all2, lat_all2, sza_all2 = \
#     read_full_day_of_crs_files(path='/Users/rcscott2/Desktop/CRS/my_output/JAN-01-2019/',
#                                file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.20190101',
#                                variable='Shortwave flux - downward - total sky',
#                                lev_arg=5)
#
#
# # CERES LW TOA flux - upwards
# # CERES downward LW surface flux - Model B
# var_all1, lon_all1, lat_all1, sza_all1 =\
#     read_full_day_of_ssf_files(path='/Users/rcscott2/Desktop/CRS/ASDC_archive/',
#                                file_struc='CER_SSF_Terra-FM1-MODIS_Edition4A_404405.20190101',
#                                variable='CERES downward SW surface flux - Model B')


#
# diff = ceres.swath_difference(field2=var_all2, field1=var_all1, day_only=True, sza=sza_all2)
#
# #var_all2[sza_all2 > 90] = np.nan
#
# lw_colormap = ceres.set_colormap(cmap_name=Balance_20, typ_arg=0)
# ceres.plot_swath(lon=lon_all2, lat=lat_all2, field=diff, nrows=1, ncols=1, cen_lon=0,
#                   varname='Shortwave flux - upward - total sky', levname='TOA', varunits='Watts per square meter',
#                   cmap=lw_colormap, cmap_lims=(-100, 100), date='20190101', date_str='01/01/2019:00-23h',
#                   nightshade=False, title_str='Terra FM-1 CERES CRS Ed4 (computed) minus CERES SSF Ed4A (observed)')


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
# # LW colormap - nice options = Deep_20_r, Thermal_20, Bilbao_20
# lw_colormap = ceres.set_colormap(cmap_name=Nuuk_20, typ_arg=0)
# ceres.plot_swath(lon=lon_all2, lat=lat_all2, field=var_all2, nrows=1, ncols=1, cen_lon=0,
#                  varname='Shortwave flux - downward - total sky', levname='TOA', varunits='Watts per square meter',
#                  cmap=lw_colormap, cmap_lims=(0, 1000), date='20190101', date_str='01/01/2019:00-23h',
#                  nightshade=False, title_str='CERES Terra-FM1 CRS Edition 4')

# # LW colormap - nice options = Deep_20_r, Thermal_20, Bilbao_20
# lw_colormap = ceres.set_colormap(cmap_name=Nuuk_20, typ_arg=0)
# ceres.plot_swath(lon=lon_all1, lat=lat_all1, field=var_all1, nrows=1, ncols=1, cen_lon=0,
#                  varname='Shortwave flux - downward - total sky', levname='TOA', varunits='Watts per square meter',
#                  cmap=lw_colormap, cmap_lims=(0, 1000), date='20190101', date_str='01/01/2019:00-23h',
#                  nightshade=False, title_str='CERES Terra-FM1 SSF Ed4A parameterized')