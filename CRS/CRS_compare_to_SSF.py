# ==========================================================================
# Author: Ryan C. Scott, ryan.c.scott@nasa.gov
# This script compares CERES footprint-level quantities (CRS, SSF, etc.)
# ==========================================================================

import sys
import numpy as np
import cerestools as ceres

# cmocean
from palettable.cmocean.sequential import Deep_20_r
from palettable.cmocean.diverging import Balance_20
from palettable.cmocean.sequential import Thermal_20
from palettable.cmocean.sequential import Dense_20
# scientific
from palettable.scientific.sequential import LaJolla_20
from palettable.scientific.sequential import LaPaz_20
from palettable.scientific.sequential import Nuuk_20
from palettable.scientific.sequential import Acton_20
from palettable.scientific.sequential import Bilbao_20
from palettable.scientific.sequential import Devon_20
from palettable.scientific.sequential import Davos_3_r
# colorbrewer
from palettable.colorbrewer.sequential import BuPu_9_r

# ==========================================================================
#
# SIMPLE PLOTS OF THE DATA
#
# ==========================================================================
#
# Plots from the TOA to the Surface

# levs = ['TOA', '70 mb', '200 mb', '500 mb', ' 850 mb', 'Surface']
#
# for i, lev in enumerate(levs):
#
#     olr, lon, lat, sza = \
#         ceres.read_day_of_crs_files(path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
#                                     file_struc='CER_CRS4_Aqua-FM3-MODIS_GH4_1111TH.20190101',
#                                     variable='Shortwave flux - upward - total sky',
#                                     lev_arg=i,
#                                     fill=True)
#
#     olr, lon, lat, sza = ceres.swath_daytime_only(olr, lon, lat, sza, 86.5)
#
#     olr_colormap = ceres.set_colormap(cmap_name=Devon_20, typ_arg=0)
#
#     ceres.plot_swath(lon=lon, lat=lat, field=olr, nrows=1, ncols=1, cen_lon=0,
#                      varname='Upward SW Radiation', levname=lev, varunits='Watts per square meter',
#                      cmap=olr_colormap, cmap_lims=(0, 1000),
#                      date='20190101', date_str='01/01/2019:00-23h',
#                      nightshade=False, title_str='CERES Aqua FM3 CRS',
#                      marker='o', msize=0.5)


# ==========================================================================
#
# Daily plots of OLR or LWd
#
# for day in range(6):
#
#     olr, lon, lat, sza = \
#         ceres.read_day_of_crs_files(path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
#                                     file_struc='CER_CRS4_Aqua-FM3-MODIS_GH4_1111TH.2019010'+str(day+1),
#                                     variable='Longwave flux - downward - total sky',
#                                     lev_arg=5,
#                                     fill=True)
#
#     olr, lon, lat, _ = ceres.swath_daytime_only(olr, lon, lat, sza, 90)
#
#     olr_colormap = ceres.set_colormap(cmap_name=LaJolla_20, typ_arg=0)
#     ceres.plot_swath(lon=lon, lat=lat, field=olr, nrows=1, ncols=1, cen_lon=0,
#                      varname='Incoming Longwave Radiation', levname='Surface', varunits='Watts per square meter',
#                      cmap=olr_colormap, cmap_lims=(100, 450),
#                      date='2019010'+str(day+1), date_str='01/0'+str(day+1)+'/2019:00-23h',
#                      nightshade=False, title_str='CERES Aqua FM3 CRS',
#                      marker='o', msize=0.5)

# ==========================================================================
#
# Daily plots of RSW or SWd

# for day in range(6):
#
#     rsw, lon, lat, sza = \
#         ceres.read_day_of_crs_files(path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
#                                     file_struc='CER_CRS4_Aqua-FM3-MODIS_GH4_1111TH.2019010'+str(day+1),
#                                     variable='Shortwave flux - downward - total sky',
#                                     lev_arg=5,
#                                     fill=True)
#
#     rsw, lon, lat, _ = ceres.swath_daytime_only(rsw, lon, lat, sza, 86.5)
#
#     rsw_colormap = ceres.set_colormap(cmap_name=LaPaz_20, typ_arg=0)
#     ceres.plot_swath(lon=lon, lat=lat, field=rsw, nrows=1, ncols=1, cen_lon=0,
#                      varname='Incoming Shortwave Radiation', levname='Surface', varunits='Watts per square meter',
#                      cmap=rsw_colormap, cmap_lims=(0, 1000),
#                      date='2019010'+str(day+1), date_str='01/0'+str(day+1)+'/2019:00-23h',
#                      nightshade=False, title_str='CERES Aqua FM3 CRS',
#                      marker='o', msize=0.5)

# ==========================================================================
#
# COMPARE CRS to CERES SSF
#
# ==========================================================================

for day in range(1):

    var_all2, lon_all2, lat_all2, sza_all2 = \
        ceres.read_day_of_crs_files(
            path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
            file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.2019010'+
                       str(day+1),
            variable='Longwave flux - upward - total sky',
            lev_arg=0, fill=True,
            dev=True)

    cf_lo, _, _, _ = \
        ceres.read_day_of_crs_files(
            path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
            file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.2019010'+
                       str(day+1),
            variable='Cloud fraction',
            lev_arg=0, fill=True,
            dev=True)

    cf_up, _, _, _ = \
        ceres.read_day_of_crs_files(
            path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
            file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.2019010'+
                       str(day+1),
            variable='Cloud fraction',
            lev_arg=1, fill=True,
            dev=True)

    ct_lo, _, _, _ = \
        ceres.read_day_of_crs_files(
            path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
            file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.2019010'+
                       str(day+1),
            variable='Cloud temperature',
            lev_arg=0, fill=True,
            dev=True)

    #ct_lo[ct_lo > 1.e6] = np.nan

    ct_up, _, _, _ = \
        ceres.read_day_of_crs_files(
            path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
            file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.2019010'+
                       str(day+1),
            variable='Cloud temperature',
            lev_arg=1, fill=True,
            dev=True)

    tot_cf = np.add(cf_lo, cf_up)

# CERES SW TOA flux - upwards
# CERES downward SW surface flux - Model B
# CERES LW TOA flux - upwards
# CERES downward LW surface flux - Model B

    var_all1, lon_all1, lat_all1, sza_all1 = \
        ceres.read_day_of_ssf_files(
            path='/Users/rcscott2/Desktop/CRS/ASDC_archive/SSF_Ed4A/',
            file_struc='CER_SSF_Terra-FM1-MODIS_Edition4A_404405.2019010' +
                        str(day+1),
            variable='CERES LW TOA flux - upwards',
            index=-1, fill=True)

    # computes var_all2 minus var_all1
    diff, fov_wgt = ceres.swath_difference(field2=var_all2,
                                           field1=var_all1,
                                           day=True,
                                           sza=sza_all2,
                                           lat=lat_all2)

    diff_colormap = ceres.set_colormap(cmap_name=Balance_20, typ_arg=0)

    cf_colormap = ceres.set_colormap(cmap_name=Davos_3_r, typ_arg=0)
    # ceres.plot_swath(
    #     lon=lon_all2, lat=lat_all2, field=fov_wgt,
    #     nrows=1, ncols=1, cen_lon=0,
    #     varname='FOV cos(latitude)',
    #     levname='', varunits='',
    #     cmap=diff_colormap, cmap_lims=(0, 1),
    #     date='2019010' + str(day + 1), date_str='01/0' + str(day + 1) + '/2019:00-23h',
    #     nightshade=False,
    #     title_str='',
    #     marker='o', msize=1)

    # ceres.plot_swath(
    #     lon=lon_all2, lat=lat_all2, field=cf_lo,
    #     nrows=1, ncols=1, cen_lon=0,
    #     varname='Cloud layer 1 cloud fraction',
    #     levname='', varunits='',
    #     cmap=cf_colormap, cmap_lims=(0, 1),
    #     date='2019010' + str(day + 1), date_str='01/0' + str(day + 1) + '/2019:00-23h',
    #     nightshade=False,
    #     title_str='',
    #     marker='o', msize=1)

    # ceres.plot_swath(
    #     lon=lon_all2, lat=lat_all2, field=cf_up,
    #     nrows=1, ncols=1, cen_lon=0,
    #     varname='Cloud layer 2 cloud fraction',
    #     levname='', varunits='',
    #     cmap=cf_colormap, cmap_lims=(0, 1),
    #     date='2019010' + str(day + 1), date_str='01/0' + str(day + 1) + '/2019:00-23h',
    #     nightshade=False,
    #     title_str='',
    #     marker='o', msize=1)
    #
    # ceres.plot_swath(
    #     lon=lon_all2, lat=lat_all2, field=tot_cf,
    #     nrows=1, ncols=1, cen_lon=0,
    #     varname='Total cloud fraction',
    #     levname='', varunits='',
    #     cmap=cf_colormap, cmap_lims=(0, 1),
    #     date='2019010' + str(day + 1), date_str='01/0' + str(day + 1) + '/2019:00-23h',
    #     nightshade=False,
    #     title_str='',
    #     marker='o', msize=1)
    #
    # ceres.plot_swath(
    #     lon=lon_all2, lat=lat_all2, field=ct_lo,
    #     nrows=1, ncols=1, cen_lon=0,
    #     varname='Cloud layer 1 temperature',
    #     levname='', varunits='K',
    #     cmap=diff_colormap, cmap_lims=(220, 320),
    #     date='2019010' + str(day + 1), date_str='01/0' + str(day + 1) + '/2019:00-23h',
    #     nightshade=False,
    #     title_str='',
    #     marker='o', msize=1)
    #
    # ceres.plot_swath(
    #     lon=lon_all2, lat=lat_all2, field=ct_up,
    #     nrows=1, ncols=1, cen_lon=0,
    #     varname='Cloud layer 2 temperature',
    #     levname='', varunits='K',
    #     cmap=diff_colormap, cmap_lims=(220, 320),
    #     date='2019010' + str(day + 1), date_str='01/0' + str(day + 1) + '/2019:00-23h',
    #     nightshade=False,
    #     title_str='',
    #     marker='o', msize=1)

    # ceres.plot_swath(
    #     lon=lon_all2, lat=lat_all2, field=mean_cloud_temp,
    #     nrows=1, ncols=1, cen_lon=0,
    #     varname='CF weighted mean lower/upper cloud temperature',
    #     levname='', varunits='K',
    #     cmap=diff_colormap, cmap_lims=(220, 320),
    #     date='2019010' + str(day + 1), date_str='01/0' + str(day + 1) + '/2019:00-23h',
    #     nightshade=False,
    #     title_str='',
    #     marker='o', msize=1)

    # diff[ct_lo > 240.] = np.nan
    #
    # ceres.plot_swath(
    #     lon=lon_all2, lat=lat_all2, field=diff,
    #     nrows=1, ncols=1, cen_lon=0,
    #     varname='Outgoing Longwave Radiation',
    #     levname='TOA', varunits='Watts per square meter',
    #     cmap=diff_colormap, cmap_lims=(-40, 40),
    #     date='2019010'+str(day+1), date_str='01/0'+str(day+1)+'/2019:00-23h',
    #     nightshade=False,
    #     title_str='Daytime \n Terra FM1 CRS (computed) minus SSF Ed4A (parameterized)',
    #     marker='o', msize=1)

    # ceres.swath_histogram_scatterplot(
    #     field2=var_all2,
    #     field1=var_all1,
    #     sza=sza_all2,
    #     lat=lat_all2,
    #     var_name='Longwave flux - upward - total sky',
    #     lev_name='TOA',
    #     f1_str='SSF observed',
    #     f2_str='CRS computed',
    #     main_title=' CRS vs SSF Ed4A',
    #     date_str='1-1-2019:00-23h',
    #     platform='Terra FM1',
    #     day=True)

# ==========================================================================
#
# COMPARE CRS to CRS - i.e., fluxes for different sky conditions
#
# ==========================================================================

# var_all3, lon_all3, lat_all3, sza_all3 = ceres.read_day_of_crs_files(
#                                path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
#                                file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20190103',
#                                variable='Shortwave flux - upward - clear sky',
#                                lev_arg=0,
#                                fill=True)
#
# var_all4, lon_all4, lat_all4, sza_all4 = ceres.read_day_of_crs_files(
#                                path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
#                                file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20190103',
#                                variable='Shortwave flux - upward - pristine sky',
#                                lev_arg=0,
#                                fill=True)
#
# aod_all, _, _, _ = ceres.read_day_of_crs_files(
#                                path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
#                                file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20190103',
#                                variable='MATCH AOT',
#                                lev_arg=-1,
#                                fill=True)
#
# # computes var_all3 minus var_all4
# diff = ceres.swath_difference(field2=var_all3, field1=var_all4, day_only=True, sza=sza_all3)
#
# # set the colormap
# diff_colormap = ceres.set_colormap(cmap_name=Balance_20, typ_arg=0)
# aod_colormap = ceres.set_colormap(cmap_name=Devon_20, typ_arg=0)
#
# # plot the swath
# ceres.plot_swath(lon=lon_all3, lat=lat_all3, field=diff, nrows=1, ncols=1, cen_lon=0,
#                  varname='Shortwave flux - upward', levname='TOA', varunits='Watts per square meter',
#                  cmap=diff_colormap, cmap_lims=(-100, 100), date='20190101', date_str='01/01/2019:00-23h',
#                  nightshade=False, title_str='Terra FM1 CERES CRS Clear-sky minus Pristine-sky',
#                  marker='o', msize=1)
#
# # Daytime only
# aod_all[sza_all3 > 90] = np.nan
#
# ceres.plot_swath(lon=lon_all3, lat=lat_all3, field=aod_all, nrows=1, ncols=1, cen_lon=0,
#                  varname='MATCH AOD', levname='', varunits='',
#                  cmap=aod_colormap, cmap_lims=(0, 1), date='20190103', date_str='01/03/2019:00-23h',
#                  nightshade=False, title_str='Terra-FM1 CERES CRS', marker='o', msize=1)


# ==========================================================================
#
# COMPARE CRS Ed X (new) to CRS Ed 2G
#
# ==========================================================================

#
# var_all3, lon_all3, lat_all3, sza_all3 = ceres.read_day_of_crs_files(
#                                path='/Users/rcscott2/Desktop/CRS_both/',
#                                file_struc='CER_CRS_Terra-FM1-MODIS_Edition2G_023034.20100601',
#                                variable='Longwave flux - upward - total',
#                                lev_arg=0,
#                                fill=True,
#                                dev=False)
#
# var_all4, lon_all4, lat_all4, sza_all4 = ceres.read_day_of_crs_files(
#                                path='/Users/rcscott2/Desktop/CRS_both/',
#                                file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20100601',
#                                variable='Longwave flux - upward - total sky',
#                                lev_arg=0,
#                                fill=True,
#                                dev=True)
#
# # computes var_all3 minus var_all4
# diff = ceres.swath_difference(field2=var_all4, field1=var_all3, day=False, sza=sza_all3)
# # set the colormap
# diff_colormap = ceres.set_colormap(cmap_name=Balance_20, typ_arg=0)
# ceres.plot_swath(lon=lon_all3, lat=lat_all3, field=diff, nrows=1, ncols=1, cen_lon=0,
#                  varname='Longwave flux - upward', levname='TOA', varunits='Watts per square meter',
#                  cmap=diff_colormap, cmap_lims=(-100, 100), date='20100601', date_str='06/01/2010:00-23h',
#                  nightshade=False, title_str='Terra FM1 CERES CRS New minus Ed 2G',
#                  marker='o', msize=1)


# ==========================================================================
#
# READ A FULL MONTH OF DATA
#
# ==========================================================================

# var_all2, lon_all2, lat_all2, sza_all2 = \
#     ceres.read_day_of_crs_files(path='/Users/rcscott2/Desktop/CRS/my_output/JAN-01-2019/',
#                                 file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.20190101',
#                                 variable='Longwave flux - downward - total sky',
#                                 lev_arg=5,
#                                 fill=True)
#
#
# # CERES LW TOA flux - upwards
# # CERES downward LW surface flux - Model B
# var_all1, lon_all1, lat_all1, sza_all1 =\
#     ceres.read_day_of_ssf_files(path='/Users/rcscott2/Desktop/CRS/ASDC_archive/',
#                                 file_struc='CER_SSF_Terra-FM1-MODIS_Edition4A_404405.20190101',
#                                 variable='CERES downward SW surface flux - Model B',
#                                 index=-1, fill=True)
#
#
# diff = ceres.swath_difference(field2=var_all4, field1=var_all3, day_only=True, sza=sza_all2)
# #
# var_all2[sza_all2 > 90] = np.nan
# #
# lw_colormap = ceres.set_colormap(cmap_name=Balance_20, typ_arg=0)
# ceres.plot_swath(lon=lon_all3, lat=lat_all3, field=diff, nrows=1, ncols=1, cen_lon=0,
#                  varname='Longwave flux - downward - total sky', levname='surface', varunits='Watts per square meter',
#                  cmap=lw_colormap, cmap_lims=(-50, 50), date='20190101', date_str='01/01/2019:00-23h',
#                  nightshade=False, title_str='Terra FM1 CERES CRS Ed4 (computed) minus CERES SSF Ed4A (observed)',
#                  marker='o', msize=1)
#
# ceres.swath_histogram_scatterplot(field2=var_all2, field1=var_all1,
#                                   var_name='Shortwave flux - downward',
#                                   lev_name='surface',
#                                   ti_str1='CRS Computed',
#                                   ti_str2='SSF Model B',
#                                   time_info='01/01/2019:00-23h',
#                                   platform='Terra FM-1',
#                                   day_only=True, sza=sza_all2)
# #
# var_all2[sza_all2 > 90] = np.nan
# var_all1[sza_all1 > 90] = np.nan
#
# # LW colormap - nice options = Deep_20_r, Thermal_20, Bilbao_20
# lw_colormap = ceres.set_colormap(cmap_name=Nuuk_20, typ_arg=0)
# ceres.plot_swath(lon=lon_all2, lat=lat_all2, field=var_all2, nrows=1, ncols=1, cen_lon=0,
#                  varname='Shortwave flux - downward - total sky', levname='TOA', varunits='Watts per square meter',
#                  cmap=lw_colormap, cmap_lims=(0, 1000), date='20190101', date_str='01/01/2019:00-23h',
#                  nightshade=False, title_str='CERES Terra-FM1 CRS', marker='o', msize=1)
#
# # LW colormap - nice options = Deep_20_r, Thermal_20, Bilbao_20
# lw_colormap = ceres.set_colormap(cmap_name=Nuuk_20, typ_arg=0)
# ceres.plot_swath(lon=lon_all1, lat=lat_all1, field=var_all1, nrows=1, ncols=1, cen_lon=0,
#                  varname='Shortwave flux - downward - total sky', levname='TOA', varunits='Watts per square meter',
#                  cmap=lw_colormap, cmap_lims=(0, 1000), date='20190101', date_str='01/01/2019:00-23h',
#                  nightshade=False, title_str='CERES Terra-FM1 SSF Ed4A parameterized', marker='o', msize=1)

# ==========================================================================
#
# Histograms and scatterplots based on the full month of data
# #
# var_all2, lon_all2, lat_all2, sza_all2 = \
#     ceres.read_month_of_crs_files(path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
#                                   file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.201901',
#                                   variable='Shortwave flux - downward - total sky',
#                                   lev_arg=5,
#                                   fill=True)
#
# var_all1, lon_all1, lat_all1, sza_all1 = \
#     ceres.read_month_of_ssf_files(path='/Users/rcscott2/Desktop/CRS/ASDC_archive/SSF_Ed4A/',
#                                   file_struc='CER_SSF_Terra-FM1-MODIS_Edition4A_404405.201901',
#                                   variable='CERES downward SW surface flux - Model B',
#                                   index=-1,
#                                   fill=True)
#
# ceres.swath_histogram_scatterplot(
#     field2=var_all2,
#     field1=var_all1,
#     sza=sza_all2,
#     lat=lat_all2,
#     var_name='Shortwave flux - downward - total sky',
#     lev_name='surface',
#     f1_str='SSF parameterized',
#     f2_str='CRS computed',
#     main_title=' CRS vs SSF Ed4A',
#     date_str='JAN 2019',
#     platform='Terra FM1',
#     day=True)