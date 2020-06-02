# ========================================================================
# Developing code to manipulate/analyze gridded CERES data
# ========================================================================

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


print('====================================')
print('\t\t\tEBAF file...\t\t\t')
print('====================================')


path = '/Users/rcscott2/Desktop/CERES/EBAF/'
file = 'CERES_EBAF-TOA_Ed4.1_Subset_200003-201910.nc'
file_path = path + file


# ========================================================================

# -------------------------------------------------------
# print information about the nc file
# -------------------------------------------------------
ceres.print_nc_file_info(file_path)

# -------------------------------------------------------
# set the colormap
# -------------------------------------------------------
cmap = ceres.set_colormap(Balance_20, 0)

# -------------------------------------------------------
# read latitude and longitude information from file
# -------------------------------------------------------

latitude, longitude, lat, lon, time = ceres.read_ebaf_geolocation(file_path=file_path)

start_mo, start_yr = ceres.read_ebaf_start_month_year(file_path=file_path)


# -------------------------------------------------------
# read variables including their name and units
# -------------------------------------------------------

lwf, lwf_name, lwf_units = ceres.read_ebaf_var(file_path=file_path, variable='toa_lw_all_mon')
swf, swf_name, swf_units = ceres.read_ebaf_var(file_path=file_path, variable='toa_sw_all_mon')
tcf, tcf_name, tcf_units = ceres.read_ebaf_var(file_path=file_path, variable='cldarea_total_daynight_mon')
cep, cep_name, cep_units = ceres.read_ebaf_var(file_path=file_path, variable='cldpress_total_daynight_mon')
cet, cet_name, cet_units = ceres.read_ebaf_var(file_path=file_path, variable='cldtemp_total_daynight_mon')
tau, tau_name, tau_units = ceres.read_ebaf_var(file_path=file_path, variable='cldtau_total_day_mon')


# -------------------------------------------------------
# compute the long-term mean and standard deviation
# -------------------------------------------------------

mean_lwf, sigma_lwf = ceres.compute_annual_climatology(lwf)
mean_tcf, sigma_tcf = ceres.compute_annual_climatology(tcf)
mean_cep, sigma_cep = ceres.compute_annual_climatology(cep)
mean_cet, sigma_cet = ceres.compute_annual_climatology(cet)
mean_tau, sigma_tau = ceres.compute_annual_climatology(tau)


# -------------------------------------------------------
# plot map of climatology
# -------------------------------------------------------

# ceres.global_map(lon=lon, lat=lat, field=sigma_lwf,
#                  varname=lwf_name, varunits=r'W m$^{-2}$', nrows=1, ncols=1, cen_lon=180,
#                  cmap=cmap, cmap_lims=(5, 40), ti_str=r'CERES-EBAF Ed4.1 St. Dev. (1$\sigma$)')

# ceres.global_map(lon=lon, lat=lat, field=mean_tcf,
#                  varname=tcf_name, varunits=tcf_units, nrows=1, ncols=1, cen_lon=180,
#                  cmap=cmap, cmap_lims=(0, 100), ti_str='CERES-EBAF Ed4.1')
# #
# ceres.global_map(lon=lon, lat=lat, field=mean_tau,
#                   varname=tau_name, varunits=tau_units, nrows=1, ncols=1, cen_lon=180,
#                   cmap=cmap, cmap_lims=(0, 20), ti_str='CERES-EBAF Ed4.1')


# -------------------------------------------------------
# compute weights for area averaging
# -------------------------------------------------------

weights = ceres.cos_lat_weight(latitude)


# -------------------------------------------------------
# compute area-weighted averages of the long-term mean
# -------------------------------------------------------

#ceres.compute_regional_averages(mean_field, latitude=latitude, w=w)

# compute anomalies by removing long-term monthly means
lwf_anomalies, lwf_seasonal_cycle = ceres.compute_monthly_anomalies(lwf, lwf_name)
tcf_anomalies, tcf_seasonal_cycle = ceres.compute_monthly_anomalies(tcf, tcf_name)
cep_anomalies, cep_seasonal_cycle = ceres.compute_monthly_anomalies(cep, cep_name)
cet_anomalies, cet_seasonal_cycle = ceres.compute_monthly_anomalies(cet, cet_name)
tau_anomalies, tau_seasonal_cycle = ceres.compute_monthly_anomalies(tau, tau_name)


# -------------------------------------------------------
# regress a field onto another
# -------------------------------------------------------

# lwf_tcf_coefs = ceres.simple_regression(tcf_anomalies, lwf_anomalies)
#
# ceres.global_map(lon=lon, lat=lat, field=lwf_tcf_coefs,
#                 varname='', varunits=r'$W$ $m^{-2}$ $\%^{-1}$', nrows=1, ncols=1, cen_lon=180,
#                 cmap=cmap, cmap_lims=(-2, 2), ti_str='TOA LW flux regressed on total CF')


# -------------------------------------------------------
# regress a field onto many other fields
# -------------------------------------------------------

# coefficients = ceres.multiple_regression(lwf_anomalies, tcf_anomalies, cet_anomalies)
#
# ceres.global_map(lon=lon, lat=lat, field=coefficients[0,:,:],
#                   varname='', varunits=r'W m$^{-2}$ %$^{-1}$', nrows=1, ncols=1, cen_lon=180,
#                   cmap=cmap, cmap_lims=(-2, 2), ti_str=r'$\partial$LW/$\partial$TCF')
#
# ceres.global_map(lon=lon, lat=lat, field=coefficients[1,:,:],
#                   varname='', varunits=r'W m$^{-2}$ K$^{-1}$', nrows=1, ncols=1, cen_lon=180,
#                   cmap=cmap, cmap_lims=(-5, 5), ti_str=r'$\partial$LW/$\partial$CET')

# -------------------------------------------------------
# compute global mean time series
# -------------------------------------------------------

global_mean_lwf = ceres.global_mean_time_series(field=lwf, weight=weights)

global_mean_lwfa = ceres.global_mean_time_series(field=lwf_anomalies, weight=weights)

# -------------------------------------------------------
# plot time series
# -------------------------------------------------------

ceres.plot_time_series(var=global_mean_lwf, name=lwf_name, units=r'W m$^{-2}$', start_mo=3, start_yr=2000)

ceres.plot_time_series(var=global_mean_lwfa, name=lwf_name, units=r'W m$^{-2}$', start_mo=3, start_yr=2000)
