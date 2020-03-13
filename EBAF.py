# ========================================================================
# Developing essential tools to manipulate/analyze CERES gridded data
# ========================================================================


import cerestools as ceres
from palettable.cmocean.diverging import Balance_19


print('====================================')
print('\t\t\tEBAF file...\t\t\t')
print('====================================')


path = '/Users/rcscott2/Desktop/CERES/EBAF/'
file = 'CERES_EBAF-TOA_Ed4.1_Subset_200003-201910.nc'
file_path = path + file


# ========================================================================


# print information about the nc file
ceres.print_nc_file_info(file_path)

# set the colormap
cmap = ceres.set_colormap(Balance_19, 0)

# read latitude and longitude information from file
latitude, longitude, lat, lon = ceres.read_ebaf_geolocation(filepath=file_path)

# read variable including its name and units
# field, name, units = ceres.read_ebaf_var(filepath=file_path, variable='cldarea_total_daynight_mon')
lwf, lwf_name, lwf_units = ceres.read_ebaf_var(filepath=file_path, variable='toa_lw_all_mon')
tcf, tcf_name, tcf_units = ceres.read_ebaf_var(filepath=file_path, variable='cldarea_total_daynight_mon')
cep, cep_name, cep_units = ceres.read_ebaf_var(filepath=file_path, variable='cldpress_total_daynight_mon')
cet, cet_name, cet_units = ceres.read_ebaf_var(filepath=file_path, variable='cldtemp_total_daynight_mon')
tau, tau_name, tau_units = ceres.read_ebaf_var(filepath=file_path, variable='cldtau_total_day_mon')

import numpy as np
tau = np.log(tau)

# compute the long-term mean and standard deviation
mean_lwf, sigma_lwf = ceres.compute_annual_climatology(lwf)
mean_tcf, sigma_tcf = ceres.compute_annual_climatology(tcf)
mean_cep, sigma_cep = ceres.compute_annual_climatology(cep)
mean_cet, sigma_cet = ceres.compute_annual_climatology(cet)
mean_tau, sigma_tau = ceres.compute_annual_climatology(tau)


# ceres.global_map(lon=lon, lat=lat, field=mean_lwf,
#                  varname=lwf_name, varunits=lwf_units, nrows=1, ncols=1, cen_lon=180,
#                  cmap=cmap, cmap_lims=(0, 300), ti_str='CERES-EBAF Ed4.1')
#
# ceres.global_map(lon=lon, lat=lat, field=mean_tcf,
#                  varname=tcf_name, varunits=tcf_units, nrows=1, ncols=1, cen_lon=180,
#                  cmap=cmap, cmap_lims=(0, 100), ti_str='CERES-EBAF Ed4.1')
#
# ceres.global_map(lon=lon, lat=lat, field=mean_tau,
#                  varname=tau_name, varunits=tau_units, nrows=1, ncols=1, cen_lon=180,
#                  cmap=cmap, cmap_lims=(0, 20), ti_str='CERES-EBAF Ed4.1')


# compute weights for area averaging
#w = ceres.cos_lat_weight(latitude)


# compute area-weighted averages of the long-term mean
#ceres.compute_regional_averages(mean_field, latitude=latitude, w=w)


# compute anomalies by removing long-term monthly means
lwf_anomalies, lwf_seasonal_cycle = ceres.compute_monthly_anomalies(lwf, lwf_name)
tcf_anomalies, tcf_seasonal_cycle = ceres.compute_monthly_anomalies(tcf, tcf_name)
cep_anomalies, cep_seasonal_cycle = ceres.compute_monthly_anomalies(cep, cep_name)
cet_anomalies, cet_seasonal_cycle = ceres.compute_monthly_anomalies(cet, cet_name)
tau_anomalies, tau_seasonal_cycle = ceres.compute_monthly_anomalies(tau, tau_name)

#lwf_tcf_coefs = ceres.simple_regression(tcf_anomalies, lwf_anomalies)
#ceres.global_map(lon=lon, lat=lat, field=lwf_tcf_coefs,
#                 varname='', varunits=r'$W$ $m^{-2}$ $\%^{-1}$', nrows=1, ncols=1, cen_lon=180,
#                 cmap=cmap, cmap_lims=(-2, 2), ti_str='LW flux regressed on total CF')


print('Performing multiple linear regression of TOA LW flux on \\'
      'cloud fraction, effective pressure, effective temperature...')

coefficients = ceres.multiple_regression(lwf_anomalies, tcf_anomalies, cep_anomalies, cet_anomalies)

ceres.global_map(lon=lon, lat=lat, field=coefficients[0,:,:],
                  varname='', varunits=r'W m$^{-2}$ %$^{-1}$', nrows=1, ncols=1, cen_lon=180,
                  cmap=cmap, cmap_lims=(-2, 2), ti_str=r'$\partial$LW/$\partial$TCF')

ceres.global_map(lon=lon, lat=lat, field=coefficients[1,:,:],
                  varname='', varunits=r'W m$^{-2}$ hPa$^{-1}$', nrows=1, ncols=1, cen_lon=180,
                  cmap=cmap, cmap_lims=(-2, 2), ti_str=r'$\partial$LW/$\partial$CEP')

ceres.global_map(lon=lon, lat=lat, field=coefficients[2,:,:],
                  varname='', varunits=r'W m$^{-2}$ K$^{-1}$', nrows=1, ncols=1, cen_lon=180,
                  cmap=cmap, cmap_lims=(-2, 2), ti_str=r'$\partial$LW/$\partial$CET')

