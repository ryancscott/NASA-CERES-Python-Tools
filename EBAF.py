

import cerestools as ceres
from palettable.cmocean.diverging import Balance_19

print('====================================')
print('\t\t\tEBAF file...\t\t\t')
print('====================================')


path = '/Users/rcscott2/Desktop/CERES/EBAF/'
file = 'CERES_EBAF-TOA_Ed4.1_Subset_200003-201910.nc'
file_path = path + file


# ========================================================================


# print information about the files
ceres.print_nc_file_info(file_path)

# read latitude and longitude information from file
latitude, longitude, lat, lon = ceres.read_ebaf_geolocation(filepath=file_path)

# compute weights for area averaging
w = ceres.cos_lat_weight(latitude)

# read variable including its name and units
# field, name, units = ceres.read_ebaf_var(filepath=file_path, variable='cldarea_total_daynight_mon')
tcf, tcf_name, tcf_units = ceres.read_ebaf_var(filepath=file_path, variable='cldarea_total_daynight_mon')
toalw, toalw_name, toalw_units = ceres.read_ebaf_var(filepath=file_path, variable='toa_lw_all_mon')

# compute the long-term mean and standard deviation
#mean_field, sigma_field = ceres.compute_climatology(field)
mean_tcf, sigma_tcf = ceres.compute_climatology(tcf)
mean_toalw, sigma_toalw = ceres.compute_climatology(toalw)

# compute area-weighted averages of the long-term mean
#ceres.compute_regional_averages(mean_field, latitude=latitude, w=w)

# compute anomalies by removing long-term monthly means
# monthly_anomalies, seasonal_cycle = ceres.compute_monthly_anomalies(field)
tcf_anomalies, tcf_seasonal_cycle = ceres.compute_monthly_anomalies(tcf, tcf_name)
toalw_anomalies, toalw_seasonal_cycle = ceres.compute_monthly_anomalies(toalw, toalw_name)

toalw_tcf_coefs = ceres.regress_fields(tcf_anomalies, toalw_anomalies)

# set the colormap
cmap = ceres.set_colormap(Balance_19, 0)


# ceres.global_map(lon=lon, lat=lat, field=mean_tcf,
#                  varname=tcf_name, varunits=tcf_units, nrows=1, ncols=1, cen_lon=180,
#                  cmap=cmap, cmap_lims=(0, 100), ti_str='CERES-EBAF Ed4.1')
#
# ceres.global_map(lon=lon, lat=lat, field=mean_toalw,
#                  varname=toalw_name, varunits=toalw_units, nrows=1, ncols=1, cen_lon=180,
#                  cmap=cmap, cmap_lims=(0, 300), ti_str='CERES-EBAF Ed4.1')


ceres.global_map(lon=lon, lat=lat, field=toalw_tcf_coefs,
                 varname='', varunits=r'$W$ $m^{-2}$ $\%^{-1}$', nrows=1, ncols=1, cen_lon=180,
                 cmap=cmap, cmap_lims=(-2, 2), ti_str='TOA LW flux regressed on total CF')



