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
field, name, units = ceres.read_ebaf_var(filepath=file_path, variable='cldarea_total_daynight_mon')

# compute the long-term mean and standard deviation
mean_field, sigma_field = ceres.compute_climatology(field)

# compute area-weighted averages of the long-term mean
ceres.compute_regional_averages(mean_field, latitude=latitude, w=w)

# compute anomalies by removing long-term monthly means
monthly_anomalies, seasonal_cycle = ceres.compute_monthly_anomalies(field)

# set the colormap
cmap = ceres.set_colormap(Balance_19, 0)

# plot global map
ceres.global_map(lon=lon, lat=lat, field=mean_field,
                 varname=name, varunits=units, nrows=1, ncols=1, cen_lon=180,
                 cmap=cmap, cmap_lims=(0, 100), ti_str='CERES-EBAF Ed4.1')

