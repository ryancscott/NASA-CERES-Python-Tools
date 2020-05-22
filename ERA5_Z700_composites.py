

import numpy as np
import era5tools as era5
from palettable.cmocean.diverging import Balance_19

cmap = Balance_19.mpl_colormap


path = '/Users/rcscott2/Desktop/SIO/West_Antarctica/Ghiz_etal/'
file1 = 'ERA5_Z700_FEB2013.nc'
file_path1 = path + file1


file2 = 'ERA5_Z700_MonthlyMeans.nc'
file_path2 = path + file2


era5.print_nc_file_info(file_path1)
era5.print_nc_file_info(file_path2)


# field for days under concern
z700, _, _ = era5.read_era5_var(file_path1, 'z')
z700 = np.true_divide(z700, 9.81)


# monthly means 1980-2010
monthly_z700, _, _ = era5.read_era5_var(file_path2, 'z')
monthly_z700 = np.true_divide(monthly_z700, 9.81)


# lat/lon vectors and grids
latitude, longitude, lat, lon, time = era5.read_era5_geolocation(file_path1)


print(monthly_z700.shape)


D_z700 = monthly_z700[2:93:3][:][:]
J_z700 = monthly_z700[3:93:3][:][:]
F_z700 = monthly_z700[4:93:3][:][:]


D_z700_mean = np.nanmean(D_z700, axis=0)
J_z700_mean = np.nanmean(J_z700, axis=0)
F_z700_mean = np.nanmean(F_z700, axis=0)

# remove long-term monthly mean for appropriate month
z700_anom = z700-F_z700_mean

# shape of each monthly mean Z700 field
print(D_z700_mean.shape)
print(J_z700_mean.shape)
print(F_z700_mean.shape)


z700_anom_composite = np.nanmean(z700_anom[:, :, :], axis=0)


ti_str = 'ERA5 Z$_{700}$ composite anomaly, FEB 2013'


era5.polar_map(lon=lon, lat=lat, field=z700_anom_composite,
                varname='', varunits='m', nrows=1, ncols=1, cen_lon=0,
                cmap=cmap, cmap_lims=(-200, 200), ti_str=ti_str)

