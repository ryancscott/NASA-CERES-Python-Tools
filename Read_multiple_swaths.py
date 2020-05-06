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

import cerestools as ceres
import numpy as np


print('============================================')
print('\tReading CRS Files...\t\t\t')
print('============================================')

len_tot = []
lat_all = np.empty([])
lon_all = np.empty([])
var_all = np.empty([])

for k in range(24):
    if k < 10:
        k = '0' + str(k)

    path = '/Users/rcscott2/Desktop/CRS/my_output/JAN-01-2019/'
    file = 'CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.20190101' + str(k)

    file_path = path+file
    print(file_path)

    lat, lon, pres_levs, obs_tim, sza = ceres.read_crs_geolocation_dev(file_path)

    # Date information
    date, date_str = ceres.get_date(file=file_path)

    swd_tot, swd_tot_name, swd_tot_units, swd_tot_lev = \
        ceres.read_crs_var_dev(file_path=file_path,
                               var_name='Longwave flux - downward - total sky',
                               lev_arg=5, fill=1)

    print(lat.shape[0])

    len_tot.append(lat.shape[0])

    lat_all = np.concatenate((lat_all, lat), axis=None)
    lon_all = np.concatenate((lon_all, lon), axis=None)
    var_all = np.concatenate((var_all, swd_tot), axis=None)

print(len_tot)

print(lat_all.shape)
print(lon_all.shape)
print(var_all.shape)


# SW colormap - nice options = LaJolla_20_r, Nuuk_20
sw_colormap = ceres.set_colormap(cmap_name=Nuuk_20, typ_arg=0)

# LW colormap - nice options = Deep_20_r, Thermal_20, Bilbao_20
lw_colormap = ceres.set_colormap(cmap_name=Deep_20_r, typ_arg=0)


ceres.plot_swath(lon=lon_all, lat=lat_all, field=var_all, nrows=1, ncols=1, cen_lon=0,
                 varname=swd_tot_name, levname=swd_tot_lev, varunits=swd_tot_units,
                 cmap=lw_colormap, cmap_lims=(150, 450), date=date, date_str='01/01/2019:00-23h',
                 nightshade=False, title_str='CERES Terra CRS Edition 4')