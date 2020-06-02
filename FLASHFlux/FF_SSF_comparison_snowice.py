
import cerestools as ceres
from palettable.cmocean.diverging import Balance_20
from palettable.scientific.sequential import Acton_20
import numpy as np

path1 = '/Users/rcscott2/Desktop/CRS/ASDC_archive/FLASHFlux/'
path2 = '/Users/rcscott2/Desktop/CRS/ASDC_archive/SSF_Ed4A/'

file1 = 'FLASH_SSF_Terra-FM1-MODIS_Version4A_400400.2019010112'
file2 = 'CER_SSF_Terra-FM1-MODIS_Edition4A_404405.2019010112'

file_path1 = path1 + file1
file_path2 = path2 + file2

var_all1, lon_all1, lat_all1, sza_all1 = ceres.read_day_of_ssf_files(
                                        path=path1,
                                        file_struc='FLASH_SSF_Terra-FM1-MODIS_Version4A_400400.20190101',
                                        variable='Snow/Ice percent coverage clear-sky overhead-sun vis albedo',#'Surface type percent coverage',
                                        index=-1,
                                        fill=False)

var_all2, lon_all2, lat_all2, sza_all2 = ceres.read_day_of_ssf_files(
                                        path=path2,
                                        file_struc='CER_SSF_Terra-FM1-MODIS_Edition4A_404405.20190101',
                                        variable='Snow/Ice percent coverage clear-sky overhead-sun vis albedo',
                                        index=-1,
                                        fill=False)


diff = var_all1 - var_all2

# Day time only
# (SZA 86.5 deg cutoff used since ADMs can't be
# applied beyond this SZA)
for i in range(len(sza_all1)):
    if sza_all1[i] >= 86.5:
        sza_all1[i] = np.nan
        diff[i] = np.nan
        lat_all1[i] = np.nan
        lon_all1[i] = np.nan
        lat_all2[i] = np.nan
        lon_all2[i] = np.nan
        var_all1[i] = np.nan
        var_all2[i] = np.nan

# ignore/remove NaNs
bad_indices = np.isnan(diff)
good_indices = ~bad_indices
lat_all1 = lat_all1[good_indices]
lon_all1 = lon_all1[good_indices]
lat_all2 = lat_all2[good_indices]
lon_all2 = lon_all2[good_indices]
var_all1 = var_all1[good_indices]
var_all2 = var_all2[good_indices]
sza_all1 = sza_all1[good_indices]
diff = diff[good_indices]

diff_colormap = ceres.set_colormap(cmap_name=Acton_20, typ_arg=1)

ceres.plot_swath(lon=lon_all2, lat=lat_all2, field=var_all2,
                 nrows=1, ncols=1, cen_lon=0,
                 varname='Snow/Ice percent coverage clear-sky overhead-sun vis albedo', levname='',
                 varunits='%',
                 cmap=diff_colormap, cmap_lims=(0, 100), date='20190101', date_str='01/01/2019:00-22h',
                 nightshade=False, title_str='FLASHFlux Terra FM1 V4A')
