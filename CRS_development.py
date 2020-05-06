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
print('\tCRS File...\t\t\t')
print('============================================')


# Aqua-FM3
# file = 'CER_CRS4_Aqua-FM3-MODIS_GH4_2222TH.2019010100'

# Terra-FM1
# file = 'CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.2010062023'
# file = 'CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.2019010100'
# file = 'CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.2019011200'
# file = 'CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.2019011212'
file = 'CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.2019010123'

path = '/Users/rcscott2/Desktop/CRS/my_output/'

file_path = path + file

print(file_path)

print('============================================')
print('\tReading Data...\t\t\t')
print('============================================')

# Geolocation
lat, lon, pres_levs, obs_tim, sza = ceres.read_crs_geolocation_dev(file_path)

# Date information
date, date_str = ceres.get_date(file=file)

# Platform = satellite + flight model
platform = ceres.get_platform_dev(file=file)

# ----------------------
# Clouds
# ----------------------
cld_frac, cld_frac_name, cld_frac_units, _ =\
    ceres.read_crs_var_dev(file_path=file_path,
                           var_name='Cloud fraction',
                           lev_arg=-1, fill=1)

cld_frac_tot = np.sum(cld_frac, axis=1)

# ----------------------
# Surface fluxes
# ----------------------
swd_tot, swd_tot_name, swd_tot_units, swd_tot_lev =\
    ceres.read_crs_var_dev(file_path=file_path,
                           var_name='Shortwave flux - downward - total sky',
                           lev_arg=0, fill=1)

swd_clr, swd_clr_name, swd_clr_units, swd_clr_lev = \
    ceres.read_crs_var_dev(file_path=file_path,
                           var_name='Shortwave flux - downward - clear sky',
                           lev_arg=5, fill=1)

swd_prs, swd_prs_name, swd_prs_units, swd_prs_lev = \
    ceres.read_crs_var_dev(file_path=file_path,
                           var_name='Shortwave flux - downward - pristine sky',
                           lev_arg=5, fill=1)
swd_tna, swd_tna_name, swd_tna_units, swd_tna_lev = \
    ceres.read_crs_var_dev(file_path=file_path,
                           var_name='Shortwave flux - downward - total sky no aerosol',
                           lev_arg=5, fill=1)

lwd_tot, lwd_tot_name, lwd_tot_units, lwd_tot_lev = \
    ceres.read_crs_var_dev(file_path=file_path,
                           var_name='Longwave flux - downward - total sky',
                           lev_arg=5, fill=1)

lwd_clr, lwd_clr_name, lwd_clr_units, lwd_clr_lev = \
    ceres.read_crs_var_dev(file_path=file_path,
                           var_name='Longwave flux - downward - clear sky',
                           lev_arg=5, fill=1)

lwd_prs, lwd_prs_name, lwd_prs_units, lwd_prs_lev = \
    ceres.read_crs_var_dev(file_path=file_path,
                           var_name='Longwave flux - downward - pristine sky',
                           lev_arg=5, fill=1)

lwd_tna, lwd_tna_name, lwd_tna_units, lwd_tna_lev = \
    ceres.read_crs_var_dev(file_path=file_path,
                           var_name='Longwave flux - downward - total sky no aerosol',
                           lev_arg=5, fill=1)

# Single level OR column integrated parameters
# aot, aot_name, aot_units, aot_lev = ceres.read_crs_var_dev(file_path=file_path, vararg=16, levarg=-1, fill=1)
# hgts, hgts_name, hgts_units, hgts_lev = ceres.read_crs_var_dev(file_path=file_path, vararg=17, levarg=-1, fill=1)
# sza, sza_name, sza_units, sza_lev = ceres.read_crs_var_dev(file_path=file_path, vararg=18, levarg=-1, fill=1)

sw_dif_tot, sw_dif_tot_name, sw_dif_tot_units, sw_dif_tot_lev =\
    ceres.read_crs_var_dev(file_path=file_path,
                           var_name='Shortwave diffuse surface flux - total sky',
                           lev_arg=-1, fill=1)

sw_dir_tot, sw_dir_tot_name, sw_dir_tot_units, sw_dir_tot_lev =\
    ceres.read_crs_var_dev(file_path=file_path,
                           var_name='Shortwave direct surface flux - total sky',
                           lev_arg=-1, fill=1)

spec_sw_tot, spec_sw_tot_name, spec_sw_tot_units, spec_sw_tot_lev, =\
    ceres.read_crs_var_dev(file_path=file_path,
                           var_name='Spectral SW surface flux - downward - total sky',
                           lev_arg=-1, fill=1)

print(spec_sw_tot.shape)

spec_lw_tot, spec_lw_tot_name, spec_lw_tot_units, spec_lw_tot_lev, =\
    ceres.read_crs_var_dev(file_path=file_path,
                           var_name='Spectral LW surface flux - downward - total sky',
                           lev_arg=-1, fill=1)

print(spec_lw_tot.shape)


print('============================================')
print('\tPlotting Data...\t\t\t')
print('============================================')

# Cloud fraction colormap
cf_colormap = ceres.set_colormap(cmap_name=Devon_20, typ_arg=0)

# SW colormap - nice options = LaJolla_20_r, Nuuk_20
sw_colormap = ceres.set_colormap(cmap_name=Nuuk_20, typ_arg=0)
# LW colormap - nice options = Deep_20_r, Thermal_20, Bilbao_20
lw_colormap = ceres.set_colormap(cmap_name=Bilbao_20, typ_arg=0)
# Flux difference colormap
fd_colormap = ceres.set_colormap(cmap_name=Balance_20, typ_arg=0)

title = platform + ' - Cloud Radiative Swath (CRS) Ed4 Development'

# Clouds

ceres.plot_swath(lon=lon, lat=lat, field=cld_frac_tot, nrows=1, ncols=1, cen_lon=0,
                 varname=cld_frac_name, levname='', varunits=cld_frac_units,
                 cmap=cf_colormap, cmap_lims=(0, 1), date=date, date_str=date_str,
                 nightshade=True, title_str=title)


# Shortwave plots

ceres.plot_swath(lon=lon, lat=lat, field=swd_tot, nrows=1, ncols=1, cen_lon=0,
                 varname=swd_tot_name, levname=swd_tot_lev, varunits=swd_tot_units,
                 cmap=lw_colormap, cmap_lims=(0, 1000), date=date, date_str=date_str,
                 nightshade=True, title_str=title)

ceres.plot_swath(lon=lon, lat=lat, field=swd_clr, nrows=1, ncols=1, cen_lon=0,
                 varname=swd_clr_name, levname=swd_clr_lev, varunits=swd_clr_units,
                 cmap=sw_colormap, cmap_lims=(0, 1000), date=date, date_str=date_str,
                 nightshade=1, title_str=title)

ceres.plot_swath(lon=lon, lat=lat, field=swd_prs, nrows=1, ncols=1, cen_lon=0,
                 varname=swd_prs_name, levname=swd_prs_lev, varunits=swd_prs_units,
                 cmap=sw_colormap, cmap_lims=(0, 1000), date=date, date_str=date_str,
                 nightshade=1, title_str=title)

ceres.plot_swath(lon=lon, lat=lat, field=swd_tna, nrows=1, ncols=1, cen_lon=0,
                 varname=swd_tna_name, levname=swd_tna_lev, varunits=swd_tna_units,
                 cmap=sw_colormap, cmap_lims=(0, 1000), date=date, date_str=date_str,
                 nightshade=1, title_str=title)

ceres.plot_swath(lon=lon, lat=lat, field=sw_dif_tot, nrows=1, ncols=1, cen_lon=0,
                 varname=sw_dif_tot_name, levname=sw_dif_tot_lev, varunits=sw_dif_tot_units,
                 cmap=sw_colormap, cmap_lims=(0, 1000), date=date, date_str=date_str,
                 nightshade=1, title_str=title)

ceres.plot_swath(lon=lon, lat=lat, field=sw_dir_tot, nrows=1, ncols=1, cen_lon=0,
                 varname=sw_dir_tot_name, levname=sw_dir_tot_lev, varunits=sw_dir_tot_units,
                 cmap=sw_colormap, cmap_lims=(0, 1000), date=date, date_str=date_str,
                 nightshade=1, title_str=title)

# band = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
# cmap_max = [5, 30, 80, 100, 150, 150, 150, 150, 100, 100, 60, 40, 30, 30]
#
# for k in range(14):
#     ceres.plot_swath(lon=lon, lat=lat, field=spec_sw_tot[:, k], nrows=1, ncols=1, cen_lon=0,
#                  varname='Narrowband SW flux' + ' - band ' + str(band[k+1]),
#                  levname='surface',
#                  varunits=spec_sw_tot_units,
#                  cmap=sw_colormap, cmap_lims=(0, cmap_max[k]),
#                  date=date, date_str=date_str,
#                  nightshade=1, title_str=title)

# Longwave plots

ceres.plot_swath(lon=lon, lat=lat, field=lwd_clr, nrows=1, ncols=1, cen_lon=0,
                 varname=lwd_clr_name, levname=lwd_clr_lev, varunits=lwd_clr_units,
                 cmap=lw_colormap, cmap_lims=(100, 400), date=date, date_str=date_str,
                 nightshade=1, title_str=title)

ceres.plot_swath(lon=lon, lat=lat, field=lwd_prs, nrows=1, ncols=1, cen_lon=0,
                 varname=lwd_prs_name, levname=lwd_prs_lev, varunits=lwd_prs_units,
                 cmap=lw_colormap, cmap_lims=(100, 400), date=date, date_str=date_str,
                 nightshade=1, title_str=title)

ceres.plot_swath(lon=lon, lat=lat, field=lwd_tna, nrows=1, ncols=1, cen_lon=0,
                 varname=lwd_tna_name, levname=lwd_tna_lev, varunits=lwd_tna_units,
                 cmap=lw_colormap, cmap_lims=(100, 400), date=date, date_str=date_str,
                 nightshade=1, title_str=title)

ceres.plot_swath(lon=lon, lat=lat, field=lwd_tot, nrows=1, ncols=1, cen_lon=0,
                 varname=lwd_tot_name, levname=lwd_tot_lev, varunits=lwd_tot_units,
                 cmap=lw_colormap, cmap_lims=(100, 400), date=date, date_str=date_str,
                 nightshade=1, title_str=title)

# Flux difference plots

# swd_clr_minus_prs = ceres.swath_difference(swd_clr, swd_prs, day_only=True, sza=sza)
#
# ceres.plot_swath(lon=lon, lat=lat, field=swd_clr_minus_prs, nrows=1, ncols=1, cen_lon=0,
#                  varname='Shortwave flux down  -  clear minus pristine',
#                  levname='at the surface', varunits='Watts per square meter',
#                  cmap=fd_colormap, cmap_lims=(-25, 25), date=date, date_str=date_str,
#                  nightshade=1, title_str=title)

# lwd_clr_minus_prs = ceres.swath_difference(lwd_clr, lwd_prs, day_only=False, sza=sza)
#
# ceres.plot_swath(lon=lon, lat=lat, field=lwd_clr_minus_prs, nrows=1, ncols=1, cen_lon=0,
#                  varname='Longwave flux down  -  clear minus pristine',
#                  levname='at the surface', varunits='Watts per square meter',
#                  cmap=fd_colormap, cmap_lims=(-25, 25), date=date, date_str=date_str,
#                  nightshade=1, title_str=title)



def plot_swath_grid(lon, lat, field,
               varname, levname, varunits,
               nrows, ncols, cen_lon,
               date, nightshade,
               date_str, title_str):
    """
    ----------------------------------------------------------------------------
    This function plots a swath of footprint-level data
    FLASHFlux, SSF, or CRS
    ----------------------------------------------------------------------------
    :param lon: FOV longitude
    :param lat: FOV latitude
    :param field: variable
    :param varname: variable name
    :param levname: level name
    :param varunits: variable units
    :param nrows: number of rows
    :param ncols: number of columns
    :param cen_lon: central longitude
    :param cmap: colormap
    :param cmap_lims: colormap limits
    :param date: date info for nightshade
    :param nightshade: whether to use nighshade feature
    :param date_str: date string
    :param title_str: title string
    ----------------------------------------------------------------------------
    :return: plot of the data
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature
    from mpl_toolkits.axes_grid1 import AxesGrid
    from cartopy.mpl.geoaxes import GeoAxes
    from cartopy.feature.nightshade import Nightshade

    # Map projection
    projection = ccrs.PlateCarree(central_longitude=cen_lon)

    # Axis class
    axes_class = (GeoAxes, dict(map_projection=projection))

    # Create figure
    fig = plt.figure(figsize=(10, 8))
    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(nrows, ncols),
                    axes_pad=(0.4, 0.4),
                    share_all=True,
                    cbar_location='right',
                    cbar_mode='each',
                    cbar_pad=0.1,
                    cbar_size='5%',
                    label_mode=1)

    # Loop over axes
    for i, ax in enumerate(axgr):
        ax.add_feature(cartopy.feature.LAND, zorder=1, facecolor='none',
                       edgecolor='darkgrey')
        ax.gridlines(color='grey', linestyle='--')
        ax.set_title(r'LW$\downarrow$ surface flux in band ' + str(i + 1) + ' - total sky', fontsize=8)
        ax.set_extent([-180, 180, -90, 90], projection)
        ax.text(0.5, -0.1, varname + ' - ' + levname + '\n' + varunits,
                va='bottom', ha='center',
                rotation='horizontal', rotation_mode='anchor',
                transform=ax.transAxes, fontsize=10)

        if nightshade == 1:
            ax.add_feature(Nightshade(date, alpha=0.1))

    # To use a different colorbar range each time, use a tuple of tuples
    # limits = ((0, 120), (0, 120), (0, 120), (0, 120), (0, 120), (0, 120), (0, 120),
    #           (0, 120), (0, 120), (0, 120), (0, 120), (0, 120), (0, 120), (0, 120))
    limits = (0, 100)
    for i in range((nrows * ncols)):
        im = axgr[i].scatter(lon, lat, c=field[:, i], s=1,
                             transform=ccrs.PlateCarree(),
                             vmin=limits[0], vmax=limits[1])
        axgr.cbar_axes[i].colorbar(im)

    for i, cax in enumerate(axgr.cbar_axes):
        if i == 15:
            break
        cax.set_yticks(np.linspace(0, 100, 5))
        cax.set_yticklabels(np.linspace(0, 100, 5), fontsize=6)

    plt.tight_layout()
    plt.show()
    return


plot_swath_grid(lon=lon, lat=lat, field=spec_lw_tot,
                varname='', levname='', varunits='',
                nrows=4, ncols=3, cen_lon=0,
                date=date, nightshade=False,
                date_str='', title_str='')
