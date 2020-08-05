# ==============================================================================
# This script is for validating footprint-level surface radiative fluxes against
# surface radiation measurements from CERES CAVE (ARM, BSRN, SURFRAD, etc.).
#
# This script reads files produced by running CRS_validation2.py and creates
# scatterplots & histograms of the difference between computed and observed
# fluxes. The results are stratified by day/night and as a function of site
# type.
#
# Site Types: 1 = CST, 2 = DES, 3 = ISL, 4 = CON, 5 = POL, 6 = BUO
#             0 = All, and change line 74 == to !=
#
# Author: Ryan Scott, SSAI
#         ryan.c.scott@nasa.gov
# ==============================================================================


import os
import sys
import numpy as np
import numpy.ma as ma
import cerestools as ceres
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from palettable.cartocolors.qualitative import Bold_6


day = False
prod = 'CRS4'
satellite = 'Terra'
flight_model = 'FM1'
yr_mon = 'JAN 2019'
site_type = 0
site_desc = 'All Sites'

# path to files containing matched FOVs & OBS
crs_obs_path = '/Users/rcscott2/Desktop/CRS/CRS_validation/_'+prod+'/' \
               'FOVs_combined_with_OBS/' + satellite + '/'

# count the number of sites w/ "site_type"
count = 0

# lists to store all of the data for "site_type" sites
all_data = []
site_ids = []
site_lats = []
site_lons = []

# loop over every file in directory containing matched FOVs & OBS
for site_file in os.listdir(crs_obs_path):

    if site_file.endswith(".txt"):
        print(site_file)

    with open(crs_obs_path+site_file, mode='r', encoding='latin1') as file:

        # read the first line to check the site type
        site_info = file.readline().strip().split(',')
        print(site_info)

        site_name = site_info[0]
        site_lat = site_info[1]
        site_lon = site_info[2]
        site_num = int(site_info[3])

        site_ids.append(site_num)
        site_lats.append(site_lat)
        site_lons.append(site_lon)

        # if the site is of the correct type, append the
        # data contents to the list containing "all_data"...
        # to consider all sites together, set the
        # site_type = 0 and use != on the next line...
        if site_num != site_type:
            count += 1

            # skip the second line
            next(file)
            for line in file:
                stripped_line = line.strip()
                line_list = stripped_line.split(',')
                all_data.append(line_list)

    file.close()

print('Number of Type-{} Sites: {}'.format(site_type, count))

site_lats = np.asarray(site_lats, dtype=np.float)
site_lons = np.asarray(site_lons, dtype=np.float)
site_ids = np.asarray(site_ids, dtype=np.int)

print('Site Latitudes:', site_lats)
print('Site Longitudes:', site_lons)
print('Site IDs:', site_ids)

cmap = ceres.set_colormap(Bold_6, typ_arg=1)
ceres.validation_sites(lon=site_lons,
                       lat=site_lats,
                       field=site_ids,
                       title_str='Surface Validation Sites',
                       date_str='January 2019',
                       cmap=cmap)

# sys.exit()


# extract fov_* to their own list
fov_year = []
fov_mon = []
fov_day = []
fov_hr = []
fov_min = []
fov_dist = []
fov_lat = []
fov_lon = []
fov_sza = []
fov_swd = []
obs_swd = []
fov_lwd = []
obs_lwd = []
fov_cf1 = []
fov_cf2 = []
# skip first two lines since they are
# the i) site info (name, lat, lon) and
# ii) header info describing each column
for el in all_data:
    fov_year.append(el[0])
    fov_mon.append(el[1])
    fov_day.append(el[2])
    fov_hr.append(el[3])
    fov_min.append(el[4])
    fov_dist.append(el[5])
    fov_lat.append(el[6])
    fov_lon.append(el[7])
    fov_sza.append(el[8])
    fov_swd.append(el[9])
    obs_swd.append(el[10])
    fov_lwd.append(el[11])
    obs_lwd.append(el[12])
    fov_cf1.append(el[13])
    fov_cf2.append(el[14])

# cast as appropriate data types
fov_year = np.asarray(fov_year, dtype=np.int)
fov_mon = np.asarray(fov_mon, dtype=np.int)
fov_day = np.asarray(fov_day, dtype=np.int)
fov_hr = np.asarray(fov_hr, dtype=np.int)
fov_min = np.asarray(fov_min, dtype=np.int)
fov_dist = np.asarray(fov_dist, dtype=np.float)
fov_lat = np.asarray(fov_lat, dtype=np.float)
fov_lon = np.asarray(fov_lon, dtype=np.float)
fov_sza = np.asarray(fov_sza, dtype=np.float)
fov_swd = np.asarray(fov_swd, dtype=np.float)
obs_swd = np.asarray(obs_swd, dtype=np.float)
fov_lwd = np.asarray(fov_lwd, dtype=np.float)
obs_lwd = np.asarray(obs_lwd, dtype=np.float)
fov_cf1 = np.asarray(fov_cf1, dtype=np.float)
fov_cf2 = np.asarray(fov_cf2, dtype=np.float)

# FOV total cloud fraction
fov_tcf = fov_cf1 + fov_cf2

# ================================

# Clear-sky FOVs only
# fov_lwd[fov_tcf > 0] = np.nan
# fov_swd[fov_tcf > 0] = np.nan

# ================================

# consider daytime (T) or nighttime (F) only
if day is True:
    day_str = 'Daytime'
    fov_lwd[fov_sza > 90] = np.nan
    fov_swd[fov_sza > 90] = np.nan
elif day is False:
    day_str = 'Nighttime'
    fov_lwd[fov_sza < 90] = np.nan
    fov_swd[fov_sza < 90] = np.nan

# number of non-NaN values
num_lw_comparisons = sum(~np.isnan(fov_lwd))
num_sw_comparisons = sum(~np.isnan(fov_swd))
print('Number of LW comparisons:', num_lw_comparisons)
print('Number of SW comparisons:', num_sw_comparisons)


# calculate validation statistics
lw_diff = fov_lwd - obs_lwd
mean_fov_lwd = np.nanmean(fov_lwd)
mean_obs_lwd = np.nanmean(obs_lwd)
lw_mean_bias = np.nanmean(lw_diff)
lw_med_bias = np.nanmedian(lw_diff)
lw_rms_diff = np.sqrt(np.nanmean(lw_diff**2))
lw_corr = ma.corrcoef(ma.masked_invalid(fov_lwd),
                      ma.masked_invalid(obs_lwd))[0][1]

sw_diff = fov_swd - obs_swd
mean_fov_swd = np.nanmean(fov_swd)
mean_obs_swd = np.nanmean(obs_swd)
sw_mean_bias = np.nanmean(sw_diff)
sw_med_bias = np.nanmedian(sw_diff)
sw_rms_diff = np.sqrt(np.nanmean(sw_diff**2))
sw_corr = ma.corrcoef(ma.masked_invalid(fov_swd),
                      ma.masked_invalid(obs_swd))[0][1]

# ==============================================================================
# LW
# ==============================================================================

if site_type == 0:
    min_lwd = 0
    max_lwd = 550
elif site_type == 5:
    min_lwd = 0
    max_lwd = 400
else:
    min_lwd = 100
    max_lwd = 500


fig = make_subplots(rows=1, cols=2)

fig.add_trace(
    go.Histogram(
        x=lw_diff,
        marker_color='#A9CCE3',
        xbins=dict(
            start=-125.0,
            end=125,
            size=5
            ),
        histnorm='probability',
        showlegend=False
    ),
    row=1, col=1)

fig.update_xaxes(
    title_text='Downwelling LW Flux Difference ('+prod+' - OBS) [W m<sup>-2</sup>]',
    range=[-125, 125],
    row=1, col=1)

fig.update_yaxes(
    title_text="Fraction of Comparisons in Bin",
    range=[0, 0.13],
    row=1, col=1)

fig.add_trace(
    go.Scatter(
        x=np.linspace(min_lwd, max_lwd, 100),
        y=np.linspace(min_lwd, max_lwd, 100),
        mode='lines',
        marker_color='#A9CCE3',
        showlegend=False
        ),
    row=1, col=2)

fig.add_trace(
    go.Scatter(
        x=obs_lwd,
        y=fov_lwd,
        mode='markers',
        showlegend=False,
        marker=dict(
            symbol='circle',
            opacity=0.8,
            color='white',
            size=3,
            line=dict(width=1)
            )
        ),
    row=1, col=2)

bin_size = (max_lwd-min_lwd)/32

fig.add_trace(
    go.Histogram2d(
        x=obs_lwd,
        y=fov_lwd,
        colorscale='BuPu',
        zauto=False,
        zmax=15,
        xbins=dict(
            start=min_lwd,
            end=max_lwd,
            size=bin_size
            ),
        ybins=dict(
            start=min_lwd,
            end=max_lwd,
            size=bin_size
            ),
        ),
    row=1, col=2)

fig.add_trace(
    go.Scatter(
        x=[min_lwd+75, min_lwd+75, min_lwd+75, min_lwd+75, min_lwd+75, min_lwd+75],
        y=[max_lwd-25, max_lwd-45, max_lwd-65, max_lwd-85, max_lwd-105, max_lwd-125],
        mode="text",
        name=" ",
        text=['N = ' + str(num_lw_comparisons),
              'Bias = ' + str(np.around(lw_mean_bias, 2)),
              'RMSD = ' + str(np.around(lw_rms_diff, 2)),
              'Corr = ' + str(np.around(lw_corr, 2)),
              prod + ' Mean = ' + str(np.around(mean_fov_lwd, 2)),
              'OBS Mean = ' + str(np.around(mean_obs_lwd, 2))
              ],
        textposition="bottom center"
        ),
    row=1, col=2)

fig.update_xaxes(
    title_text="Surface Observed Downwelling LW Flux [W m<sup>-2</sup>]",
    range=[min_lwd, max_lwd],
    row=1, col=2)

fig.update_yaxes(
    title_text=prod+' Downwelling LW Flux [W m<sup>-2</sup>]',
    range=[min_lwd, max_lwd],
    row=1, col=2)

if prod == 'CRS4' or prod == 'SSF4A':
    sup_title_str = 'CERES ' + satellite + ' ' + flight_model + " " + \
                    prod + ' Surface Validation - ' + yr_mon +\
                    ' - ' + site_desc + ' - ' + day_str
elif prod == 'FF3C' or prod == 'FF4A':
    sup_title_str = 'FLASHFlux ' + 'V' + prod[-2:] + ' ' + satellite + ' ' \
                    + flight_model + ' Surface Validation - ' \
                    + yr_mon + ' - ' + site_desc + ' - ' + day_str

fig.update_layout(
    title_text=sup_title_str,
    height=700,
    plot_bgcolor='white',
    xaxis=dict(
        showline=True,
        showgrid=True,
        linecolor='rgb(0, 0, 0)'
        ),
    yaxis=dict(
        showline=True,
        showgrid=True,
        linecolor='rgb(0, 0, 0)'
        )
    )

fig.show()

fig.write_image('/Users/rcscott2/Desktop/CRS/CRS_validation/_'+prod+'/Validation_figs/'
                + satellite + '-' + flight_model +
                '-' + site_desc[0:3] + '-JAN-2019-' + day_str + '_LW.pdf',
                width=1400, height=700)


# ==============================================================================
# SW
# ==============================================================================

min_swd = 0
max_swd = 1300

fig = make_subplots(rows=1, cols=2)

fig.add_trace(
    go.Histogram(
        x=sw_diff,
        marker_color='#A9CCE3',
        xbins=dict(
            start=-125.0,
            end=125,
            size=5
            ),
        histnorm='probability',
        showlegend=False
    ), row=1, col=1)

fig.update_xaxes(
    title_text='Downwelling SW Flux Difference ('+prod+' - OBS) [W m<sup>-2</sup>]',
    range=[-125, 125],
    row=1, col=1)

fig.update_yaxes(
    title_text='Fraction of Comparisons in Bin',
    range=[0, 0.13],
    row=1, col=1)

fig.add_trace(
    go.Scatter(
        x=np.linspace(min_swd, max_swd, 100),
        y=np.linspace(min_swd, max_swd, 100),
        mode='lines',
        marker_color='#A9CCE3',
        showlegend=False
    ), row=1, col=2)

fig.add_trace(
    go.Scatter(
        x=obs_swd,
        y=fov_swd,
        mode='markers',
        showlegend=False,
        marker=dict(
            symbol='circle',
            opacity=0.8,
            color='white',
            size=3,
            line=dict(width=1)
        )
    ), row=1, col=2)

bin_size = (max_swd-min_swd)/32

fig.add_trace(
    go.Histogram2d(
        x=obs_swd,
        y=fov_swd,
        colorscale='BuPu',
        zmax=15,
        xbins=dict(
            start=min_swd,
            end=max_swd,
            size=bin_size),
        ybins=dict(
            start=min_swd,
            end=max_swd,
            size=bin_size),
        zauto=False
        ),
    row=1, col=2)

fig.add_trace(
    go.Scatter(
        x=[min_swd+250, min_swd+250, min_swd+250, min_swd+250, min_swd+250, min_swd+250],
        y=[max_swd-125, max_swd-165, max_swd-205, max_swd-245, max_swd-285, max_swd-325],
        mode="text",
        name=" ",
        text=['N = ' + str(num_sw_comparisons),
              'Bias = ' + str(np.around(sw_mean_bias, 2)),
              'RMSD = ' + str(np.around(sw_rms_diff, 2)),
              'Corr = ' + str(np.around(sw_corr, 2)),
              prod + ' Mean = ' + str(np.around(mean_fov_swd, 2)),
              'OBS Mean = ' + str(np.around(mean_obs_swd, 2))
              ],
        textposition="bottom center"
        ),
    row=1, col=2)

fig.update_xaxes(
    title_text='Surface Observed Downwelling SW Flux [W m<sup>-2</sup>]',
    range=[min_swd, max_swd],
    row=1, col=2)

fig.update_yaxes(
    title_text=prod+' Downwelling SW Flux [W m<sup>-2</sup>]',
    range=[min_swd, max_swd],
    row=1, col=2)

fig.update_layout(
    title_text=sup_title_str,
    height=700,
    plot_bgcolor='white',
    xaxis=dict(
        showline=True,
        showgrid=True,
        linecolor='rgb(0, 0, 0)'
        ),
    yaxis=dict(
        showline=True,
        showgrid=True,
        linecolor='rgb(0, 0, 0)'
        )
    )

if day is True:
    fig.show()
    fig.write_image(
        '/Users/rcscott2/Desktop/CRS/CRS_validation/_'+prod+'/Validation_figs/' +
        satellite + '-' + flight_model +
        '-' + site_desc[0:3] + '-JAN-2019-' + day_str + '_SW.pdf',
        width=1400, height=700
    )

# ===========================================================================
