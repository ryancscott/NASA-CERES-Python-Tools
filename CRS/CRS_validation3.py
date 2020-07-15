
import os
import numpy as np
import numpy.ma as ma
import plotly.graph_objects as go
from plotly.subplots import make_subplots

satellite = 'Aqua'
flight_model = 'FM3'
yr_mon = 'JAN 2019'
surface_type = 5
site_desc = 'Polar Sites'

crs_obs_path = obs_path = '/Users/rcscott2/Desktop/CRS/CRS_validation/' \
                          'FOVs_combined_with_OBS/' + satellite + '/'

# count the number of sites of "surface_type"
count = 0
# list to store all of the data for "surface_type" sites
all_data = []

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
        site_type = int(site_info[3])

        # if the site is of the correct type, append the
        # data contents to the list containing "all_data"
        if site_type == surface_type:
            count += 1

            next(file)
            for line in file:
                stripped_line = line.strip()
                line_list = stripped_line.split(',')
                all_data.append(line_list)

    file.close()

print('Number of type-{} sites: {}'.format(surface_type, count))

print(all_data)

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


# calculate statistics
lw_diff = fov_lwd - obs_lwd
lw_mean_bias = np.nanmean(lw_diff)
lw_rms_diff = np.sqrt(np.nanmean(lw_diff**2))
lw_corr = ma.corrcoef(ma.masked_invalid(fov_lwd), ma.masked_invalid(obs_lwd))[0][1]

sw_diff = fov_swd - obs_swd
sw_mean_bias = np.nanmean(sw_diff)
sw_rms_diff = np.sqrt(np.nanmean(sw_diff**2))
sw_corr = ma.corrcoef(ma.masked_invalid(fov_swd), ma.masked_invalid(obs_swd))[0][1]

# min_lwd = 100
# max_lwd = 500
# # LW validation figure
#
# fig = go.Figure()
#
# fig.add_trace(go.Scatter(
#     x=np.linspace(min_lwd, max_lwd, 100),
#     y=np.linspace(min_lwd, max_lwd, 100),
#     showlegend=False
# ))
#
# fig.add_trace(go.Scatter(
#     x=obs_lwd,
#     y=fov_lwd,
#     mode='markers',
#     showlegend=False,
#     marker=dict(
#         symbol='circle',
#         opacity=0.8,
#         color='white',
#         size=4,
#         line=dict(width=1)
#     )
# ))
#
# bin_size = (max_lwd-min_lwd)/32
#
# fig.add_trace(go.Histogram2d(
#     x=obs_lwd,
#     y=fov_lwd,
#     colorscale='BuPu',
#     zmax=6,
#     xbins=dict(start=min_lwd, end=max_lwd, size=bin_size),
#     ybins=dict(start=min_lwd, end=max_lwd, size=bin_size),
#     zauto=False,
# ))
#
# fig.add_trace(go.Scatter(
#     x=[min_lwd+75, min_lwd+75, min_lwd+75, min_lwd+75],
#     y=[max_lwd-25, max_lwd-45, max_lwd-65, max_lwd-85],
#     mode="text",
#     name=" ",
#     text=['N = ' + str(len(all_data)),
#           'Bias = ' + str(np.around(lw_mean_bias, 2)),
#           'RMSD = ' + str(np.around(lw_rms_diff, 3)),
#           'Corr = ' + str(np.around(lw_corr, 3))],
#     textposition="bottom center"
# ))
#
# fig.update_layout(
#     title=satellite + ' ' + flight_model + " CRS Ed4 Surface Validation - " + yr_mon + ' - ' + site_desc,
#     xaxis_title="Surface Observed Downwelling LW Flux ",
#     yaxis_title="CRS Computed Downwelling LW Flux",
#     xaxis=dict(ticks='', showgrid=False, zeroline=False, nticks=20, range=[min_lwd, max_lwd]),
#     yaxis=dict(ticks='', showgrid=False, zeroline=False, nticks=20, range=[min_lwd, max_lwd]),
#     autosize=False,
#     height=650,
#     width=650,
#     hovermode='closest',
# )

#fig.show()
# fig.write_image('/Users/rcscott2/Desktop/'+satellite+'-'+flight_model+'-JAN-2019.pdf')


min_lwd = 0
max_lwd = 400

fig = make_subplots(rows=1, cols=2)

fig.add_trace(go.Histogram(x=lw_diff,
                           marker_color='#A9CCE3',
                           xbins=dict(
                                start=-125.0,
                                end=125,
                                size=5
                           ),
                           histnorm='probability',
                           showlegend=False), row=1, col=1)
fig.update_xaxes(title_text="Downwelling LW Flux Difference (CRS - OBS) [W m<sup>-2</sup>]",
                 range=[-125, 125],
                 row=1, col=1)
fig.update_yaxes(title_text="Fraction of Comparisons in Bin",
                 row=1, col=1)

fig.add_trace(go.Scatter(
    x=np.linspace(min_lwd, max_lwd, 100),
    y=np.linspace(min_lwd, max_lwd, 100),
    mode='lines',
    marker_color='#A9CCE3',
    showlegend=False
), row=1, col=2)

fig.add_trace(go.Scatter(
    x=obs_lwd,
    y=fov_lwd,
    mode='markers',
    showlegend=False,
    marker=dict(
        symbol='circle',
        opacity=0.8,
        color='white',
        size=4,
        line=dict(width=1)
    )
), row=1, col=2)

bin_size = (max_lwd-min_lwd)/32

fig.add_trace(go.Histogram2d(
    x=obs_lwd,
    y=fov_lwd,
    colorscale='BuPu',
    zmax=15,
    xbins=dict(start=min_lwd, end=max_lwd, size=bin_size),
    ybins=dict(start=min_lwd, end=max_lwd, size=bin_size),
    zauto=False
), row=1, col=2)

fig.add_trace(go.Scatter(
    x=[min_lwd+75, min_lwd+75, min_lwd+75, min_lwd+75],
    y=[max_lwd-25, max_lwd-45, max_lwd-65, max_lwd-85],
    mode="text",
    name=" ",
    text=['N = ' + str(len(all_data)),
          'Bias = ' + str(np.around(lw_mean_bias, 2)),
          'RMSD = ' + str(np.around(lw_rms_diff, 3)),
          'Corr = ' + str(np.around(lw_corr, 3))],
    textposition="bottom center"
), row=1, col=2)

fig.update_xaxes(title_text="Surface Observed Downwelling LW Flux [W m<sup>-2</sup>]",
                 range=[min_lwd, max_lwd],
                 row=1, col=2)

fig.update_yaxes(title_text="CRS Computed Downwelling LW Flux [W m<sup>-2</sup>]",
                 range=[min_lwd, max_lwd],
                 row=1, col=2)

fig.update_layout(title_text='CERES ' + satellite + ' ' + flight_model +
                             " CRS Ed4 Surface Validation - " +
                             yr_mon + ' - ' + site_desc, height=700,
                  plot_bgcolor='white',
                  xaxis=dict(
                      showline=True,
                      showgrid=True,
                      linecolor='rgb(0, 0, 0)'
                  ),
                  yaxis=dict(
                      showline=True,
                      showgrid=True,
                      linecolor='rgb(0, 0, 0)'))

fig.show()

fig.write_image('/Users/rcscott2/Desktop/' + satellite + '-' + flight_model +
                '-' + site_desc[0:5] + '-JAN-2019.pdf',
                width=1400, height=700)