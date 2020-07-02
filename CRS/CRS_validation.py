# ==============================================================================
# This script is used for the validation of the CRS computed surface fluxes.
#
# Author: Ryan C. Scott
#         ryan.c.scott@nasa.gov
# ==============================================================================

import numpy as np
import jul2greg
import cerestools as ceres

# path to text file containing information about surface flux validation sites
sites_path = '/Users/rcscott2/Desktop/FLASHFlux/FOR_Ryan_SSF/site.txt'

# satellite and date information for output files
satellite = 'Terra-FM1'
date = 'JAN-2019'

# open the validation site text file, loop over and read each line
list_of_lists = []
with open(sites_path, mode='r', encoding='latin1') as file:
    for line in file:
        stripped_line = line.strip()
        line_list = stripped_line.split()
        list_of_lists.append(line_list)

file.close()

# for el in list_of_lists:
#     print(el)

# extract lat, lon, etc. to their own list
val_site_lats = []
val_site_lons = []
val_site_name = []
val_site_desc = []
for el in list_of_lists:
    val_site_lats.append(el[0])
    val_site_lons.append(el[1])
    val_site_name.append(el[2])
    val_site_desc.append(el[3:])

# convert list(str) of lat/lon to np arrays
val_site_lats = np.asarray(val_site_lats, dtype=np.float)
val_site_lons = np.asarray(val_site_lons, dtype=np.float)

# create a list of tuples containing site lat, lon, & name
sites = list(zip(val_site_lats, val_site_lons, val_site_name))
# convert to list of lists
sites = [list(el) for el in sites]

# construct an output file for each site, add to "sites" list of lists
for i, site in enumerate(sites):
    file = site[2] + '_' + satellite + '_' + date + '.txt'
    sites[i].append(file)

# for footprints near the surface site, output the
# year, month, day, hour, minute, second (this info
# is needed to match w/ the surface observations) of
# the FOV by calling jul2greg.daycnv, along with
# relevant fluxes and other parameters

yr, mn, day, hr, mi, sec, ms = jul2greg.daycnv(xjd=2458484.583331018, mode="dtlist")

print(yr, mn, day, hr, mi, sec, ms)

# # read in full month of CRS files
# var, lon, lat, sza, = ceres.read_day_of_crs_files(
#     path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019/',
#     file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.20190100',
#     variable='Julian Time',
#     lev_arg=-1,
#     fill=True)

fov_var, fov_lon, fov_lat, fov_sza = ceres.read_day_of_crs_files(
                               path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
                               file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.20190101',
                               variable='Julian Time',
                               lev_arg=-1,
                               fill=True)


# print out information from FOV extracted over site
header = ['year', 'month', 'day', 'hour', 'min', 'sec', 'lat', 'lon', 'sza', 'etc.']

# function calculate distances between CERES FOVs & validation site
def haversine(cer_lat, cer_lon, site_lat, site_lon):
    import math

    lon1, lon2 = cer_lon, site_lon
    lat1, lat2 = cer_lat, site_lat

    R = 6371000.  # radius of Earth in meters
    phi_1 = math.radians(cer_lat)
    phi_2 = math.radians(site_lat)

    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lon2 - lon1)

    a = math.sin(delta_phi / 2.0) ** 2 + math.cos(phi_1) * math.cos(phi_2) * math.sin(delta_lambda / 2.0) ** 2

    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    meters = R * c        # output distance in meters
    km = meters / 1000.0  # output distance in kilometers

    meters = round(meters, 3)
    km = round(km, 3)

    print(f"Distance: {km} km")

    return km


for site in sites:                            # for each site...
    for i, fov in enumerate(list(zip(fov_lat, fov_lon))):   # loop over CERES FOVs...
        print('Computing distance between FOV', fov, 'and', site[2], 'at', str(site[0]), str(site[1]))
        dist = haversine(fov[0], fov[1], site[0], site[1])
        if dist <= 10:
            break









