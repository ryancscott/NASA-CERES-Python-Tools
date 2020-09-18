# ==============================================================================
# This script is used for validating CRS computed surface fluxes against surface
# radiation measurements from ARM, BSRN, SURFRAD, Buoys, etc.
#
# This script isolates all CERES FOVs located within 10 km of the surface
# validation sites listed in sites.txt, and output the data to text files.
# The output files are used by CRS_validation2.py to extract the surface
# observations at the time of the CERES FOV.
#
# Author: Ryan Scott, SSAI
#         ryan.c.scott@nasa.gov
# ==============================================================================

import sys
import numpy as np
import cerestools as ceres
from palettable.cartocolors.qualitative import Bold_6

# path to text file containing information about surface validation sites
sites_path = '/Users/rcscott2/Desktop/CRS/CRS_validation/sites_01-2019.txt'

# satellite and date information for output files
satellite = 'Terra-FM1'
date = 'JAN-2019'
prod = 'CRS'

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
val_site_type = []
val_site_desc = []
for el in list_of_lists:
    val_site_lats.append(el[0])
    val_site_lons.append(el[1])
    val_site_name.append(el[2])
    val_site_type.append(el[3])
    val_site_desc.append(el[4:])

# convert list(str) of lat/lon to np arrays
val_site_lats = np.asarray(val_site_lats, dtype=np.float)
val_site_lons = np.asarray(val_site_lons, dtype=np.float)
val_site_type = np.asarray(val_site_type, dtype=np.int)

# create a list of tuples containing site lat, lon, & name
sites = list(zip(val_site_name, val_site_lats, val_site_lons, val_site_type))
# convert to list of lists
sites = [list(el) for el in sites]

# construct an output file for each site, add to "sites" list of lists
for i, site in enumerate(sites):
    file = site[0] + '_' + prod + '_' + satellite + '_' + date + '.txt'
    sites[i].append(file)

cmap = ceres.set_colormap(Bold_6, typ_arg=1)
ceres.validation_sites(lon=val_site_lons,
                       lat=val_site_lats,
                       type_ids=val_site_type,
                       title_str='CERES CAVE Surface Validation Sites',
                       date_str='',
                       cmap=cmap)

# sys.exit()

# read in full month of CRS data
fov_lon, fov_lat, fov_tim, fov_sza, fov_swd, fov_lwd, fov_cf1, fov_cf2 = \
    ceres.read_crs_files_validation(
        path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
        file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.201901')


# ==============================================================================
#
# For each surface validation site, loop over CERES footprints & compute the
# distance to each site. Extract FOVs within 10 km of the site and write out
# relevant data to a .txt file.
#
# ==============================================================================


for site in sites:

    # header describing data in each file
    header = str(site[0]) + ', ' + str(site[1]) + ', ' + str(site[2]) + ', ' \
             + str(site[3]) + '\n' + \
             'jultim  year  mon day hour min  sec    dist    fov_lat     fov_lon  ' \
             '   fov_sza      fov_swd     fov_lwd    fov_cf1    fov_cf2 '

    # open site output file
    file = open('/Users/rcscott2/Desktop/CRS/CRS_validation/_CRS/Extracted_FOVs/'+site[4], 'w')
    file.write(header)
    file.write('\n')

    for i, fov in enumerate(list(zip(fov_lat, fov_lon))):

        # code runs faster if print statement is suppressed
        # print('Distance between footprint', i, fov, 'and',
        #      site[2], 'at', str(site[0]), str(site[1]))

        dist = ceres.haversine(fov[0], fov[1], site[1], site[2])

        if dist <= 10:

            # FOV julian time conversion to gregorian date
            yr, mn, day, hr, mi, sec, _ = \
                ceres.jul2greg(xjd=fov_tim[i], mode="dtlist")

            # data to write to file
            data = [fov_tim[i], yr, mn, day, hr, mi, sec, dist,
                    fov[0], fov[1], fov_sza[i], fov_swd[i],
                    fov_lwd[i], fov_cf1[i], fov_cf2[i]]

            # output data to the file
            for el in data:
                file.write(str(el))
                file.write(',  ')

            file.write('\n')


# ==============================================================================

