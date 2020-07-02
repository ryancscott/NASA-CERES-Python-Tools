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

# read in full month of CRS files
# fov_tim, fov_lon, fov_lat, fov_sza, = \
#     ceres.read_month_of_crs_files(
#         path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019/',
#         file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_2222TH.20190100',
#         variable='Julian Time',
#         lev_arg=-1,
#         fill=True)



def read_month_of_crs_files_validation(path, file_struc):
    """
    ----------------------------------------------------------------------------
    This function loops over and reads an entire month of CRS data
    ----------------------------------------------------------------------------
    :param path: path to files
    :param file_struc: file structure (without day & hour portion at the end)
    :param variable: variable to be read from file
    :param lev_arg: level argument (0 = TOA, 5 = sfc)
    :return: variable, lat, lon, sza
    ----------------------------------------------------------------------------
    """
    import numpy as np

    print('============================================')
    print('\tReading CRS Files...\t\t\t')
    print('============================================')

    len_tot = []
    sza_all = np.empty([])
    lat_all = np.empty([])
    lon_all = np.empty([])
    tim_all = np.empty([])
    swd_all = np.empty([])
    lwd_all = np.empty([])

    for d in range(1, 8, 1):
        if d < 10:
            d = '0' + str(d)

        for k in range(24):
            if k < 10:
                k = '0' + str(k)

            file = file_struc + str(d) + str(k)

            file_path = path + file
            print(file_path)

            lat, lon, pres_levs, obs_tim, sfc_ind, sza \
                = ceres.read_crs_geolocation_dev(file_path)

            swd, _, _, _ = \
                ceres.read_crs_var_dev(
                    file_path=file_path,
                    var_name='Shortwave flux - downward - total sky',
                    lev_arg=5,
                    fill=False)

            lwd, _, _, _ = \
                ceres.read_crs_var_dev(
                    file_path=file_path,
                    var_name='Longwave flux - downward - total sky',
                    lev_arg=5,
                    fill=False)

            len_tot.append(lat.shape[0])
            sza_all = np.concatenate((sza_all, sza), axis=None)
            lat_all = np.concatenate((lat_all, lat), axis=None)
            lon_all = np.concatenate((lon_all, lon), axis=None)
            tim_all = np.concatenate((tim_all, obs_tim), axis=None)
            swd_all = np.concatenate((tim_all, obs_tim), axis=None)
            lwd_all = np.concatenate((tim_all, obs_tim), axis=None)

    print(len_tot)
    print(lat_all.shape)

    return lon_all, lat_all, tim_all, sza_all, swd_all, lwd_all


# read in full day of CRS files
fov_lon, fov_lat, fov_tim, fov_sza, fov_swd, fov_lwd = \
    read_month_of_crs_files_validation(
        path='/Users/rcscott2/Desktop/CRS/my_output/JAN-2019_/',
        file_struc='CER_CRS4_Terra-FM1-MODIS_GH4_1111TH.201901')


def haversine(cer_lat, cer_lon, site_lat, site_lon):
    """
    -------------------------------------------------------------------
    Function calculates the distance between CERES FOVs and
    the surface validation site using the haversine formula
    -------------------------------------------------------------------
    """
    import math

    lon1, lon2 = cer_lon, site_lon
    lat1, lat2 = cer_lat, site_lat

    r = 6371000.  # radius of Earth in meters
    phi_1 = math.radians(cer_lat)
    phi_2 = math.radians(site_lat)

    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lon2 - lon1)

    a = math.sin(delta_phi / 2.0) ** 2 + \
        math.cos(phi_1) * math.cos(phi_2) * math.sin(delta_lambda / 2.0) ** 2

    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    dist_km = (r * c) / 1000.0  # output distance in kilometers
    dist_km = round(dist_km, 3)

    # print(f"Distance: {km} km")

    return dist_km

# ==============================================================================
#
# For each surface validation site, loop over and compute the distance of all
# CERES footprints to the site. Extract FOVs within 10 km of site and write out
# relevant data to a text file.
#
# ==============================================================================


for site in sites:

    # header describing data in each file
    header = str(site[2]) + ': ' + str(site[0]) + ', ' + str(site[1]) + '\n' + \
             'year  mon day hour min  sec    dist    fov_lat     fov_lon '

    # open site output file
    file = open('/Users/rcscott2/Desktop/CRS/Validation_files/'+site[3], 'w')
    file.write(header)
    file.write('\n')

    for i, fov in enumerate(list(zip(fov_lat, fov_lon))):

        # code runs faster if print statement is suppressed
        # print('Distance between footprint', i, fov, 'and',
        #      site[2], 'at', str(site[0]), str(site[1]))

        dist = haversine(fov[0], fov[1], site[0], site[1])

        if dist <= 10:

            # FOV julian time conversion to gregorian date
            yr, mn, day, hr, mi, sec, _ = \
                jul2greg.daycnv(xjd=fov_tim[i], mode="dtlist")

            # data to write to file
            data = [yr, mn, day, hr, mi, sec, dist,
                    fov[0], fov[1], fov_sza[i], fov_swd[i],
                    fov_lwd[i]]

            # output data to the file
            for el in data:
                file.write(str(el))
                file.write(',  ')

            file.write('\n')

















