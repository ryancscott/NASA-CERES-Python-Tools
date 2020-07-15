

import os
import numpy as np
import cerestools as ceres
import matplotlib.pyplot as plt

satellite = 'Aqua'
flight_model = 'FM3'
yr_mon = 'JAN-2019'

obs_path = '/Users/rcscott2/Desktop/CRS/CRS_validation/' \
           'Surface_radiation_data/JAN-2019/'
crs_path = '/Users/rcscott2/Desktop/CRS/CRS_validation/' \
           'Extracted_FOVs/' + satellite + '/'


def plot_ts_scatplots():
    # scatter plots
    plt.figure(figsize=(12, 7))
    plt.subplot(1, 2, 1)
    plt.scatter(obs_lwd, fov_lwd)
    x = np.linspace(100, 500, 100)
    y = x
    plt.plot(x, y)
    plt.xlim([100, 500])
    plt.ylim([100, 500])
    plt.axis('square')
    plt.title(r'LW$\downarrow$ Flux [W m$^{-2}$]')
    plt.xlabel('Surface Observed')
    plt.ylabel('CRS Computed')
    plt.grid()

    plt.subplot(1, 2, 2)
    plt.scatter(obs_swd, fov_swd)
    x2 = np.linspace(0, 1400, 100)
    y2 = x2
    plt.plot(x2, y2)
    plt.xlim([0, 1400])
    plt.ylim([0, 1400])
    plt.axis('square')
    plt.title(r'SW$\downarrow$ Flux [W m$^{-2}$]')
    plt.xlabel('Surface Observed')
    plt.ylabel('CRS Computed')
    plt.grid()
    plt.show()

    # LW time series
    plt.figure(figsize=(12, 7))
    plt.subplot(2, 1, 1)
    plt.plot(site_lwd)
    plt.plot(inst_ind, fov_lwd, 'o')
    xticklabels = [i+1 for i in range(31)]
    plt.xticks(ticks=range(0, 44640, 1440), labels=xticklabels, fontsize=10)
    plt.grid()
    plt.ylabel(r'LW$\downarrow$ Flux [W m$^{-2}$]')
    plt.title('Surface Validation Site: ' + site_name)
    plt.legend(['obs', 'CRS'])

    # SW time series
    plt.subplot(2, 1, 2)
    plt.plot(site_swd)
    plt.plot(inst_ind, fov_swd, 'o')
    xticklabels = [i + 1 for i in range(31)]
    plt.xticks(ticks=range(0, 44640, 1440), labels=xticklabels, fontsize=10)
    plt.grid()
    plt.xlabel('Day of Month')
    plt.ylabel(r'SW$\downarrow$ Flux [W m$^{-2}$]')
    plt.show()

# ==============================================================================
# For every site with surface radiation measurements during "yr_mon"
# i)   read the binary surface radiation data,
# ii)  read in the .txt file of the CRS FOVs near the site,
# iii) extract surface data coincident in time with the CRS FOVs
# iv)  write to new file the matched data
# ==============================================================================


# count the # of files available for yr_mon
count = 0

# loop over files in the directory
for obs_file in os.listdir(obs_path):
    if obs_file.endswith("01"):
        print(os.path.join(obs_path, obs_file))
        count += 1

        # 3-letter site identifier
        site_name = obs_file[4:7]

        # read surface radiation data
        site_csza, site_lwu, site_lwd, \
            site_swdif, site_swu, site_swdir, site_swd = \
            ceres.read1min_binary_cave_rad_obs(file_path=
                                               os.path.join(obs_path, obs_file))

        # crs file
        crs_file = site_name+'_'+satellite+'-'+flight_model+'_'+yr_mon+'.txt'

        # path to CRS FOVs extracted over surface sites
        crs_file_path = crs_path + crs_file

        list_of_lists = []
        with open(crs_file_path, mode='r', encoding='latin1') as file:
            for line in file:
                stripped_line = line.strip()
                line_list = stripped_line.split(',')
                list_of_lists.append(line_list)

        file.close()

        # read first line to get site info
        site_info = list_of_lists[0]
        print('Site name, lat, lon, type:', site_info)
        site_name = site_info[0]
        site_lat = site_info[1]
        site_lon = site_info[2]
        site_type = site_info[3]
        print(site_name)
        print(site_lat)
        print(site_lon)
        print(site_type)

        # extract fov_* to their own list
        fov_year = []
        fov_mon = []
        fov_day = []
        fov_hr = []
        fov_min = []
        fov_sec = []
        fov_dist = []
        fov_lat = []
        fov_lon = []
        fov_sza = []
        fov_swd = []
        fov_lwd = []
        # skip first two lines since they are
        # the i) site info (name, lat, lon) and
        # ii) header info describing each column
        for el in list_of_lists[2:]:
            fov_year.append(el[0])
            fov_mon.append(el[1])
            fov_day.append(el[2])
            fov_hr.append(el[3])
            fov_min.append(el[4])
            fov_sec.append(el[5])
            fov_dist.append(el[6])
            fov_lat.append(el[7])
            fov_lon.append(el[8])
            fov_sza.append(el[9])
            fov_swd.append(el[10])
            fov_lwd.append(el[11])

        # cast as appropriate data types
        fov_year = np.asarray(fov_year, dtype=np.int)
        fov_mon = np.asarray(fov_mon, dtype=np.int)
        fov_day = np.asarray(fov_day, dtype=np.int)
        fov_hr = np.asarray(fov_hr, dtype=np.int)
        fov_min = np.asarray(fov_min, dtype=np.int)
        fov_lat = np.asarray(fov_lat, dtype=np.float)
        fov_lon = np.asarray(fov_lon, dtype=np.float)
        fov_dist = np.asarray(fov_dist, dtype=np.float)
        fov_sza = np.asarray(fov_sza, dtype=np.float)
        fov_swd = np.asarray(fov_swd, dtype=np.float)
        fov_lwd = np.asarray(fov_lwd, dtype=np.float)

        # FOV cos(SZA)
        fov_mu0 = np.cos(fov_sza * (np.pi / 180))

        # next we need to get the index of the
        # surface observations at the FOV time
        day_ind = 1440 * (fov_day - 1)  # 1440 min/d * num days passed - 1
        hr_ind = 60 * (fov_hr - 1)  # 60 min/hr  * num hrs passed - 1

        # calculate instantaneous index to get obs @ FOV time
        inst_ind = day_ind + hr_ind + fov_min

        # subtract one since Python uses 0-based indexing
        inst_ind = inst_ind - 1

        # extract surface LW obs coincident with CERES FOV
        obs_lwd = site_lwd[inst_ind]

        # range of indices centered on the instantaneous match
        # for extracting/averaging the SW data
        # print('Instantaneous - 15:', inst_ind - 15)
        # print('Instantaneous + 15:', inst_ind + 15)

        ind_range = list(zip(inst_ind - 10, inst_ind + 10))
        # print(ind_range)

        sw_ind = []
        obs_swd = np.empty([len(obs_lwd)])
        obs_mu0 = np.empty([len(obs_lwd)])
        for i, el in enumerate(ind_range):
            sw_ind[:] = range(el[0], el[1], 1)
            # print(sw_ind)
            obs_swd[i] = np.nanmean(site_swd[sw_ind])
            obs_mu0[i] = np.nanmean(site_csza[sw_ind])

        # convert CRS SWd using Dave Rutan's normalization formula
        fov_swd = fov_swd * (obs_mu0 / fov_mu0)

        # info about the FOV and surface obs time, flux, SZA
        print('FOV day:\n', fov_day)
        print('FOV hour:\n', fov_hr)
        print('FOV min:\n', fov_min)
        print('OBS instantaneous index:\n', inst_ind)
        print('FOV LWd:\n', fov_lwd)
        print('OBS LWd:\n', obs_lwd)
        print('FOV SWd:\n', fov_swd)
        print('OBS SWd:\n', obs_swd)
        print('OBS cos(SZA):\n', obs_mu0)
        print('FOV cos(SZA):\n', fov_mu0)

        # plot time series and scatter plots for each site
        # plot_ts_scatplots()

        # output FOV and site information to file
        new_file = site_name+'_'+'CRS+OBS_' + satellite + '-' + flight_model + '_' + yr_mon + '.txt'

        # output new file for each site
        header = str(site_name) + ', ' + site_lat + ', ' + site_lon + ', ' + site_type + '\n' + \
            'year  mon day hour min   dist    fov_lat     fov_lon  ' \
            '   fov_sza     fov_swd     obs_swd     fov_lwd     obs_lwd'

        file = open('/Users/rcscott2/Desktop/CRS/CRS_validation/FOVs_combined_with_OBS/'
                    + satellite + '/' + new_file, 'w')
        file.write(header)
        file.write('\n')

        # data to write to file
        data = [fov_year, fov_mon, fov_day, fov_hr, fov_min, fov_dist,
                fov_lat, fov_lon, fov_sza, fov_swd, obs_swd, fov_lwd, obs_lwd]

        # output data to the file
        for i in range(len(fov_year)):
            for el in data:
                file.write(str(el[i]))
                file.write(',  ')
            file.write('\n')


print('There are {} sites with data during {}'.format(count, yr_mon))