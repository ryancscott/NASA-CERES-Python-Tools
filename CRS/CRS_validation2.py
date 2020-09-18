# ==============================================================================
# This script is used for validating CRS-computed surface fluxes against surface
# radiation measurements from ARM, BSRN, SURFRAD, etc.
#
# This script isolates and extracts surface radiation measurements at various
# surface sites coincident in time with the CERES FOV and outputs the combined
# data to text files.
#
# Author: Ryan Scott, SSAI
#         ryan.c.scott@nasa.gov
# ==============================================================================


import os
import numpy as np
import cerestools as ceres
import matplotlib.pyplot as plt

satellite = 'Aqua'
flight_model = 'FM3'
yr_mon = 'JAN-2019'
prod = 'CRS'

obs_path = '/Users/rcscott2/Desktop/CRS/CRS_validation/' \
           'Surface_radiation_data/JAN-2019/'
crs_path = '/Users/rcscott2/Desktop/CRS/CRS_validation/_' + prod + '/' \
           'Extracted_FOVs/' + yr_mon + '/' + satellite + '/'


def site_time_series_scatterplots(location):
    """
    ----------------------------------------------------------------------------
    This function plots two figures per site showing (i) scatterplots of the
    observed versus computed downward LW and SW fluxes and (ii) time series
    of the data and the corresponding fluxes from CRS.
    ----------------------------------------------------------------------------
    :param location: 3-letter site identifier [str]
    :return: plots
    ----------------------------------------------------------------------------
    """
    import numpy.ma as ma

    # LW daytime only
    fov_lwd_day = fov_lwd[fov_sza < 90]
    obs_lwd_day = obs_lwd[fov_sza < 90]

    # LW nighttime only
    fov_lwd_ngt = fov_lwd[fov_sza > 90]
    obs_lwd_ngt = obs_lwd[fov_sza > 90]

    print('FOV_LWD shape:', fov_lwd.shape)
    print('FOV_LWD day shape:', fov_lwd_day.shape)
    print('FOV_LWD ngt shape:', fov_lwd_ngt.shape)

    # LW stats - daytime
    lw_num_day = sum(~np.isnan(fov_lwd_day))
    lw_diff_day = fov_lwd_day - obs_lwd_day
    lw_mean_bias_day = np.nanmean(lw_diff_day)
    lw_rms_diff_day = np.sqrt(np.nanmean(lw_diff_day**2))
    lw_corr_day = ma.corrcoef(ma.masked_invalid(fov_lwd_day),
                              ma.masked_invalid(obs_lwd_day))[0][1]
    # LW stats - nighttime
    lw_num_ngt = sum(~np.isnan(fov_lwd_ngt))
    lw_diff_ngt = fov_lwd_ngt - obs_lwd_ngt
    lw_mean_bias_ngt = np.nanmean(lw_diff_ngt)
    lw_rms_diff_ngt = np.sqrt(np.nanmean(lw_diff_ngt ** 2))
    lw_corr_ngt = ma.corrcoef(ma.masked_invalid(fov_lwd_ngt),
                              ma.masked_invalid(obs_lwd_ngt))[0][1]

    # SW - daytime
    fov_swd[fov_sza > 90] = np.nan
    obs_swd[fov_sza > 90] = np.nan

    # SW stats - daytime
    sw_num = sum(~np.isnan(fov_swd))
    sw_diff = fov_swd - obs_swd
    sw_mean_bias = np.nanmean(sw_diff)
    sw_rms_diff = np.sqrt(np.nanmean(sw_diff**2))
    sw_corr = ma.corrcoef(ma.masked_invalid(fov_swd),
                          ma.masked_invalid(obs_swd))[0][1]

    # scatter plots
    plt.figure(figsize=(9, 8))
    plt.rcParams['axes.axisbelow'] = True
    plt.suptitle(yr_mon[:3] + ' ' + yr_mon[-4:] + '\n' + prod + ' ' + satellite + ' ' + flight_model +
                 '\n Surface Validation Site: ' + location)

    plt.subplot(2, 2, 1)
    plt.grid(zorder=0)
    plt.scatter(obs_lwd_day, fov_lwd_day, c=obs_lwd_day, cmap='cividis', s=12, zorder=3)
    x = np.linspace(100, 500, 100)
    plt.plot(x, x)
    plt.xlim([100, 500])
    plt.ylim([100, 500])
    plt.axis('square')
    plt.ylabel(prod + r' Computed LW$\downarrow$ Flux [W m$^{-2}$]', fontsize=9)
    plt.xlabel(r'Surface Measured LW$\downarrow$ Flux [W m$^{-2}$]', fontsize=9)
    # show basic descriptive statistics
    lw_text_str_day = "Daytime \n" + "N = " + str(lw_num_day) + "\n" + \
                      r"Bias ($\bar{\Delta}$) = " + str(np.around(lw_mean_bias_day, 2)) + "\n" + \
                      "RMSD = " + str(np.around(lw_rms_diff_day, 2)) + "\n" + \
                      "Corr = " + str(np.around(lw_corr_day, 2))
    props = dict(facecolor='white', alpha=0.85)
    plt.text(100, 500, lw_text_str_day, fontsize=8, verticalalignment='top', bbox=props)

    plt.subplot(2, 2, 2)
    plt.scatter(obs_swd, fov_swd, c=obs_swd, cmap='plasma', s=12, zorder=3)
    x2 = np.linspace(0, 1200, 100)
    plt.plot(x2, x2)
    plt.xlim([0, 1200])
    plt.ylim([0, 1200])
    plt.axis('square')
    plt.xlabel(r'Surface Measured SW$\downarrow$ Flux [W m$^{-2}$]', fontsize=9)
    plt.ylabel(prod+r' Computed SW$\downarrow$ Flux [W m$^{-2}$]', fontsize=9)
    plt.grid(zorder=0)
    # show basic descriptive statistics
    sw_text_str = "Daytime \n" + "N = " + str(sw_num) + "\n" + \
                  r"Bias ($\bar{\Delta}$) = " + str(np.around(sw_mean_bias, 2)) + "\n" + \
                  "RMSD = " + str(np.around(sw_rms_diff, 2)) + "\n" + \
                  "Corr = " + str(np.around(sw_corr, 2))
    props = dict(facecolor='white', alpha=0.85)
    plt.text(10, 1200, sw_text_str, fontsize=8, verticalalignment='top', bbox=props)

    plt.subplot(2, 2, 3)
    plt.scatter(obs_lwd_ngt, fov_lwd_ngt, c=obs_lwd_ngt, cmap='cividis', s=12, zorder=3)
    x3 = np.linspace(100, 500, 10)
    plt.plot(x3, x3)
    plt.xlim([100, 500])
    plt.ylim([100, 500])
    plt.axis('square')
    plt.xlabel(r'Surface Measured LW$\downarrow$ Flux [W m$^{-2}$]', fontsize=9)
    plt.ylabel(prod + r' Computed LW$\downarrow$ Flux [W m$^{-2}$]', fontsize=9)
    plt.grid(zorder=0)
    # show basic descriptive statistics
    lw_text_str_ngt = "Nighttime \n" + "N = " + str(lw_num_ngt) + "\n" + \
                      r"Bias ($\bar{\Delta}$) = " + str(np.around(lw_mean_bias_ngt, 2)) + "\n" + \
                      "RMSD = " + str(np.around(lw_rms_diff_ngt, 2)) + "\n" + \
                      "Corr = " + str(np.around(lw_corr_ngt, 2))
    props = dict(facecolor='white', alpha=0.85)
    plt.text(100, 500, lw_text_str_ngt, fontsize=8, verticalalignment='top', bbox=props)

    plt.savefig(
        '/Users/rcscott2/Desktop/CRS/CRS_validation/_' + prod + '/Validation_figs/Site_figs/' +
        satellite + '-' + flight_model +
        '-' + location + '-JAN-2019-' + '1' + '.pdf',
        width=1400, height=700)

    #plt.show()

    # =================================================

    # LWd time series
    plt.figure(figsize=(12, 7))
    plt.subplot(2, 1, 1)
    plt.plot(site_lwd)
    plt.plot(inst_ind, fov_lwd, 'o')
    xticklabels = [i+1 for i in range(31)]
    plt.xticks(ticks=range(0, 44640, 1440), labels=xticklabels, fontsize=10)
    plt.grid()
    plt.ylabel(r'LW$\downarrow$ Flux [W m$^{-2}$]')
    plt.title(yr_mon[:3] + ' ' + yr_mon[-4:] + '\n' + prod + ' ' + satellite + ' ' + flight_model +
                 '\n Surface Validation Site: ' + location)
    plt.legend(['obs', prod])

    # SWd time series
    plt.subplot(2, 1, 2)
    plt.plot(site_swd)
    plt.plot(inst_ind, fov_swd, 'o')
    xticklabels = [i + 1 for i in range(31)]
    plt.xticks(ticks=range(0, 44640, 1440), labels=xticklabels, fontsize=10)
    plt.grid()
    plt.xlabel('Day of Month')
    plt.ylabel(r'SW$\downarrow$ Flux [W m$^{-2}$]')

    plt.savefig(
        '/Users/rcscott2/Desktop/CRS/CRS_validation/_' + prod + '/Validation_figs/Site_figs/' +
        satellite + '-' + flight_model +
        '-' + location + '-JAN-2019-' + '2' + '.pdf',
        width=1400, height=700)

    #plt.show()

# ==============================================================================
# For every site with surface radiation measurements during "yr_mon"
# i)   read the binary surface radiation data
# ii)  read the .txt file of the CRS FOVs near the site
# iii) extract surface obs coincident in time with the CRS FOVs
# iv)  write the matched data to a new file
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
            ceres.read1min_binary_rad_obs(file_path=
                                          os.path.join(obs_path, obs_file))

        # crs file
        crs_file = site_name+'_'+prod+'_'+satellite+'-'+flight_model+'_'+yr_mon+'.txt'

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
        fov_jtim = []
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
        fov_cf1 = []
        fov_cf2 = []
        # skip first two lines since they are
        # the i) site info (name, lat, lon) and
        # ii) header info describing each column
        for el in list_of_lists[2:]:
            fov_jtim.append(el[0])
            fov_year.append(el[1])
            fov_mon.append(el[2])
            fov_day.append(el[3])
            fov_hr.append(el[4])
            fov_min.append(el[5])
            fov_sec.append(el[6])
            fov_dist.append(el[7])
            fov_lat.append(el[8])
            fov_lon.append(el[9])
            fov_sza.append(el[10])
            fov_swd.append(el[11])
            fov_lwd.append(el[12])
            fov_cf1.append(el[13])
            fov_cf2.append(el[14])

        # cast as appropriate data types
        fov_jtim = np.asarray(fov_jtim, dtype=np.float)
        fov_year = np.asarray(fov_year, dtype=np.int)
        fov_mon = np.asarray(fov_mon, dtype=np.int)
        fov_day = np.asarray(fov_day, dtype=np.int)
        fov_hr = np.asarray(fov_hr, dtype=np.int)
        fov_min = np.asarray(fov_min, dtype=np.int)
        fov_sec = np.asarray(fov_sec, dtype=np.int)
        fov_lat = np.asarray(fov_lat, dtype=np.float)
        fov_lon = np.asarray(fov_lon, dtype=np.float)
        fov_dist = np.asarray(fov_dist, dtype=np.float)
        fov_sza = np.asarray(fov_sza, dtype=np.float)
        fov_swd = np.asarray(fov_swd, dtype=np.float)
        fov_lwd = np.asarray(fov_lwd, dtype=np.float)
        fov_cf1 = np.asarray(fov_cf1, dtype=np.float)
        fov_cf2 = np.asarray(fov_cf2, dtype=np.float)

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
        minutes = 7
        ind_range = list(zip(inst_ind - minutes, inst_ind + minutes))
        # print(ind_range)

        sw_ind = []
        obs_swd = np.empty([len(obs_lwd)])
        obs_mu0 = np.empty([len(obs_lwd)])
        for i, el in enumerate(ind_range):
            sw_ind[:] = range(el[0], el[1], 1)
            print('SW indices used:', sw_ind)
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
        site_time_series_scatterplots(site_name)

        # output FOV and site information to file
        new_file = site_name+'_'+prod+'+OBS_' + satellite + '-' + flight_model + '_' + yr_mon + '.txt'

        # output new file for each site
        header = str(site_name) + ', ' + site_lat + ', ' + site_lon + ', ' + site_type + '\n' + \
            'julian_time   year  mon day hour min  sec  dist    fov_lat     fov_lon  ' \
            '   fov_sza     fov_swd     obs_swd     fov_lwd     obs_lwd' \
            '   fov_cf1    fov_cf2'

        file = open('/Users/rcscott2/Desktop/CRS/CRS_validation/_'+prod+'/FOVs_combined_with_OBS/'
                    + yr_mon + '/' + satellite + '/' + new_file, 'w')
        file.write(header)
        file.write('\n')

        # data to write to file
        data = [fov_jtim, fov_year, fov_mon, fov_day, fov_hr, fov_min, fov_sec, fov_dist,
                fov_lat, fov_lon, fov_sza, fov_swd, obs_swd, fov_lwd,
                obs_lwd, fov_cf1, fov_cf2]

        # output data to the file
        for i in range(len(fov_year)):
            for el in data:
                file.write(str(el[i]))
                file.write(',  ')
            file.write('\n')


print('There are {} sites with data during {}'.format(count, yr_mon))
