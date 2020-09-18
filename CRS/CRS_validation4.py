# ==============================================================================
# This script is for validating footprint-level surface radiative fluxes against
# surface radiation measurements from CERES CAVE (ARM, BSRN, SURFRAD, etc.).
#
# This script matches FOVs from multiple products: CRS4, SSF4A, FF3C (FF4A, etc.)
#
# Author: Ryan Scott, SSAI
#         ryan.c.scott@nasa.gov
# ==============================================================================

import os
import sys
import itertools
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

# modify before each run
satellite = 'Terra'
flight_model = 'FM1'
yr_mon = 'JAN-2019'
day = False

# loop over files from each product
product = ['CRS', 'SSF4A', 'FF3C']# , 'FF4A_test']

crs_data = []
ssf_data = []
ff3_data = []
ff4_data = []

# =============================================================================


def read_validation_text_file(path_to_file):

    # list of lists to be populated and returned
    data = []

    with open(path_to_file, mode='r', encoding='latin1') as file:

        # read the first line to check the site type
        site_info = file.readline().strip().split(',')
        print(site_info)

        # count the number of lines in file
        num_lines = sum(1 for line in open(path_to_file, mode='r', encoding='latin1'))
        print('# of FOVs at', site_info[0], ': ', num_lines - 2)

        # skip the second line
        next(file)
        # read rest of file line-by-line
        for line in file:
            stripped_line = line.strip()
            line_list = stripped_line.split(',')
            data.append(line_list)

    file.close()

    return data


# =============================================================================


crs_path = '/Users/rcscott2/Desktop/CRS/CRS_validation/_CRS/' \
           'FOVs_combined_with_OBS/' + yr_mon + '/' + satellite + '/'
ssf_path = '/Users/rcscott2/Desktop/CRS/CRS_validation/_SSF4A/' \
           'FOVs_combined_with_OBS/' + yr_mon + '/' + satellite + '/'
ff3_path = '/Users/rcscott2/Desktop/CRS/CRS_validation/_FF3C/' \
           'FOVs_combined_with_OBS/' + yr_mon + '/' + satellite + '/'
ff4_path = '/Users/rcscott2/Desktop/CRS/CRS_validation/_FF4A_test/' \
           'FOVs_combined_with_OBS/' + yr_mon + '/' + satellite + '/'


# obtain a list of validation site names (3 letter IDs)
site_names = []
site_type = []

for site_file in os.listdir(crs_path):

    with open(crs_path+site_file, mode='r', encoding='latin1') as file:
        # read the first line to get the site name
        site_info = file.readline().strip().split(',')
        site_names.append(site_info[0])
        site_type.append(site_info[3])
    file.close()

print('Available sites: \n', site_names)

# number of common FOVs
count = 0

# collect footprints common to all products
crs_data2 = []
ssf_data2 = []
ff3_data2 = []
ff4_data2 = []

# loop over sites
for site in site_names:

    prod_file = []
    for prod in product:
        prod_file.append(site + '_' + prod + '+OBS_' + satellite + '-' +
                         flight_model + '_' + yr_mon + '.txt')

    crs_file = prod_file[0]
    ssf_file = prod_file[1]
    ff3_file = prod_file[2]
    # ff4_file = prod_file[3]

    print(crs_file)
    print(ssf_file)
    print(ff3_file)
    # print(ff4_file)

    crs_data = read_validation_text_file(crs_path + crs_file)
    ssf_data = read_validation_text_file(ssf_path + ssf_file)
    ff3_data = read_validation_text_file(ff3_path + ff3_file)
    # ff4_data = read_validation_text_file(ff4_path + ff4_file)

    # print(crs_data)
    # print(ssf_data)
    # print(ff3_data)
    # print(ff4_data)

    # isolate FOV distance-to-site info...
    # perhaps replace with Julian time?
    crs_jtim = []
    ssf_jtim = []
    ff3_jtim = []
   # ff4_jtim = []
    for el in crs_data:
        crs_jtim.append(el[0])
    for el in ssf_data:
        ssf_jtim.append(el[0])
    for el in ff3_data:
        ff3_jtim.append(el[0])
    # for el in ff4_data:
    #     ff4_jtim.append(el[0])

    # only use footprints when all 3 products have the
    # same computed dist to the site, i.e., footprints are the same
    intersection = set(crs_jtim).intersection(ssf_jtim, ff3_jtim)#, ff4_dist)
    num_common_fovs = len(intersection)
    print('# of common FOVs at site:', num_common_fovs)

    count += num_common_fovs

    # for el in crs_data:
    #     if el[5] in list(intersection):
    #         crs_data2.append(el)
    # for el in ssf_data:
    #     if el[5] in list(intersection):
    #         ssf_data2.append(el)
    # for el in ff3_data:
    #     if el[5] in list(intersection):
    #         ff3_data2.append(el)

    for crs_fov, ssf_fov, ff3_fov, inter in \
            list(zip(crs_data, ssf_data, ff3_data, list(intersection))):
        if crs_fov[0] and ssf_fov[0] and ff3_fov[0] in list(intersection):
            print(crs_fov, ssf_fov, ff3_fov, inter)
            crs_data2.append(crs_fov)
            ssf_data2.append(ssf_fov)
            ff3_data2.append(ff3_fov)
            # ff4_data2.append(ff4_fov)


print(count, 'common FOVs...')
print('Number of matched FOVs, CRS:', len(crs_data2))
print('Number of matched FOVs, SSF:', len(ssf_data2))
print('Number of matched FOVs, FF3:', len(ff3_data2))
# print('Number of matched FOVs, FF4:', len(ff4_data2))


# print(crs_data2)
# print(ssf_data2)
# print(ff3_data2)
# print(ff4_data2)

# FOV info - common to all (if matched up correctly! verify!)
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
# SW
crs_swd = []
ssf_swd = []
ff3_swd = []
ff4_swd = []
obs_swd = []
# LW
crs_lwd = []
ssf_lwd = []
ff3_lwd = []
ff4_lwd = []
obs_lwd = []
# crs_cf1 = []
# crs_cf2 = []
# ssf_cf1 = []
# ssf_cf2 = []
# ff3_cf1 = []
# ff3_cf2 = []

# skip first two lines since they are
# the i) site info (name, lat, lon) and
# ii) header info describing each column
for el in crs_data2:
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

    crs_swd.append(el[11])
    obs_swd.append(el[12])
    crs_lwd.append(el[13])
    obs_lwd.append(el[14])
    # crs_cf1.append(el[13])
    # crs_cf2.append(el[14])

for el in ssf_data2:
    ssf_swd.append(el[11])
    ssf_lwd.append(el[13])
    # ssf_cf1.append(el[13])
    # ssf_cf2.append(el[14])

for el in ff3_data2:
    ff3_swd.append(el[11])
    ff3_lwd.append(el[13])
    # ff3_cf1.append(el[13])
    # ff3_cf2.append(el[14])

# for el in ff4_data2:
#     ff4_swd.append(el[10])
#     ff4_lwd.append(el[12])
#     # ff3_cf1.append(el[13])
#     # ff3_cf2.append(el[14])

# cast as appropriate data types
fov_jtim = np.asarray(fov_jtim, dtype=np.float)
fov_year = np.asarray(fov_year, dtype=np.int)
fov_mon = np.asarray(fov_mon, dtype=np.int)
fov_day = np.asarray(fov_day, dtype=np.int)
fov_hr = np.asarray(fov_hr, dtype=np.int)
fov_min = np.asarray(fov_min, dtype=np.int)
fov_dist = np.asarray(fov_dist, dtype=np.float)
fov_lat = np.asarray(fov_lat, dtype=np.float)
fov_lon = np.asarray(fov_lon, dtype=np.float)
fov_sza = np.asarray(fov_sza, dtype=np.float)
# SW
crs_swd = np.asarray(crs_swd, dtype=np.float)
ssf_swd = np.asarray(ssf_swd, dtype=np.float)
ff3_swd = np.asarray(ff3_swd, dtype=np.float)
# ff4_swd = np.asarray(ff4_swd, dtype=np.float)
obs_swd = np.asarray(obs_swd, dtype=np.float)
# LW
crs_lwd = np.asarray(crs_lwd, dtype=np.float)
ssf_lwd = np.asarray(ssf_lwd, dtype=np.float)
ff3_lwd = np.asarray(ff3_lwd, dtype=np.float)
# ff4_lwd = np.asarray(ff4_lwd, dtype=np.float)
obs_lwd = np.asarray(obs_lwd, dtype=np.float)


# =====================================================


# consider daytime (T) or nighttime (F) only
if day is True:
    day_str = 'Daytime'
    crs_lwd[fov_sza > 90] = np.nan
    crs_swd[fov_sza > 90] = np.nan
    ssf_lwd[fov_sza > 90] = np.nan
    ssf_swd[fov_sza > 90] = np.nan
    ff3_lwd[fov_sza > 90] = np.nan
    ff3_swd[fov_sza > 90] = np.nan
    # ff4_lwd[fov_sza > 90] = np.nan
    # ff4_swd[fov_sza > 90] = np.nan
    obs_swd[fov_sza > 90] = np.nan
    obs_lwd[fov_sza > 90] = np.nan

elif day is False:
    day_str = 'Nighttime'
    crs_lwd[fov_sza < 90] = np.nan
    crs_swd[fov_sza < 90] = np.nan
    ssf_lwd[fov_sza < 90] = np.nan
    ssf_swd[fov_sza < 90] = np.nan
    ff3_lwd[fov_sza < 90] = np.nan
    ff3_swd[fov_sza < 90] = np.nan
    # ff4_lwd[fov_sza < 90] = np.nan
    # ff4_swd[fov_sza < 90] = np.nan
    obs_swd[fov_sza < 90] = np.nan
    obs_lwd[fov_sza < 90] = np.nan


# =====================================================

# ignore footprint if one of the FOV fluxes or the OBS are NaN
for i, (crs, ssf, ff3, obs) in enumerate(list(zip(crs_lwd, ssf_lwd, ff3_lwd, obs_lwd))):
    if np.isnan(crs) or np.isnan(ssf) or np.isnan(ff3) or np.isnan(obs):
        crs_lwd[i] = np.nan
        ssf_lwd[i] = np.nan
        ff3_lwd[i] = np.nan
        # ff4_lwd[i] = np.nan
        obs_lwd[i] = np.nan
        # print(i, crs_lwd[i], ssf_lwd[i], ff3_lwd[i], obs_lwd[i])
    print(i, crs_lwd[i], ssf_lwd[i], ff3_lwd[i], obs_lwd[i])

# remove NaNs so that only valid FOVs remain
crs_lwd = crs_lwd[~np.isnan(crs_lwd)]
ssf_lwd = ssf_lwd[~np.isnan(ssf_lwd)]
ff3_lwd = ff3_lwd[~np.isnan(ff3_lwd)]
# ff4_lwd = ff4_lwd[~np.isnan(ff4_lwd)]
obs_lwd = obs_lwd[~np.isnan(obs_lwd)]

# print(crs_lwd)
# print(ssf_lwd)
# print(ff3_lwd)

N_crs_lwd = sum(~np.isnan(crs_lwd))
N_ssf_lwd = sum(~np.isnan(ssf_lwd))
N_ff3_lwd = sum(~np.isnan(ff3_lwd))
# N_ff4_lwd = sum(~np.isnan(ff4_lwd))
N_obs_lwd = sum(~np.isnan(obs_lwd))

print('=====================')
print('N_{CRS} = ', N_crs_lwd)
print('N_{SSF} = ', N_ssf_lwd)
print('N_{FF3} = ', N_ff3_lwd)
# print('N_{FF4} = ', N_ff4_lwd)
# print('N_{OBS} = ', N_obs_lwd)

# =====================================================


crs_lw_diff = crs_lwd - obs_lwd
crs_mean_fov_lwd = np.nanmean(crs_lwd)
crs_mean_obs_lwd = np.nanmean(crs_lwd)
crs_lw_mean_bias = np.nanmean(crs_lw_diff)
crs_lw_med_bias = np.nanmedian(crs_lw_diff)
crs_lw_rms_diff = np.sqrt(np.nanmean(crs_lw_diff**2))
crs_lw_corr = ma.corrcoef(ma.masked_invalid(crs_lwd),
                          ma.masked_invalid(obs_lwd))[0][1]


ssf_lw_diff = ssf_lwd - obs_lwd
ssf_mean_fov_lwd = np.nanmean(ssf_lwd)
ssf_mean_obs_lwd = np.nanmean(ssf_lwd)
ssf_lw_mean_bias = np.nanmean(ssf_lw_diff)
ssf_lw_med_bias = np.nanmedian(ssf_lw_diff)
ssf_lw_rms_diff = np.sqrt(np.nanmean(ssf_lw_diff**2))
ssf_lw_corr = ma.corrcoef(ma.masked_invalid(ssf_lwd),
                          ma.masked_invalid(obs_lwd))[0][1]


ff3_lw_diff = ff3_lwd - obs_lwd
ff3_mean_fov_lwd = np.nanmean(ff3_lwd)
ff3_mean_obs_lwd = np.nanmean(ff3_lwd)
ff3_lw_mean_bias = np.nanmean(ff3_lw_diff)
ff3_lw_med_bias = np.nanmedian(ff3_lw_diff)
ff3_lw_rms_diff = np.sqrt(np.nanmean(ff3_lw_diff**2))
ff3_lw_corr = ma.corrcoef(ma.masked_invalid(ff3_lwd),
                          ma.masked_invalid(obs_lwd))[0][1]

# ff4_lw_diff = ff4_lwd - obs_lwd
# ff4_mean_fov_lwd = np.nanmean(ff4_lwd)
# ff4_mean_obs_lwd = np.nanmean(ff4_lwd)
# ff4_lw_mean_bias = np.nanmean(ff4_lw_diff)
# ff4_lw_med_bias = np.nanmedian(ff4_lw_diff)
# ff4_lw_rms_diff = np.sqrt(np.nanmean(ff4_lw_diff**2))
# ff4_lw_corr = ma.corrcoef(ma.masked_invalid(ff4_lwd),
#                           ma.masked_invalid(obs_lwd))[0][1]

print('============================================================')
print('LW Bias:', crs_lw_mean_bias, ssf_lw_mean_bias, ff3_lw_mean_bias)
print('LW RMSD:', crs_lw_rms_diff, ssf_lw_rms_diff, ff3_lw_rms_diff)
print('LW Corr:', crs_lw_corr, ssf_lw_corr, ff3_lw_corr)
print('============================================================')
# =====================================================

crs_lw_stats = "N = " + str(sum(~np.isnan(crs_lwd))) + "\n" + \
               r"Bias ($\bar{\Delta}$) = " + str(np.around(crs_lw_mean_bias, 2)) + "\n" + \
               "RMSD = " + str(np.around(crs_lw_rms_diff, 2)) + "\n" + \
               "Corr = " + str(np.around(crs_lw_corr, 2))

ssf_lw_stats = "N = " + str(sum(~np.isnan(ssf_lwd))) + "\n" + \
               r"Bias ($\bar{\Delta}$) = " + str(np.around(ssf_lw_mean_bias, 2)) + "\n" + \
               "RMSD = " + str(np.around(ssf_lw_rms_diff, 2)) + "\n" + \
               "Corr = " + str(np.around(ssf_lw_corr, 2))

ff3_lw_stats = "N = " + str(sum(~np.isnan(ff3_lwd))) + "\n" + \
               r"Bias ($\bar{\Delta}$) = " + str(np.around(ff3_lw_mean_bias, 2)) + "\n" + \
               "RMSD = " + str(np.around(ff3_lw_rms_diff, 2)) + "\n" + \
               "Corr = " + str(np.around(ff3_lw_corr, 2))

fig, axs = plt.subplots(2, 3, figsize=(14, 8))
fig.suptitle(r'Surface Longwave (LW$\downarrow$)'+' Flux Validation \n '
             'Comparison of CERES CRS, SSF Ed4A, and FF v3C \n' +
              satellite + ' ' + flight_model + ' - ' +
              yr_mon[:3] + ' ' + yr_mon[-4:] + ' - ' + day_str + ' Only')

min_lwd = 50
max_lwd = 500
x = np.linspace(min_lwd, max_lwd, 100)

# 2d hist bins
bins = np.arange(min_lwd, max_lwd+50, 12.5)

axs[0, 0].plot(x, x, 'k')
_, _, _, im = axs[0, 0].hist2d(obs_lwd, crs_lwd, bins=[bins, bins], cmap='BuPu', vmin=0, vmax=10)
cb_ax = fig.add_axes([0.905, 0.53, 0.015, 0.35])
cb = plt.colorbar(im, cax=cb_ax)
# axs[0, 0].plot(obs_lwd, crs_lwd, 'o', color='grey', markersize=0.4)
axs[0, 0].scatter(obs_lwd, crs_lwd, c=obs_lwd, cmap='plasma', s=0.3)
axs[0, 0].set_xlim(min_lwd, max_lwd)
axs[0, 0].set_ylim(min_lwd, max_lwd)
axs[0, 0].set_ylabel(r'CRS Computed LW$\downarrow$ [W m$^{-2}$]')
axs[0, 0].set_xlabel(r'Measured LW$\downarrow$ [W m$^{-2}$]')
axs[0, 0].grid()
axs[0, 0].set_axisbelow(True)
axs[0, 0].text(75, 370, crs_lw_stats)

axs[0, 1].plot(x, x, 'k')
axs[0, 1].hist2d(obs_lwd, ssf_lwd, bins=[bins, bins], cmap='BuPu', vmin=0, vmax=10)
# axs[0, 1].plot(obs_lwd, ssf_lwd, 'o', color='grey', markersize=0.4)
axs[0, 1].scatter(obs_lwd, ssf_lwd, c=obs_lwd, cmap='plasma', s=0.3)
axs[0, 1].set_xlim(min_lwd, max_lwd)
axs[0, 1].set_ylim(min_lwd, max_lwd)
axs[0, 1].set_ylabel(r'SSF Ed4A (Model B) LW$\downarrow$ [W m$^{-2}$]')
axs[0, 1].set_xlabel(r'Measured LW$\downarrow$ [W m$^{-2}$]')
axs[0, 1].grid()
axs[0, 1].set_axisbelow(True)
axs[0, 1].text(75, 370, ssf_lw_stats)

axs[0, 2].plot(x, x, 'k')
axs[0, 2].hist2d(obs_lwd, ff3_lwd, bins=[bins, bins], cmap='BuPu', vmin=0, vmax=10)
# axs[0, 2].plot(obs_lwd, ff3_lwd, 'o', color='grey', markersize=0.4)
axs[0, 2].scatter(obs_lwd, ff3_lwd, c=obs_lwd, cmap='plasma', s=0.3)
axs[0, 2].set_xlim(min_lwd, max_lwd)
axs[0, 2].set_ylim(min_lwd, max_lwd)
axs[0, 2].set_ylabel(r'FF v3C (Model B) LW$\downarrow$ [W m$^{-2}$]')
axs[0, 2].set_xlabel(r'Measured LW$\downarrow$ [W m$^{-2}$]')
axs[0, 2].grid()
axs[0, 2].set_axisbelow(True)
axs[0, 2].text(75, 370, ff3_lw_stats)

axs[1, 0].hist(crs_lwd-obs_lwd, bins=np.arange(-102.5, 105, 5), color='#A9CCE3')
axs[1, 0].set_xlim(-125, 125)
axs[1, 0].set_ylim(0, 100)
axs[1, 0].set_xlabel(r'Flux Difference [W m$^{-2}$]'+'\n'+'(CRS - OBS)')
axs[1, 0].set_ylabel('# of CERES FOVs per Bin')
axs[1, 0].grid()
axs[1, 0].set_axisbelow(True)

axs[1, 1].hist(ssf_lwd-obs_lwd, bins=np.arange(-102.5, 105, 5), color='#A9CCE3')
axs[1, 1].set_xlim(-125, 125)
axs[1, 1].set_ylim(0, 100)
axs[1, 1].set_xlabel(r'Flux Difference [W m$^{-2}$]'+'\n'+'(SSF4A - OBS)')
axs[1, 1].grid()
axs[1, 1].set_axisbelow(True)

axs[1, 2].hist(ff3_lwd-obs_lwd, bins=np.arange(-102.5, 105, 5), color='#A9CCE3')
axs[1, 2].set_xlim(-125, 125)
axs[1, 2].set_ylim(0, 100)
axs[1, 2].set_xlabel(r'Flux Difference [W m$^{-2}$]'+'\n'+'(FFv3C - OBS)')
axs[1, 2].grid()
axs[1, 2].set_axisbelow(True)

plt.show()


# =====================================================

for i, (crs, ssf, ff3, obs) in enumerate(list(zip(crs_swd, ssf_swd, ff3_swd, obs_swd))):
    if np.isnan(crs) or np.isnan(ssf) or np.isnan(ff3) or np.isnan(obs):
        crs_swd[i] = np.nan
        ssf_swd[i] = np.nan
        ff3_swd[i] = np.nan
        # ff4_swd[i] = np.nan
        obs_swd[i] = np.nan
        # print(i, crs_swd[i], ssf_swd[i], ff3_swd[i], obs_swd[i])
    print(i, crs_swd[i], ssf_swd[i], ff3_swd[i], obs_swd[i])


crs_swd = crs_swd[~np.isnan(crs_swd)]
ssf_swd = ssf_swd[~np.isnan(ssf_swd)]
ff3_swd = ff3_swd[~np.isnan(ff3_swd)]
# ff4_swd = ff4_swd[~np.isnan(ff4_swd)]
obs_swd = obs_swd[~np.isnan(obs_swd)]

# print(crs_swd)
# print(ssf_swd)
# print(ff3_swd)

print('N_{CRS} = ', sum(~np.isnan(crs_swd)))
print('N_{SSF} = ', sum(~np.isnan(ssf_swd)))
print('N_{FF3} = ', sum(~np.isnan(ff3_swd)))
print('N_{OBS} = ', sum(~np.isnan(obs_swd)))


crs_sw_diff = crs_swd - obs_swd
crs_mean_fov_swd = np.nanmean(crs_swd)
crs_mean_obs_swd = np.nanmean(crs_swd)
crs_sw_mean_bias = np.nanmean(crs_sw_diff)
crs_sw_med_bias = np.nanmedian(crs_sw_diff)
crs_sw_rms_diff = np.sqrt(np.nanmean(crs_sw_diff**2))
crs_sw_corr = ma.corrcoef(ma.masked_invalid(crs_swd),
                          ma.masked_invalid(obs_swd))[0][1]

ssf_sw_diff = ssf_swd - obs_swd
ssf_mean_fov_swd = np.nanmean(ssf_swd)
ssf_mean_obs_swd = np.nanmean(ssf_swd)
ssf_sw_mean_bias = np.nanmean(ssf_sw_diff)
ssf_sw_med_bias = np.nanmedian(ssf_sw_diff)
ssf_sw_rms_diff = np.sqrt(np.nanmean(ssf_sw_diff**2))
ssf_sw_corr = ma.corrcoef(ma.masked_invalid(ssf_swd),
                          ma.masked_invalid(obs_swd))[0][1]

ff3_sw_diff = ff3_swd - obs_swd
ff3_mean_fov_swd = np.nanmean(ff3_swd)
ff3_mean_obs_swd = np.nanmean(ff3_swd)
ff3_sw_mean_bias = np.nanmean(ff3_sw_diff)
ff3_sw_med_bias = np.nanmedian(ff3_sw_diff)
ff3_sw_rms_diff = np.sqrt(np.nanmean(ff3_sw_diff**2))
ff3_sw_corr = ma.corrcoef(ma.masked_invalid(ff3_swd),
                          ma.masked_invalid(obs_swd))[0][1]

#
# ff4_sw_diff = ff4_swd - obs_swd
# ff4_mean_fov_swd = np.nanmean(ff4_swd)
# ff4_mean_obs_swd = np.nanmean(ff4_swd)
# ff4_sw_mean_bias = np.nanmean(ff4_sw_diff)
# ff4_sw_med_bias = np.nanmedian(ff4_sw_diff)
# ff4_sw_rms_diff = np.sqrt(np.nanmean(ff4_sw_diff**2))
# ff4_sw_corr = ma.corrcoef(ma.masked_invalid(ff4_swd),
#                           ma.masked_invalid(obs_swd))[0][1]

print('SW Bias:', crs_sw_mean_bias, ssf_sw_mean_bias, ff3_sw_mean_bias)# , ff4_sw_mean_bias)
print('SW RMSD:', crs_sw_rms_diff, ssf_sw_rms_diff, ff3_sw_rms_diff)# , ff4_sw_rms_diff)
print('SW Corr:', crs_sw_corr, ssf_sw_corr, ff3_sw_corr)# , ff4_sw_corr)

crs_sw_stats = "N = " + str(sum(~np.isnan(crs_swd))) + "\n" + \
               r"Bias ($\bar{\Delta}$) = " + str(np.around(crs_sw_mean_bias, 2)) + "\n" + \
               "RMSD = " + str(np.around(crs_sw_rms_diff, 2)) + "\n" + \
               "Corr = " + str(np.around(crs_sw_corr, 2))

ssf_sw_stats = "N = " + str(sum(~np.isnan(ssf_swd))) + "\n" + \
               r"Bias ($\bar{\Delta}$) = " + str(np.around(ssf_sw_mean_bias, 2)) + "\n" + \
               "RMSD = " + str(np.around(ssf_sw_rms_diff, 2)) + "\n" + \
               "Corr = " + str(np.around(ssf_sw_corr, 2))

ff3_sw_stats = "N = " + str(sum(~np.isnan(ff3_swd))) + "\n" + \
               r"Bias ($\bar{\Delta}$) = " + str(np.around(ff3_sw_mean_bias, 2)) + "\n" + \
               "RMSD = " + str(np.around(ff3_sw_rms_diff, 2)) + "\n" + \
               "Corr = " + str(np.around(ff3_sw_corr, 2))

# ff4_sw_stats = "N = " + str(sum(~np.isnan(ff4_swd))) + "\n" + \
#                r"Bias ($\bar{\Delta}$) = " + str(np.around(ff4_sw_mean_bias, 2)) + "\n" + \
#                "RMSD = " + str(np.around(ff4_sw_rms_diff, 2)) + "\n" + \
#                "Corr = " + str(np.around(ff4_sw_corr, 2))


# SW figure
fig, axs = plt.subplots(2, 3, figsize=(14, 8))
fig.suptitle(r'Surface Shortwave (SW$\downarrow$)'+' Flux Validation \n '
             'Comparison of CERES CRS, SSF Ed4A, and FF v3C \n' +
             satellite + ' ' + flight_model + ' - ' +
             yr_mon[:3] + ' ' + yr_mon[-4:] + ' - ' + day_str + ' Only')

min_swd = 0
max_swd = 1100
x = np.linspace(min_swd, max_swd, 100)

# 2d hist bins
bins = np.arange(min_swd, max_swd+50, 25)

axs[0, 0].plot(x, x, 'k')
axs[0, 0].hist2d(obs_swd, crs_swd, bins=[bins, bins], cmap='BuPu', vmin=0, vmax=10)
# axs[0, 0].plot(obs_swd, crs_swd, 'o', color='grey', markersize=0.4)
axs[0, 0].scatter(obs_swd, crs_swd, c=crs_swd, cmap='plasma', s=0.3)
axs[0, 0].set_xlim(min_swd, max_swd)
axs[0, 0].set_ylim(min_swd, max_swd)
axs[0, 0].set_ylabel(r'CRS Computed SW$\downarrow$ [W m$^{-2}$]')
axs[0, 0].set_xlabel(r'Measured SW$\downarrow$ [W m$^{-2}$]')
axs[0, 0].grid()
axs[0, 0].set_axisbelow(True)
axs[0, 0].text(80, 800, crs_sw_stats)

axs[0, 1].plot(x, x, 'k')
axs[0, 1].hist2d(obs_swd, ssf_swd, bins=[bins, bins], cmap='BuPu', vmin=0, vmax=10)
# axs[0, 1].plot(obs_swd, ssf_swd, 'o', color='grey', markersize=0.4)
axs[0, 1].scatter(obs_swd, ssf_swd, c=ssf_swd, cmap='plasma', s=0.3)
axs[0, 1].set_xlim(min_swd, max_swd)
axs[0, 1].set_ylim(min_swd, max_swd)
axs[0, 1].set_ylabel(r'SSF Ed4A (Model B) SW$\downarrow$ [W m$^{-2}$]', labelpad=-0.25)
axs[0, 1].set_xlabel(r'Measured SW$\downarrow$ [W m$^{-2}$]')
axs[0, 1].grid()
axs[0, 1].set_axisbelow(True)
axs[0, 1].text(80, 800, ssf_sw_stats)

axs[0, 2].plot(x, x, 'k')
_, _, _, im = axs[0, 2].hist2d(obs_swd, ff3_swd, bins=[bins, bins], cmap='BuPu', vmin=0, vmax=10)
cb_ax = fig.add_axes([0.905, 0.53, 0.015, 0.35])
cb = plt.colorbar(im, cax=cb_ax)
# axs[0, 2].plot(obs_swd, ff3_swd, 'o', color='grey', markersize=0.4)
axs[0, 2].scatter(obs_swd, ff3_swd, c=ff3_swd, cmap='plasma', s=0.3)
axs[0, 2].set_xlim(min_swd, max_swd)
axs[0, 2].set_ylim(min_swd, max_swd)
# axs[0, 2].set_ylabel(r'FF v4A$_{\beta}$ (Model B) SW$\downarrow$ [W m$^{-2}$]', labelpad=-0.25)
axs[0, 2].set_ylabel(r'FF v3C (Model B) SW$\downarrow$ [W m$^{-2}$]', labelpad=-0.25)
axs[0, 2].set_xlabel(r'Measured SW$\downarrow$ [W m$^{-2}$]')
axs[0, 2].grid()
axs[0, 2].set_axisbelow(True)
axs[0, 2].text(80, 800, ff3_sw_stats)


axs[1, 0].hist(crs_swd-obs_swd, bins=np.arange(-152.5, 155, 5), color='#A9CCE3')
axs[1, 0].set_xlim(-150, 150)
axs[1, 0].set_ylim(0, 50)
axs[1, 0].set_xlabel(r'Flux Difference [W m$^{-2}$]'+'\n'+'(CRS - OBS)')
axs[1, 0].set_ylabel('# of CERES FOVs per Bin')
axs[1, 0].grid()
axs[1, 0].set_axisbelow(True)

axs[1, 1].hist(ssf_swd-obs_swd, bins=np.arange(-152.5, 155, 5), color='#A9CCE3')
axs[1, 1].set_xlim(-150, 150)
axs[1, 1].set_ylim(0, 50)
axs[1, 1].set_xlabel(r'Flux Difference [W m$^{-2}$]'+'\n'+'(SSF4A - OBS)')
axs[1, 1].grid()
axs[1, 1].set_axisbelow(True)

axs[1, 2].hist(ff3_swd-obs_swd, bins=np.arange(-152.5, 155, 5), color='#A9CCE3')
axs[1, 2].set_xlim(-150, 150)
axs[1, 2].set_ylim(0, 50)
axs[1, 2].set_xlabel(r'Flux Difference [W m$^{-2}$]'+'\n'+r'(FFv3C - OBS)')
axs[1, 2].grid()
axs[1, 2].set_axisbelow(True)

plt.show()
