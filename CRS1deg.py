import numpy as np
import cerestools as ceres

path = '/Users/rcscott2/Desktop/CERES/EBAF/'
file = 'CERES_EBAF-TOA_Ed4.1_Subset_200003-201910.nc'
file_path = path + file

ceres.print_nc_file_info(file_path)

lat, lat_name, lat_units = ceres.read_ebaf_var(file_path=file_path, variable='lat')
lon, lon_name, lon_units = ceres.read_ebaf_var(file_path=file_path, variable='lon')

# lat = np.array(lat)
# lat_bins = lat - 0.5
# print(lat_bins)
lat_bins = np.arange(-90, 91)
lon_bins = np.arange(0, 361)
print(lat_bins)
print(lon_bins)

# fake FOV lat data
lat_data = np.random.random_sample(20)
lat_data = 180*lat_data-90
print(lat_data)

# fake FOV lon data
lon_data = np.random.random_sample(20)
lon_data = 360*lon_data
print(lon_data)

lat_ind = np.digitize(lat_data, lat_bins)
lon_ind = np.digitize(lon_data, lon_bins)


print(lat_ind)
print(lon_ind)

for n in range(lat_data.size):
    print(lat_bins[lat_ind[n]-1], "<=", lat_data[n], "<", lat_bins[lat_ind[n]], '...', lat_ind[n])
    print(lon_bins[lon_ind[n]-1], "<=", lon_data[n], "<", lon_bins[lon_ind[n]], '...', lon_ind[n])


# gridded = np.empty([180, 360])
# for i in range(181):
#     for j in range(361):
#         for n in range(20):
#             condition = (i == lat_ind[n] and j == lon_ind[n])
#             gridded[i, j] = np.extract(condition, lat_data)


print(np.sort(lat_ind))
print(np.sort(lat_data))

print(np.sort(lon_ind))
print(np.sort(lon_data))

# x = np.array([0.2, 6.4, 3.0, 1.6])
# bins = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])
# inds = np.digitize(x, bins)
#
# for n in range(x.size):
#     print(bins[inds[n]-1], "<=", x[n], "<", bins[inds[n]])
#
#
# # lat
# # lat_bins
# # ind
