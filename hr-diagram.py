import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from astropy.io import fits
from scipy.stats import gaussian_kde

filter1 = "g"
filter2 = "i"

data_loc = "analyse/hr-data/star_data_"

def dist(x,y):
    return np.sqrt(x**2 + y**2)

def small_dist(mag,mag_combo,mag_iso,mag_combo_iso):
    distance = []   
    for i in range(len(mag[filter2])):
        temp = []
        for j in range(len(mag_iso[filter2])):
            temp.append(dist(mag_combo_iso[j]- mag_combo[i],mag_iso[filter2][j]-mag[filter2][i]))
        distance.append(min(temp))
    return distance

plt.figure(figsize=(7, 7))
mag_iso = {
    "g" : [],
    "r" : [],
    "i" : []
}

with open("analyse/isochrone/isochrone.dat.txt", "r") as f:
    for lines in f:
        temp = lines.split() 
        if len(temp) == 33:
            mag_iso["i"].append(float(temp[31])+6.27)
            mag_iso["r"].append(float(temp[30])+6.27)
            mag_iso["g"].append(float(temp[29])+6.27)

mag_combo_iso = []
for i in range(len(mag_iso[filter2])):
    mag_combo_iso.append(mag_iso[filter1][i] - mag_iso[filter2][i])
#plt.scatter(mag_combo_iso, mag_iso[filter2], s=10, color='blue')
plt.plot(mag_combo_iso,mag_iso[filter2],'-',color="blue")

mag = {
    "i" : [],
    "r" : [],
    "g" : []
}

id = []
x = []
y = []


# data_file_i_1 = np.asarray(fits.getdata('analyse/calibrated_img/1_60s_m44_i.fit'), dtype=np.float32)
# data_file_i_2 = np.asarray(fits.getdata('analyse/calibrated_img/2_60s_m44_i.fit'), dtype=np.float32)

# plt.figure(figsize=(8, 8))
# plt.imshow(data_file_i_1, cmap='gray', norm=LogNorm())

with open(f"{data_loc}fov_1.csv", 'r') as input_bestand:
    for regel in input_bestand:
        data_opgeknipt = regel.split(',')
        if data_opgeknipt[0] == "ster_idx": 
            continue
        id.append(int(data_opgeknipt[0]))
        x.append(float(data_opgeknipt[1]))
        y.append(float(data_opgeknipt[2]))
        mag['i'].append(float(data_opgeknipt[6])+27.9635)
        mag['g'].append(float(data_opgeknipt[8])+28.9266)
        mag['r'].append(float(data_opgeknipt[7])+28.6672)

# mag_combo = []
# for i in range(len(mag[filter2])):
#     mag_combo.append(mag[filter1][i] - mag[filter2][i])

# x_cl = []
# y_cl = []
# combo_in_cl = []
# mag_cl = {
#     "g" : [],
#     "r" : [],
#     "i": []
# }

# distance = small_dist(mag,mag_combo,mag_iso,mag_combo_iso)

# for i in range(len(distance)):
#     if distance[i] < .8:
#         mag_cl['g'].append(mag['g'][i])
#         mag_cl["i"].append(mag['i'][i])
#         mag_cl['r'].append(mag['r'][i])
#         combo_in_cl.append(mag_combo[i])
#         x_cl.append(x[i])
#         y_cl.append(y[i])

# plt.title("FoV 1")
# plt.plot(x, y,'o', color='blue', markersize=5)
# plt.plot(x_cl, y_cl,'o', color='red', markersize=5)
# plt.savefig('fov_1_clusters.png')
# plt.show()
#plt.scatter(combo_in_cl, mag_cl[filter2], s=10, color='red')

#plt.scatter(mag_combo, mag[filter2], s=10, color='black')

# for j in mag:
#    mag[j].clear()


#x = []
#y = []
#id = []

# plt.figure(figsize=(8, 8))
# plt.imshow(data_file_i_2, cmap='gray', norm=LogNorm())

with open(f"{data_loc}fov_2.csv", 'r') as input_bestand:
    for regel in input_bestand:
        data_opgeknipt = regel.split(',')
        if data_opgeknipt[0] == "ster_idx": 
            continue
        id.append(int(data_opgeknipt[0])+229)
        x.append(float(data_opgeknipt[1]))
        y.append(float(data_opgeknipt[2]))
        mag['i'].append(float(data_opgeknipt[6])+27.9635)
        mag['g'].append(float(data_opgeknipt[7])+28.9266)
        mag['r'].append(float(data_opgeknipt[8])+28.6672)

mag_combo = []
for i in range(len(mag[filter2])):
    mag_combo.append(mag[filter1][i] - mag[filter2][i])

plt.scatter(mag_combo, mag[filter2], s=10, color='black')

# x_cl = []
# y_cl = []

# distance = small_dist(mag,mag_combo,mag_iso,mag_combo_iso)

# for i in range(len(distance)):
#     if distance[i] < .8:
#         mag_cl['g'].append(mag['g'][i])
#         mag_cl["i"].append(mag['i'][i])
#         mag_cl['r'].append(mag['r'][i])
#         combo_in_cl.append(mag_combo[i])
#         x_cl.append(x[i])
#         y_cl.append(y[i])

# plt.title("FoV 2")
# plt.plot(x, y,'o', color='blue', markersize=5)
# plt.plot(x_cl, y_cl,'o', color='red', markersize=5)
# plt.savefig('fov_2_clusters.png')
# plt.show()


# for i in range(1,3):
#     for j in mag:
#         mag[j].clear()

#     with open(f"{data_loc}fov_{i}_b.csv", 'r') as input_bestand:
#         for regel in input_bestand:
#             data_opgeknipt = regel.split(',')
#             if data_opgeknipt[0] == "ster_idx": 
#                 continue
#             mag['i'].append(float(data_opgeknipt[6])+26.1890)
#             mag['g'].append(float(data_opgeknipt[7])+27.4417)
#             mag['r'].append(float(data_opgeknipt[8])+26.3697)

#     mag_combo = []
#     for j in range(len(mag[filter2])):
#         mag_combo.append(mag[filter1][j] - mag[filter2][j])
#     if i == 1:
#         plt.scatter(mag_combo, mag[filter2], s=10, color='yellow')
#     else: 
#         plt.scatter(mag_combo, mag[filter2], s=10, color='green')

combo_in_cl = []
mag_cl = {
    "g" : [],
    "r" : [],
    "i": []
}

distance = small_dist(mag,mag_combo,mag_iso,mag_combo_iso)
lim_cl = .77

for i in range(len(distance)):
    if distance[i] <= lim_cl:
        mag_cl['g'].append(mag['g'][i])
        mag_cl["i"].append(mag['i'][i])
        mag_cl['r'].append(mag['r'][i])
        combo_in_cl.append(mag_combo[i])
        #x_cl.append(x[i])
        #y_cl.append(y[i])

plt.scatter(combo_in_cl, mag_cl[filter2], s=10, color='red')

plt.xlabel(f'{filter1} − {filter2}  (kleur)')
plt.ylabel(f'Instrumentele magnitude ({filter2})')
plt.ylim(5,20)
plt.xlim(right=5)
plt.gca().invert_yaxis()
plt.savefig('final_hr_iso.png')
plt.clf()

fd_bins = np.histogram_bin_edges(distance, bins='fd')
n_bins = int(len(fd_bins) * 1.5)

plt.hist(distance, bins=n_bins, density=True, edgecolor='black')

kde = gaussian_kde(distance)
x = np.linspace(min(distance), max(distance), 500)

plt.plot(x, kde(x), linewidth=2)
plt.plot([lim_cl,lim_cl],[0,1],linewidth=2,color="red")
plt.ylim(0,.7)
plt.savefig('histo_iso.png')
plt.show()