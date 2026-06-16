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

#plt.figure(figsize=(7, 7))
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
#plt.plot(mag_combo_iso,mag_iso[filter2],'-',color="blue")

mag = {
    "i" : [],
    "r" : [],
    "g" : []
}

id = []
x = []
y = []


data_file_i_1 = np.asarray(fits.getdata('analyse/calibrated_img/1_60s_m44_i.fit'), dtype=np.float32)
data_file_i_2 = np.asarray(fits.getdata('analyse/calibrated_img/2_60s_m44_i.fit'), dtype=np.float32)

plt.figure(figsize=(8, 8))
plt.imshow(data_file_i_1, cmap='gray', norm=LogNorm())

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

mag_combo = []
for i in range(len(mag[filter2])):
    mag_combo.append(mag[filter1][i] - mag[filter2][i])

x_cl = []
y_cl = []
combo_in_cl = []
mag_cl = {
    "g" : [],
    "r" : [],
    "i": []
}

distance = small_dist(mag,mag_combo,mag_iso,mag_combo_iso)

for i in range(len(distance)):
    if distance[i] < .8:
        mag_cl['g'].append(mag['g'][i])
        mag_cl["i"].append(mag['i'][i])
        mag_cl['r'].append(mag['r'][i])
        combo_in_cl.append(mag_combo[i])
        x_cl.append(x[i])
        y_cl.append(y[i])

plt.title("FoV 1")
plt.plot(x, y,'o', color='blue', markersize=5)
plt.plot(x_cl, y_cl,'o', color='red', markersize=5)
plt.plot([457,360,463,581,488,2669,2115,3146,2941,245], [943,1134,1184,1851,2144,630,1831,2406,3709,3104],'o', color='red', markersize=5)
plt.savefig('fov_1_clusters.png')
plt.show()
#plt.scatter(combo_in_cl, mag_cl[filter2], s=10, color='red')

#plt.scatter(mag_combo, mag[filter2], s=10, color='black')

for j in mag:
   mag[j].clear()


x = []
y = []
id = []

plt.figure(figsize=(8, 8))
plt.imshow(data_file_i_2, cmap='gray', norm=LogNorm())

with open(f"{data_loc}fov_2.csv", 'r') as input_bestand:
    for regel in input_bestand:
        data_opgeknipt = regel.split(',')
        if data_opgeknipt[0] == "ster_idx": 
            continue
        id.append(int(data_opgeknipt[0]))
        x.append(float(data_opgeknipt[1]))
        y.append(float(data_opgeknipt[2]))
        mag['i'].append(float(data_opgeknipt[6])+27.9635)
        mag['g'].append(float(data_opgeknipt[7])+28.9266)
        mag['r'].append(float(data_opgeknipt[8])+28.6672)

mag_combo = []
for i in range(len(mag[filter2])):
    mag_combo.append(mag[filter1][i] - mag[filter2][i])

plt.scatter(mag_combo, mag[filter2], s=10, color='black')

x_cl = []
y_cl = []

distance = small_dist(mag,mag_combo,mag_iso,mag_combo_iso)

for i in range(len(distance)):
    if distance[i] < .8:
        mag_cl['g'].append(mag['g'][i])
        mag_cl["i"].append(mag['i'][i])
        mag_cl['r'].append(mag['r'][i])
        combo_in_cl.append(mag_combo[i])
        x_cl.append(x[i])
        y_cl.append(y[i])

plt.title("FoV 2")
plt.plot(x, y,'o', color='blue', markersize=5)
plt.plot(x_cl, y_cl,'o', color='red', markersize=5)
plt.plot([67,1250,879,1847,3413,3720,1010], [1108,1478,2193,3443,2408,2555,2222],'o', color='red', markersize=5)
#plt.savefig('fov_2_clusters.png')
plt.show()