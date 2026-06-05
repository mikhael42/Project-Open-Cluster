from astropy.io import fits
from astropy.visualization import ZScaleInterval, PowerStretch, ImageNormalize
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm
import numpy as np
import photutils
import photutils.detection
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry, ApertureStats
# %% meta-data

informatie = fits.getheader('2_60s_m44_i.fit') #voer hier de naam van jouw bestand in 
print(informatie)

# %% data inladen en zichtbaar maken
data_file = fits.getdata('2_60s_m44_i.fit') #ook hier dezelfde naam van jouw bestand invoeren

lijst_sterren_x = [165, 485, 386, 1430, 1662, 2041, 2126, 2339, 2777, 3016, 1430, 1809, 1364, 154, 268, 2534, 2541, 3023, 3038, 2935, 529, 742, 849, 1750, 2155, 2269, 3718, 3612, 1717, 2799, 2994, 2942, 555, 2357, 908, 1000, 1864, 1806, 2574, 3123, 3597, 485, 787, 636, 908, 43, 143, 764, 349, 779, 1857, 1875, 1953, 2048, 2192, 2615, 2254, 2968, 3248, 3251, 3630, 3336, 3711, 3937, 1405, 1943, 2936, 3495, 3782, 3863]
lijst_sterren_y = [413, 307, 292, 399, 237, 553, 527, 560, 546, 513, 741, 829, 1024, 686, 719, 829, 1057, 932, 1013, 1234, 1094, 1344, 1602, 1373, 1186, 1318, 1373, 1212, 1705, 1554, 1631, 1763, 1848, 1808, 2150, 2047, 2282, 2433, 2455, 2260, 2267, 2617, 2738, 2529, 2944, 2977, 2955, 3529, 3717, 3846, 2731, 2786, 2904, 2834, 3029, 3165, 3286, 2735, 3033, 3430, 3404, 3713, 3860, 2780, 2736, 3869, 2304, 501, 434, 285]
#plt.imshow(data_file)
plt.imshow (data_file,
            cmap='grey',
            norm=LogNorm())

#plt.scatter(lijst_sterren_x, lijst_sterren_y)

#for i, (x_val, y_val) in enumerate(zip(lijst_sterren_x, lijst_sterren_y), start=1):
    #plt.annotate(
        #str(i),              # tekst: 1, 2, 3, ...
        #(x_val, y_val),      # positie van het punt
        #xytext=(5, 5),       # verschuiving t.o.v. het punt
        #textcoords="offset points"
    #)



#voor z-scale 
#norm = ImageNormalize(
    #data_file,
    #interval=ZScaleInterval(),
    #stretch=PowerStretch(0.5)
#)

#plt.figure(figsize = (8,8))
iets = photutils.detection.DAOStarFinder(1700, 9.0, ratio=1.0, theta=0.0, sigma_radius=1.2, exclude_border=False, brightest=None, peakmax=None, xycoords=None, min_separation=0.0)

ans = iets.find_stars(data_file, mask=0)
print(ans)
#plt.plot(lijst_sterren_x, lijst_sterren_y, 'o', color='red',)
plt.plot(ans['xcentroid'], ans['ycentroid'], 'o', markersize=5, label='Detected Stars')
plt.show()

#voorbeeld coordinaten ster (x, y) = (3292, 1027)

#data_fotometrie = data_file[900:1300, 3000:3600] #pas deze waarden aan zodat je een ster in het midden hebt
#plt.imshow (data_fotometrie,
            #cmap='grey',
            #norm=LogNorm())


#plt.show()


#positions = [(293, 125)] #voer hier de coordinaten van jouw ster in. 
#aperture = CircularAperture(positions, r=6)
#annulus_aperture = CircularAnnulus(positions, r_in=15, r_out=30)

#plt.imshow (data_fotometrie,
            #cmap='grey',
            #norm=LogNorm())


#plt.show()


#positions = [(293, 125)] #voer hier de coordinaten van jouw ster in. 
#aperture = CircularAperture(positions, r=6)
#annulus_aperture = CircularAnnulus(positions, r_in=15, r_out=30)

#plt.imshow (data_fotometrie,
            #cmap='grey',
            #norm=LogNorm())


#plt.show()


#positions = [(293, 125)] #voer hier de coordinaten van jouw ster in. 
#aperture = CircularAperture(positions, r=6)
#annulus_aperture = CircularAnnulus(positions, r_in=15, r_out=30)

#plt.imshow (data_fotometrie,
            #cmap='grey',
            #norm=LogNorm())


#plt.show()


#positions = [(293, 125)] #voer hier de coordinaten van jouw ster in. 
#aperture = CircularAperture(positions, r=6)
#annulus_aperture = CircularAnnulus(positions, r_in=15, r_out=30)

#plt.imshow (data_fotometrie,
            #cmap='grey',
            #norm=LogNorm())


#plt.show()


#positions = [(293, 125)] #voer hier de coordinaten van jouw ster in. 
#aperture = CircularAperture(positions, r=6)
#annulus_aperture = CircularAnnulus(positions, r_in=15, r_out=30)

#plt.imshow (data_fotometrie,
            #cmap='grey',
            #norm=LogNorm())


#plt.show()


#positions = [(293, 125)] #voer hier de coordinaten van jouw ster in. 
#aperture = CircularAperture(positions, r=6)
#annulus_aperture = CircularAnnulus(positions, r_in=15, r_out=30)

#plt.imshow (data_fotometrie,
            #cmap='grey',
            #norm=LogNorm())


#plt.show()


#positions = [(293, 125)] #voer hier de coordinaten van jouw ster in. 
#aperture = CircularAperture(positions, r=6)
#annulus_aperture = CircularAnnulus(positions, r_in=15, r_out=30)

#plt.imshow (data_fotometrie,
            #cmap='grey',
            #norm=LogNorm())

##ap_patches = aperture.plot(color='yellow', lw=2,

                           #label='Photometry aperture')

#ann_patches = annulus_aperture.plot(color='red', lw=2,

                                    #label='Background annulus')

#handles = (ap_patches[0], ann_patches[0])

#plt.legend(loc=(0.17, 0.05), facecolor='#458989', labelcolor='white',

           #handles=handles, prop={'weight': 'bold', 'size': 6})

#phot_table = aperture_photometry(data_fotometrie, aperture)

#plt.show() 

#aperstats = ApertureStats(data_fotometrie, annulus_aperture)

#bkg_mean = aperstats.mean
#
#aperture_area = aperture.area_overlap(data_fotometrie)

#total_bkg = bkg_mean * aperture_area 

#nu nog omrekenen naar magnitude 

#calibration_constant  = 26.57

#phot_bkgsub = phot_table['aperture_sum'] - total_bkg

#phot_table['total_bkg'] = total_bkg

#phot_table['aperture_sum_bkgsub'] = phot_bkgsub

#phot_table['instrumental magnitude'] = -2.5*np.log10(phot_bkgsub)

#phot_table['calibrated magnitude'] = -2.5*np.log10(phot_bkgsub) + calibration_constant

#for col in phot_table.colnames:

    #phot_table[col].info.format = '%.8g'  # for consistent table output



#print(phot_table)

#plt.savefig('check_image_1.png') 


