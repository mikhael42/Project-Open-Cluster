from astropy.io import fits
from astropy.visualization import ZScaleInterval, PowerStretch, ImageNormalize
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm
import numpy as np
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry, ApertureStats
# %% meta-data

informatie = fits.getheader('1_60s_m44_g.fit') #voer hier de naam van jouw bestand in 
print(informatie)

# %% data inladen en zichtbaar maken
data_file = fits.getdata('1_60s_m44_g.fit') #ook hier dezelfde naam van jouw bestand invoeren

#plt.imshow(data_file)
plt.imshow (data_file,
            cmap='grey',
            norm=LogNorm())
plt.show()

#voor z-scale 
#norm = ImageNormalize(
    #data_file,
    #interval=ZScaleInterval(),
    #stretch=PowerStretch(0.5)
#)

plt.figure(figsize = (8,8))


#voorbeeld coordinaten ster (x, y) = (3292, 1027)

data_fotometrie = data_file[900:1300, 3000:3600] #pas deze waarden aan zodat je een ster in het midden hebt
plt.imshow (data_fotometrie,
            cmap='grey',
            norm=LogNorm())

plt.show()


positions = [(292, 117)] #voer hier de coordinaten van jouw ster in. 
aperture = CircularAperture(positions, r=8)
annulus_aperture = CircularAnnulus(positions, r_in=12, r_out=20)

plt.imshow (data_fotometrie,
            cmap='grey',
            norm=LogNorm())

ap_patches = aperture.plot(color='yellow', lw=2,

                           label='Photometry aperture')

ann_patches = annulus_aperture.plot(color='red', lw=2,

                                    label='Background annulus')

handles = (ap_patches[0], ann_patches[0])

plt.legend(loc=(0.17, 0.05), facecolor='#458989', labelcolor='white',

           handles=handles, prop={'weight': 'bold', 'size': 6})

phot_table = aperture_photometry(data_fotometrie, aperture)

plt.show() 

aperstats = ApertureStats(data_fotometrie, annulus_aperture)

bkg_mean = aperstats.mean

aperture_area = aperture.area_overlap(data_fotometrie)

total_bkg = bkg_mean * aperture_area 

#nu nog omrekenen naar magnitude 

calibration_constant  = 26.57

phot_bkgsub = phot_table['aperture_sum'] - total_bkg

phot_table['total_bkg'] = total_bkg

phot_table['aperture_sum_bkgsub'] = phot_bkgsub

phot_table['instrumental magnitude'] = -2.5*np.log10(phot_bkgsub)

phot_table['calibrated magnitude'] = -2.5*np.log10(phot_bkgsub) + calibration_constant

for col in phot_table.colnames:

    phot_table[col].info.format = '%.8g'  # for consistent table output



print(phot_table)

#plt.savefig('check_image_1.png') 


