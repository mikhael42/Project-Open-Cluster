from astropy.io import fits
from astropy.visualization import ZScaleInterval, PowerStretch, ImageNormalize
import matplotlib.pyplot as plt
import numpy as np


# %% meta-data

informatie = fits.getheader('2_1_m44_g.fit') #voer hier de naam van jouw bestand in 
print(informatie)

# %% data inladen en zichtbaar maken
data_file = fits.getdata('2_1_m44_g.fit') #ook hier dezelfde naam van jouw bestand invoeren

#plt.imshow(data_file)

#plt.show()#

# %% Data beter zichtbaar maken: speel een beetje met de waarden voor vmin en vmax
norm = ImageNormalize(
    data_file,
    interval=ZScaleInterval(),
    stretch=PowerStretch(0.5)  # experimenteer met deze waarde
)

plt.figure(figsize = (8,8))
plt.imshow(data_file, cmap='gray', norm=norm, origin='lower')
plt.colorbar()
plt.show()
