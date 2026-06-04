from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import astroalign as aa


bright = True
filter = "i"
darktime = "1s"
fov = "2"
begin = 1
end = 23

cal_loc = ("analyse/master_calibration")
fits_loc = (f"analyse/2_m44_chart{fov}")

flat = fits.getdata(f"{cal_loc}/2_master_flat_{filter}.fit")
dark = fits.getdata(f"{cal_loc}/2_master_dark_{darktime}.fit")
#bias = fits.getdata(f"{cal_loc}/2_master_bias.fit")


def calibration(light):
    reduced = light - dark
    return reduced

if bright:
    filter = filter + "_b"


science_list = []

for i in range(begin,end):
    light_header = fits.getheader(f"{fits_loc}/m44_{fov}-00{i+1}_{filter}.fit")
    light = fits.getdata(f"{fits_loc}/m44_{fov}-00{i+1}_{filter}.fit")
    science_list.append(calibration(light)/flat)

target = science_list[0]
aligned_img = []
for i in range(3,len(science_list)):
    #if i == 17 or i == 18:
    #    continue
    source = science_list[i]
    registered_image, footprint = aa.register(
    np.nan_to_num(source, nan=0.0, posinf=0.0, neginf=0.0),
    np.nan_to_num(target, nan=0.0, posinf=0.0, neginf=0.0))
    #plt.imshow(registered_image, vmin=np.min(registered_image), vmax=np.median(registered_image), cmap="gray")
    #plt.show()
    aligned_img.append(registered_image)

stacked_img = np.sum(aligned_img, axis=0)

calibratedfit = fits.PrimaryHDU(stacked_img)

calibratedfit.header["filter"] = filter

calibratedfit.writeto(f"analyse/calibrated_img/{fov}_{darktime}_m44_{filter}.fit", overwrite=True)
