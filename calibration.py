from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import ccdproc as ccdp
from astropy import units as u
from astropy.nddata import CCDData



cal_loc = ("analyse/master_calibration")
fits_loc = ("analyse/2_m44_chart2")

filter = "i"
darktime = "5s"
set = "2"

flat = fits.getdata(f"{cal_loc}/2_master_flat_{filter}.fit")
dark = fits.getdata(f"{cal_loc}/2_master_dark_{darktime}.fit")
#bias = fits.getdata(f"{cal_loc}/2_master_bias.fit")


def calibration(light):
    reduced = light - dark
    return reduced


stacked_img = []

for i in range(24,33):
    light = fits.getdata(f"{fits_loc}/m44_2-00{i+1}_{filter}_b.fit")
    stacked_img.append(calibration(light))

#gain_constant = 26
constant = 1

stacked_img = np.sum(stacked_img, axis=0)
stacked_img = stacked_img/flat
stacked_img = stacked_img * constant
calibratedfit = fits.PrimaryHDU(stacked_img)

calibratedfit.header["filter"] = filter

calibratedfit.writeto(f"analyse/calibrated_img/2_{set}_m44_{filter}_b.fit", overwrite=True)
