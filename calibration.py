from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import ccdproc as ccdp
from astropy import units as u
from astropy.nddata import CCDData



cal_loc = ("analyse/master_calibration")
fits_loc = ("analyse/2_m44_chart1")

flat = fits.getdata(f"{cal_loc}/2_master_flat_r.fit")
dark = fits.getdata(f"{cal_loc}/2_master_dark_60s.fit")
bias = fits.getdata(f"{cal_loc}/2_master_bias.fit")
light = fits.getdata(f"{fits_loc}/m44_1-0001_r.fit")
#light = CCDData(fits.getdata(f"{fits_loc}/m44_1-0001_g.fit"), unit=u.adu)
#bias = CCDData(fits.getdata(f"{cal_loc}/2_master_bias.fit"), unit=u.adu)
#dark = CCDData(fits.getdata(f"{cal_loc}/2_master_dark_60s.fit"), unit=u.adu)
#flat = CCDData(fits.getdata(f"{cal_loc}/2_master_flat_g.fit"), unit=u.adu)
#dark.header = {'exptime': 60.0}
#light.header = {'exptime': 60.0}
gain_constant = 100

def calibration(light):
    const_temp = dark-bias
    reduced = light - dark
    reduced = reduced/(flat-dark)
    reduced = reduced * gain_constant
    #reduced = ccdp.subtract_bias(light, bias)
    #reduced = ccdp.subtract_dark(reduced, dark, exposure_time="exptime", exposure_unit=u.second)
    #reduced = ccdp.flat_correct(reduced, flat)
    return reduced


calibratedfit = fits.PrimaryHDU(calibration(light))

calibratedfit.writeto("analyse/calibrated_img/2_m44_1-1_r.fit", overwrite=True)
>>>>>>> Stashed changes

for i in range(10):
    print("later for automation")
