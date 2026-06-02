from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import pathlib as Path

cal_loc = Path("analyse/master_calibration")


def calibration():
    flat = 
    dark =
    bias = 
    light = 
    constant_gain = 

    const = (dark+bias)
    reduced = light - const
    reduced *= constant_gain
    calibrated = reduced/(flat-const)



for i in range(10):
    print("later for automation")


