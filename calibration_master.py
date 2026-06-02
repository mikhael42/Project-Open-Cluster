from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import sigma_clip


#lists to biases on my laptop

biases_2 = ["analyse/2_calibration/bias/Calibration-0001_bias.fit",
            "analyse/2_calibration/bias/Calibration-0002_bias.fit",
            "analyse/2_calibration/bias/Calibration-0003_bias.fit",
            "analyse/2_calibration/bias/Calibration-0004_bias.fit",
            "analyse/2_calibration/bias/Calibration-0005_bias.fit",
            "analyse/2_calibration/bias/Calibration-0006_bias.fit",
            "analyse/2_calibration/bias/Calibration-0007_bias.fit",
            "analyse/2_calibration/bias/Calibration-0008_bias.fit",
            "analyse/2_calibration/bias/Calibration-0009_bias.fit",
            "analyse/2_calibration/bias/Calibration-0010_bias.fit",
            "analyse/2_calibration/bias/Calibration-0011_bias.fit",
            "analyse/2_calibration/bias/Calibration-0012_bias.fit",
            "analyse/2_calibration/bias/Calibration-0013_bias.fit",
            "analyse/2_calibration/bias/Calibration-0014_bias.fit",
            "analyse/2_calibration/bias/Calibration-0015_bias.fit",
            "analyse/2_calibration/bias/Calibration-0016_bias.fit",
            "analyse/2_calibration/bias/Calibration-0017_bias.fit",
            "analyse/2_calibration/bias/Calibration-0018_bias.fit",
            "analyse/2_calibration/bias/Calibration-0019_bias.fit",
            "analyse/2_calibration/bias/Calibration-0020_bias.fit"]

flat_i_2 = [
    "analyse/2_calibration/flat-0001_i.fit",
    "analyse/2_calibration/flat-0002_i.fit",
    "analyse/2_calibration/flat-0003_i.fit",
    "analyse/2_calibration/flat-0004_i.fit",
    "analyse/2_calibration/flat-0005_i.fit",
    "analyse/2_calibration/flat-0006_i.fit",
    "analyse/2_calibration/flat-0007_i.fit",
    "analyse/2_calibration/flat-0008_i.fit",
    "analyse/2_calibration/flat-0009_i.fit",
    "analyse/2_calibration/flat-0010_i.fit",
]
flat_r_2 = [
    "analyse/2_calibration/flat-0001_r.fit",
    "analyse/2_calibration/flat-0002_r.fit",
    "analyse/2_calibration/flat-0003_r.fit",
    "analyse/2_calibration/flat-0004_r.fit",
    "analyse/2_calibration/flat-0005_r.fit",
    "analyse/2_calibration/flat-0006_r.fit",
    "analyse/2_calibration/flat-0007_r.fit",
    "analyse/2_calibration/flat-0008_r.fit",
    "analyse/2_calibration/flat-0009_r.fit",
    "analyse/2_calibration/flat-0010_r.fit",
]
flat_g_2 = [
    "analyse/2_calibration/flat-0001_g.fit",
    "analyse/2_calibration/flat-0002_g.fit",
    "analyse/2_calibration/flat-0003_g.fit",
    "analyse/2_calibration/flat-0004_g.fit",
    "analyse/2_calibration/flat-0005_g.fit",
    "analyse/2_calibration/flat-0006_g.fit",
    "analyse/2_calibration/flat-0007_g.fit",
    "analyse/2_calibration/flat-0008_g.fit",
    "analyse/2_calibration/flat-0009_g.fit",
    "analyse/2_calibration/flat-0010_g.fit",
]
dark_1s_2 = [
    "analyse/2_calibration/Calibration-0001_1s.fit",
    "analyse/2_calibration/Calibration-0002_1s.fit",
    "analyse/2_calibration/Calibration-0003_1s.fit",
    "analyse/2_calibration/Calibration-0004_1s.fit",
    "analyse/2_calibration/Calibration-0005_1s.fit",
    "analyse/2_calibration/Calibration-0006_1s.fit",
    "analyse/2_calibration/Calibration-0007_1s.fit",
    "analyse/2_calibration/Calibration-0008_1s.fit",
    "analyse/2_calibration/Calibration-0009_1s.fit",
    "analyse/2_calibration/Calibration-0010_1s.fit",
]
dark_05s_2 = [
    "analyse/2_calibration/Calibration-0001_05s.fit",
    "analyse/2_calibration/Calibration-0002_05s.fit",
    "analyse/2_calibration/Calibration-0003_05s.fit",
    "analyse/2_calibration/Calibration-0004_05s.fit",
    "analyse/2_calibration/Calibration-0005_05s.fit",
    "analyse/2_calibration/Calibration-0006_05s.fit",
    "analyse/2_calibration/Calibration-0007_05s.fit",
    "analyse/2_calibration/Calibration-0008_05s.fit",
    "analyse/2_calibration/Calibration-0009_05s.fit",
    "analyse/2_calibration/Calibration-0010_05s.fit",
]
dark_5s_2 = [
    "analyse/2_calibration/Calibration-0001_5s.fit",
    "analyse/2_calibration/Calibration-0002_5s.fit",
    "analyse/2_calibration/Calibration-0003_5s.fit",
    "analyse/2_calibration/Calibration-0004_5s.fit",
    "analyse/2_calibration/Calibration-0005_5s.fit",
    "analyse/2_calibration/Calibration-0006_5s.fit",
    "analyse/2_calibration/Calibration-0007_5s.fit",
    "analyse/2_calibration/Calibration-0008_5s.fit",
    "analyse/2_calibration/Calibration-0009_5s.fit",
    "analyse/2_calibration/Calibration-0010_5s.fit",
]
dark_60s_2 = [
    "analyse/2_calibration/Calibration-0001_60s.fit",
    "analyse/2_calibration/Calibration-0002_60s.fit",
    "analyse/2_calibration/Calibration-0003_60s.fit",
    "analyse/2_calibration/Calibration-0004_60s.fit",
    "analyse/2_calibration/Calibration-0005_60s.fit",
    "analyse/2_calibration/Calibration-0006_60s.fit",
    "analyse/2_calibration/Calibration-0007_60s.fit",
    "analyse/2_calibration/Calibration-0008_60s.fit",
    "analyse/2_calibration/Calibration-0009_60s.fit",
    "analyse/2_calibration/Calibration-0010_60s.fit",
]

biases_1 = ["analyse/1_calibration/Calibration-0001_bias.fit",
            "analyse/1_calibration/Calibration-0002_bias.fit",
            "analyse/1_calibration/Calibration-0003_bias.fit",
            "analyse/1_calibration/Calibration-0004_bias.fit",
            "analyse/1_calibration/Calibration-0005_bias.fit",
            "analyse/1_calibration/Calibration-0006_bias.fit",
            "analyse/1_calibration/Calibration-0007_bias.fit",
            "analyse/1_calibration/Calibration-0008_bias.fit",
            "analyse/1_calibration/Calibration-0009_bias.fit",
            "analyse/1_calibration/Calibration-0010_bias.fit",
            "analyse/1_calibration/Calibration-0011_bias.fit",
            "analyse/1_calibration/Calibration-0012_bias.fit",
            "analyse/1_calibration/Calibration-0013_bias.fit",
            "analyse/1_calibration/Calibration-0014_bias.fit",
            "analyse/1_calibration/Calibration-0015_bias.fit",
            "analyse/1_calibration/Calibration-0016_bias.fit",
            "analyse/1_calibration/Calibration-0017_bias.fit",
            "analyse/1_calibration/Calibration-0018_bias.fit",
            "analyse/1_calibration/Calibration-0019_bias.fit",
            "analyse/1_calibration/Calibration-0020_bias.fit"]

flat_r_1 = [
    "analyse/1_calibration/flat-0001r.fit",
    "analyse/1_calibration/flat-0002r.fit",
    "analyse/1_calibration/flat-0003r.fit",
    "analyse/1_calibration/flat-0004r.fit",
    "analyse/1_calibration/flat-0005r.fit",
    "analyse/1_calibration/flat-0006r.fit",
    "analyse/1_calibration/flat-0007r.fit",
    "analyse/1_calibration/flat-0008r.fit",
    "analyse/1_calibration/flat-0009r.fit",
    "analyse/1_calibration/flat-0010r.fit",
]
flat_g_1 = [
    "analyse/1_calibration/flat-0001g.fit",
    "analyse/1_calibration/flat-0002g.fit",
    "analyse/1_calibration/flat-0003g.fit",
    "analyse/1_calibration/flat-0004g.fit",
    "analyse/1_calibration/flat-0005g.fit",
    "analyse/1_calibration/flat-0006g.fit",
    "analyse/1_calibration/flat-0007g.fit",
    "analyse/1_calibration/flat-0008g.fit",
    "analyse/1_calibration/flat-0009g.fit",
    "analyse/1_calibration/flat-0010g.fit",
]
dark_5s_1 = [
    "analyse/1_calibration/Calibration-0001_dark_5s.fit",
    "analyse/1_calibration/Calibration-0002_dark_5s.fit",
    "analyse/1_calibration/Calibration-0003_dark_5s.fit",
    "analyse/1_calibration/Calibration-0004_dark_5s.fit",
    "analyse/1_calibration/Calibration-0005_dark_5s.fit",
    "analyse/1_calibration/Calibration-0006_dark_5s.fit",
    "analyse/1_calibration/Calibration-0007_dark_5s.fit",
    "analyse/1_calibration/Calibration-0008_dark_5s.fit",
    "analyse/1_calibration/Calibration-0009_dark_5s.fit",
    "analyse/1_calibration/Calibration-0010_dark_5s.fit",
]

def mastercombining(list):
    allfits = [fits.getdata(i) for i in list]
    stack = np.stack(allfits, axis=0)

    clipped = sigma_clip(stack, sigma=5, axis=0, masked=False)

    masterfit = np.nanmean(clipped, axis=0) #masterfit = np.nanmean(clipped, axis=0)
    return masterfit

directory_mean2 = ["2_master_bias","2_master_flat_i","2_master_flat_r","2_master_flat_g","2_master_dark_1s","2_master_dark_0.5s","2_master_dark_5s","2_master_dark_60s", "1_master_bias", "1_master_flat_r", "1_master_flat_g", "1_master_dark_5s"]
directory_list2 = [biases_2,flat_i_2,flat_r_2,flat_g_2,dark_1s_2,dark_05s_2,dark_5s_2,dark_60s_2, biases_1, flat_r_1, flat_g_1, dark_5s_1]

#file to put in should be changed

fileofmean = f"analyse/master_calibration/{directory_mean2[11]}.fit"

masterfit = fits.PrimaryHDU(mastercombining(directory_list2[11]))

masterfit.header["TYPE"] = "dark" #"bias" "dark" "flat"
masterfit.header["exposure"] = 5 #60 #.5 1 
#masterfit.header["filter"] = "g" #"i" #"r" "g"

masterfit.writeto(fileofmean, overwrite=True)

