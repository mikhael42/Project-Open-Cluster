from astropy.io import fits
import matplotlib.pyplot as plt
import astroalign as aa
from matplotlib.colors import LogNorm
import numpy as np
import photutils
import photutils.detection
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry, ApertureStats
from astropy.modeling import models, fitting
import pandas as pd
from lmfit import models as md

#hiermee haal je de data uit de fits. Het moet met np.asarray om ze later te kunnen alignen 
data_file_i_1 = np.asarray(fits.getdata('1_60s_m44_i.fit'), dtype=np.float32)
data_file_g_1 = np.asarray(fits.getdata('1_60s_m44_g.fit'), dtype=np.float32)
data_file_r_1 = np.asarray(fits.getdata('1_60s_m44_r.fit'), dtype=np.float32)

data_file_i_2 = np.asarray(fits.getdata('2_60s_m44_i.fit'), dtype=np.float32)
data_file_g_2 = np.asarray(fits.getdata('2_60s_m44_g.fit'), dtype=np.float32)
data_file_r_2 = np.asarray(fits.getdata('2_60s_m44_r.fit'), dtype=np.float32)

#vervolgens alignen we de filters g en r naar i. Zo zorgen we ervoor dat sterren in alle drie de filter dezelfde coordinaten hebben
aligned_g_1, _ = aa.register(data_file_g_1, data_file_i_1)
aligned_r_1, _ = aa.register(data_file_r_1, data_file_i_1)

aligned_g_2, _ = aa.register(data_file_g_2, data_file_i_2)
aligned_r_2, _ = aa.register(data_file_r_2, data_file_i_2)

#Om vervolgens sterren te kunnen vinden moeten we eerst nog fwhm bepalen 
# dit voegen we samen tot een def-functie 


def vind_sterren(data_file, threshold): #input: data_file en threshold(simpelweg de achtergrond intensiteit)
    # Stap 1: fwhm bepalen
    finder_grof = photutils.detection.DAOStarFinder(threshold, fwhm=8.0, roundness_range=(-0.5, 0.5), sharpness_range=(0.3, 0.8))
    ans_grof = finder_grof.find_stars(data_file)

    #om onze fit niet te uitgebreid te maken gebruiken we alleen 10 helderste sterren 
    idx = np.argsort(ans_grof['flux'])[-10:] #hiermee maak je een tabel met gegevens van 10 helderste sterren 
    lijst_fwhm = []

    for n in idx: #vervolgens berekenen we per ster wat bijbehorende fwhm is 
        x0 = int(ans_grof['xcentroid'][n])
        y0 = int(ans_grof['ycentroid'][n])

        cutout = data_file[y0-15:y0+16, x0-15:x0+16]
        y, x = np.mgrid[:cutout.shape[0], :cutout.shape[1]]

        g_init = models.Gaussian2D(
            amplitude=cutout.max(),
            x_mean=cutout.shape[1] / 2,
            y_mean=cutout.shape[0] / 2
        )
        fit_g = fitting.LevMarLSQFitter()
        g = fit_g(g_init, x, y, cutout)

        fwhm = 2.355 * (g.x_stddev.value + g.y_stddev.value) / 2
        lijst_fwhm.append(fwhm) 
    #als output hebben we nu een lijst met verschillende fwhm's 

    gem_fwhm = np.mean(lijst_fwhm) #we berekenen de gemmidelde fwhm
    print(f"Gemiddelde fwhm: {gem_fwhm:.2f} pixels")
    #dit wordt onze fwhm bij het vinden van sterren (DAO Star Finder)

    # Stap 2: detecteren van sterren 
    #eerst vinden we alle sterren
    finder_echt = photutils.detection.DAOStarFinder(threshold, fwhm=gem_fwhm, roundness_range=(-0.5, 0.5), sharpness_range=(0.3, 0.8)) 
    #vervolgens maken we een tabel met de gegevens per ster. 
    temp = finder_echt.find_stars(data_file)

    #vervolgens slaan we deze tabel op als een csv bestand (verander de naam van het bestand)
    with open("info_sterren_i_fov1", "w") as f:
        for i in temp:
            f.writelines(str(i))
    return temp

#om nu verder te kunnen gaan we eerst alle sterren detecteren in i-filter
#we gebruiken de def met de volgende input: data file van i-filter, achtergrond intensiteit (te vinden met behulp van DS9)
ans_ref_1 = vind_sterren(data_file_i_1, threshold=1000)

ans_ref_2 = vind_sterren(data_file_i_2, threshold=1700)

#ons output is dus een tabels - "temp", met de gegevans per ster (id, x, y, sharpness, roundness1, roundness2, n_pixels, peak, flux, mag, daofind_mag)

# overbelichte sterren weghalen
#bij deze stap moeten we nog met alkaar afstemmen wat de mask-voorwaarde moet zijn
mask_1 = ans_ref_1['peak'] < 700000
mask_2 = ans_ref_2['peak'] < 700000
ans_ref_1 = ans_ref_1[mask_1]
ans_ref_2 = ans_ref_2[mask_2]
print(f"Aantal sterren na saturatiefilter FoV1: {len(ans_ref_1)}") #vervolgens printen we de totale hoeveelheid gedetecteerde sterren
print(f"Aantal sterren na saturatiefilter FoV2: {len(ans_ref_2)}") 
#nu weten we de x- en de y-coordinaten van de sterren in i-filter
#de filters zijn de aligned en dus hebben de sterren precies dezelfde coordinaten in alle drie de filters 
#we maken een 'lijst' met alle coordinaten van de gedetecteerde sterren 
positions_1 = np.column_stack((ans_ref_1['xcentroid'], ans_ref_1['ycentroid'], ))
positions_2 = np.column_stack((ans_ref_2['xcentroid'], ans_ref_2['ycentroid'], ))
#vervolgens gaan we de fluxen bepalen van de desbetreffende coordinaten, maar dan in g- en r-filter 
aperture_1 = CircularAperture(positions_1, r=10)
annulus_1  = CircularAnnulus(positions_1, r_in=15, r_out=20)

aperture_2 = CircularAperture(positions_2, r=10)
annulus_2  = CircularAnnulus(positions_2, r_in=15, r_out=20)

#fluxen berkenenen gebeurt hier 
phot_i_1 = aperture_photometry(data_file_i_1, aperture_1)
phot_g_1 = aperture_photometry(aligned_g_1,   aperture_1)
phot_r_1 = aperture_photometry(aligned_r_1,   aperture_1)

phot_i_2 = aperture_photometry(data_file_i_2, aperture_2)
phot_g_2 = aperture_photometry(aligned_g_2,   aperture_2)
phot_r_2 = aperture_photometry(aligned_r_2,   aperture_2)

#vervolgens bereken we de magnitudes per ster per filter 
for phot_1, data in zip([phot_i_1, phot_g_1, phot_r_1],
                      [data_file_i_1, aligned_g_1, aligned_r_1]):
    stats = ApertureStats(data, annulus_1)
    phot_1['flux_netto'] = phot_1['aperture_sum'] - stats.mean * aperture_1.area
    phot_1['mag'] = -2.5 * np.log10(np.abs(phot_1['flux_netto']))

for phot_2, data in zip([phot_i_2, phot_g_2, phot_r_2],
                      [data_file_i_2, aligned_g_2, aligned_r_2]):
    stats = ApertureStats(data, annulus_2)
    phot_2['flux_netto'] = phot_2['aperture_sum'] - stats.mean * aperture_2.area
    phot_2['mag'] = -2.5 * np.log10(np.abs(phot_2['flux_netto']))


lijst_mag_i_1 = np.array(phot_i_1['mag'])
lijst_mag_g_1 = np.array(phot_g_1['mag'])
lijst_mag_r_1 = np.array(phot_r_1['mag'])

lijst_mag_i_2 = np.array(phot_i_2['mag'])
lijst_mag_g_2 = np.array(phot_g_2['mag'])
lijst_mag_r_2 = np.array(phot_r_2['mag'])

ster_nrs_1 = np.array([82, 157, 126, 147, 56, 145, 209, 8, 60, 29]) - 1
ster_nrs_2 = np.array([70, 66, 52, 81, 69, 89, 93, 112, 101, 87])-1

inst_i = np.concatenate([np.array(lijst_mag_i_1[ster_nrs_1]), np.array(lijst_mag_i_2[ster_nrs_2])])
inst_g= np.concatenate([np.array(lijst_mag_g_1[ster_nrs_1]), np.array(lijst_mag_g_2[ster_nrs_2])])
inst_r = np.concatenate([np.array(lijst_mag_r_1[ster_nrs_1]), np.array(lijst_mag_r_2[ster_nrs_2])])

echte_mag_i_1 = np.array([14.301, 11.07, 11.0, 10.48, 10.042, 14.488, 10.01, 16.594, 14.238, 18.162])
echte_mag_g_1 = np.array([15.228, 11.61, 11.54, 11.07, 10.756, 14.608, 12.67, 17.468, 12.75, 17.624])
echte_mag_r_1 = np.array([14.341, 11.22, 11.15, 10.65, 10.118, 14.455, 10.07, 16.922, 14.035, 17.812])

echte_mag_i_2 = np.array([13.184, 13.99, 10.591, 12.03, 14.044, 13.29, 11.673, 14.491, 9.672, 14.6])
echte_mag_g_2 = np.array([13.643, 15.82, 11.609, 13.13, 14.103, 13.962, 12.779, 15.609, 10.014, 16.7])
echte_mag_r_2 = np.array([13.393, 14.6, 13.556, 12.36, 14.714, 13.504, 14.029, 14.746, 9.713, 15.49])

echte_mag_i = np.concatenate([echte_mag_i_1, echte_mag_i_2])
echte_mag_g = np.concatenate([echte_mag_g_1, echte_mag_g_2])
echte_mag_r = np.concatenate([echte_mag_r_1, echte_mag_r_2])

def fit_function(x, calibratie_constante):
    return x + calibratie_constante

def calibratie_drie_filters(in_inst_i, in_inst_g, in_inst_r, in_echte_mag_i, in_echte_mag_g, in_echte_mag_r):
    for inst, catalogus, naam in zip([in_inst_i, in_inst_g, in_inst_r], [in_echte_mag_i, in_echte_mag_g, in_echte_mag_r], ['i', 'g', 'r']):
        model = md.Model(fit_function, independent_vars=['x'])
        params = model.make_params(calibratie_constante = 26.57)

        result = model.fit(catalogus, x=inst, params=params)

        calibratie_constante = result.params['calibratie_constante'].value
        onzekerheid = result.params['calibratie_constante'].stderr

        print(f"\nFilter {naam}: calibratie constante = {calibratie_constante:.4f} +/- {onzekerheid:.4f}")

#70, 66, 52, 81, 69                                  89, 93, 112, 101, 87 (2)
#82, 157, 126, 147,                                                                         56, 145, 209, 8, 60, 29 (1)
lijst_inst_mag_i_2 = [phot_i_2['mag'][69], phot_i_2['mag'][65], phot_i_2['mag'][51], phot_i_2['mag'][80], phot_i_2['mag'][68], phot_i_2['mag'][88], phot_i_2['mag'][92], phot_i_2['mag'][111], phot_i_2['mag'][100], phot_i_2['mag'][86]]
lijst_inst_mag_i_1 = [phot_i_1['mag'][81], phot_i_1['mag'][156], phot_i_1['mag'][125], phot_i_1['mag'][146], phot_i_1['mag'][55], phot_i_1['mag'][144], phot_i_1['mag'][208], phot_i_1['mag'][7], phot_i_1['mag'][59], phot_i_1['mag'][28]]

lijst_inst_mag_g_2 = [phot_g_2['mag'][69], phot_g_2['mag'][65], phot_g_2['mag'][51], phot_g_2['mag'][80], phot_g_2['mag'][68], phot_g_2['mag'][88], phot_g_2['mag'][92], phot_g_2['mag'][111], phot_g_2['mag'][100], phot_g_2['mag'][86]]
lijst_inst_mag_g_1 = [phot_g_1['mag'][81], phot_g_1['mag'][156], phot_g_1['mag'][125], phot_g_1['mag'][146], phot_g_1['mag'][55], phot_g_1['mag'][144], phot_g_1['mag'][208], phot_g_1['mag'][7], phot_g_1['mag'][59], phot_g_1['mag'][28]]

lijst_inst_mag_r_2 = [phot_r_2['mag'][69], phot_r_2['mag'][65], phot_r_2['mag'][51], phot_r_2['mag'][80], phot_r_2['mag'][68], phot_r_2['mag'][88], phot_r_2['mag'][92], phot_r_2['mag'][111], phot_r_2['mag'][100], phot_r_2['mag'][86]]
lijst_inst_mag_r_1 = [phot_r_1['mag'][81], phot_r_1['mag'][156], phot_r_1['mag'][125], phot_r_1['mag'][146], phot_r_1['mag'][55], phot_r_1['mag'][144], phot_r_1['mag'][208], phot_r_1['mag'][7], phot_r_1['mag'][59], phot_r_1['mag'][28]]

#calibratie_drie_filters(inst_i, inst_g, inst_r, echte_mag_i, echte_mag_g, echte_mag_r) #twee fov samen

#calibratie_drie_filters(lijst_mag_i_1[ster_nrs_1],lijst_mag_g_1[ster_nrs_1], lijst_mag_r_1[ster_nrs_1], echte_mag_i_1, echte_mag_g_1, echte_mag_r_1) #fov_1

#calibratie_drie_filters(lijst_mag_i_2[ster_nrs_2],lijst_mag_g_2[ster_nrs_2], lijst_mag_r_2[ster_nrs_2], echte_mag_i_2, echte_mag_g_2, echte_mag_r_2) #fov_2



plt.plot(lijst_inst_mag_i_2, echte_mag_r_2, 'o', color='red')
plt.plot(lijst_inst_mag_i_1, echte_mag_r_2, 'o', color='black')
plt.xlabel('instrumentale mag (r-filter)')
plt.ylabel('echte mag (r-filter)')
plt.savefig('fit_twee_fovs_r_10.png')
