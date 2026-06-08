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

# %% Data inladen
data_file_i = np.asarray(fits.getdata('1_60s_m44_i.fit'), dtype=np.float32)
data_file_g = np.asarray(fits.getdata('1_60s_m44_g.fit'), dtype=np.float32)
data_file_r = np.asarray(fits.getdata('1_60s_m44_r.fit'), dtype=np.float32)

# %% Alignment naar i als referentie
aligned_g, _ = aa.register(data_file_g, data_file_i)
aligned_r, _ = aa.register(data_file_r, data_file_i)

# %% Functie om FWHM te meten en sterren te vinden
def vind_sterren(data_file, threshold):
    # Stap 1: grove sterdetectie om FWHM te schatten
    finder_grof = photutils.detection.DAOStarFinder(threshold, fwhm=8.0)
    ans_grof = finder_grof.find_stars(data_file)

    # Top 10 helderste sterren gebruiken voor FWHM-meting
    idx = np.argsort(ans_grof['flux'])[-10:]
    lijst_fwhm = []

    for n in idx:
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

    gem_fwhm = np.mean(lijst_fwhm)
    print(f"Gemiddelde FWHM: {gem_fwhm:.2f} pixels")

    # Stap 2: echte detectie met gemeten FWHM
    finder_echt = photutils.detection.DAOStarFinder(threshold, fwhm=gem_fwhm)
    temp = finder_echt.find_stars(data_file)

    print(temp.columns)
    
    with open("temp.csv", "w") as f:
        for i in temp:
            f.writelines(str(i))
    return temp

# %% Sterren zoeken ALLEEN in i (referentiefilter)
ans_ref = vind_sterren(data_file_i, threshold=1700)

# overbelichte sterren weghalen
mask = ans_ref['peak'] < 800000
ans_ref = ans_ref[mask]
print(f"Aantal sterren na saturatiefilter: {len(ans_ref)}")

# %% Aperture photometry op dezelfde posities in alle drie filters
positions = np.column_stack((ans_ref['xcentroid'], ans_ref['ycentroid'], ))

aperture = CircularAperture(positions, r=10)
annulus  = CircularAnnulus(positions, r_in=15, r_out=20)

phot_i = aperture_photometry(data_file_i, aperture)
phot_g = aperture_photometry(aligned_g,   aperture)
phot_r = aperture_photometry(aligned_r,   aperture)

# %% Achtergrond aftrekken en magnitude berekenen
for phot, data in zip([phot_i, phot_g, phot_r],
                      [data_file_i, aligned_g, aligned_r]):
    stats = ApertureStats(data, annulus)
    phot['flux_netto'] = phot['aperture_sum'] - stats.mean * aperture.area
    phot['mag'] = -2.5 * np.log10(np.abs(phot['flux_netto']))

# %% Plot: sterren op afbeelding
plt.figure(figsize=(8, 8))
plt.imshow(data_file_i, cmap='gray', norm=LogNorm())
#plt.plot(ans_ref['xcentroid'], ans_ref['ycentroid'],
         #'o', color='blue', markersize=5, label='Sterren (i-filter)')

for i, (x_val, y_val) in enumerate(zip(ans_ref['xcentroid'], ans_ref['ycentroid']), start=1):
        plt.annotate(
            str(i),              # tekst: 1, 2, 3, ...
            (x_val, y_val),      # positie van het punt
            xytext=(5, 5),       # verschuiving tov het punt
            textcoords="offset points"
            )

plt.legend()
plt.title('Gedetecteerde sterren — M44')
plt.savefig('gedetecteerd.png')

# %% CMD: g - i vs mag_i
kleur_ri = phot_r['mag'] - phot_i['mag']

plt.figure(figsize=(7, 7))
plt.scatter(kleur_ri, phot_r['mag'], s=10, color='black')

for i, (x_val, y_val) in enumerate(zip(kleur_ri, phot_r['mag']), start=1):
    plt.annotate(
    str(i),              # tekst: 1, 2, 3, ...
    (x_val, y_val),      # positie van het punt
    xytext=(5, 5),       # verschuiving t.o.v. het punt
    textcoords="offset points"
    )

plt.xlabel('r − i  (kleur)')
plt.ylabel('Instrumentele magnitude (r)')
plt.gca().invert_yaxis()

plt.savefig('check_hr_1.png')

#tabel = pd.DataFrame({
    #'ster_idx':   np.arange(1, len(ans_ref) + 1),
    #'x':         np.array(ans_ref['xcentroid']),
    #'y':         np.array(ans_ref['ycentroid']),
    #'flux_i':    np.array(phot_i['flux_netto']),
    #'flux_g':    np.array(phot_g['flux_netto']),
    #'flux_r':    np.array(phot_r['flux_netto']),
    #'mag_i':     np.array(phot_i['mag']),
    #'mag_g':     np.array(phot_g['mag']),
    ##'mag_r':     np.array(phot_r['mag']),
#})

#print(type.phot_g)
#print(tabel)
#tabel.to_csv('fov_2_sterren_gegevens.csv', index=False)
#misschien heeft het te maken met het feit dat sommige sterren niet zo goed te zien zijn in andere filters dan i