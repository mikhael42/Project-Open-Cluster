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

#hiermee haal je de data uit de fits. Het moet met np.asarray om ze later te kunnen alignen 
data_file_i = np.asarray(fits.getdata('2_60s_m44_i.fit'), dtype=np.float32)
data_file_g = np.asarray(fits.getdata('2_60s_m44_g.fit'), dtype=np.float32)
data_file_r = np.asarray(fits.getdata('2_60s_m44_r.fit'), dtype=np.float32)

#vervolgens alignen we de filters g en r naar i. Zo zorgen we ervoor dat sterren in alle drie de filter dezelfde coordinaten hebben
aligned_g, _ = aa.register(data_file_g, data_file_i)
aligned_r, _ = aa.register(data_file_r, data_file_i)

#Om vervolgens sterren te kunnen vinden moeten we eerst nog fwhm bepalen 
# dit voegen we samen tot een def-functie 


def vind_sterren(data_file, threshold): #input: data_file en threshold(simpelweg de achtergrond intensiteit)
    # Stap 1: fwhm bepalen
    finder_grof = photutils.detection.DAOStarFinder(threshold, fwhm=8.0)
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
    finder_echt = photutils.detection.DAOStarFinder(threshold, fwhm=gem_fwhm) 
    #vervolgens maken we een tabel met de gegevens per ster. 
    temp = finder_echt.find_stars(data_file)

    #vervolgens slaan we deze tabel op als een csv bestand (verander de naam van het bestand)
    with open("info_sterren_i_fov", "w") as f:
        for i in temp:
            f.writelines(str(i))
    return temp

#om nu verder te kunnen gaan we eerst alle sterren detecteren in i-filter
#we gebruiken de def met de volgende input: data file van i-filter, achtergrond intensiteit (te vinden met behulp van DS9)
ans_ref = vind_sterren(data_file_i, threshold=1700)
pd.DataFrame(np.array(ans_ref)).to_csv('info_sterren_i_fov', index=False)
#ons output is dus een tabels - "temp", met de gegevans per ster (id, x, y, sharpness, roundness1, roundness2, n_pixels, peak, flux, mag, daofind_mag)

# overbelichte sterren weghalen
#bij deze stap moeten we nog met alkaar afstemmen wat de mask-voorwaarde moet zijn
mask = ans_ref['peak'] < 800000
ans_ref = ans_ref[mask]
print(f"Aantal sterren na saturatiefilter: {len(ans_ref)}") #vervolgens printen we de totale hoeveelheid gedetecteerde sterren

#nu weten we de x- en de y-coordinaten van de sterren in i-filter
#de filters zijn de aligned en dus hebben de sterren precies dezelfde coordinaten in alle drie de filters 
#we maken een 'lijst' met alle coordinaten van de gedetecteerde sterren 
positions = np.column_stack((ans_ref['xcentroid'], ans_ref['ycentroid'], ))

#vervolgens gaan we de fluxen bepalen van de desbetreffende coordinaten, maar dan in g- en r-filter 
aperture = CircularAperture(positions, r=10)
annulus  = CircularAnnulus(positions, r_in=15, r_out=20)

#fluxen berkenenen gebeurt hier 
phot_i = aperture_photometry(data_file_i, aperture)
phot_g = aperture_photometry(aligned_g,   aperture)
phot_r = aperture_photometry(aligned_r,   aperture)

#vervolgens bereken we de magnitudes per ster per filter 
for phot, data in zip([phot_i, phot_g, phot_r],
                      [data_file_i, aligned_g, aligned_r]):
    stats = ApertureStats(data, annulus)
    phot['flux_netto'] = phot['aperture_sum'] - stats.mean * aperture.area
    phot['mag'] = -2.5 * np.log10(np.abs(phot['flux_netto']))

#we maken een overzicht waarbij wij de sterren in een paaltje visueel weergeven en benummeren 
plt.figure(figsize=(8, 8))
plt.imshow(data_file_i, cmap='gray', norm=LogNorm())
plt.plot(ans_ref['xcentroid'], ans_ref['ycentroid'],
         'o', color='blue', markersize=5, label='Sterren (i-filter)')

for i, (x_val, y_val) in enumerate(zip(ans_ref['xcentroid'], ans_ref['ycentroid']), start=1):
        plt.annotate(
            str(i),              # tekst: 1, 2, 3, ...
            (x_val, y_val),      # positie van het punt
            xytext=(5, 5),       # verschuiving tov het punt
            textcoords="offset points"
            )

plt.legend()
plt.title('Gedetecteerde sterren — M44')
#we slaan dit afbeelding op
plt.savefig('gedetecteerde_sterren_i_fov2.png')

#vervolgens stellen we een HR-diagram op 
#We hebben een verschil nodig tussen twee filters
kleur_ri = phot_r['mag'] - phot_i['mag']

plt.figure(figsize=(7, 7))
plt.scatter(kleur_ri, phot_r['mag'], s=10, color='black')

#vervolgens nummeren we ook in ons HR de sterren 
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

plt.savefig('voorlopig_hr_fov2.png')

#uiteindelijk maken we nog een tabel met PANDAS waarin we per ster alle magnitudes weergeven per filter 
#deze zetten we vervolgens op GitHub zodat Teo ze in de totale HR verwerkt 

tabel_gegevens_per_ster = pd.DataFrame({
    'ster_idx':   np.arange(1, len(ans_ref) + 1),
    'x':         np.array(ans_ref['xcentroid']),
    'y':         np.array(ans_ref['ycentroid']),
    'flux_i':    np.array(phot_i['flux_netto']),
    'flux_g':    np.array(phot_g['flux_netto']),
    'flux_r':    np.array(phot_r['flux_netto']),
    'mag_i':     np.array(phot_i['mag']),
    'mag_g':     np.array(phot_g['mag']),
    'mag_r':     np.array(phot_r['mag']),
})
tabel_gegevens_per_ster.to_csv('sterren_gegevens_fov2.csv', index=False)

#check jouw output
#1. waarde van gem. fwhm (moet ergens tussen 3 en 12) - als het meer is dan 20 -> waarschijnlijk een fout -> check jouw code 
#2. csv-bestand info_sterren_i_fov
#3. afbeelding met gedetecteerde sterren - "gedetecteerde_sterren_i_fov2.png"
#4. hr-diagram - "voorlopig_hr_fov2.png"
#5. csv-bestand sterren_gegevens_fov2.csv