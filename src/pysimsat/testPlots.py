import pysimsat.pysimsat.MissionClasses as mc
import matplotlib.pyplot as plt
import numpy as np
from xspec import *
import pysimsat.pysimsat.MissionClasses as mc
import requests
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import commentjson
import xraydb
import pysimsat.pysimsat.MissionClasses as MissionClasses
params = {
    "axes.labelsize": 15,
    "font.size": 15,
    "legend.fontsize": 15,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "font.family": "serif",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.top": True,
    "ytick.right": True,
    "xtick.direction": "in",
    "ytick.direction": "in",
}
plt.rcParams.update(params)

#This is for if you need to get the classes for any of the plots. You can also import from genRSP.py

'''instrument = 'ASTROSAT'

with open("instrumentCharacteristics.json") as f:
    jsons = commentjson.load(f)
chars = jsons[instrument]

rspname = chars["rsp_name"]
num_det_pixels = chars["num_det_pixels"]
exposureTime = 300 #seconds

#This is for our specific Cubesat
orb = mc.Orbit(chars["altitude"], chars['inclination'])
geo = mc.geometry(chars['config']) #can also input chars['config']. I did not want to do that.
mission1 = mc.Mission(instrument, chars['e_min'], chars['e_max'])
mask = mc.tantalum(chars['mask_thickness'])
cztDetector = mc.czt(geometry=geo, orbit=orb, mission= mission1, optics= mask, res= chars["spec_resolution"], grad=chars["spec_gradient"], low_ecut=chars["low_ecut"])
background = mc.BackgroundModel(detector=cztDetector)


print(geo.__str__())
'''

#Plot the different background components (CXB and albedo)
'''
energy = np.linspace(mission1.energymin, mission1.energymax, 500)
cxb = np.array([e*e*background.cxb(e, geo.fov_sr)/geo.fov_sr for e in energy])
albedo = np.array([e*e*background.albedo(e, geo.fov_sr)/geo.fov_sr for e in energy])

print(cxb[::5])
print("\n")
print(albedo[::5])
print("\n")

plt.figure(figsize=(8,5))
plt.scatter(energy, cxb, color="black")
plt.scatter(energy, albedo, color="red")
plt.xlabel("Energy (keV)")
plt.ylabel(r"Energy$^2$ $\times$ $\frac{dN}{dt}$   keV$^2$(photons cm$^{-2}$ s$^{-1}$ sr$^{-1}$ keV$^{-1}$)")
plt.title("CXB and albedo")
plt.savefig("outputs/sampleBKG.png", dpi=300, bbox_inches="tight")
plt.close()
'''

#See what is in your arf file
'''
with fits.open("response_files/swiftbat.arf") as hdul:
    arf = hdul["SPECRESP"].data #type: ignore

energy = (arf["ENERG_LO"] + arf["ENERG_HI"])/2
area = arf["SPECRESP"]

plt.plot(energy, area)
plt.xlabel("Energy (keV)")
plt.ylabel("Effective area (cm$^2$)")
plt.yscale("log")
plt.xscale("log")
plt.savefig("outputs/effectiveAreatest.png", dpi=300, bbox_inches="tight")
plt.show()
'''

#Plot the effective area on a given range
'''
w_Cd = 0.425
w_Zn = 0.028
w_Te = 0.547
density = 5.78 #g/cm^3

energy_vals = np.linspace(1, 1000, 990)

atten_const = np.array([(w_Cd * xraydb.mu_elam('Cd', e*1000) + w_Zn * xraydb.mu_elam('Zn', e*1000) + w_Te * xraydb.mu_elam('Te', e*1000))*density for e in energy_vals])
det_abs = (1-np.exp(-atten_const * geo.detthickness * 0.1)) #The 0.1 is to convert from mm to cm, since the thickness is in mm and the attenuation constant is in cm^-1.


effective_area = [cztDetector.effective_area(energy=e) for e in energy_vals]
optic_trans = [mask.transmission(e) for e in energy_vals]

fig, ax1 = plt.subplots(figsize=(8, 5))

# Left y-axis: Effective area
ax1.plot(energy_vals, effective_area, color="black", linewidth=2, label="Effective Area")
ax1.set_xlabel(r"Energy ($keV$)")
ax1.set_ylabel(r"Effective Area ($cm^2$)", color="black")
#ax1.set_ylim(0, 1500)
ax1.tick_params(axis='y', labelcolor='black')
ax1.set_xscale("log")

# Right y-axis: Transmission/Absorption
ax2 = ax1.twinx()
ax2.plot(energy_vals, det_abs * 100,color="tab:blue", linestyle="--", label="Detector Absorption")
ax2.plot(energy_vals, np.array(optic_trans) * 100, color="tab:red", linestyle="-.", label="Optics Transmission")
ax2.set_ylabel(r"Efficiency (%)")
ax2.set_ylim(0, 100)
ax2.tick_params(axis='y')

# Combine legends from both axes
lines = ax1.get_lines() + ax2.get_lines()
labels = [line.get_label() for line in lines]
ax1.legend(lines, labels, loc="best")

plt.title("ASTROSAT Effective Area and Efficiencies")
plt.tight_layout()
plt.savefig("outputs/EffectiveArea.png", dpi=300, bbox_inches="tight")
plt.show()

'''
#This is to plot the actual SWIFT spectrum on the range from 10-200keV
'''
with fits.open("spectrum_files/bat_total.pha") as hdul:
    spec = hdul["SPECTRUM"].data #type: ignore
    ebounds = hdul["EBOUNDS"].data#type: ignore
    rate = np.array(spec["RATE"])
    elo = np.array(ebounds["E_MIN"])
    ehi = np.array(ebounds["E_MAX"])

energy = (elo + ehi) / 2
dE = ehi - elo

rate_per_keV = rate / dE

mask = (energy >= 10) & (energy <= 200)

energy_plot = energy[mask]
rate_plot = rate_per_keV[mask]

plt.figure(figsize=(8,5))
plt.step(energy_plot, rate_plot, where="mid", color="black")
plt.xlabel("Energy (keV)")
plt.ylabel("rate (counts/s/keV)")
plt.xlim(10, 200)
plt.ylim(10**-1, 5*10**2)
plt.xscale("log")
plt.yscale("log")
plt.title("Swift/BAT Crab Spectrum with Background")
plt.savefig("outputs/real_crabbkg_spectrum.png", dpi=300, bbox_inches="tight")
plt.close()
'''

#This is code to output the SNR of the actual SWIFT data
'''
with fits.open("spectrum_files/bat_justcrab.pha") as hdul:
    spectru = hdul["SPECTRUM"]
    spec = spectru.data

    exposure = spectru.header["EXPOSURE"]
    rspfile = spectru.header["RESPFILE"]
    arffile = spectru.header["ANCRFILE"]

rate = spec["RATE"]
source_err = spec["STAT_ERR"]
err = spec["STAT_ERR"]
snr = spec["SNR"]

snr = np.sum(rate) / np.sqrt(np.sum(err**2))

print("SNR: ",snr)
#print(spec.columns)
#print(rate.shape)
print("Rate: ", np.sum(rate))
print("stdv", np.sqrt(np.sum(err**2)))
print("exposure: ",exposure)
print("counts: " ,np.sum(rate)*exposure)
'''

'''
hdu = fits.open("spectrum_files/sw01132968000bevpb_sk.img")
print(hdu[0].header["OBJECT"])
print(hdu[0].header["RA_OBJ"])
print(hdu[0].header["DEC_OBJ"])

subprocess.run([
    "batcelldetect",
    "infile=spectrum_files/sw01132968000bevpb_sk.img",
    "outfile=crab_catalog.fits"
])
'''


'''
hdul = fits.open("crab_catalog.fits")

cat = hdul["BAT_CATALOG"].data

for row in cat:
    print(
        row["NAME"],
        row["RA_OBJ"],
        row["DEC_OBJ"],
        row["SNR"],
        row["COUNTS"],
        row["EXPOSURE"]
    )
print(np.sum(cat["COUNTS"]))
'''