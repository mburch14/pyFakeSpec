import MissionClasses as mc
import matplotlib.pyplot as plt
import numpy as np
from xspec import *
import MissionClasses as mc
import requests
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import commentjson
import genRSP
'''
instrument = 'SWIFTBAT'

with open("instrumentCharacteristics.json") as f:
    jsons = commentjson.load(f)
chars = jsons[instrument]

rspname = chars["rsp_name"]
num_det_pixels = chars["num_det_pixels"]
exposureTime = 300 #seconds

#This is for our specific Cubesat
orb = mc.Orbit(chars["altitide"], chars['inclination'])
geo = mc.geometry(chars['config']) #can also input chars['config']. I did not want to do that.
mission1 = mc.Mission(instrument, chars['e_min'], chars['e_max'])
mask = mc.lead(chars['mask_thickness'])
cztDetector = mc.czt(geometry=geo, orbit=orb, mission= mission1, optics= mask, res= chars["spec_resolution"], grad=chars["spec_gradient"])
optics = mc.lead(thickness=chars["mask_thickness"])
background = mc.BackgroundModel(detector=cztDetector)
'''
'''energy = np.linspace(mission1.energymin, mission1.energymax, 500)
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
plt.close()'''

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
'''with fits.open("spectrum_files/background.pha") as hdul:
    data = hdul["SPECTRUM"].data #type: ignore
    print(data["COUNTS"][:5])

AllData.clear()
AllModels.clear()

AllData("spectrum_files/background.pha")
spec = AllData(1)

# Energy bin edges
energies = np.array(spec.energies) #type: ignore
elow = energies[:,0]
ehigh = energies[:,1]
dE = ehigh - elow
energy = (elow + ehigh) / 2
Plot.xAxis = "keV"
Plot("data")

counts = np.array(Plot.y())
countRate = ((counts / 300) / dE) / 32768 #puts the y-axis in the correct units.

print(counts[:5])
print(dE[:5])
print(countRate[:5])
print(energies[:5])'''
'''
energy_vals = np.linspace(10, 1000, 990)
effective_area = [genRSP.cztDetector.effective_area(energy=e) for e in energy_vals]

plt.figure(figsize=(8,5))
plt.plot(energy_vals, effective_area, color="black")
plt.xlabel("Energy (keV)")
plt.ylabel(r"Effective Area ($cm^2$)")
plt.xscale('log')
plt.title(f'effective area of BAT')
plt.savefig("outputs/EffectiveArea.png", dpi=300, bbox_inches="tight")
plt.close()'''

'''with fits.open("spectrum_files/bat_justcrab.pha") as hdul:
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
plt.ylabel("Count rate (counts/s/keV)")
plt.xlim(10, 200)
plt.xscale("log")
plt.yscale("log")
plt.title("Swift/BAT Crab Spectrum without Background")
plt.savefig("outputs/real_justcrab_spectrum.png", dpi=300, bbox_inches="tight")
plt.close()'''  
with fits.open("spectrum_files/bat_justcrab.pha") as hdul:
    spec = hdul["SPECTRUM"].data

rate = np.array(spec["RATE"])
err = np.array(spec["STAT_ERR"])

snr_bins = rate / err

source = np.sum(rate)
error = np.sqrt(np.sum(err**2))

snr = source / error

print("BAT SNR =", snr)