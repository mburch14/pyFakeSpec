from xspec import *
import MissionClasses as mc
import requests
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import json
instrument = 'SWIFTBAT'

with open("instrumentCharacteristics.json", "r") as file:
    jsons = json.load(file)
chars = jsons[instrument]

#This is for our specific Cubesat
orb = mc.Orbit(chars['altitude'], chars['inclination'])
geo = mc.geometry() #can also input chars['config']. I did not want to do that.
mission1 = mc.Mission(instrument, chars['e_min'], chars['e_max'])
cztDetector = mc.czt(geo, orb, mission1, res = chars["spec_resolution"], grad=chars["spec_gradient"])
background = mc.CXB(detector=cztDetector)
exposureTime = 10000 #seconds

AllData.clear()
AllModels.clear()

'''
arfname = cztDetector.gen_arf(energy_lo = cztDetector.energy_low, energy_hi = cztDetector.energy_high, arf=chars["arf_name"])
rspname =cztDetector.gen_rsp(arfname, rsp = chars["rsp_name"])
backgroundname = background.gen_spectrum_table()

#turns the ASCII file into an xspec model. 
subprocess.run(["flx2tab", "cxb_background.dat", "cxb", "cxb.mod"])
m = Model("atable{cxb.mod}")

#run xspec on the model using the response files. 
fake = FakeitSettings(response= 'cubesat.rsp', exposure= exposureTime, fileName="simulated.pha")
AllData.fakeit(1, fake)

AllData("simulated.pha")
Plot.xAxis = "keV"
Plot("data")
energy = Plot.x()
counts = Plot.y()
errors = Plot.yErr()
numcounts = np.array(counts) / exposureTime
plt.figure(figsize=(8,5))
plt.step(energy, numcounts, where="mid", color="black")
plt.errorbar(energy, numcounts, yerr=np.array(errors)/exposureTime, fmt="none", color="black", alpha=0.5)
plt.xlabel("Energy (keV)")
plt.ylabel("Count rate (counts/s/bin)")
plt.yscale("log")
plt.xscale('log')
plt.title("Simulated CXB Spectrum")
plt.savefig("simulated_spectrum.png", dpi=300, bbox_inches="tight")
plt.close()
'''

'''
#playing around and graphing the x-ray background.
CXB_model = CXB(mission=mission1)
Albedo_model = Albedo(mission=mission1)
energies = np.logspace(.5, 2.5, 300)  # Energy range from ~3 keV to ~300 keV
plt.figure(figsize=(10, 6))
plt.plot(energies, [geo.fov_sr * CXB_model.photonIntensity(energy) for energy in energies], color='green', label='CXB Background')
plt.plot(energies, [geo.fov_sr * Albedo_model.photonIntensity(energy) for energy in energies], color='blue', label='Albedo Background')
plt.plot(energies, [geo.fov_sr * (CXB_model.photonIntensity(energy) + Albedo_model.photonIntensity(energy)) for energy in energies], color='red', label='Total Background')
plt.xlabel('Energy (keV)')
plt.ylabel(r"$flux$ ( photons keV$^{-1}$ cm$^{-2}$ s$^{-1}$)")
plt.title('Background vs Energy')
plt.xscale('log')
plt.yscale('log')
plt.xlim(1, 1000)
plt.legend()
plt.savefig("cxb.png", dpi=300)
'''

#This is for testing the ARF file from the HEASARC calibration database.
'''
arf = fits.open("cubesat.arf")

#arf.info() #This is just to check the contents of the ARF file.
aeff = arf[1].data["SPECRESP"]
#print(aeff[-5:]) #Check the last 5 values of the effective area.

plt.plot(arf[1].data["ENERG_LO"], aeff, color='purple', marker='o')
plt.xlabel('Energy (keV)')
plt.ylabel('Effective Area (cm²)')
plt.xscale('log')
plt.title('Effective Area vs Energy for CdZnTe Detector')
plt.savefig("effective_area.png", dpi=300)
'''

#This is for testing the effective area of the CdZnTe detector.
'''
det = czt(geo, orb, mission1)
energies = np.linspace(10e3, 500e3, 100)
area = det.effective_area(energies)
plt.plot(energies/1e3, area, color='purple', marker='o')
plt.xlabel('Energy (keV)')
plt.ylabel('Effective Area (cm²)')
plt.xscale('log')
plt.title('Effective Area vs Energy for CdZnTe Detector')
plt.savefig("effective_area.png", dpi=300)
'''

#This is for testing the RMF file from the HEASARC calibration database.
'''
url = "https://swift.gsfc.nasa.gov/proposals/bat.rsp"
r = requests.get(url)
r.raise_for_status()
with open("bat.rsp", "wb") as f:
    f.write(r.content)
rmf = fits.open("bat.rsp")
matrix = rmf["EBOUNDS"].data
print(matrix[:5])
'''
