from xspec import *
import matplotlib.pyplot as plt
import numpy as np
import genRSP as genRSP
from astropy.io import fits

exposureTime = genRSP.exposureTime
rspname = genRSP.rspname
num_det_pixels = genRSP.num_det_pixels

AllData.clear()
AllModels.clear()

#this is the source that we are doing. This is not right. We are not considering the field of view of the simulation. We need to fix this. It should be the exact same as the background.
crab = Model("powerlaw")
crab.powerlaw.PhoIndex = 2.15 #type: ignore
crab.powerlaw.norm = 10.17 #type: ignore

#run xspec on the model using the response files. 
fake = FakeitSettings(response= rspname, exposure= exposureTime, fileName="spectrum_files/observation.pha", background = 'spectrum_files/background.pha')  #type: ignore
AllData.fakeit(1, fake)

AllData.clear()
AllData("spectrum_files/observation.pha")
spec = AllData(1)

# Energy bin edges
energies = np.array(spec.energies) #type: ignore
elow = energies[:,0]
ehigh = energies[:,1]
dE = ehigh - elow
energy = (elow + ehigh) / 2

with fits.open("spectrum_files/observation.pha") as hdul:
    pha = hdul["SPECTRUM"].data #type: ignore
    counts = np.array(pha["COUNTS"])
countRate = ((counts / exposureTime) / dE) /num_det_pixels #puts the y-axis in the correct units.

plt.figure(figsize=(8,5))
plt.step(energy, countRate, where="mid", color="black")
plt.xlabel("Energy (keV)")
plt.ylabel("Count rate (counts/s/keV)")
plt.yscale("log")
plt.xscale('log')
plt.title("Simulated Crab Spectrum")
plt.savefig("outputs/Crab_spectrum.png", dpi=300, bbox_inches="tight")
plt.close()