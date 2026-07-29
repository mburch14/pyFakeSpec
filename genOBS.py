from xspec import *
import matplotlib.pyplot as plt
import numpy as np
import genRSP as genRSP
from astropy.io import fits

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

exposureTime = genRSP.exposureTime
rspname = genRSP.rspname
num_det_pixels = genRSP.num_det_pixels

AllData.clear()
AllModels.clear()

#this is the source that we are doing.
crab = Model("powerlaw")
crab.PhoIndex = 2.15
crab.norm = 10.17

#run xspec on the model using the response files. 
fake = FakeitSettings(response= rspname, exposure= str(exposureTime), fileName="spectrum_files/observation.pha", background = 'spectrum_files/background.pha')
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
countRate = ((counts / exposureTime) / dE) #puts the y-axis in the correct units.

plt.figure(figsize=(8,5))
plt.step(energy, countRate, where="mid", color="black")
plt.xlabel("Energy (keV)")
plt.ylabel("rate (counts/s/keV)")
plt.yscale("log")
plt.xscale('log')
plt.title("Simulated Crab Spectrum with Background")
plt.savefig("outputs/Sim_Crab_spectrum.png", dpi=300, bbox_inches="tight")
plt.close()