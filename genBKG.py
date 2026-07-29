from xspec import *
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import genRSP
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
background = genRSP.background
rspname = genRSP.rspname
num_det_pixels = genRSP.num_det_pixels

backgroundname = background.gen_spectrum_table(output = 'spectrum_files/background.dat', albedo= False, particle= False)

AllData.clear()
AllModels.clear()

#turns the ASCII file into an xspec model. 
subprocess.run(["flx2tab", "spectrum_files/background.dat", "bkg", "bkg.mod"])
background = Model("atable{bkg.mod}")

#run xspec on the model using the response files. 
fake = FakeitSettings(response= rspname, exposure= str(exposureTime), fileName="spectrum_files/background.pha")
AllData.fakeit(1, fake)

AllData.clear()
AllData("spectrum_files/background.pha")

# Energy bin edges
spec = AllData(1)
energies = np.array(spec.energies) #type: ignore
elow = energies[:,0]
ehigh = energies[:,1]
dE = ehigh - elow
energy = (elow + ehigh) / 2

with fits.open("spectrum_files/background.pha") as hdul:
    pha = hdul["SPECTRUM"].data #type: ignore
    counts = np.array(pha["COUNTS"])
countRate = ((counts / exposureTime) / dE)/num_det_pixels #puts the y-axis in the correct units.

plt.figure(figsize=(8,5))
plt.step(energy, countRate, where="mid", color="black")
plt.xlabel("Energy (keV)")
plt.ylabel("Count rate (counts/s/keV/detector)")
plt.yscale("log")
plt.xscale('log')
plt.title("Simulated background Spectrum")
plt.savefig("outputs/sim_background_spectrum.png", dpi=300, bbox_inches="tight")
plt.close()