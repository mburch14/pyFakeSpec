from xspec import *
import MissionClasses as mc
import requests
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import genRSPfile

exposureTime = genRSPfile.exposureTime
background = genRSPfile.background
rspname = genRSPfile.rspname

backgroundname = background.gen_spectrum_table()

#turns the ASCII file into an xspec model. 
subprocess.run(["flx2tab", "background.dat", "bkg", "bkg.mod"])
background = Model("atable{bkg.mod}")

#run xspec on the model using the response files. 
fake = FakeitSettings(response= rspname, exposure= str(exposureTime), fileName="background.pha")
AllData.fakeit(1, fake)

AllData("background.pha")
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
countRate = (counts / exposureTime) / dE #puts the y-axis in the correct units.
errors = Plot.yErr()

plt.figure(figsize=(8,5))
plt.step(energy, countRate, where="mid", color="black")
plt.errorbar(energy, countRate, yerr=np.array(errors)/exposureTime/dE, fmt="none", color="black", alpha=0.5)
plt.xlabel("Energy (keV)")
plt.ylabel("Count rate (counts/s/keV)")
plt.yscale("log")
plt.xscale('log')
plt.title("Simulated background Spectrum")
plt.savefig("background_spectrum.png", dpi=300, bbox_inches="tight")
plt.close()