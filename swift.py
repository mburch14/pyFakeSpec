from xspec import *
import MissionClasses as mc
import requests
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import commentjson
import genRSPfile

exposureTime = genRSPfile.exposureTime
rspname = genRSPfile.rspname

AllData.clear()
AllModels.clear()

#this is the source that we are doing. This is not right. We are not considering the field of view of the simulation. We need to fix this. It should be the exact same as the background.
crab = Model("powerlaw")
crab.powerlaw.PhoIndex = 2.15 #type: ignore
crab.powerlaw.norm = 10.17 #type: ignore

#run xspec on the model using the response files. 
fake = FakeitSettings(response= rspname, exposure= exposureTime, fileName="source.pha")  #type: ignore
AllData.fakeit(1, fake)

AllData("source.pha")
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
plt.title("Simulated Crab Spectrum")
plt.savefig("Crab_spectrum.png", dpi=300, bbox_inches="tight")
plt.close()