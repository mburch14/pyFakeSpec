from xspec import *
import MissionClasses as mc
import requests
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import commentjson
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

instrument = 'SWIFTBAT'

with open("instrumentCharacteristics.json") as f:
    jsons = commentjson.load(f)
chars = jsons[instrument]

rspname = chars["rsp_name"]
num_det_pixels = chars["num_det_pixels"]
exposureTime = 2655 #seconds

#This is for our specific Cubesat
orb = mc.Orbit(chars["altitide"], chars['inclination'])
geo = mc.geometry(chars['config']) #can also input chars['config']. I did not want to do that.
mission1 = mc.Mission(instrument, chars['e_min'], chars['e_max'])
mask = mc.lead(chars['mask_thickness'])
cztDetector = mc.czt(geometry=geo, orbit=orb, mission= mission1, optics= mask, res= chars["spec_resolution"], grad=chars["spec_gradient"], low_ecut=chars{"low_ecut"})
optics = mc.lead(thickness=chars["mask_thickness"])
background = mc.BackgroundModel(detector=cztDetector)

def main():
    print(f'field of view: {geo.fov_sr} steradians')

    #Generates the .arf and the .rsp file to be used by xspec. Then, generates a ASCII file for the background spectrum.
    arfname = cztDetector.gen_arf(energy_lo = cztDetector.energy_low, energy_hi = cztDetector.energy_high, arf=chars["arf_name"])
    cztDetector.gen_rsp(arfname, rsp = rspname)

    energy_vals = np.linspace(mission1.energymin, mission1.energymax, mission1.energymax - mission1.energymin + 1)
    effective_area = [cztDetector.effective_area(energy=e) for e in energy_vals]

    plt.figure(figsize=(8,5))
    plt.plot(energy_vals, effective_area, color="black")
    plt.xlabel("Energy (keV)")
    plt.ylabel(r"Effective Area ($cm^2$)")
    plt.xscale('log')
    plt.title(f'effective area of {instrument}')
    plt.savefig("outputs/EffectiveArea.png", dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()