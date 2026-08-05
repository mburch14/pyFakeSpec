from xspec import *
import matplotlib.pyplot as plt
import numpy as np
import subprocess
from astropy.io import fits
import os

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

def generate_background_spectrum(background, mission, exposureTime, spec_dir, resp_dir, albedo = True, particle = True):


    #remove existing background models if they exist.
    for f in ["bkg.mod", "bkg.mod.gz"]:
        if os.path.exists(f):
            os.remove(f)

    #remove existing background spectrum if it exists.
    if os.path.exists(f"{spec_dir}/background.pha"):
        os.remove(f"{spec_dir}/background.pha")
    if os.path.exists(f"{spec_dir}/background.dat"):
        os.remove(f"{spec_dir}/background.dat")

    backgroundname = background.gen_spectrum_table(output = f'{spec_dir}/background.dat', albedo= albedo, particle= particle)

    AllData.clear()
    AllModels.clear()

    #turns the ASCII file into an xspec model. 
    subprocess.run(["flx2tab", backgroundname, "bkg", "bkg.mod"], check = True)
    bkg_model = Model("atable{bkg.mod}")

    #run xspec on the model using the response files. 
    fake = FakeitSettings(response= f"{resp_dir}/{mission.name}.rsp", exposure= str(exposureTime), fileName=f"{spec_dir}/background.pha")
    AllData.fakeit(1, fake)


def plot_background_spectrum(num_det_pixels, exposureTime, output_dir, spec_dir):

    AllData.clear()
    AllData(f"{spec_dir}/background.pha")
    # Energy bin edges
    spec = AllData(1)
    energies = np.array(spec.energies) #type: ignore
    elow = energies[:,0]
    ehigh = energies[:,1]
    dE = ehigh - elow
    energy = (elow + ehigh) / 2

    with fits.open(f"{spec_dir}/background.pha") as hdul:
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
    plt.savefig(f"{output_dir}/sim_background_spectrum.png", dpi=300, bbox_inches="tight")
    plt.close()


def gen_background(background, mission, exposureTime, num_det_pixels, output_dir, spec_dir, resp_dir):
    generate_background_spectrum(background, mission, exposureTime, spec_dir, resp_dir)
    plot_background_spectrum(num_det_pixels, exposureTime, output_dir, spec_dir)