from xspec import *
import matplotlib.pyplot as plt
import numpy as np
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


def generate_observation_spectrum(mission, sourcechars, spec_dir, resp_dir):


    #remove existing observation spectrum if it exists.
    if os.path.exists(f"{spec_dir}/observation.pha"):
        os.remove(f"{spec_dir}/observation.pha")
    if os.path.exists(f"{spec_dir}/observation_bkg.pha"):
        os.remove(f"{spec_dir}/observation_bkg.pha")

    AllData.clear()
    AllModels.clear()

    #this is the source that we are doing.
    sourcemodel = Model(sourcechars["sourceShape"])
    sourcemodel.powerlaw.PhoIndex = sourcechars["phoIndex"] #type: ignore
    sourcemodel.powerlaw.norm = sourcechars["normalization"] #type: ignore

    #run xspec on the model using the response files. 
    fake = FakeitSettings(response= f"{resp_dir}/{mission.name}.rsp", exposure= sourcechars["exposure"], fileName=f"{spec_dir}/observation.pha", background = f"{spec_dir}/background.pha")
    AllData.fakeit(1, fake)


def plot_observation_spectrum(sourcename, exposureTime, output_dir, spec_dir):
    AllData.clear()
    AllData(f"{spec_dir}/observation.pha")
    spec = AllData(1)

    # Energy bin edges
    energies = np.array(spec.energies) #type: ignore
    elow = energies[:,0]
    ehigh = energies[:,1]
    dE = ehigh - elow
    energy = (elow + ehigh) / 2

    with fits.open(f"{spec_dir}/observation.pha") as hdul:
        pha = hdul["SPECTRUM"].data #type: ignore
        counts = np.array(pha["COUNTS"])
    countRate = ((counts / exposureTime) / dE) #puts the y-axis in the correct units.

    plt.figure(figsize=(8,5))
    plt.step(energy, countRate, where="mid", color="black")
    plt.xlabel("Energy (keV)")
    plt.ylabel("rate (counts/s/keV)")
    plt.yscale("log")
    plt.xscale('log')
    plt.title(f"Simulated {sourcename} Spectrum with Background")
    plt.savefig(f"{output_dir}/Sim_{sourcename}_spectrum.png", dpi=300, bbox_inches="tight")
    plt.close()
    

def calculate_snr(spec_dir):
    obs = fits.getdata(f"{spec_dir}/observation.pha", 1)["COUNTS"].sum() #type: ignore
    bkg = fits.getdata(f"{spec_dir}/observation_bkg.pha", 1)["COUNTS"].sum() #type: ignore

    print("Source counts:", obs-bkg)
    print("Background counts:", bkg)
    print('SNR:', (obs-bkg)/np.sqrt(bkg))


def gen_observation(mission, sourcechars, output_dir, spec_dir, resp_dir):
    generate_observation_spectrum(mission, sourcechars, spec_dir, resp_dir)
    plot_observation_spectrum(sourcechars["name"], sourcechars["exposure"], output_dir, spec_dir)
    calculate_snr(spec_dir)