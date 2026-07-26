from astropy.io import fits
import numpy as np

obs = fits.getdata("spectrum_files/observation.pha", 1)["COUNTS"].sum() #type: ignore
bkg = fits.getdata("spectrum_files/observation_bkg.pha", 1)["COUNTS"].sum() #type:ignore

print("Observation counts:", obs)
print("Background counts:", bkg)

print('SNR:', obs/np.sqrt(bkg))


