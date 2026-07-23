from xspec import *
import MissionClasses as mc
import requests
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import commentjson

instrument = 'SWIFTBAT'

with open("instrumentCharacteristics.json") as f:
    jsons = commentjson.load(f)
chars = jsons[instrument]

rspname = chars["rsp_name"]
exposureTime = 10000 #seconds

#This is for our specific Cubesat
orb = mc.Orbit(chars["altitide"], chars['inclination'])
geo = mc.geometry() #can also input chars['config']. I did not want to do that.
mission1 = mc.Mission(instrument, chars['e_min'], chars['e_max'])
mask = mc.lead(chars['mask_thickness'])
cztDetector = mc.czt(geometry=geo, orbit=orb, mission= mission1, optics= mask, res= chars["spec_resolution"], grad=chars["spec_gradient"])
optics = mc.lead(thickness=chars["mask_thickness"])
background = mc.CXB(detector=cztDetector)

#Generates the .arf and the .rsp file to be used by xspec. Then, generates a ASCII file for the background spectrum.
arfname = cztDetector.gen_arf(energy_lo = cztDetector.energy_low, energy_hi = cztDetector.energy_high, arf=chars["arf_name"])
rspname =cztDetector.gen_rsp(arfname, rsp = rspname)