#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 12:17:36 2026

@author: lilianaflores
"""

import numpy as np
import pandas as pd
import os 
import sys
from astropy.io import fits
import matplotlib.pyplot as plt
sys.path.insert(0, os.getcwd() + '/../')
from abs_function_module import wavelength_to_velocity



sdss_file = os.getcwd() + '/SDSS_fullspecs/spec-allepoch-59187-4382190169.fits'

'''
When using SDSS files:
    - Must correct for redshift
    - using full spec more info: https://data.sdss.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/full/FIELD/MJD/specFull.html
      for 249 cases of the parent sample in xray with z > 1.9 and SNR > 10
    - For the case 59187-4382190169 the field number is 0 so full spec cannot be downloaded. I was able to download an allepoch
      version that still extends past the data limit where spectra is cutoff on the files provided by Hiremath.
      - https://skyserver.sdss.org/dr19/VisualTools/explore/summary?sId=7084023500000005918702060103
      - download by pasting https://data.sdss.org/sas/dr19/spectro/boss/redux/v6_1_3/spectra/lite/allepoch/59187/spec-allepoch-59187-4382190169.fits
      - use spec-allepoch-59187-4382190169.fits for this case

    - more info on specFULL vs allepoch
'''


sdss = fits.open(sdss_file)
sdss_tab = sdss[1].data
sdss_tab.columns

sdss_wavelength = 10**sdss_tab['LOGLAM'] # log10(Angs)
sdss_flux = sdss_tab['FLUX'] # 10^-17 ergs/s/cm^2/Angs
sdss_ivar = sdss_tab['IVAR']
sdss_rest_lam = sdss_tab['WRESL'] # ang   initially thought could be ret wavelength but I don't think so

sdss_noise = np.zeros_like(sdss_flux)
sdss_redshift = np.zeros_like(sdss_flux)
    
sdss.close()

wave = sdss_wavelength
#wave = sdss_rest_lam

z = 2.606959 #real z
#z = 0 

beta2 = wavelength_to_velocity(z, wave)


plt.plot(beta2, sdss_flux, color = 'red', label='sdss full')
plt.xlim(-70000,0)
plt.ylim(0,60)
plt.show()

'''
Scale of flux seems off here when comparing to the un normalized spectra of the files Hiremath provided.

'''
