#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 14:05:05 2022

@author: mikelcharles
"""


## include SDSS name, plate, mjd, fiber, spec name, z, SNR (calculated SNR)

#%%
import os
import sys
from astropy.io import fits
sys.path.insert(0, os.getcwd() + '/../' + 'DR16Q/') # changes the directory to the DR16Q --> all paths after this will need to be written as if this was in the top level of the DR16Q
from utility_functions import read_file, clear_file, append_row_to_csv#, print_to_file

fits_direc = os.getcwd() + '/../DR16Q_DATA/DR16Q_v4.fits'

ehvo_info_file = os.getcwd() + '/DR16_ehvo_info_file.csv'

#-- dr16
ehvo_file = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + '/DR16Q_EHVO/DR16_EHVO_sorted_norm.csv'
zem, snr, ehvo_spectra_list = read_file(ehvo_file)


#%%
##------ open fits table
hdu = fits.open(fits_direc)
data = hdu[1].data

fits_SDSS_NAME = data['SDSS_NAME'].astype(str)

fits_PLATE = data['PLATE'].astype(str) 
fits_MJD = data['MJD'].astype(str)
fits_FIBER = data['FIBERID'].astype(str)

fits_ZEM = data['Z'].astype(str)

hdu.close()


#%%
##------ clear csv files
clear_file(ehvo_info_file)

fits_names = []

ehvo_SDSS_NAME = []
ehvo_PLATE = []
ehvo_MJD = []
ehvo_FIBER = []
ehvo_SPEC_NAME = []
ehvo_ZEM = []
ehvo_SNR = []



for i in range(len(fits_PLATE)):
    ehvo_SDSS_NAME = []
    ehvo_PLATE = []
    ehvo_MJD = []
    ehvo_FIBER = []
    ehvo_SPEC_NAME = []
    ehvo_ZEM = []
    ehvo_SNR = []
    
    fits_names.append('spec-' + fits_PLATE[i].zfill(4) + '-' + fits_MJD[i].zfill(4) + '-' + fits_FIBER[i].zfill(4))
    
    for ehvo in range(len(ehvo_spectra_list)):
        ehvo_spec = ehvo_spectra_list[ehvo][:-11]
        
        if ehvo_spec == fits_names[i]:
            ehvo_PLATE.append(fits_PLATE[i].zfill(4))
            ehvo_MJD.append(fits_MJD[i].zfill(4))
            ehvo_FIBER.append(fits_FIBER[i].zfill(4))
            
            ehvo_SDSS_NAME.append(fits_SDSS_NAME[i])
            ehvo_SPEC_NAME.append('spec-' + ehvo_PLATE[0] + '-' + ehvo_MJD[0] + '-' + ehvo_FIBER[0] + '-dered.dr16')
            
            ehvo_ZEM.append(fits_ZEM[i])
            ehvo_SNR.append(snr[ehvo]) 
            
            ehvo_fields = [ehvo_SDSS_NAME, ehvo_PLATE, ehvo_MJD, ehvo_FIBER, ehvo_SPEC_NAME, ehvo_ZEM, ehvo_SNR]
            append_row_to_csv(ehvo_info_file, ehvo_fields)
            
            
            
            
            