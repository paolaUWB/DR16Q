# -*- coding: utf-8 -*-
"""
move_files.py
=============

This code sorts spectra files of duplicate observations of objects to directories (organized by object). 

@author Mikel Charles

"""

#%%
import os
import sys
import numpy as np
from astropy.io import fits
import shutil
# from VARIABILITY.EHVO_files import FITS_DUPLICATES_FILE
sys.path.insert(0, os.getcwd() + '/../' + 'DR16Q') # changes the directory to the DR16Q --> all paths after this will need to be written as if this was in the top level of the DR16Q
from utility_functions import print_to_file, read_file, clear_file

DR = '16' # which data release are you working with? DRX (i.e. '9' or '16') 

#-- fits table 
fits_direc = os.getcwd() + '/../DR16Q_DATA/DR16Q_v4.fits'

#-- location of object directories
make_direc = os.getcwd() + '/VARIABILITY/DATA_VARIABILITY/'

#-- files to read
data_req_file = os.getcwd() + '/VARIABILITY/data_request.csv'
ehvo_dup_file = os.getcwd() + '/VARIABILITY/ehvo_duplicate_file.csv'

##------ read EHVOs
#-- dr9
if DR == '9':
    ehvo_file = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + '/DR9Q_EHVO/DR9_EHVO_sorted_norm.csv'
    zem, snr, ehvo_spectra_list = read_file(ehvo_file)
    
#-- dr16
if DR == '16':
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

fits_nspec = data['NSPEC'].astype(int)

fits_duplicate_PLATE = data['PLATE_DUPLICATE'].astype(str)
fits_duplicate_MJD = data['MJD_DUPLICATE'].astype(str)
fits_duplicate_FIBER = data['FIBERID_DUPLICATE'].astype(str)

hdu.close()


#%%
##------ move ehvo files to proper directory
ehvo_plate = []
ehvo_mjd = []
ehvo_fiber = []

with open(ehvo_dup_file) as e:
    for line in e: 
        ehvo_row = line.split(",")
        ehvo_plate.append(ehvo_row[1])
        ehvo_mjd.append(ehvo_row[2])
        ehvo_fiber.append(ehvo_row[3])

for file in os.listdir(os.getcwd() + '/VARIABILITY/'):
    
    
# for k in range(len(fits_SDSS_NAME)):
#     for l in range(len(ehvo_plate)):
#         if (fits_PLATE[k] == ehvo_plate[l]) and (fits_MJD[k] == ehvo_mjd[l]) and (fits_FIBER[k] == ehvo_fiber[l]):
            

# if DR == '16':
#     shutil.copyfile(SPEC_DIREC + EHVO_spectra_list[k], MAKE_DIREC + foldernames[i] + '/' + EHVO_spectra_list[k])
# if DR == '9': 
#     EHVO_spectra_list_norm = EHVO_spectra_list[k][:-4] + 'norm.dr9'
#     shutil.copyfile(SPEC_DIREC + '/' + EHVO_spectra_list_norm, os.getcwd() + '/VARIABILITY/DATA_VARIABILITY/' + foldernames[i] + '/' + EHVO_spectra_list_norm)




#%%
duplicate_plate = []
duplicate_mjd = []
duplicate_fiber = []


with open(data_req_file) as f:    
    for line in f: 
        row = line.split(",")
        duplicate_plate.append(row[1])
        duplicate_mjd.append(row[2])
        duplicate_fiber.append(row[3])

for i in range(len(fits_SDSS_NAME)):
    
    for j in range(len(duplicate_plate) - 1):
        print(duplicate_plate[i+1])
        if fits_SDSS_NAME[i] == duplicate_plate[i+1]:
            
    
        


