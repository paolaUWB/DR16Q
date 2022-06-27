# -*- coding: utf-8 -*-
"""
Created on Sun May  8 12:14:27 2022

@author: Dakota
"""

import numpy as np
import os
from astropy.io import fits
import shutil

directory = os.getcwd()
fits_file = os.getcwd() + '/../ + VARIABILITY_DATA/DR16Q_v4.fits' ## may change this?? 

hdu = fits.open(fits_file)
data = hdu[1].data

fits_plate = data['PLATE   '].astype(str)
fits_mjd = data['MJD     '].astype(str)
fits_fiber = data['FIBERID '].astype(str)
fits_objid = data['OBJID   ']
fits_duplicate_plate = data['PLATE_DUPLICATE'].astype(str)
fits_duplicate_mjd = data['MJD_DUPLICATE'].astype(str)
fits_duplicate_fiber = data['FIBERID_DUPLICATE'].astype(str)
hdu.close()

fits_names = []
fits_duplicates = []
duplicates = []

for i in range(len(fits_plate)):
    duplicates = []
    fits_names.append('spec-' + fits_plate[i].zfill(4) + '-' + fits_mjd[i] + \
                      '-' + fits_fiber[i].zfill(4))
    if (i % 50000) == 0:
        print("You have loaded", i, "spectra!")
    if np.any(np.array(fits_duplicate_plate[i]) != '-1'):
        for j in range(74):
            if fits_duplicate_plate[i][j] != '-1':
                duplicates.append('spec-' + fits_duplicate_plate[i][j].zfill(4) + '-' + \
                                  fits_duplicate_mjd[i][j] + '-' + fits_duplicate_fiber[i][j].zfill(4))
    if duplicates != []:
        fits_duplicates.append(duplicates)

    else:
        fits_duplicates.append('-1')
#%%    
filenames = []
foldernames = []
indexs = []
for file in os.listdir(directory):
    loc = 0
    if file.startswith('spec'):
        filenames.append(file)
    if file.endswith('dered.txt'):
        try:
            loc = np.where(file[:-10] == np.array(fits_names))[0][0]
        except IndexError:
            pass
        if loc == 0:
            for i in range(len(fits_duplicates)):
                if type(fits_duplicates[i]) == list:
                    for j in range(len(fits_duplicates[i])):
                        if file[:-10] == fits_duplicates[i][j]:
                            loc = i
    if loc != 0:
        indexs.append(loc)

    if file.startswith('spec'):
        foldernames.append('J' + fits_objid[loc])

og_names = np.array(fits_names)[indexs]
data_dir = os.cwd() + '/DATA/DR9Q_SNR10' #possibly: /../ + /DATA/DR9Q_SNR10 #'C:/Users/Dakota/Documents/GitHub/DR16Q/DATA/DR9Q_SNR10'

for i in range(len(filenames)):
    try:
        os.mkdir(directory + '/DATA_VARIABILITY/' + foldernames[i])
    except FileExistsError:
        pass
    os.replace(directory + '/' + filenames[i], directory + \
               '/DATA_VARIABILITY/' + foldernames[i] + '/' + filenames[i])
#%%
for i in range(len(og_names)):
    try:
        loc = np.where(og_names[i] + '-dered.dr16' == np.array(os.listdir(data_dir)))[0][0]
        print('yes')
        
    except IndexError:
        pass