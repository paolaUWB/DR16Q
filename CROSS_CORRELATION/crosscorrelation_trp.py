#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 16:05:09 2022

@author: Tzitzi
"""

import sys
import numpy as np
import scipy.stats as stat
from pylab import*
from sympy import sympify
from scipy import*
from astropy import*
from scipy.stats import ks_2samp
from astropy import constants as const
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from astropy.io import fits
from utility_functions import read_spectra, append_row_to_csv, clear_file, print_to_file, read_list_spectra
import os
from os.path import exists


#specnumDR16=XXX
#specnumDR14=XXX
#specnumEHVO=98 

#The inputs in this program should be:
infoDR16 = os.getcwd() + "/DR16Q_v4 .fits" #DR16 fits file
infoDR14 = os.getcwd() + "/dr14q_spec_prop.fits"  #DR14Q table fits file from Rakshit+2020
infoEHVO = os.getcwd() + "/DR16Q/DR16Q_EHVO/good_fit_EHVO.csv" #the list of EHVOs (?? with more info or just the list)
infoPARENT = os.getcwd() + "/DR16Q/DR16_parent_sample.csv" #the list of DR16 parent sample


#%%
#-----reads csv + fits files

SPECTRA_FILE_NAME, REDSHIFT, CALCULATED_SNR = read_list_spectra(infoEHVO, ['SPECTRA FILE NAME', 'REDSHIFT', 'CALCULATED SNR'])

parent_spectra, parent_z, parent_snr = read_list_spectra(infoPARENT, ['SPECTRA', 'z', 'SNR'])

hdu_16 = fits.open(infoDR16)
data_16 = hdu_16[1].data
SDSS_name_16 = data_16['SDSS_NAME'].astype(str)
plate_16 = data_16['PLATE  '].astype(str)
mjd_16 = data_16['MJD     '].astype(str)
fiber_16 = data_16['FIBERID '].astype(str)
SNR_16 = data_16['SN_MEDIAN_ALL'].astype(str)
redshift_16 = data_16['z'].astype(str)
fits_16__duplicate_PLATE = data_16['PLATE_DUPLICATE'].astype(str)
fits_16__duplicate_MJD = data_16['MJD_DUPLICATE'].astype(str)
fits_16__duplicate_FIBER = data_16['FIBERID_DUPLICATE'].astype(str)
hdu_16.close()



hdu_14 = fits.open(infoDR14)
data_14 = hdu_14[1].data
SDSS_name_14 = data_16['SDSS_NAME'].astype(str)
plate_14 = data_14['PLATE  '].astype(str)
mjd_14 = data_14['MJD     '].astype(str)
fiber_14 = data_14['FIBERID '].astype(str)
log_mbh_14 = data_14['LOG_MBH '].astype(str)
log_mbh_err_14 = data_14['LOG_MBH_ERR'].astype(str)
q_mbh_14 = data_14['QUALITY_MBH'].astype(str)
log_lbol_14 = data_14['LOG_LBOL'].astype(str)
q_lbol_14 = data_14['QUALITY_LBOL'].astype(str)
log_redd_14 = data_14['LOG_REDD'].astype(str)
q_redd_14 = data_14['QUALITY_REDD'].astype(str)
bi_civ_14 = data_14['BI_CIV  '].astype(str)
bi_civ_err_14 = data_14['ERR_BI_CIV'].astype(str)
bal_flag_14 = data_14['BAL_FLAG'].astype(str)
hdu_14.close()

#%% --- NAMING SPECTRA IN FITS FILES 
#for i in range(len(plate_16)):
 #   spectra_name_16 = "spec-" + str(plate_16[i]).zfill(4) + "-" + str(mjd_16[i]) + "-" + str(fiber_16[i]).zfill(4) + "-dered.dr16"
  #  print(spectra_name_16)
#print_to_file(spectra_name_a + ' ' + str(BI_CIV_a[i]) + ' '+ str(BI_CIV_err_a[i]) + ' '+ str(BI_ratio_a[i]), BAL_BI_FILE_A)  *IGNORE*

# for i in range(len(plate_14)):
#     spectra_name_14 = "spec-" + str(plate_14[i]).zfill(4) + "-" + str(mjd_14[i]) + "-" + str(fiber_14[i]).zfill(4) + "-dered.dr16"
#     print(spectra_name_14)
#%%
for i in range(len(parent_spectra)):
    SDSS_name_ALL_16 = []
    if spectra_name_16 == parent_spectra[i]:
        SDSS_name_ALL_16.append(SDSS_name_16)
print(SDSS_name_ALL_16)

#%%
dups_file = os.getcwd() + '/CROSS_CORRELATION/duplicates_with_info.csv'

##------ clear csv files
clear_file(dups_file)

##------ create directories
fits_names = []
fits_duplicates = []

ehvo_names = []
ehvo_duplicates = []

foldernames = []

ehvo_fields = ['sdssName', 'ehvoPlate', 'ehvoMJD', 'ehvoFiber', 'duplicate_spectra']
append_row_to_csv(ehvo_dup_file, ehvo_fields)

data_req_fields = ['sdssName', 'Plate', 'MJD', 'Fiber']
append_row_to_csv(data_req_file, data_req_fields)


#%%
for i in range(len(fits_duplicate_PLATE)):
    duplicates = []
    test_duplicates = []
    
    fits_names.append('spec-' + fits_PLATE[i].zfill(4) + '-' + fits_MJD[i] + '-' + fits_FIBER[i].zfill(4))
    foldernames.append('J' + fits_SDSS_NAME[i]) ## move this lower (only needs to be created for ehvo duplicates)
    
    if (i % 50000) == 0:
        print('You have loaded ', i, ' spectra!')

    if fits_duplicate_PLATE[i][0] != '-1':
        
        for ehvo in range(len(ehvo_spectra_list)):
            if DR == '16':
                ehvo_spec = ehvo_spectra_list[ehvo][:-11]
            if DR == '9': 
                ehvo_spec = ehvo_spectra_list[ehvo][:-4]


            if ehvo_spec == fits_names[i]:
                ehvo_PLATE = fits_PLATE[i].zfill(4)
                ehvo_MJD = fits_MJD[i].zfill(4)
                ehvo_FIBER = fits_FIBER[i].zfill(4)
                
                for j in range(fits_nspec[i]):
                    dup_spec_name = 'spec-' + fits_duplicate_PLATE[i][j].zfill(4) + '-' + fits_duplicate_MJD[i][j].zfill(4) + '-' + fits_duplicate_FIBER[i][j].zfill(4)
                    duplicates.append(dup_spec_name)
                    data_req_fields = [fits_SDSS_NAME[i], fits_duplicate_PLATE[i][j].zfill(4), fits_duplicate_MJD[i][j].zfill(4), fits_duplicate_FIBER[i][j].zfill(4)]
                    append_row_to_csv(data_req_file, data_req_fields)
                
                ehvo_fields = [fits_SDSS_NAME[i], ehvo_PLATE, ehvo_MJD, ehvo_FIBER, duplicates]
                append_row_to_csv(ehvo_dup_file, ehvo_fields)
                
#%%
for i in range(len(fits_16__duplicate_PLATE)):
    print("pls run im gonna cry if it doesnt")
    for j in range(len(plate_14)):
        print("if it doesnt get here im gonna drink bleach")
        if ((fits_16__duplicate_PLATE[j] == plate_14[i]) & (fits_16__duplicate_MJD[j] == mjd_14[i]) & (fits_16__duplicate_FIBER[j] == fiber_14[i])):
            spec_name = "spec-" + fits_16__duplicate_PLATE[j].zfill(4) + '-' + fits_16__duplicate_MJD[j].zfill(4) + '-' + fits_16__duplicate_FIBER[j].zfill(4) + '-dered.txt'
            print(spec_name)
            #fields = [spec_name, fits_Z[i], fits_SNR[i]]
            #append_row_to_csv(FILE, fields)