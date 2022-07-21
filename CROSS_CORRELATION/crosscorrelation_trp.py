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
infoDR16 = "/Volumes/MyPassport/Fits_Files/DR16Q_v4 .fits" #DR16 fits file
infoDR14 = "/Volumes/MyPassport/Fits_Files/dr14q_spec_prop.fits"  #DR14Q table fits file from Rakshit+2020

infoEHVO = '/Volumes/MyPassport/DR16Q/DR16Q_EHVO' #the list of EHVOs (?? with more info or just the list)
SPECTRA_FILE_NAME, REDSHIFT, CALCULATED_SNR = read_list_spectra(infoEHVO + '/good_fit_EHVO.csv', ['SPECTRA FILE NAME', 'REDSHIFT', 'CALCULATED SNR'])

infoparent = '/Volumes/MyPassport/DR16Q' #the list of DR16 parent sample
parent_spectra, parent_z, parent_snr = read_list_spectra(infoparent + '/DR16_parent_sample.csv', ['SPECTRA', 'z', 'SNR'])


hdu_16 = fits.open(infoDR16)
data_16 = hdu_16[1].data
SDSS_name_16 = data_16['SDSS_NAME']
plate_16 = data_16['PLATE  ']
mjd_16 = data_16['MJD     ']
fiber_16 = data_16['FIBERID ']
SNR_16 = data_16['SN_MEDIAN_ALL']
redshift_16 = data_16['z']
hdu_16.close()

hdu_14 = fits.open(infoDR14)
data_14 = hdu_14[1].data
SDSS_name_14 = data_16['SDSS_NAME']
plate_14 = data_14['PLATE  ']
mjd_14 = data_14['MJD     ']
fiber_14 = data_14['FIBERID ']
log_mbh_14 = data_14['LOG_MBH ']
log_mbh_err_14 = data_14['LOG_MBH_ERR']
q_mbh_14 = data_14['QUALITY_MBH']
log_lbol_14 = data_14['LOG_LBOL']
q_lbol_14 = data_14['QUALITY_LBOL']
log_redd_14 = data_14['LOG_REDD']
q_redd_14 = data_14['QUALITY_REDD']
bi_civ_14 = data_14['BI_CIV  ']
bi_civ_err_14 = data_14['ERR_BI_CIV']
bal_flag_14 = data_14['BAL_FLAG']
hdu_14.close()



# DR16_spec_name = OUT_DIREC + "/" + "BAL_BI_A.txt" *IGNORE*
# clear_file(BAL_BI_FILE_A)  *IGNORE*
# BAL_check_a_PDF = PdfPages('BAL_check_a.pdf')  *IGNORE*
for i in range(len(plate_16)):
    spectra_name_16 = "spec-" + str(plate_16[i]).zfill(4) + "-" + str(mjd_16[i]) + "-" + str(fiber_16[i]).zfill(4) + "-dered.dr16"
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
        