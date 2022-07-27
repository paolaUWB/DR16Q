#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 16:05:09 2022

@author: Tzitzi
"""

import sys
import numpy as np
from astropy import*
from astropy import constants as const
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from astropy.io import fits
from utility_functions import read_spectra, append_row_to_csv, clear_file, print_to_file, read_list_spectra
import os
from os.path import exists
import pandas as pd
import csv

#specnumDR16=XXX
#specnumDR14=XXX
#specnumEHVO=98 
specnumPARENT=int(18165)

#The inputs in this program should be:
infoDR16 = os.getcwd() + "/DR16Q_v4 .fits" #DR16 fits file
infoDR14 = os.getcwd() + "/dr14q_spec_prop.fits"  #DR14Q table fits file from Rakshit+2020
infoEHVO = os.getcwd() + "/DR16Q/DR16Q_EHVO/good_fit_EHVO.csv" #the list of EHVOs (?? with more info or just the list)
infoPARENT = os.getcwd() + "/DR16Q/DR16_parent_sample.csv" #the list of DR16 parent sample


#%%
#-----reads csv + fits files

#SPECTRA_FILE_NAME, REDSHIFT, CALCULATED_SNR = read_list_spectra(infoEHVO, ['SPECTRA FILE NAME', 'REDSHIFT', 'CALCULATED SNR'])

#parent_spectra, parent_z, parent_snr = read_list_spectra(infoPARENT, ['SPECTRA', 'z', 'SNR'])  #just 


hdu_16 = fits.open(infoDR16)
data_16 = hdu_16[1].data
#print(data_16.columns)
SDSS_name_16 = data_16['SDSS_NAME'].astype(str)
plate_16 = data_16['PLATE  ']
mjd_16 = data_16['MJD']
fiber_16 = data_16['FIBERID ']
SNR_16 = data_16['SN_MEDIAN_ALL']
redshift_16 = data_16['z']
fits_16__duplicate_PLATE = data_16['PLATE_DUPLICATE']
fits_16__duplicate_MJD = data_16['MJD_DUPLICATE']
fits_16__duplicate_FIBER = data_16['FIBERID_DUPLICATE']
hdu_16.close()


#%%

hdu_14 = fits.open(infoDR14)
data_14 = hdu_14[1].data
SDSS_name_14 = data_16['SDSS_NAME'].astype(str)
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

#%%
df = pd.read_csv(infoPARENT, header=None)
parent_first_col = df[df.columns[0]].to_numpy()
spec_name_PARENT = parent_first_col.tolist()
#print(spec_name_PARENT)

FILE = os.getcwd() + "/DR16Q/CROSS_CORRELATION/DR16PA_DR14INFO.csv"
clear_file(FILE)

MBH_parent=[]
errMBH_parent=[]

for ii in range(18165):
    aa=spec_name_PARENT[ii]
    #print(aa)
    bb=aa.split('-')
    #print(bb)
    cc=bb[3].split('.')
    #print(cc)
    plate_PARENT = int (bb[1])
   # print(plate_PARENT)
    mjd_PARENT = int (bb[2])
    fiber_PARENT = int (cc[0])
    #print("pls run im gonna cry if it doesnt")
    vv, = np.where((plate_14[:] == plate_PARENT) & (mjd_14[:] == mjd_PARENT) & (fiber_14[:] == fiber_PARENT))
    #print(plate_14[:])
    if (len(vv,) != 0):
        MBH_parent.append = log_mbh_14[vv,]
        #print(vv)
        print(SDSS_name_14[vv,],plate_14[vv,],mjd_14[vv,],fiber_14[vv,],log_mbh_14[vv,],log_lbol_14[vv,], q_lbol_14[vv,], log_redd_14[vv,], q_redd_14[vv,])
        fields = [SDSS_name_14[vv,], plate_14[vv,], mjd_14[vv,], fiber_14[vv,], log_mbh_14[vv,],log_mbh_err_14[vv,], log_lbol_14[vv,], q_lbol_14[vv,], log_redd_14[vv,], q_redd_14[vv,]]
        append_row_to_csv(FILE, fields)
            
       # outdata.write((str(SDSS_name_14[vv,])),plate_14[vv,],mjd_14[vv,],fiber_14[vv,],log_mbh_14[vv,],log_lbol_14[vv,], log_redd_14[vv,])
        #rows = [{'SDSS_Names':str(SDSS_name_14[vv,])}]
        #with open('DR16PA_DR14INFO.csv', 'w', encoding='UTF8', newline='') as f:
         #   writer = csv.DicWriter(f, fieldnames=fieldnames)
          #  writer.writeheader()
           # writer.writerows(rows)
           #def append_row_to_csv(file_name: str, fields: list):
               #with open(file_name, 'a') as f:
                #   writer = csv.writer(f)
                 #  writer.writerow(fields)
