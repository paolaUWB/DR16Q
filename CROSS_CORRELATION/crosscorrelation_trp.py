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
<<<<<<< Updated upstream
=======
import pandas as pd
import csv
>>>>>>> Stashed changes


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
<<<<<<< Updated upstream
=======
fits_16__duplicate_PLATE = data_16['PLATE_DUPLICATE']
fits_16__duplicate_MJD = data_16['MJD_DUPLICATE']
fits_16__duplicate_FIBER = data_16['FIBERID_DUPLICATE']
>>>>>>> Stashed changes
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

<<<<<<< Updated upstream


# DR16_spec_name = OUT_DIREC + "/" + "BAL_BI_A.txt" *IGNORE*
# clear_file(BAL_BI_FILE_A)  *IGNORE*
# BAL_check_a_PDF = PdfPages('BAL_check_a.pdf')  *IGNORE*
for i in range(len(plate_16)):
    spectra_name_16 = "spec-" + str(plate_16[i]).zfill(4) + "-" + str(mjd_16[i]) + "-" + str(fiber_16[i]).zfill(4) + "-dered.dr16"
=======
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
           
#%% --- NAMING SPECTRA IN FITS FILES 
#for i in range(len(plate_16)):
 #   spectra_name_16 = "spec-" + str(plate_16[i]).zfill(4) + "-" + str(mjd_16[i]) + "-" + str(fiber_16[i]).zfill(4) + "-dered.dr16"
>>>>>>> Stashed changes
  #  print(spectra_name_16)
#print_to_file(spectra_name_a + ' ' + str(BI_CIV_a[i]) + ' '+ str(BI_CIV_err_a[i]) + ' '+ str(BI_ratio_a[i]), BAL_BI_FILE_A)  *IGNORE*

# for i in range(len(plate_14)):
#     spectra_name_14 = "spec-" + str(plate_14[i]).zfill(4) + "-" + str(mjd_14[i]) + "-" + str(fiber_14[i]).zfill(4) + "-dered.dr16"
#     print(spectra_name_14)
#%%
<<<<<<< Updated upstream
for i in range(len(parent_spectra)):
    SDSS_name_ALL_16 = []
    if spectra_name_16 == parent_spectra[i]:
        SDSS_name_ALL_16.append(SDSS_name_16)
print(SDSS_name_ALL_16)
        
=======
# =============================================================================
# for i in range(len(fits_duplicate_PLATE)):
#     duplicates = []
#     test_duplicates = []
#     
#     fits_names.append('spec-' + fits_PLATE[i].zfill(4) + '-' + fits_MJD[i] + '-' + fits_FIBER[i].zfill(4))
#     foldernames.append('J' + fits_SDSS_NAME[i]) ## move this lower (only needs to be created for ehvo duplicates)
#     
#     if (i % 5000) == 0:
#         print('You have loaded ', i, ' spectra!')
# 
#     if fits_duplicate_PLATE[i][0] != '-1':
#         
#         for ehvo in range(len(ehvo_spectra_list)):
#             if DR == '16':
#                 ehvo_spec = ehvo_spectra_list[ehvo][:-11]
#             if DR == '9': 
#                 ehvo_spec = ehvo_spectra_list[ehvo][:-4]
# 
# 
#             if ehvo_spec == fits_names[i]:
#                 ehvo_PLATE = fits_PLATE[i].zfill(4)
#                 ehvo_MJD = fits_MJD[i].zfill(4)
#                 ehvo_FIBER = fits_FIBER[i].zfill(4)
#                 
#                 for j in range(fits_nspec[i]):
#                     dup_spec_name = 'spec-' + fits_duplicate_PLATE[i][j].zfill(4) + '-' + fits_duplicate_MJD[i][j].zfill(4) + '-' + fits_duplicate_FIBER[i][j].zfill(4)
#                     duplicates.append(dup_spec_name)
#                     data_req_fields = [fits_SDSS_NAME[i], fits_duplicate_PLATE[i][j].zfill(4), fits_duplicate_MJD[i][j].zfill(4), fits_duplicate_FIBER[i][j].zfill(4)]
#                     append_row_to_csv(data_req_file, data_req_fields)
#                 
#                 ehvo_fields = [fits_SDSS_NAME[i], ehvo_PLATE, ehvo_MJD, ehvo_FIBER, duplicates]
#                 append_row_to_csv(ehvo_dup_file, ehvo_fields)
# =============================================================================
                
# =============================================================================
# #%%
# for i in range(len(fits_16__duplicate_PLATE)):
#     print("pls run im gonna cry if it doesnt")
#     for j in range(len(plate_14)):
#         print("if it doesnt get here im gonna drink bleach")
#         if ((fits_16__duplicate_PLATE[j] == plate_14[i]) & (fits_16__duplicate_MJD[j] == mjd_14[i]) & (fits_16__duplicate_FIBER[j] == fiber_14[i])):
#             spec_name_duplicate_dr16 = "spec-" + fits_16__duplicate_PLATE[j].zfill(4) + '-' + fits_16__duplicate_MJD[j].zfill(4) + '-' + fits_16__duplicate_FIBER[j].zfill(4) + '-dered.txt'
#             print(spec_name_duplicate_dr16)
#             #fields = [spec_name, fits_Z[i], fits_SNR[i]]
#             #append_row_to_csv(FILE, fields)
# =============================================================================

#%%
# =============================================================================
#  
#     Read the plate in this line [i] (plate_parentDR16)
#    Read the mjd
#   
#      vv,=where((PLATE_in14[:] == plate_parentDR16) & (MJD_DR7_in7[:] == ff2) & (FIBERID_DR7_in7[:] == ff3))
#       if (len(vv) != 0):
# 	MBH_in9[i]=MBHcorr_in7[vv]
# 
# aa='spec-000-1111-222.dr16'
# 
# bb=aa.split('-')
# 
# print(bb)
# ['spec', '000', '1111', '222.dr16']
# 
# bb[1]
# Out[65]: '000'
# 
# plate_parentD16 = int(bb[1])
# 
# print(plate_parentD16)
# 
# Output from spyder call 'get_namespace_view':
# 0
# 
# bb[2]
# Out[68]: '1111'
# 
# mjd_parentD16 = int(bb[2])
# 
# print(mjd_parentD16)
# 
# Output from spyder call 'get_var_properties':
# 1111
# 
# =============================================================================
>>>>>>> Stashed changes
