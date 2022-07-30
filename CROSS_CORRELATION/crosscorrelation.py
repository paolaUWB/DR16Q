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

<<<<<<< Updated upstream:CROSS_CORRELATION/crosscorrelation_trp.py


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
        
=======


#%%
df = pd.read_csv(infoPARENT, header=None)
parent_first_col = df[df.columns[0]].to_numpy()
spec_name_PARENT = parent_first_col.tolist()
#print(spec_name_PARENT)

FILE = os.getcwd() + "/CROSS_CORRELATION/DR16PA_DR14INFO.csv"
clear_file(FILE)

# save_new_output_file = 'yes' ## DO YOU WANT TO SAVE TO THE OUTPUT FILES? 'yes'/'no'


# if save_new_output_file =='yes':
#     DR16_DR14file = os.getcwd() + "/CROSS_CORRELATION/" + "DR1614infotest" + ".csv"

# else:
#     FILE = os.getcwd() + "/CROSS_CORRELATION/DR16PA_DR14INFO.csv"
#     clear_file(FILE)

MBH_parent=[]
errMBH_parent=[]
PAheaders = ['SDSS_Names_14', 'Plate_14', 'Mjd_14', 'Fiber_14', 'Mbh_14','Mbh_err_14', 'Lbol_14', 'Quality_Lbol_14', 'Redd_14', 'Quality_Redd_14']

with open(FILE, 'w') as file:
    dw = csv.DictWriter(file, fieldnames = PAheaders)
    dw.writeheader() 


    
# reader = csv.reader(file)
# writer = csv.writer(csvfile, fmtparams)

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
        MBH_parent.append(log_mbh_14[vv,])
        # MBH_parent.append = log_mbh_14[vv,]
        #print(vv)
        #print(SDSS_name_14[vv,],plate_14[vv,],mjd_14[vv,],fiber_14[vv,],log_mbh_14[vv,],log_mbh_err_14[vv,],log_lbol_14[vv,], q_lbol_14[vv,], log_redd_14[vv,], q_redd_14[vv,])
        SDSS_name0 = SDSS_name_14[vv,]
        fields = [SDSS_name0[0], int(plate_14[vv,]), int(mjd_14[vv,]), int(fiber_14[vv,]), float(log_mbh_14[vv,]),float(log_mbh_err_14[vv,]), float(log_lbol_14[vv,]), int(q_lbol_14[vv,]), float(log_redd_14[vv,]), int(q_redd_14[vv,])]
        append_row_to_csv(FILE, fields)
        
        


#%%
#Histogram codes:
MBH_parent=np.hstack(MBH_parent)
MBH_parent_toplot=MBH_parent[np.hstack(np.where(MBH_parent != 0.))]
 
fig=plt.figure(1)
 
iqr = np.subtract(*np.percentile(MBH_parent_toplot, [75, 25]))
nhist=(max(MBH_parent_toplot)-min(MBH_parent_toplot))/(2*iqr*(len(MBH_parent_toplot)**(-1/3)))
 
bins=np.linspace(min(MBH_parent_toplot),max(MBH_parent_toplot),int(nhist))
 
plt.hist([MBH_parent_toplot],bins,color=['black'],label=['parent'],histtype='step')
# This for later when you have the 3 samples: plt.hist([MBH,MBH_BAL_toplot,MBH_EHVO_toplot],bins,color=['black','blue','red'],label=['parent','BAL','EHVO'],histtype='step')
 
plt.xlabel(r'log($M_{\mathrm{BH}}/M_{\odot})$')
#plt.ylabel('$n/N_{tot}$')
plt.ylabel('Number')
plt.legend(loc='upper left')


