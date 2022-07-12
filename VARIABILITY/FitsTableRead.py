#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 13:47:52 2020

File to open table.fits files and save them as a table

@author: Paola
"""

from astropy.io import fits

import matplotlib.pyplot as plt 
import numpy as np

#Variables: 

duplicates = 'no'#/'yes' 


file = fits.open('/Users/Paola/QUASAR/Work_EHVO/DATA/DR16Q_v4.fits', memmap=True)



data = file[1].data

#print(file[1].columns)
columns = ['SDSS_NAME','PLATE', 'MJD', 'FIBERID', 'Z', 'ZWARNING', 'Z_ERR', 'SOURCE_Z','BAL_PROB','BI_CIV','ERR_BI_CIV','AI_CIV','ERR_AI_CIV','BI_SiIV','ERR_BI_SiIV','AI_SIIV','ERR_AI_SiIV','NSPEC','SN_MEDIAN_ALL', 'M_I','EXTINCTION','MJD_DUPLICATE']

zem=data['z']

snr=data['SN_MEDIAN_ALL']

mjd = data['MJD']

aa=np.where((snr > 10.) & (zem > 1.9) & (mjd > 55752))
#aa=np.where((zem > 1.8) & (mjd > 55752))
print(len(aa[0]))


cutoffdata=data[aa]

print(cutoffdata['SDSS_NAME'])

NBINS = 200

z_hist = plt.hist(cutoffdata['Z'], NBINS)

# 9999/spec-9999-56789-0123.fits


#file_name='/Users/Paola/QUASAR/Work_EHVO/ROUTINES/list_DR16Q_cutoff_z_and_mjd_and_SN_MEDIAN_ALL_test2.txt'
file_name='/Users/Paola/QUASAR/Work_EHVO/ROUTINES/lala.txt'

f = open(file_name, 'w')


for qq in range(0,len(aa[0])):
    cutoffdata_now=cutoffdata[qq]
    sdssnamenow=cutoffdata_now['SDSS_NAME']
    znow=cutoffdata_now['Z']
    platenow=cutoffdata_now['PLATE']
    mjdnow=cutoffdata_now['MJD']
    fiberidnow=cutoffdata_now['FIBERID']
    nspecnow=cutoffdata_now['NSPEC']
    nspecbossnow=cutoffdata_now['NSPEC_BOSS']
    plateduplicatenow=cutoffdata_now['PLATE_DUPLICATE']
    mjdduplicatenow=cutoffdata_now['MJD_DUPLICATE']
    fiberidduplicatenow=cutoffdata_now['FIBERID_DUPLICATE']
    
    nn=str(platenow).zfill(4)+'/spec-'+str(platenow).zfill(4)+'-'+str(mjdnow).zfill(5)+'-'+str(fiberidnow).zfill(4)+'.fits'
    nn=str(platenow).zfill(4)+'  '+str(mjdnow).zfill(5)+'  '+str(fiberidnow).zfill(4)+'  '+ str(sdssnamenow)+'  '+str(round(znow,3))
    
    f.write(nn)
    f.write("\n")
    
    print('here: '+str(qq)+' '+str(nspecbossnow))
    
    if duplicates == 'yes':
        if nspecbossnow > 0:
            for qq2 in range(0,nspecbossnow):
                platenow=plateduplicatenow[qq2]
                mjdnow=mjdduplicatenow[qq2]
                fiberidnow=fiberidduplicatenow[qq2]
                if mjdnow > 55752:
                    nn=str(platenow).zfill(4)+'  '+str(mjdnow).zfill(5)+'  '+str(fiberidnow).zfill(4)+'  '+ str(sdssnamenow)+'  '+str(round(znow,3))
                   # nn=str(platenow).zfill(4)+'/spec-'+str(platenow).zfill(4)+'-'+str(mjdnow).zfill(5)+'-'+str(fiberidnow).zfill(4)+'.fits'
                    f.write(nn)
                    f.write("\n")
                    print('did it! '+str(qq2)+' '+str(nspecbossnow))
                    

        
f.close()
