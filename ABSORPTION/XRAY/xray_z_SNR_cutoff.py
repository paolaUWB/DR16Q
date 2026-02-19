#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 12:32:28 2026

@author: lilianaflores
"""

import sys
import os
from astropy.io import fits
import numpy as np
import pandas as pd
from useful_wavelength_flux_error_modules import calculate_snr
from data_types import Range
from utility_functions import clear_file, read_list_spectra, read_spectra, append_row_to_csv

################### Changable Variables #############################################################################################################

file = os.getcwd() + '/BI_identified_sources_Hiremath2025.fits'

NORM_DIREC = os.getcwd() + "/Recons_Hiremath2025/"

OUT_DIREC = os.getcwd() + "/spec_lists/"

z_cutoff_value = 1.9

WAVELENGTH_FOR_SNR = Range(1276., 1400.) #tighter wavelength range for quasar sample after z cutoff needed due to wavelength min of some quasars

snr_cutoff_value = 10

csv_file_name = '/xray_list_SNR10_z1.9.csv'

#############################################################################################################

#Reading in data file and extracting data columns ........................................................................................................
data = fits.open(file)
#print(data[1].columns) #run to print column names
data_tab = data[1].data


#not used in this program..
BI = data_tab['BI']
Vmax = data_tab['VMAX']
Vmin = data_tab['VMIN']
EW = data_tab['EW']
MJD = data_tab['MJD']


specnames = data_tab['spectralist']
SNR = data_tab['SN_MEDIAN_ALL']
z = data_tab['zfinal']

#cutting folder name from specnames
for i in range(np.size(specnames)):
    specnames[i] = specnames[i][20:]
    

#cutoff with redshift z>1.9 .............................................................................
z_cut = np.where(z>z_cutoff_value)

data_tab_z_cut = data_tab.copy()[z_cut]
specnames_z = data_tab_z_cut['spectralist']

z_med = data_tab_z_cut['zfinal']
#print(np.size(z_med)) #1067 quasars after z>1.9 cutoff

#calculateing SNR values of parent sample with redshift cutoff applied. 1067 quasars ....................................................
calc_snr = []
min_wave = []
for i in np.arange(np.size(specnames_z)):
    specname_z = specnames_z[i]

    File = fits.open(NORM_DIREC + str(specname_z))
    dat = File[1].data   
    
    
    wavelength = dat['Wave']
    flux = dat["Flux"]
    normalized_error = dat["Noise"] #is this normalized as it is though? Does the error array need to be divided by recon?
    recon = dat["Recon"]
    normalized_flux = flux/recon
    min_wave.append(np.min(wavelength))
    
    normalized_error[normalized_error==0] = np.nan
    
    z = 0 #Only because the spectra is already in restframe.

    #calculate_snr(wavelength, z: float, WAVELENGTH_FOR_SNR: range, flux, error)
    snr_mean = calculate_snr(wavelength, z, WAVELENGTH_FOR_SNR, normalized_flux, normalized_error)
    #          calculate_snr(wavelength, z, WAVELENGTH_FOR_SNR, flux, error)

    calc_snr.append(snr_mean)

#print(max(min_wave)) #After the z>1.9 cutoff the maximum minimum wavelength value is 1275.2637064597302 so the range for the SNR calc must be greater

calc_snr = np.array(calc_snr)
print(f'# of incorrect SNR calculations: {np.sum(calc_snr==np.inf)}')

#applying SNR>10 cutoff ........................................................................................................
SNR_mask = np.where(calc_snr>snr_cutoff_value)

specnames = specnames_z[SNR_mask]
index = np.array(np.arange(1,np.size(specnames)+1))
zeros = np.zeros_like(index, dtype=np.int32)

#creating dataframe with columns needed
da = {'SPECTRA INDEX':index,	
       'SPECTRA FILE NAME': specnames , 
       'NORM SPECTRA FILE NAME': specnames,	
       #'REDSHIFT': zeros,	#redshift columns of zeros because spectra are in restframe
       'REDSHIFT': z_med[SNR_mask],	#redshift columns of zeros because spectra are in restframe
       'CALCULATED SNR': calc_snr[SNR_mask]}


df = pd.DataFrame(da)
print(np.max(df['SPECTRA INDEX'])) #Result is 250 with SNR cutoff of 10 and redshift cutoff of 1.9

df.to_csv(OUT_DIREC + csv_file_name,index=False)


# Pick only the rows that passed both cutoffs
filtered_data = data_tab_z_cut[SNR_mask]
#print(filtered_data.columns)
#print(len(filtered_data))


la = {'MJD':filtered_data['MJD'], 
      'RUN2D':filtered_data['RUN2D'], 
      'CATALOGID':filtered_data['CATALOGID'], 
      'SDSS_ID':filtered_data['SDSS_ID'],
      'Z':filtered_data['zfinal']}


# Cutting z and SNR filtered data to the columns needed for cross search with SDSS to find field vals
df_filtered_cut = pd.DataFrame(la,index=None)
print(len(df_filtered_cut))


# Saving to csv
csv_file_name2 = '/xray_SDSSCross_SNR10_z1.9.csv'
df_filtered_cut.to_csv(OUT_DIREC + csv_file_name2, index=False)












