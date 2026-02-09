#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  6 15:25:15 2026

@author: lilianaflores
"""

from matplotlib.backends.backend_pdf import PdfPages
from abs_plot_module import draw_abs_figure
from utility_functions import read_list_spectra
from astropy.io import fits
from abs_function
import os
import sys

#defining the config file
CONFIG_FILE1 = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()+"/spec_lists/xray_parent_list.csv" #2281 full parent sample
CONFIG_FILE2 = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()+"/spec_lists/xray_list_SNR10_z1.9.csv" #250 with SNR>10 & z>1.9
CONFIG_FILE3 = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()+"/spec_lists/xray_list_SNR10.csv" #268 with SNR>10

CONFIG_FILE = CONFIG_FILE2

# range of spectra you are working with from the good_fit.csv file
if CONFIG_FILE == CONFIG_FILE1:
    STARTS_FROM, ENDS_AT = 1, 2281
elif CONFIG_FILE == CONFIG_FILE2:
    STARTS_FROM, ENDS_AT = 1, 250
elif CONFIG_FILE == CONFIG_FILE3:
    STARTS_FROM, ENDS_AT = 1, 268
#STARTS_FROM, ENDS_AT = 1, 268 #uncomment to override if statement auto selection of range



norm_spectra_list, redshift_list, calc_snr_list = read_list_spectra(CONFIG_FILE, ["NORM SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR"]) 



vlast = []
# whether abs_count or all_count is used is based on the value of all_plot_and_text
abs_count = 0 # counter for amount of spectra that have absorption when all_plot_and_text = no and for text files
all_count = 0 # counter for all spectra ran when all_plot_and_text = yes

#keeping track of how many spectra do not have data points within the velocity limits
not_in_range = []

# loops over each spectra from a specified starting and ending point
for spectra_index in range(STARTS_FROM, ENDS_AT + 1): 
    #keeping track of how many spectra do not have data points within the velocity limits
    
    # rounding the numbers of the redshift, calculated snr and setting the norm file name to the current file name from the csv
    #z = round(redshift_list[spectra_index - 1], 5)
    z = 0  #EDITED to be able to include correct z value in plot but not correct for redshift when plotting in velocity
    z_final = round(redshift_list[spectra_index - 1], 2)

    calc_snr = round(calc_snr_list[spectra_index - 1], 5)
    norm_spectrum_file_name = norm_spectra_list[spectra_index - 1]
    
    #EDITED: To accomodate xray spectra in .fits format....................................................
    
    # from the norm spectra name retrieving it's wavelength, normalized flux, and normalized error (in this case from NORM_DRXQ)
   
    # setting a variable for each of those values from the spectra
    #wavelength, normalized_flux, normalized_error = read_spectra(norm_spectra_data) vvvvvv
    #STOP Liliana come fix this before uploading to Github
    File = fits.open(os.getcwd() + "/Recons_Hiremath2025/" + str(norm_spectrum_file_name))
    data = File[1].data   
    
    
    wavelength = data['Wave']
    flux = data["Flux"]
    normalized_error = data["Noise"] #is this normalized as it is though? Does the error array need to be divided by recon?
    recon = data["Recon"]
    normalized_flux = flux/recon
    
    #EDITED: This section so that the velocity search limit would adjust if the spectra does not reach the full range
    beta_test = wavelength_to_velocity(z, wavelength)
    
    
    VELOCITY_LIMIT = Range(-30000, -60000)

    min_beta = np.min(beta_test)
    
    if min_beta < VELOCITY_LIMIT.start:
        # Compute the new end but clamp it at -60000
        new_end = min(-min_beta, 60000)
        VELOCITY_LIMIT = Range(-30000, -new_end)
    else:
        VELOCITY_LIMIT = VELOCITY_LIMIT
