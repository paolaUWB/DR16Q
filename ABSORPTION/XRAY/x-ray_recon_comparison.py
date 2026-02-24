#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  6 15:25:15 2026

@author: lilianaflores
"""

import os
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages
from spectra_comparison import Plot_spec_compare_morphed, Plot_spec_compare_og
sys.path.insert(0, os.getcwd()+'/../')
from utility_functions import read_list_spectra
sys.path.insert(0, os.getcwd()+'/../')

"""
Great work on the x-ray_recon_comparison.py program! Some ideas for improvements to the program and the plotting functions called:
The function already allows for providing xlimits. We should run with x limits -70,000 to 0 in order to see the range we are looking for EHVOs in.
Automatic y limits: Could be by having the function find the max y value of the spectra within the x limits.
Clean any excess unused code or comments
Makes things more compact (optional): Combine the functions for og vs recon and morphed vs recon. Could be with if statements and a function parameter.
"""
#defining the config file
CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()+"/spec_lists/xray_list_SNR10_z1.9.csv" #250 with SNR>10 & z>1.9

# range of spectra you are working with from the good_fit.csv file
STARTS_FROM, ENDS_AT = 1, 250 # eventually 1 -> 250

norm_spectra_list, redshift_list, calc_snr_list = read_list_spectra(CONFIG_FILE, ["NORM SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR"]) 


vlast = []

# whether abs_count or all_count is used is based on the value of all_plot_and_text|
abs_count = 0 # counter for amount of spectra that have absorption when all_plot_and_text = no and for text files
all_count = 0 # counter for all spectra ran when all_plot_and_text = yes

#keeping track of how many spectra do not have data points within the velocity limits
not_in_range = []

# loops over each spectra from a specified starting and ending point
output_folder = os.path.join(os.getcwd(), "OUTPUT_FILES_PDF")

output_morphed_pdf = os.path.join(output_folder, "spectra_output_morphed.pdf")
output_og_pdf = os.path.join(output_folder, "spectra_output_og.pdf")
xlims = -70000, 0

with PdfPages(output_morphed_pdf) as pdf:
    for spectra_index in range(STARTS_FROM, ENDS_AT + 1):
        #closes the global storage of figures
        plt.close('all')

        plt.figure()
        
        #keeping track of how many spectra do not have data points within the velocity limits
        norm_spectrum_file_name = norm_spectra_list[spectra_index - 1]
        #gets the cwd and and names each files corectly
        File = os.getcwd() + "/Recons_Hiremath2025/" + str(norm_spectrum_file_name)

        #Plot of x-ray selected quasar spectra comparing the Morphed, Normalized, and Reconstruction spectra. Found in spectra_comparison.py
        fig = Plot_spec_compare_morphed(File, xlims)
        
        #Saves the current figure into a page of the pdf
        pdf.savefig(fig)
        plt.show(fig)
        plt.close(fig)



with PdfPages(output_og_pdf) as pdf:     
    for spectra_index in range(STARTS_FROM, ENDS_AT + 1): 
        #closes the global storage of figures
        plt.close('all')
        
        plt.figure()
        #keeping track of how many spectra do not have data points within the velocity limits
        norm_spectrum_file_name = norm_spectra_list[spectra_index - 1]
        File = (os.getcwd() + "/Recons_Hiremath2025/" + str(norm_spectrum_file_name))
        
        #Plot of x-ray selected quasar spectra comparing the original spectrum and reconstruction. Found in spectra_comparison.py
        fig = Plot_spec_compare_og(File, xlims)
        
        #creates a fig and uses the fucntion found in spectra_comparison.py
        pdf.savefig(fig)
        plt.show(fig)
        plt.close(fig)
