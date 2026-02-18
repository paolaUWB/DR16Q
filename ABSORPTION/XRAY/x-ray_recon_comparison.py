#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  6 15:25:15 2026

@author: lilianaflores
"""

import os
import matplotlib.pyplot as plt
import sys
import numpy as np
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
from spectra_comparison import Plot_spec_compare_morphed, Plot_spec_compare_og
sys.path.insert(0, os.getcwd()+'/../')
from abs_plot_module import draw_abs_figure
from utility_functions import read_list_spectra
sys.path.insert(0, os.getcwd()+'/../')
from data_types import Range
from abs_function_module import smooth, abs_parameters_plot_optional, wavelength_to_velocity


#defining the config file

CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()+"/spec_lists/xray_list_SNR10_z1.9.csv" #250 with SNR>10 & z>1.9

# range of spectra you are working with from the good_fit.csv file

STARTS_FROM, ENDS_AT = 1, 25 # eventually 1 -> 250

norm_spectra_list, redshift_list, calc_snr_list = read_list_spectra(CONFIG_FILE, ["NORM SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR"]) 



vlast = []
# whether abs_count or all_count is used is based on the value of all_plot_and_text|
abs_count = 0 # counter for amount of spectra that have absorption when all_plot_and_text = no and for text files
all_count = 0 # counter for all spectra ran when all_plot_and_text = yes

#keeping track of how many spectra do not have data points within the velocity limits
not_in_range = []

# loops over each spectra from a specified starting and ending point
output_morphed_pdf = "spectra_output_morphed.pdf"
output_og_pdf = "spectra_output_og.pdf"

with PdfPages(output_morphed_pdf) as pdf:
    for spectra_index in range(STARTS_FROM, ENDS_AT + 1):

        norm_spectrum_file_name = norm_spectra_list[spectra_index - 1]
        File = os.getcwd() + "/Recons_Hiremath2025/" + str(norm_spectrum_file_name)

        fig = Plot_spec_compare_morphed(File)

        pdf.savefig(fig)
        plt.show()
        plt.close(fig)

#for spectra_index in range(STARTS_FROM, ENDS_AT + 1): 
#    #keeping track of how many spectra do not have data points within the velocity limits
#    
#    norm_spectrum_file_name = norm_spectra_list[spectra_index - 1]
#    File = (os.getcwd() + "/Recons_Hiremath2025/" + str(norm_spectrum_file_name))
#    fig = Plot_spec_compare_morphed(File)
#    PdfPages.savefig(fig)
#    fig.save()
#    pyplot.gcf()
#    plt.show()
#    plt.close()
    
    '''
    data = File[1].data   
    
    
    we want to compare, original, reconstruction and error in terms of velocity (check line 75)
    look at spectra comparison program, examples of what the flux needs to be
    
    can feed the reconstruction and the original through the function draw_abs_figure; inside abs_plot_module
    alternatively, call the function defined within the spectra_comparison that can make the plots displayed, fairly automated
    
    whichever one was chosen, have it produce into a pdf, we found a function in absorption.py; check line 13
    
    '''
    #draw_abs_figure(spectra_count_abs, spectra_index, velocity, flux_normalized, error, savefile_name, spectra_name, redshift, snr, max_peak):
    
    
    '''
    RECONSTRUCTION_MORPHED_OUTPUT_PLOT_PDF = PdfPages(OUT_DIREC + 'BI' + str(BALNICITY_INDEX_LIMIT) + '.pdf')
    RECONSTRUCTION_MORPHED_OUTPUT_PLOT_PDF = PdfPages(OUT_DIREC + 'MORPHED' + str(File))
    
    
    # creates directory for output files
    OUT_DIREC = os.getcwd() + "/OUTPUT_FILES_PDF/"
    '''
with PdfPages(output_og_pdf) as pdf:     
    for spectra_index in range(STARTS_FROM, ENDS_AT + 1): 
        #keeping track of how many spectra do not have data points within the velocity limits
        norm_spectrum_file_name = norm_spectra_list[spectra_index - 1]
        File = (os.getcwd() + "/Recons_Hiremath2025/" + str(norm_spectrum_file_name))
        
        #draw_abs_figure(spectra_count_abs, spectra_index, velocity, flux_normalized, error, savefile_name, spectra_name, redshift, snr, max_peak):
        fig = Plot_spec_compare_og(File)
        
        pdf.savefig(fig)
        plt.show()
        plt.close(fig)
