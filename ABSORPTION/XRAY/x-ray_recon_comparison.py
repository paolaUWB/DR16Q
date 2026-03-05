#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  6 15:25:15 2026

@author: lilianaflores and elijahfacklam
"""

import os
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages
from spectra_comparison import Plot_spec_compare_morphed, Plot_spec_compare_og
sys.path.insert(0, os.getcwd()+'/../')
from utility_functions import read_list_spectra
sys.path.insert(0, os.getcwd()+'/../')


######################################## PATH FINDING ########################################

#defining the config file

'''
# 2nd csv with all file names of the SDSS spectrum
# Give dummy zeros for redshift and SNR 
# Make sure the names are the same and match match middle set of numbers
# Make sure to loop them in order together
# Make a folder if you dont have one from the output_folder
# Make a changable variable for when you want to show fig
'''
CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()+"/spec_lists/xray_list_SNR10_z1.9.csv" #250 with SNR>10 & z>1.9

# range of spectra you are working with from the good_fit.csv file
STARTS_FROM, ENDS_AT = 1, 250 # eventually 1 -> 250

norm_spectra_list, redshift_list, calc_snr_list = read_list_spectra(CONFIG_FILE, ["NORM SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR"]) 

########################################## VARIABLES ##########################################

# loops over each spectra from a specified starting and ending point
output_folder = os.path.join(os.getcwd(), "OUTPUT_FILES_PDF")

output_morphed_pdf = os.path.join(output_folder, "spectra_output_morphed.pdf")
output_og_pdf = os.path.join(output_folder, "spectra_output_og.pdf")
xlims = -70000, 0

# Modes can be "morphed", "og", or "both"
MODE = "both"
show_plot = True

########################################## FUNCTIONS ##########################################

def Spectra_Comparison_Generate_PDF(output_path, plot_function, xlims, show_plot = False):
    
    with PdfPages(output_path) as pdf:
        for spectra_index in range(STARTS_FROM, ENDS_AT + 1):
            #closes the global storage of figures
            plt.close('all')
    
            plt.figure()
            
            #keeping track of how many spectra do not have data points within the velocity limits
            norm_spectrum_file_name = norm_spectra_list[spectra_index - 1]
            #gets the cwd and and names each files corectly
            File = os.getcwd() + "/Recons_Hiremath2025/" + str(norm_spectrum_file_name)
    
            #Plot of x-ray selected quasar spectra comparing the Morphed, Normalized, and Reconstruction spectra. Found in spectra_comparison.py
            fig = plot_function(File, xlims)
            
            #Saves the current figure into a page of the pdf
            pdf.savefig(fig)
            if show_plot == True:
                plt.show(fig)
            plt.close(fig)
            
if MODE == "morphed":
    Spectra_Comparison_Generate_PDF(output_morphed_pdf, Plot_spec_compare_morphed, xlims, show_plot)
    
elif MODE == "og":
    Spectra_Comparison_Generate_PDF(output_og_pdf, Plot_spec_compare_og, xlims, show_plot)

elif MODE == "both":
    Spectra_Comparison_Generate_PDF(output_morphed_pdf, Plot_spec_compare_morphed, xlims, show_plot)
    Spectra_Comparison_Generate_PDF(output_og_pdf, Plot_spec_compare_og, xlims, show_plot)