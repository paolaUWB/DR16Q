#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 11:39:37 2026

@author: lilianaflores
"""


###############################################################################################################################
########################################## IMPORTS ############################################################################
import os
import sys
import numpy as np 
import math
#from numpy.lib.function_base import append #Remove: Unused in program and incompatible with updated version of numpy
from matplotlib.backends.backend_pdf import PdfPages
sys.path.insert(0, os.getcwd() + '/../' ) # changes the directory to the DR16Q --> all paths after this will need to be written as if this was in the top level of the DR16Q
from utility_functions import clear_file, read_list_spectra, read_spectra, append_row_to_csv
from data_types import Range
from abs_function_module import smooth, abs_parameters_plot_optional
from abs_plot_module import draw_abs_figure
from abs_function_module import wavelength_to_velocity
import matplotlib.py as plt

'''
###############################################################################################################################
################################ IGONORE: TESTING OUTPUT WITH DR9Q FILES ######################################################
# defining the config file
CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/test_absorption/EHVOcases_updatedredshift.csv" # testing

# sets the directory to find the normalized data files
SPEC_DIREC = os.getcwd() + "/test_absorption/EHVOnorm/" # testing

#BI_INDEX_LIMIT should be 1000 to get accurate results for testing

# be sure to uncomment this and comment out CONFIG_FILE and SPEC_DIREC
'''

###############################################################################################################################
############################## CHANGEABLE VARIABLES ###########################################################################

#defining the config file
CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/../DR16Q_EHVO/good_fit_EHVO.csv" #"/OUTPUT_FILES/NORMALIZATION/good_fit_EHVO.csv" #good_fit_EHVO.csv" ##_newSNR_flagged_but_ok.csv #_EHVO.csv" 

# directory of where normalized data files are
# data NOT on github but local computer
#NORM_DIREC = os.getcwd() + '/../' + "NORM_DR16Q/"

NORM_DIREC = os.getcwd() + "/../" + "/DR16Q_EHVO/NORM_DR16Q_EHVO/"

# creates directory for output files
OUT_DIREC = os.getcwd() + "/OUTPUT_FILES/"

# do you want to use smoothed norm flux/error
# boxcar_size must always be an odd integer
want_to_smooth = 'no' 
boxcar_size = 11

# plot all cases or only those with absorption
# and provide text file for all cases or only those with absorption 
# yes for everything, no for only absorption
all_plot_and_text = 'yes'

# lower limit of absorption width to be flagged 
BALNICITY_INDEX_LIMIT = 2000

#xlimits for plotting
xlow = None #No value -> VELOCITY_LIMIT range is used +10,000 on either end to plot xlimits
xhigh = None

#xlow = -80000
#xhigh = -10000

# limits on velocity     min,   max
VELOCITY_LIMIT = Range(-30000, -60000.)

# range of spectra you are working with from the good_fit.csv file
STARTS_FROM, ENDS_AT = 1, 98 #Note that the end is inclusive

# what percentage value you want to go below the continuum
percent = 0.9

# whether you want to output a csv table of your run
want_csv = 'yes'

# Do you want to use a specific reference wavelength?
# data from Verner table
wavelength_CIV_emit1 = 1548.1950
wavelength_CIV_emit2 = 1550.7700
wavelength_SiIV_emit1 = 1393.755
wavelength_SiIV_emit2 = 1402.770
avr_CIV_doublet = 1549.0524 # weighted average
avr_SiIV_doublet = 1396.747 # weighted average

ref_wavelength = avr_CIV_doublet
#ref_wavelength = wavelength_CIV_emit1

###############################################################################################################################
######################################## OUTPUT FILES #########################################################################

# set name of output .txt file with absorption values
#ABSORPTION_VALUES = OUT_DIREC + "/" + 'BI' + str(BALNICITY_INDEX_LIMIT) + '.txt'

# set name of output pdf with plots 
ABSORPTION_OUTPUT_PLOT_zoom = PdfPages(OUT_DIREC + 'BI' + str(BALNICITY_INDEX_LIMIT) + '.pdf') 

ABSORPTION_TABLE = OUT_DIREC + 'absorption_mask_ranges.csv'

###############################################################################################################################
######################################### MAIN CODE ###########################################################################

# clear files
if __name__ == "__main__":
    if (want_csv == 'yes'):
        clear_file(ABSORPTION_TABLE)

# read list of normalized spectra, zem, and calculated snr from csv file (in this case good_normalization.csv)
# and set variable name to each value
norm_spectra_list, redshift_list, calc_snr_list, depth_flag = read_list_spectra(CONFIG_FILE, ["NORM SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR", "NEEDS RECALCULATION"]) 
vlast = []
# whether abs_count or all_count is used is based on the value of all_plot_and_text
abs_count = 0 # counter for amount of spectra that have absorption when all_plot_and_text = no and for text files
all_count = 0 # counter for all spectra ran when all_plot_and_text = yes

if (want_csv == 'yes'):
    field = ['NORM SPECTRA FILE NAME',"REDSHIFT", "CALCULATED SNR", "NEEDS RECALCULATION", "MASK RANGES"]
    append_row_to_csv(ABSORPTION_TABLE, field)

# loops over each spectra from a specified starting and ending point
for spectra_index in range(STARTS_FROM, ENDS_AT + 1):

    # rounding the numbers of the redshift, calculated snr and setting the norm file name to the current file name from the csv
    z = round(redshift_list[spectra_index - 1], 5)
    calc_snr = round(calc_snr_list[spectra_index - 1], 5)
    norm_spectrum_file_name = norm_spectra_list[spectra_index - 1]

    # from the norm spectra name retrieving it's wavelength, normalized flux, and normalized error (in this case from NORM_DRXQ)
    print(str(spectra_index), "current spectra file name:", norm_spectrum_file_name)
    norm_spectra_data = np.loadtxt(NORM_DIREC + norm_spectrum_file_name)

    # setting a variable for each of those values from the spectra
    wavelength, normalized_flux, normalized_error = read_spectra(norm_spectra_data)

    # smoothing the flux and error based on what the user wants (yes or no)
    if want_to_smooth == 'yes':
        normalized_flux = smooth(normalized_flux, boxcar_size)
        normalized_error = smooth(normalized_error, boxcar_size) / math.sqrt(boxcar_size)
        
    
    beta = wavelength_to_velocity(z, wavelength)
    
    if depth_flag == 'Y':
        
        plt.plot(beta, wavelength)
        plt.show()
        plt.close()
        
        
    # then some sort of user input function for identifying masked regions