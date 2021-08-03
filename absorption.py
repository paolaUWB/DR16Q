"""
=============
absorption.py
=============

@author Wendy Garcia Naranjo, Mikel Charles, Nathnael Kahassai, Michael Parker
based on code prepared by Abdul Khatri and Paola Rodriguez Hidalgo

Creates figure to visually inspect the absorption. Calculates absorption parameters 
(BALnicity Index [BI], vmin, vmax, equivalent width and depth) for a list of spectra.

Notes
-----
Usage:
    python absorption_code.py spectra_data_list.csv (the default is ...).

Input file:
    This program takes a CSV file with the format ``spectrum_name``, ``z``, ``snr``.

Parameters
----------
    Spectra base path (path to where spectra are stored on disk).
"""

###############################################################################################################################
########################################## IMPORTS ############################################################################
import os
import sys
import numpy as np 
import math
from numpy.lib.function_base import append
from matplotlib.backends.backend_pdf import PdfPages
from utility_functions import clear_file, read_list_spectra, read_spectra
from data_types import Range
from abs_plot import draw_abs_figure
from basic_absorption_parameters import smooth, absorption_parameters_with_plot

###############################################################################################################################
################################ IGONORE: TESTING OUTPUT WITH DR9Q FILES ######################################################
# defining the config file
CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/test_absorption/EHVOcases_updatedredshift.csv" # testing

# sets the directory to find the normalized data files
SPEC_DIREC = os.getcwd() + "/test_absorption/EHVOnorm/" # testing

#BI_INDEX_LIMIT should be 1000 to get accurate results for testing

# be sure to uncomment this and comment out CONFIG_FILE and SPEC_DIREC

###############################################################################################################################
############################## CHANGEABLE VARIABLES ###########################################################################

# input which data release you are working with [input the number as a string i.e. '9']
DR = '16'
'''
#defining the config file
CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/OUTPUT_FILES/NORMALIZATION/good_fit.csv" 

# sets the directory to find the normalized data files
SPEC_DIREC = os.getcwd() + "/DATA/NORM_DR" + DR + "Q/" 
'''
# testing SPEC_DIREC
#SPEC_DIREC = "/Users/wendygarcia/Documents/GitHub/DATA/NORM_DR16Q_HIDE/"

# creates directory for output files
OUT_DIREC = os.getcwd() + "/OUTPUT_FILES/ABSORPTION/"

# do you want to use smoothed norm flux/error
# boxcar_size must always be an odd integer
want_to_smooth = 'no' 
boxcar_size = 3

# plot all cases or only those with absorption
# and provide text file for all cases or only those with absorption 
# yes for everything, no for only absorption
all_plot_and_text = 'yes'

# lower limit of absorption width to be flagged 
BALNICITY_INDEX_LIMIT = 1000

# limits on velocity     min,   max
VELOCITY_LIMIT = Range(-30000, -60000.)

# range of spectra you are working with from the good_normalization.csv file
STARTS_FROM, ENDS_AT = 1, 10

# what percentage value you want to go below the continuum
percent = 0.9

###############################################################################################################################
######################################## OUTPUT FILES #########################################################################

# set name of output .txt file with absorption values
ABSORPTION_VALUES = OUT_DIREC + "/" + "absorption_measurements_test.txt"

# set name of output pdf with plots 
ABSORPTION_OUTPUT_PLOT_PDF = PdfPages('absorption_BI' + str(BALNICITY_INDEX_LIMIT) + '_test.pdf') 

###############################################################################################################################
######################################### MAIN CODE ###########################################################################

# clear files
if __name__ == "__main__":
    clear_file(ABSORPTION_VALUES)

# read list of normalized spectra, zem, and calculated snr from csv file (in this case good_normalization.csv)
# and set variable name to each value
norm_spectra_list, redshift_list, calc_snr_list = read_list_spectra(CONFIG_FILE, ["NORM SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR"]) 
vlast = []
abs_count = 0
abs_count_all = 0
all_count = 0

# loops over each spectra from a specified starting and ending point
for spectra_index in range(STARTS_FROM, ENDS_AT + 1):

    # rounding the numbers of the redshift, calculated snr and setting the norm file name to the current file name from the csv
    z = round(redshift_list[spectra_index - 1], 5)
    calc_snr = round(calc_snr_list[spectra_index - 1], 5)
    current_spectrum_file_name = norm_spectra_list[spectra_index - 1]

    # from the norm spectra name retrieving it's wavelength, normalized flux, and normalized error (in this case from NORM_DRXQ)
    print(str(spectra_index), "current spectra file name:", current_spectrum_file_name)
    current_spectra_data = np.loadtxt(SPEC_DIREC + current_spectrum_file_name)

    # setting a variable for each of those values
    wavelength, normalized_flux, normalized_error = read_spectra(current_spectra_data)

    # smoothing the spectra based on whether the user wants it or not
    if want_to_smooth == 'yes':
        normalized_flux = smooth(normalized_flux, boxcar_size)
        normalized_error = smooth(normalized_error, boxcar_size) / math.sqrt(boxcar_size)

    # getting various BI-related values from the absorption_parameters_with_plot function
    BI_total, BI_individual, BI_all, vmins, vmaxs, EW_individual, final_depth_individual, final_depth_all_individual, beta = absorption_parameters_with_plot(z, wavelength, normalized_flux, BALNICITY_INDEX_LIMIT, VELOCITY_LIMIT, percent)

    ############################# putting things into a text file or plot #######################################
    if (all_plot_and_text == 'yes'):
        all_count += 1
        if (len(vmaxs) != 0):
            text = [f"{str(all_count)}: {current_spectrum_file_name}",
                    f"BI ({VELOCITY_LIMIT.start} > v > {VELOCITY_LIMIT.end}): {BI_total}",
                    f"vmins: {vmins}",
                    f"vmaxs: {vmaxs}",
                    f"BI_individual: {BI_individual}",
                    f"EW_individual: {EW_individual}",
                    f"Depth: {final_depth_individual}"]
            vlast.extend(['\n'.join(text), '\n'])

        draw_abs_figure(all_count, beta, normalized_flux, normalized_error, ABSORPTION_OUTPUT_PLOT_PDF, current_spectrum_file_name, z, calc_snr)
    
    else: 
        if (len(vmaxs) != 0):
            abs_count += 1
            text = [f"{str(abs_count)}: {current_spectrum_file_name}",
                    f"BI ({VELOCITY_LIMIT.start} > v > {VELOCITY_LIMIT.end}): {BI_total}",
                    f"vmins: {vmins}",
                    f"vmaxs: {vmaxs}",
                    f"BI_individual: {BI_individual}",
                    f"EW_individual: {EW_individual}",
                    f"Depth: {final_depth_individual}"]
            vlast.extend(['\n'.join(text), '\n'])
            draw_abs_figure(abs_count, beta, normalized_flux, normalized_error, ABSORPTION_OUTPUT_PLOT_PDF, current_spectrum_file_name, z, calc_snr)

    #####################################################################################################################
    
    final_depth_all_individual.append(final_depth_individual)

    #if (len(vmaxs) != 0) or (all_plot_and_text == 'yes'):
        #vmins_all.append(vmins)
        #vmaxs_all.append(vmaxs)

BI_all= np.array(BI_all)

vmins = np.array(vmins)
vmaxs = np.array(vmaxs)

ABSORPTION_OUTPUT_PLOT_PDF.close()

vmins_final, vmaxs_final = [], []

'''
# creating list of all vmaxs
for loop in range(0, len(vmaxs_all)):
    vmaxs_final.append(str(vmaxs_all[loop])+ ',' )

# creating list of all vins
for loop2 in range(0, len(vmins_all)):
    vmins_final.append(str(vmins_all[loop2])+ ',' ) 
                    
vmaxs_final = np.array(vmaxs_final)
vmins_final = np.array(vmins_final)
'''

np.savetxt(ABSORPTION_VALUES, vlast, fmt='%s')