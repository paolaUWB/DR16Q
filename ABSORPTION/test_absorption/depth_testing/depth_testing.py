"""
=============
depth_testing.py
=============

@author Wendy Garcia Naranjo

Short description:
    Testing depth value equation through using various boxcar size values and plotting them. Uses absorption.py as a base.

Input file:
    This program takes a CSV file with the format ``spectrum_name``, ``z``, ``snr``. In this paticular case it is depth_test.csv. 
    From those spectra names it reads in the wavelength, flux, and error from NORM_DR16Q.
"""

###############################################################################################################################
########################################## IMPORTS ############################################################################
import os
import sys
import numpy as np 
from matplotlib.backends.backend_pdf import PdfPages
from utility_functions import clear_file, read_list_spectra, read_spectra
from data_types import Range
from abs_function_module import abs_parameters_plot_optional
from depth_test_calc_graph import depth_figure, depth_testing_calculation, new_smooth_normalzied_flux_value
###############################################################################################################################
############################## CHANGEABLE VARIABLES ###########################################################################

#defining the config file
CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/test_absorption/depth_testing/depth_test.csv"

# directory of where normalized data files are
# data NOT on github but local computer
NORM_DIREC = os.getcwd() + '/../' + "NORM_DR16Q/"

# creates directory for output files
OUT_DIREC = os.getcwd() + "/test_absorption/depth_testing"

# boxcar_size must always be an odd integer
# orginal spectra graphs by default
# the same spectra with boxcar_size_first value is graphed with smooth function
# same thing for boxcar_size_second
boxcar_size_all_values = [21, 35, 51]

# lower limit of absorption width to be flagged 
BALNICITY_INDEX_LIMIT = 2000

# limits on velocity     min,   max
VELOCITY_LIMIT = Range(-30000, -60000.)

# range of spectra you are working with from the depth_test.csv file
STARTS_FROM, ENDS_AT = 1, 10

# what percentage value you want to go below the continuum
percent = 0.9

###############################################################################################################################
######################################## OUTPUT FILES #########################################################################

# set name of output .txt file with absorption values
ABSORPTION_VALUES = OUT_DIREC + "/" + "depth_test.txt"

# set name of output pdf with plots 
ABSORPTION_OUTPUT_PLOT_PDF = PdfPages('depth_test.pdf') 

###############################################################################################################################
######################################### MAIN CODE ###########################################################################

# clear files
if __name__ == "__main__":
    clear_file(ABSORPTION_VALUES)

# read list of normalized spectra, zem, and calculated snr from csv file (in this case good_normalization.csv)
# and set variable name to each value
norm_spectra_list, redshift_list, calc_snr_list = read_list_spectra(CONFIG_FILE, ["NORM SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR"]) 
vlast = []

all_count = 0 # counter for all spectra ran

# loops over each spectra from a specified starting and ending point
for spectra_index in range(STARTS_FROM, ENDS_AT + 1):

    # rounding the numbers of the redshift, calculated snr and setting the norm file name to the current file name from the csv
    z = round(redshift_list[spectra_index - 1], 5)
    calc_snr = round(calc_snr_list[spectra_index - 1], 5)
    norm_spectrum_file_name = norm_spectra_list[spectra_index - 1]

    # from the norm spectra name retrieving it's wavelength, normalized flux, and normalized error (in this case from NORM_DR16Q)
    print(str(spectra_index), "current spectra file name:", norm_spectrum_file_name)
    norm_spectra_data = np.loadtxt(NORM_DIREC + norm_spectrum_file_name)

    # setting a variable for each of those values from the spectra
    wavelength, normalized_flux, normalized_error = read_spectra(norm_spectra_data)

    # getting various BI-related values from the absorption_parameters_with_plot function
    BI_total, BI_individual, BI_all, vmins, vmaxs, EW_individual, final_depth_individual, final_depth_all_individual, beta, vminindex_for_range, vmaxindex_for_range = abs_parameters_plot_optional(
        z, wavelength, normalized_flux, BALNICITY_INDEX_LIMIT, VELOCITY_LIMIT, percent)

    # getting the smoothed normalized flux value based on the snr ratio
    smooth_norm_flux_value, number = new_smooth_normalzied_flux_value(normalized_flux, calc_snr, boxcar_size_all_values)

    ################################################# SMOOTH SPECTRA ##########################################
    BI_total_sm, BI_individual_sm, BI_all_sm, vmins_sm, vmaxs_sm, EW_individual_sm, final_depth_individual_sm, final_depth_all_individual_sm, beta_sm, vminindex_for_range_sm, vmaxindex_for_range_sm = abs_parameters_plot_optional(
        z, wavelength, smooth_norm_flux_value, BALNICITY_INDEX_LIMIT, VELOCITY_LIMIT, percent, plots = "no")
    ##########################################################################################################

    # calculating seven point avergae using smoothed normalized flux value chosen from new_smooth_normalzied_flux_value function
    final_depth_seven_point_average = depth_testing_calculation(smooth_norm_flux_value, vminindex_for_range_sm, vmaxindex_for_range_sm,)

    max_peak = np.max(normalized_flux[vmaxindex_for_range + 1 : vminindex_for_range + 1])

    ############################# putting things into a text file or plot #######################################
    all_count += 1
    text = [f"{str(all_count)} tot",
            f"{norm_spectrum_file_name}",
            f"BI ({VELOCITY_LIMIT.start} > v > {VELOCITY_LIMIT.end}): {BI_total}",
            f"vmins: {vmins}",
            f"vmaxs: {vmaxs}",
            f"BI_individual: {BI_individual}",
            f"EW_individual: {EW_individual}",
            f"Original Graph Depth: {final_depth_individual}",
            f"Smooth Graph Depth: {final_depth_individual_sm}",
            f"Depth Seven point average: {final_depth_seven_point_average}"]
    vlast.extend(['\n'.join(text), '\n'])
    depth_figure(all_count, beta, normalized_flux, normalized_error, ABSORPTION_OUTPUT_PLOT_PDF, norm_spectrum_file_name, z, calc_snr, max_peak, smooth_norm_flux_value, number)

    #####################################################################################################################

ABSORPTION_OUTPUT_PLOT_PDF.close()
np.savetxt(ABSORPTION_VALUES, vlast, fmt='%s')