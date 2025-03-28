"""
absorption
==========

@author Wendy Garcia Naranjo, Mikel Charles, Nathnael Kahassai, Michael Parker 
based on code prepared by Abdul Khatri and Paola Rodriguez Hidalgo

Short description:
    Creates text file and plot for BI.

Extended description:
    Loops through a list of spectra and uses absorption_parameters_with_plot from basic_absorption_parameters.py to receive: 
    BI_total, vmins, vmaxs, BI_individual, EW_individual, final_depth_individual values and a plot is created to show where 
    CIV, CII, and OI would be *if* the EHVO absorption found was due to SiIV. From those values a plot and text file 
    are created and saved.

Input file:
    This program takes a CSV file with the format ``spectrum_name``, ``z``, ``snr``. In this paticular case it is good_fit.csv. 
    From those spectra names it reads in the wavelength, flux, and error from NORM_DR16Q.
"""

###############################################################################################################################
########################################## IMPORTS ############################################################################
import os
import sys
import numpy as np 
import math
from numpy.lib.function_base import append
from matplotlib.backends.backend_pdf import PdfPages
sys.path.insert(0, os.path.dirname(os.getcwd())) # changes the directory to start at DR16Q --> all paths after this will need to be written as if this was in the top level of the DR16Q
print(os.getcwd())
from utility_functions import clear_file, read_list_spectra, read_spectra, append_row_to_csv
from data_types import Range
from abs_function_module import smooth, abs_parameters_plot_optional
from abs_plot_module import draw_abs_figure

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

# limits on velocity     min,   max
VELOCITY_LIMIT = Range(-30000, -60000.)

# range of spectra you are working with from the good_fit.csv file
STARTS_FROM, ENDS_AT = 1, 98

# what percentage value you want to go below the continuum
percent = 0.9

# whether you want to output a csv table of your run
want_csv = 'yes'

###############################################################################################################################
######################################## OUTPUT FILES #########################################################################

# set name of output .txt file with absorption values
ABSORPTION_VALUES = OUT_DIREC + "/" + 'BI' + str(BALNICITY_INDEX_LIMIT) + '.txt'

# set name of output pdf with plots 
ABSORPTION_OUTPUT_PLOT_PDF = PdfPages(OUT_DIREC + 'BI' + str(BALNICITY_INDEX_LIMIT) + '.pdf') 

ABSORPTION_TABLE = OUT_DIREC + 'absorption_table.csv'

###############################################################################################################################
######################################### MAIN CODE ###########################################################################

# clear files
if __name__ == "__main__":
    clear_file(ABSORPTION_VALUES)
    if (want_csv == 'yes'):
        clear_file(ABSORPTION_TABLE)

# read list of normalized spectra, zem, and calculated snr from csv file (in this case good_normalization.csv)
# and set variable name to each value
norm_spectra_list, redshift_list, calc_snr_list = read_list_spectra(CONFIG_FILE, ["NORM SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR"]) 
vlast = []
# whether abs_count or all_count is used is based on the value of all_plot_and_text
abs_count = 0 # counter for amount of spectra that have absorption when all_plot_and_text = no and for text files
all_count = 0 # counter for all spectra ran when all_plot_and_text = yes

if (want_csv == 'yes'):
    field = ['NORM SPECTRA FILE NAME','BI TOTAL','BI INDIVIDUAL','VMINS', 'VMAXS', 'EW INDIVIDUAL', 'DEPTH']
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

    # getting various BI-related values from the absorption_parameters_with_plot function
    BI_total, BI_individual, BI_all, vmins, vmaxs, EW_individual, final_depth_individual, final_depth_all_individual, beta, vminindex_for_range, vmaxindex_for_range = abs_parameters_plot_optional(
        z, wavelength, normalized_flux, BALNICITY_INDEX_LIMIT, VELOCITY_LIMIT, percent)

    max_peak = np.max(normalized_flux[vmaxindex_for_range + 1 : vminindex_for_range + 1])

    ############################# putting things into a text file or plot #######################################
    fields = [norm_spectrum_file_name, BI_total, BI_individual, vmins, vmaxs, EW_individual, final_depth_individual]
    
    if (all_plot_and_text == 'yes'): # plot all is yes, graph everything but only text file for when abs is found
        all_count += 1
        if (len(vmaxs) != 0): # text file created only when absorption is found
            abs_count += 1
            text = [f"{str(abs_count)} abs | {str(all_count)} tot",
                    f"{norm_spectrum_file_name}",
                    f"BI ({VELOCITY_LIMIT.start} > v > {VELOCITY_LIMIT.end}): {BI_total}",
                    f"vmins: {vmins}",
                    f"vmaxs: {vmaxs}",
                    f"BI_individual: {BI_individual}",
                    f"EW_individual: {EW_individual}",
                    f"Depth: {final_depth_individual}"]
            vlast.extend(['\n'.join(text), '\n'])
            abs = abs_count
        else: # create graph no matter what
            abs = 'no'
        draw_abs_figure(
            abs, all_count, beta, normalized_flux, normalized_error, ABSORPTION_OUTPUT_PLOT_PDF, norm_spectrum_file_name, z, calc_snr, max_peak)
        # whether you want to create a master csv table or not
        if (want_csv == 'yes'):
            append_row_to_csv(ABSORPTION_TABLE, fields)  
    else: # plot all is no and only create text file and graph of cases where absorption is found
        all_count += 1
        if (len(vmaxs) != 0):
            abs_count += 1
            text = [f"{str(abs_count)} abs | {str(all_count)} tot",
                    f"{norm_spectrum_file_name}",
                    f"BI ({VELOCITY_LIMIT.start} > v > {VELOCITY_LIMIT.end}): {BI_total}",
                    f"vmins: {vmins}",
                    f"vmaxs: {vmaxs}",
                    f"BI_individual: {BI_individual}",
                    f"EW_individual: {EW_individual}",
                    f"Depth: {final_depth_individual}"]
            vlast.extend(['\n'.join(text), '\n'])
            abs = abs_count
            draw_abs_figure(
                abs_count, all_count, beta, normalized_flux, normalized_error, ABSORPTION_OUTPUT_PLOT_PDF, norm_spectrum_file_name, z, calc_snr, max_peak)
        
        # whether you want to create a master csv table or not
        if (want_csv == 'yes'):
            append_row_to_csv(ABSORPTION_TABLE, fields)  
    #####################################################################################################################
    
    final_depth_all_individual.append(final_depth_individual)

    # testing
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

