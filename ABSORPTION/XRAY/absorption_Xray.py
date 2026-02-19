"""
absorption
==========

@author Wendy Garcia Naranjo, Mikel Charles, Nathnael Kahassai, Michael Parker 
based on code prepared by Abdul Khatri and Paola Rodriguez Hidalgo

This is an EDITED version of absorption.py to accomodate xray project data.
Search "EDITED" to see changes made - LEF

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
#from numpy.lib.function_base import append
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
sys.path.insert(0, os.getcwd() + '/../../') # changes the directory to the DR16Q --> all paths after this will need to be written as if this was in the top level of the DR16Q
from utility_functions import clear_file, read_list_spectra, read_spectra, append_row_to_csv
from data_types import Range
sys.path.insert(0, os.getcwd() + '/../') # changes the directory to the DR16Q --> all paths after this will need to be written as if this was in the top level of the DR16Q
from abs_function_module import smooth, abs_parameters_plot_optional, wavelength_to_velocity
from abs_plot_module import draw_abs_figure
import matplotlib.pyplot as plt


###############################################################################################################################
############################## CHANGEABLE VARIABLES ###########################################################################

#defining the config file
CONFIG_FILE1 = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()+"/spec_lists/xray_parent_list.csv" #2281 full parent sample
CONFIG_FILE2 = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()+"/spec_lists/xray_list_SNR10_z1.9.csv" #250 with SNR>10 & z>1.9
CONFIG_FILE3 = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()+"/spec_lists/xray_list_SNR10.csv" #268 with SNR>10

CONFIG_FILE = CONFIG_FILE2


# directory of where normalized data files are
# data NOT on github but local computer
NORM_DIREC = os.getcwd() + "/Recons_Hiremath2025/"


# creates directory for output files
OUT_DIREC = os.getcwd() + "/OUTPUT_FILES/"



# do you want to use smoothed norm flux/error
# boxcar_size must always be an odd integer
want_to_smooth = 'no' 
boxcar_size = 5

# plot all cases or only those with absorption
# and provide text file for all cases or only those with absorption 
# yes for everything, no for only absorption
all_plot_and_text = 'no'

# lower limit of absorption width to be flagged 
BALNICITY_INDEX_LIMIT = 2000

xlow = -70000
xhigh = -20000

# limits on velocity
VELOCITY_LIMIT = Range(-30000, -60000)

# range of spectra you are working with from the good_fit.csv file
if CONFIG_FILE == CONFIG_FILE1:
    STARTS_FROM, ENDS_AT = 1, 2281
elif CONFIG_FILE == CONFIG_FILE2:
    STARTS_FROM, ENDS_AT = 1, 250
elif CONFIG_FILE == CONFIG_FILE3:
    STARTS_FROM, ENDS_AT = 1, 268
#STARTS_FROM, ENDS_AT = 1, 268 #uncomment to override if statement auto selection of range

# what percentage value you want to go below the continuum
percent = 0.9

# whether you want to output a csv table of your run
want_csv = 'yes'

###############################################################################################################################
######################################## OUTPUT FILES #########################################################################
if CONFIG_FILE == CONFIG_FILE1:
    xtra_name = 'parent'
elif CONFIG_FILE == CONFIG_FILE2:
    xtra_name = 'z_SNR'
elif CONFIG_FILE == CONFIG_FILE3:
    xtra_name = 'SNR'
    
if want_to_smooth == 'yes':
    smooth_name = '_smooth'
elif want_to_smooth == 'no':
    smooth_name = ''

# set name of output .txt file with absorption values
ABSORPTION_VALUES = OUT_DIREC + "/" + 'BI' + str(BALNICITY_INDEX_LIMIT) + '_P0.9_' + xtra_name + '.txt'

# set name of output pdf with plots 
ABSORPTION_OUTPUT_PLOT_PDF = PdfPages(OUT_DIREC + 'BI' + str(BALNICITY_INDEX_LIMIT) + '_P0.9_' + xtra_name + smooth_name + '.pdf') 

ABSORPTION_TABLE = OUT_DIREC + 'absorption_table_P0.9_' + xtra_name + '.csv'

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
    #........................................................................................................

    # smoothing the flux and error based on what the user wants (yes or no)
    if want_to_smooth == 'yes':
        normalized_flux = smooth(normalized_flux, boxcar_size)
        normalized_error = smooth(normalized_error, boxcar_size) / math.sqrt(boxcar_size)

    # getting various BI-related values from the absorption_parameters_with_plot function
    BI_total, BI_individual, BI_all, vmins, vmaxs, EW_individual, final_depth_individual, final_depth_all_individual, beta, vminindex_for_range, vmaxindex_for_range = abs_parameters_plot_optional(
        z, wavelength, normalized_flux, BALNICITY_INDEX_LIMIT, VELOCITY_LIMIT, percent=percent)

    #........................................................................................................
    # EDITED: I rewrote this section so that if there is no spectra in between the searching VELOCITY_LIMITS there is still a max_peak value to pass through
    #max_peak = np.max(normalized_flux[vmaxindex_for_range + 1 : vminindex_for_range + 1])
    
    
    if np.size(normalized_flux[vmaxindex_for_range + 1 : vminindex_for_range + 1]) > 0:
        max_peak = np.max(normalized_flux[vmaxindex_for_range + 1 : vminindex_for_range + 1])
        
    else:
        max_peak = 3 #this is just used in plotting the y limits if there is no spectra in the limits 
        not_in_range.append(spectra_index)
    #........................................................................................................
    
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
            #EDITED replaced z to z_final to include correct z value in plot but not correct for redshift when plotting
        draw_abs_figure(
            abs, all_count, beta, normalized_flux, normalized_error, ABSORPTION_OUTPUT_PLOT_PDF, norm_spectrum_file_name[:-5], z_final, calc_snr, max_peak, VELOCITY_LIMIT, percent, xlow, xhigh)
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
            #EDITED replaced z to z_final to include correct z value in plot but not correct for redshift when plotting
            draw_abs_figure(
                abs_count, all_count, beta, normalized_flux, normalized_error, ABSORPTION_OUTPUT_PLOT_PDF, norm_spectrum_file_name[:-5], z_final, calc_snr, max_peak, VELOCITY_LIMIT, percent, xlow, xhigh)
        
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


np.savetxt(ABSORPTION_VALUES, vlast, fmt='%s')

