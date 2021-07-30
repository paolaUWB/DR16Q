"""
==============================
basic_absorption_parameters.py
==============================

@author Wendy Garcia Naranjo, Mikel Charles, Nathnael Kahassai, Michael Parker
based on code prepared by Abdul Khatri and Paola Rodriguez Hidalgo

Calculates absorption parameters 
(BALnicity Index BI, vmin and vmax, equivalent width, and depth) for a list of spectra.

Notes
-----
Input file:
    This program takes a CSV file with the format ``spectrum_name``, ``z``, ``snr``.

 Outputs:    
    In case the naming gets confusing, these are list outputs:
    BI_all 
        Is the array of all the total BALnicity indices for all the spectra (not to be confused with BI_total, which doesn't contain *all* the BI values)
    BI_all_individual
        Is the list containing all BI values with some spectra having sublists of multiple BI (due to multiple BALs in a single spectrum)
    vmins_final and vmaxs_final 
        Are the arrays of lower and upper bounds for each trough (same number of elements and subelements as BI_all_individual, EW_all individual, and final_depth_all_individual)
    EW_all_individual
        Is the list of all equivalent widths for each BAL (sublists for multiple BALs)
    final_depth_all_individual 
        Is the list of all the depths (and sublists for multiple BALs)
    

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
from matplotlib import pyplot as plt
from scipy import signal
from numpy.lib.function_base import append
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
from utility_functions import print_to_file, clear_file, read_list_spectra, read_spectra, wavelength_to_velocity
from data_types import Range


######################################## TESTING OUTPUT WITH DR9Q FILES #######################################################
# defining the config file
CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/test_absorption/EHVOcases_updatedredshift.csv" # testing

# sets the directory to find the normalized data files
SPEC_DIREC = os.getcwd() + "/test_absorption/EHVOnorm/" # testing

#BI_INDEX_LIMIT should be 1000 to get accurate results for testing

###############################################################################################################################
############################## CHANGEABLE VARIABLES ###########################################################################

# input which data release you are working with [input the number as a string i.e. '9']
DR = '16'

'''
#defining the config file
CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/OUTPUT_FILES/NORMALIZATION/good_normalization.csv" 

# sets the directory to find the normalized data files
SPEC_DIREC = os.getcwd() + "/DATA/NORM_DR" + DR + "Q/" 
'''

# creates directory for output files
OUT_DIREC = os.getcwd() + "/OUTPUT_FILES/ABSORPTION/"

# do you want to use smoothed norm flux/error
# boxcar_size must always be an odd integer
want_to_smooth = 'no' 
boxcar_size = 101 

# plot all cases or only those with absorption
# and provide text file for all cases or only those with absorption 
all_plot_and_text = 'yes'

# lower limit of absorption width to be flagged 
BALNICITY_INDEX_LIMIT = 1000 

# limits on velocity     min,   max
VELOCITY_LIMIT = Range(-30000, -60000.)

# range of spectra you are working with from the good_normalization.csv file
STARTS_FROM, ENDS_AT = 1, 40

###############################################################################################################################
######################################## OUTPUT FILES #########################################################################

## No longer outputs text file, rather this module is meant to calculate the absorption
## parameters into arrays for use elsewhere. 

#________________________________________________________
# set name of output .txt file with absorption values
ABSORPTION_VALUES = OUT_DIREC + "/" + "absorption_measurements_test.txt"


###############################################################################################################################
######################################### FUNCTION(S) #########################################################################

def smooth(norm_flux, box_size):   
    """Function: 

    Parameters:
    -----------
    norm_flux : 
        Normalized flux to be smoothed.

    box_size: int
        This is the number of points that are smoothed into one. Always be sure to use an odd 
        number, because we need the same amount of points on each side of the data point to be
        smoothed.

    Returns:
    --------
    y_smooth
    """    
    y_smooth = signal.savgol_filter(norm_flux,box_size,2)
    return y_smooth

######################################### VARIABLES ###########################################################################

brac_all = []
vmins, vmaxs, vmins_all, vmaxs_all, delta_v_all, vlast = [], [], [], [], [], [] # v = velocity
final_depth_individual, final_depth_all_individual = [], []
BI_all, BI_total, BI_ind_sum, BI_individual, BI_all_individual, BI_ind, BI_mid = [], [], [], [], [], [], []
EW_individual, EW_ind, EW_all_individual = [], [], [] #EW = equivalent width

non_trough_count = 999 # arbitrary large number that we will never reach

delta_v = 0 #change in velocity
sum_of_deltas = 0        

count_v = 0 # variable initialization to get into vmin/vmax loop

###############################################################################################################################
######################################### MAIN CODE ###########################################################################

# clear files
if __name__ == "__main__":
    clear_file(ABSORPTION_VALUES)

# read list of normalized spectra, zem, and calculated snr from csv file (in this case good_normalization.csv)
# and set variable name to each value
#norm_spectra_list, redshift_list, calc_snr_list = read_list_spectra(CONFIG_FILE, ["NORM SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR"])
norm_spectra_list, redshift_list, calc_snr_list = read_list_spectra(CONFIG_FILE, ["NORM SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR"]) # testing

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

    # re-nitialize all variables for each spectrum #####################
    vmins, vmaxs = [], []
    BI_mid, BI_individual = [], []
    EW_individual = []
    beta = []
    index_depth_final, flux_depth, final_depth_individual = [], [], []

    non_trough_count = 100
  
    delta_v = 0 #change in velocity
    sum_of_deltas = 0               
    count_v = 0   # variable initialization to get into vmin/vmax loop
    ####################################################################

    # transform the wavelength array to velocity (called "beta") based on the CIV doublet: 
    beta = wavelength_to_velocity(z, wavelength)

    # finding and labeling index of beta that we will be looping through
                                                        # start,  end
                                                        #   min,  max
    if any(beta): # for reference VELOCITY_LIMIT = Range(-30000, -60000.))
        try:
            vmaxindex_for_range = np.max(np.where(beta <= VELOCITY_LIMIT.end)) # index value of VELOCITY_LIMIT.end or closest value
        except:
            vmaxindex_for_range = 0  
    try:
        vminindex_for_range = np.min(np.where(beta >= VELOCITY_LIMIT.start)) # index value of VELOCITY_LIMIT.start or closest value
     
    except:
        vminindex_for_range = np.where(beta == np.min(beta)) 
       
    velocity_range_index = np.arange(vmaxindex_for_range, vminindex_for_range) # from left to right
    velocity_range_index  = np.array(velocity_range_index[::-1])   # from right to left (reversed list)
                                                                   # ^^^^^^^^ 0 to -60000

    # looping through the velocity ranges
    for current_velocity_index in velocity_range_index:
        C = 0 # C will be 0 or 1 and is the C used in the integral for the calculation of BI
        # ([1 - f(v)/0.9] = bracket) > 0 when there is an absorption feature 
        # bracket is the things inside the bracket from the BI integral calculation 
        bracket = (1. - (normalized_flux[current_velocity_index] / 0.9))
        
        # Handle 3-point spike limit ###################################################
        if bracket > 0:
            non_trough_count = 0    
        else:
            non_trough_count += 1
            bracket = 0
            
        if((bracket > 0) or (non_trough_count <= 3)):
            delta_v = beta[current_velocity_index] - beta[current_velocity_index - 1]
            sum_of_deltas += delta_v
            brac_all.append(bracket)
            delta_v_all.append(delta_v)

            EW = bracket * delta_v
            EW = np.round(EW, 5)
            EW_ind.append(EW)   
        ################################################################################## 

            # BI calculation #################################################################################################
            if sum_of_deltas >= BALNICITY_INDEX_LIMIT: # passing the BALNICITY_INDEX_LIMIT (in this case 2,000 km/s) threshold
                C = 1  #set to 1 only if square bracket is continuously positive over a velocity interval            
                BI = (bracket * C) * (delta_v) #Calculate BAL for this delta_v
                BI_mid.append(np.round(BI, 5)) #Append to intermediate results
                BI_ind.append(np.round(BI, 5)) 

                 
                ############################# V MIN CALCULATION ##################################################
                if count_v == 0 and non_trough_count == 0:  
                    vmins_index = np.min(np.where(beta >= (beta[current_velocity_index] + BALNICITY_INDEX_LIMIT))) # vmins occurs current beta plus BALNICITY_INDEX_LIMIT
                    vmins.append(np.round(beta[vmins_index], 5))                    


                ############################################################################################

                    count_v = 1
                
                bracket_1 = (1. - (normalized_flux[current_velocity_index - 1] / 0.9))
                bracket_2 = (1. - (normalized_flux[current_velocity_index - 2] / 0.9))
                bracket_3 = (1. - (normalized_flux[current_velocity_index - 3] / 0.9))
                bracket_4 = (1. - (normalized_flux[current_velocity_index - 4] / 0.9))

                ############################# V MAX CALCULATION ######################################################
                if (((bracket > 0 and bracket_1 < 0 and bracket_2 < 0 and bracket_3 < 0 and bracket_4 < 0 and count_v == 1)) or (current_velocity_index == vmaxindex_for_range)):  
                    vmaxs_index = np.min(np.where (beta >= beta[current_velocity_index]))
                    vmaxs.append(np.round(beta[current_velocity_index], 4))
                 
                    
                
                    ############################################################################################

                    BI_ind_sum = np.round(sum(BI_ind), 2)
                    BI_individual.append(BI_ind_sum) # this array contains one single BI value of each absorption feature in a single spectrum
                    BI_ind = []
                    
                    EW_ind_sum = np.round(sum(EW_ind), 2)
                    EW_individual.append(EW_ind_sum)
                    EW_ind = []
                    
                    # calculating the depth of each individual absorption trough
                    final_depth = np.round((1. - np.min(normalized_flux[vmaxs_index:vmins_index])), 2)
                    final_depth_individual.append(final_depth)
                    
                    count_v = 0 
        
        else: #if the bracket value is not more than zero (so if we don't have absorption feature)
            sum_of_deltas = 0 # this is b/c we do not want to keep counting the width of the absorption feature if it is not wider than the BALnicity_index_limit
            count_v = 0 # this is b/c if the code encounters another absorption feature which is wider than 600km/s, the code is going to go through the if statement on line 242
            EW_ind = []
          
        if current_velocity_index == vmaxindex_for_range:
            BI_total = np.round(sum(BI_mid), 2)         
            BI_all.append(BI_total)    
            BI_all_individual.append(BI_individual)
            EW_all_individual.append(EW_individual)
"""
    ################################################ putting the information into a text file #######################################
    if (len(vmaxs) != 0) or (all_plot_and_text == 'yes'):
        text = [f"{spectra_index}: {current_spectrum_file_name}",
                f"BI ({VELOCITY_LIMIT.start} > v > {VELOCITY_LIMIT.end}): {BI_total}",
                f"vmins: {vmins}",
                f"vmaxs: {vmaxs}",
                f"BI_individual: {BI_individual}",
                f"EW_individual: {EW_individual}",
                f"Depth: {final_depth_individual}"]
        vlast.extend(['\n'.join(text), '\n'])
   ###################################################################################################################################
   Line: final_depth_all_individual.append(final_depth_individual)
    and the if statement after will need to be indented to run the textfile output above
"""
final_depth_all_individual.append(final_depth_individual)

if (len(vmaxs) != 0) or (all_plot_and_text == 'yes'):
    vmins_all.append(vmins)
    vmaxs_all.append(vmaxs)
      
BI_all= np.array(BI_all)

vmins = np.array(vmins)
vmaxs = np.array(vmaxs)

#ABSORPTION_OUTPUT_PLOT_PDF.close()

vmins_final, vmaxs_final = [], []

for loop in range(0, len (vmaxs_all)):
    vmaxs_final.append(str(vmaxs_all[loop])+ ',' )

for loop2 in range(0, len(vmins_all)):
    vmins_final.append(str(vmins_all[loop2])+ ',' ) 
                    
vmaxs_final = np.array(vmaxs_final)
vmins_final = np.array(vmins_final)


#print("BI all individual is ", type(BI_all_individual))
#print("BI all is", type(BI_all))
#print("vmins final is ", type(vmins_final))
#print("vmaxs_final is ", type(vmaxs_final))
#print("EW is ", type(EW_all_individual))
#print("final depth all individual is", type(final_depth_all_individual))

np.savetxt(ABSORPTION_VALUES, vlast, fmt='%s')