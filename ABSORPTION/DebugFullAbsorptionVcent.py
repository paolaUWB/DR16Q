# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 14:05:26 2024

@author: phyfr
"""

import os
import sys
import numpy as np 
import math
from numpy.lib.function_base import append
from matplotlib.backends.backend_pdf import PdfPages

if os.getcwd()[-5:] != 'DR16Q': # Sets directory to the default location for DR16Q
    user = os.getlogin()
    drive = os.getenv('SYSTEMDRIVE')
    path = f"{drive}/Users/{user}"

    os.chdir(path + '/Documents/Github/DR16Q')

sys.path.append(os.getcwd() + '\ABSORPTION')

from utility_functions import clear_file, read_list_spectra, read_spectra, append_row_to_csv
from data_types import Range
from abs_function_moduleModular import smooth, abs_parameters_plot_optional, wavelength_to_velocity, black_line, vmin_line, vmin_plot_IF, span_vmin_vmax, vmax_plot_span_IF
from abs_plot_module import draw_abs_figure
import matplotlib.pyplot as plt
from time import process_time
import warnings

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
CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/VARIABILITY/good_fit_EHVO_variability_for_absorption_code_updatedRankineRedshifts.csv" #"/OUTPUT_FILES/NORMALIZATION/good_fit_EHVO.csv" #good_fit_EHVO.csv" ##_newSNR_flagged_but_ok.csv #_EHVO.csv" 

# directory of where normalized data files are
# data NOT on github but local computer
#NORM_DIREC = os.getcwd() + '/../' + "NORM_DR16Q/"

NORM_DIREC = os.getcwd()  + "/VARIABILITY/Variability_Norms_Unsmoothed/"

# creates directory for output files
OUT_DIREC = os.getcwd() + "/ABSORPTION/OUTPUT_FILES/"

# do you want to use smoothed norm flux/error
# boxcar_size must always be an odd integer
want_to_smooth = 'no' 
boxcar_size = 11

# plot all cases or only those with absorption
# and provide text file for all cases or only those with absorption 
# yes for everything, no for only absorption
all_plot_and_text = 'yes'

# lower limit of absorption width to be flagged 
BALNICITY_INDEX_LIMIT = 500

# limits on velocity     min,   max
VELOCITY_LIMIT = Range(-30000, -60000.)

# range of spectra you are working with from the good_fit.csv file
STARTS_FROM, ENDS_AT = 46,46


trough_detection_percent = 0.9 # what percentage value you want to go below the continuum to detect a trough
limits_percent = 0.9 # What percentage do you want to identifty the limits of the troughs (Vmin/Vmax)


# whether you want to output a csv table of your run
want_csv = 'yes'

###############################################################################################################################
######################################## OUTPUT FILES #########################################################################

# set name of output .txt file with absorption values
ABSORPTION_VALUES = OUT_DIREC + "/" + 'BI' + str(BALNICITY_INDEX_LIMIT) + 'BugVCent.txt'

# set name of output pdf with plots 
ABSORPTION_OUTPUT_PLOT_PDF = PdfPages(OUT_DIREC + '1-108BIVCentBug' + str(BALNICITY_INDEX_LIMIT) + 'VCent.pdf') 

ABSORPTION_TABLE = OUT_DIREC + 'BI2000VCentBug.csv'

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
    field = ['NORM SPECTRA FILE NAME','BI TOTAL','BI INDIVIDUAL','VMINS', 'VMAXS', 'EW INDIVIDUAL', 'DEPTH', 'VCentral']
    append_row_to_csv(ABSORPTION_TABLE, field)
    
    # rounding the numbers of the redshift, calculated snr and setting the norm file name to the current file name from the csv
z = round(redshift_list[STARTS_FROM - 1], 5)
calc_snr = round(calc_snr_list[STARTS_FROM - 1], 5)
norm_spectrum_file_name = norm_spectra_list[STARTS_FROM - 1]

# from the norm spectra name retrieving it's wavelength, normalized flux, and normalized error (in this case from NORM_DRXQ)
print(str(STARTS_FROM), "current spectra file name:", norm_spectrum_file_name)
norm_spectra_data = np.loadtxt(NORM_DIREC + norm_spectrum_file_name)

# setting a variable for each of those values from the spectra
wavelength, normalized_flux, normalized_error = read_spectra(norm_spectra_data)

# smoothing the flux and error based on what the user wants (yes or no)
if want_to_smooth == 'yes':
    normalized_flux = smooth(normalized_flux, boxcar_size)
    normalized_error = smooth(normalized_error, boxcar_size) / math.sqrt(boxcar_size)
    
    
    
#Compatability variabiles
velocity_limits = VELOCITY_LIMIT
EW_percent = 1.0
BI_percent = 0.9
absorption_cutoff = 4
plots = 'yes'




 # variables #########################################################################################################
vmins, vmaxs = [], [] # v = velocity
final_depth_individual, final_depth_all_individual, avg_depth = [], [], []
BI_all, BI_total, BI_ind_sum, BI_individual, BI_all_individual, BI_ind, BI_mid = [], [], [], [], [], [], []
EW_individual, EW_ind, EW_all_individual = [], [], [] #EW = equivalent width
non_trough_count = 9999 # arbitrary large number that we will never reach
delta_v = 0 #change in velocity
sum_of_deltas = 0        
count_v = 0 # variable initialization to get into vmin/vmax loop
V_calc, V_weight = 0,0
v_cent = []
vmaxs_index, vmins_index = [],[]
# transform the wavelength array to velocity (called "beta") based on the CIV doublet: 
beta = wavelength_to_velocity(z, wavelength)

# finding and labeling index of beta that we will be looping through ################################################
                                                    # start,  end
                                                    #   min,  max
              # for reference VELOCITY_LIMIT = Range(-30000, -60000.))
if any(beta):
    try:
        vmaxindex_for_range = np.max(np.where(beta <= velocity_limits.end)) # index value of velocity_limits.end or closest value
    except:
        vmaxindex_for_range = 0  
try:
    vminindex_for_range = np.min(np.where(beta >= velocity_limits.start)) # index value of velocity_limits.start or closest value

except:
    vminindex_for_range = np.where(beta == np.min(beta)) 

velocity_range_index = np.arange(vmaxindex_for_range, vminindex_for_range) # from left to right
velocity_range_index  = np.array(velocity_range_index[::-1])   # from right to left (reversed list)
                                                            # ^^^^^^^^ 0 to -60000
# Warnings + Errors
if trough_detection_percent >= EW_percent:
    raise ValueError('Trough detection percent must be smaller than EW percent')
    
if absorption_cutoff < 4:
    warning_message = 'Low absorption cutoffs will cause exponentially more trough detections, in some case saturating the entire search area'
    warnings.warn(warning_message, UserWarning)
    
# looping through the velocity ranges ##############################################################################
for current_velocity_index in velocity_range_index:
    C = 0 # C will be 0 or 1 and is the C used in the integral for the calculation of BI
    # ([1 - f(v)/0.9] = bracket) > 0 when there is an absorption feature 
    # bracket is the things inside the bracket from the BI integral calculation 
    bracket = (1. - (normalized_flux[current_velocity_index] / trough_detection_percent)) # Determines when we have a trough
    EW_bracket = (1. - (normalized_flux[current_velocity_index] / EW_percent))
    BI_bracket = (1. - (normalized_flux[current_velocity_index]/ BI_percent)) # BI calculation. BI is by default defined as below 0.9
    

    # handle 3-point spike limit ###################################################################################
    if bracket > 0:
        non_trough_count = 0    
    else:
        non_trough_count += 1
        bracket = 0
        
    if EW_bracket < 0: # This prevents negative values being calcuated when bracket > 0 but bracket_two is < 0
        EW_bracket = 0
        
    if BI_bracket < 0 :
        BI_bracket = 0
        
    if((bracket > 0) or (non_trough_count <= 3)):
        delta_v = beta[current_velocity_index] - beta[current_velocity_index - 1]
        sum_of_deltas += delta_v

        EW = EW_bracket * delta_v # EW calculation
        EW = round(EW, 5) 
        EW_ind.append(EW)   
        
        print("beta and 1-normalized flux",beta[current_velocity_index], (1-normalized_flux[current_velocity_index]))
        # A weighted average is the sum(number * weights)/ sum(weights), in this case the number is our velocity
        V_calc_delta = ((1-normalized_flux[current_velocity_index])*beta[current_velocity_index]) 
        if V_calc_delta > 0:
            V_calc_delta = 0
        V_calc += V_calc_delta # Summing (velocity * weights)
        V_weight_delta = (1-normalized_flux[current_velocity_index])
        if V_weight_delta < 0:
            V_weight_delta = 0
        V_weight += V_weight_delta # Summing weights
        print('Vcalc and Vweight', V_calc, V_weight)
        
        
        # BI calculation ###########################################################################################
        if sum_of_deltas >= BALNICITY_INDEX_LIMIT: # passing the BALNICITY_INDEX_LIMIT (in this case 2,000 km/s) threshold
            C = 1  #set to 1 only if square bracket is continuously positive over a velocity interval            
            BI = (BI_bracket * C) * (delta_v) #Calculate BAL for this delta_v. Note BI is DEFINED as below 0.9, do not change this!
            BI_mid.append(round(BI, 5)) #Append to intermediate results
            BI_ind.append(round(BI, 5)) 
            
            # vMIN calculation + plotting ###########################ssdssssss####################################################
            if count_v == 0 and non_trough_count == 0: 
                
                vmins_index = np.min(np.where(beta >= (beta[current_velocity_index] + BALNICITY_INDEX_LIMIT))) #Basically find the next point
                vmins.append(round(beta[vmins_index], 5))
                plot_index = current_velocity_index
                if trough_detection_percent != limits_percent: # We don't need to calculate new vmin/vmax, added BI and EW if trough_detection_percent and limits_percent are the same
                    
                    vmin_searcher = velocity_range_index[velocity_range_index >= vmins_index] #Should be taking all of the points after vmin, this is used later to search for the limits_percent trough
                    vmin_searcher = vmin_searcher[::-1] # Reversing order, again...
                    
                    for current_vmin_index in vmin_searcher: # Search again for limits trough min
                        num_above_percentage = 0
                        
                        BI_bracket_min = (1. - (normalized_flux[current_vmin_index] / BI_percent))
                        if BI_bracket_min < 0:
                            BI_bracket_min = 0
                            
                        EW_bracket_min = (1. - (normalized_flux[current_vmin_index] / EW_percent))
                        if EW_bracket_min < 0: # We might get spikes above limits_percent, but we don't want to add their negative EW
                            EW_bracket_min = 0
                            
                        delta_v = beta[current_vmin_index] - beta[current_vmin_index - 1]
                        
                        EW = EW_bracket_min * delta_v # Calculates the EW gained from searching for the limits_trough minimum
                        EW = round(EW,5)
                        EW_ind.append(EW)
                        
                        BI_min = (BI_bracket_min* delta_v) # Calculates the BI gained from searching for the limits_trough minimum
                        BI_mid.append(round(BI_min,5))
                        BI_ind.append(round(BI_min, 5))
                        
                        for i in range(1, absorption_cutoff + 1): # Checks the next points, if they are all above limits_percent then we have our new min
                            if normalized_flux[current_vmin_index + i] >= limits_percent:
                                num_above_percentage = num_above_percentage + 1
                                
                            else:
                                break
                        if num_above_percentage >= absorption_cutoff: # If => our cutoff or we reach the end of the range, we save the vmin and then break out
                            vmins[-1] = (round(beta[current_vmin_index], 5))
                            vmins_index = current_vmin_index
                            plot_index = beta[current_vmin_index] - BALNICITY_INDEX_LIMIT
                            plot_index = np.min(np.where(beta >= plot_index)) # This index is used for calculating the wavelist below
                            
                            break
                    
                
                if plots == 'yes':
                   
                    # plotting notable vertical line of v min occurance in absorption found
                    vmin_line(beta, vmins_index)

                    # plotting notable vertical line of v min occurance of where CIV would be *if* the EHVO 
                    # absorption found was due to SiIV
                    wavelist = vmin_plot_IF(beta, wavelength, plot_index, BALNICITY_INDEX_LIMIT)
                    carbon_iv = wavelist[0] # the vmin value of where CIV would be *if* the EHVO absorption found was due to SiIV
                    carbon_ii = wavelist[1] # the vmin value of where CII should be *if* the EHVO absorption found was due to SiIV
                    oxygen_i = wavelist[2] # the vmin value of where OI should be *if* the EHVO absorption found was due to SiIV               
                else: 
                    pass
                
                count_v = 1 # Variable used to prevent calculating multiple vmins per trough
            
            num_above_percentage = 0
            for i in range(1, absorption_cutoff + 1): # Finding max for trough_detection_percent, better than old method
                if normalized_flux[current_velocity_index - i] >= trough_detection_percent:
                    num_above_percentage = num_above_percentage + 1
                else:
                    break

            # vMAX calculation + plotting #############################################################################
            if (((bracket > 0 and num_above_percentage >= absorption_cutoff and count_v == 1)) or 
                (current_velocity_index == vmaxindex_for_range)): 

                vmaxs_index = current_velocity_index # Should be faster than before
                vmaxs.append(round(beta[current_velocity_index], 5))
                v_cent.append(V_calc/V_weight) # Dividing the sums
                avg_depth.append(V_weight/(vmins_index-vmaxs_index)) # V_weight is the sum of all the depths in the trough, the avg_depth can be easily calculated by dividing by the number of pixels
                
                print('HERE')
                if trough_detection_percent != limits_percent: # Same for Vmax
                    vmax_searcher = velocity_range_index[velocity_range_index <= vmaxs_index] # list to search for new limits_percent max
                    
                    for current_vmax_index in vmax_searcher:
                        num_above_percentage = 0
                        
                        BI_bracket_max = (1. - (normalized_flux[current_vmax_index] / BI_percent))
                        if BI_bracket_max < 0:
                            BI_bracket_max = 0
                            
                        EW_bracket_max = (1. - (normalized_flux[current_vmax_index] / EW_percent))
                        if EW_bracket_max < 0:
                            EW_bracket_max = 0
                        
                        delta_v = beta[current_vmax_index] - beta[current_vmax_index - 1]
                        
                        EW = EW_bracket_max * delta_v #Calculates the EW gained from searching for the limits_percent maximum
                        EW = round(EW,5)
                        EW_ind.append(EW)
                        
                        BI_max = (BI_bracket_max * delta_v) # Calculates the BI gained from searching for the limits_percent maximum
                        BI_mid.append(round(BI_max,5))
                        BI_ind.append(round(BI_max, 5))

                        for i in range(1, absorption_cutoff + 1):
                            if normalized_flux[current_vmax_index - i] >= limits_percent:
                                num_above_percentage = num_above_percentage + 1
                            else:
                                break
                        if num_above_percentage >= absorption_cutoff or current_vmax_index == vmaxindex_for_range: # If => our cutoff or we reach the end of the range, we save the vmax and then break out
                            vmaxs[-1] = (round(beta[current_vmax_index], 5))
                            vmaxs_index = current_vmax_index
                            break

                if plots == 'yes':
                    # plotting the red bar between vmins and vmaxs index
                    span_vmin_vmax(beta, vmins_index, vmaxs_index)
                    # plotting different colored bars based on where CIV, CII, and OI would be *if* the 
                    # EHVO absorption found was due to SiIV
                    vmax_plot_span_IF(beta, wavelength, vmaxs_index, carbon_iv, carbon_ii, oxygen_i)
                    plot_min = 1- beta[vmins_index]/-70000
                    plot_max = 1 -beta[vmaxs_index]/-70000
                    plt.axhline(1.37, plot_min, plot_max, color = 'k')

                else: 
                    pass

                BI_ind_sum = round(sum(BI_ind), 2)
                BI_individual.append(BI_ind_sum) # this array contains one single BI value of each absorption feature in a single spectrum
                BI_ind = []
                
                EW_ind_sum = round(sum(EW_ind), 2)
                EW_individual.append(EW_ind_sum)
                EW_ind = []
                
                # depth calculation ##################################################################################
                final_depth = round(((1. - np.min(normalized_flux[vmaxs_index:vmins_index]))),2) 
               
                final_depth_individual.append(final_depth)
                final_depth = []
                
                V_calc = 0
                V_weight = 0
                sum_of_deltas = 0 # Added over old code. Prevents bug with vmin calculations/double counting
                count_v = 0 
    #if the bracket value is not more than zero (so if we don't have absorption feature)
    else: 
        sum_of_deltas = 0 # to reset counting the width of the absorption feature if it is not wider than the BI_index_limit
        count_v = 0 # reset count in case there is another absorption feature that is wider than 2,000km/s
        EW_ind = []
        V_calc = 0
        V_weight = 0
    
    if current_velocity_index == vmaxindex_for_range:
        BI_total = round(sum(BI_mid), 2)         
        BI_all.append(BI_total)    
        BI_all_individual.append(BI_individual)
        EW_all_individual.append(EW_individual)

final_depth_all_individual.append(final_depth_individual)   
BI_all= np.array(BI_all)
vmins = np.array(vmins)
vmaxs = np.array(vmaxs)
BI_all= np.array(BI_all)

vmins = np.array(vmins)
vmaxs = np.array(vmaxs)
v_cent = np.array(v_cent)