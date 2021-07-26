#Add keyword argument detailing the specific range over which BI is calculated. Outputs will need to be total BI_EHVO, BI per abs feature,
#vmax, and vmin per absorption feature, and depth of absorption feature
 
  # Initialize all variables for each spectrum (again, clean so it is not so many lines)
from absorption import BALNICITY_INDEX_LIMIT, BI_all, VELOCITY_LIMIT
import os
import sys
import numpy as np 
from matplotlib import pyplot as plt
from scipy import signal
from numpy.lib.function_base import append
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
from utility_functions import print_to_file, clear_file, read_list_spectra, read_spectra, wavelength_to_velocity
from data_types import Range, RangesData, FigureData, FigureDataOriginal, FlaggedSNRData, DataNormalized 
from abs_plot import draw_abs_figure 


vmins, vmaxs=[]
BI_total, EW_individual=[] #EW is equivalent width = (vmax-vmin)
BI_individual=[] #BALnicity index will be single-valued per spectrum in terms of km/s

# verner table data
wavelength_CIV_emit1=  1550.7700
wavelength_CIV_emit2 = 1548.1950
avr_CIV_doublet = 1549.0524 #weighted average
avr_SiIV_doublet = 1396.747 # weighted average; individuals: 1402.770, 1393.755
CII_emitted = 1335.313 # (weighted average); individuals:
OI_emitted = 1303.4951 # weighted average; individuals pag 20 in Verner Table
avr_NV_doublet = 1240.15 # weighted average; individuals: 1242.80, 1238.82
avr_OVI_doublet = 1033.8160 # weighted average; individuals: 1037.6167, 1031. 9261


######################################### VARIABLES ###########################################################################

brac_all, delta_v_all = [], []
absspeccount = 0
count = 0
BI = 0
vmins, vmaxs, vmins_all, vmaxs_all = [], [], [], [] # v = velocity
final_depth_individual, final_depth_all_individual = [], []
BI_total, BI_ind_sum, BI_individual, BI_all_individual, BI_mid, BI_ind, BI_all = [], [], [], [], [], [], []
EW_individual, EW_all_individual, vlast, EW_ind = [], [], [], [] #EW = equivalent width

###############################################################################################################################
######################################### MAIN CODE ###########################################################################

    ################################# INITIALIZING  VARIABLES IN LOOP ##############################################################
    index_depth_final, flux_depth, final_depth_individual = [], [], []
    non_trough_count = 999 # arbitrary large number that we will never reach

    delta_v = 0 #change in velocity
    sum_of_deltas = 0
    bb = -1 #???
    BALNICITY_INDEX_LIMIT =2000
    
    count_v = 0 # variable initialization to get into vmin/vmax loop
    ###############################################################################################################################

    # Calculate BI, v_min and v_max by looping through the beta array in the velocity limits
    # Calculate depth of each individual absorption trough
    # VVVVVVVV add these things later VVVVVVVV once the module is made
    # BI_ehvo, BI_abs, v_min, v_max, EW, depth = basic_absorption_parameters(wavelength, normalized_flux, z, VELOCITY_LIMIT.end, VELOCITY_LIMIT.start)

                                                        #   min,  max
    if beta.any(): # for reference VELOCITY_LIMIT = Range(-30000, -60000.))
        try:
            vmaxindex_for_range = np.max(np.where(beta <= VELOCITY_LIMIT.end)) #index value of the starting point (on the very left) -- index value of VELOCITY_LIMIT.end
        except:
            vmaxindex_for_range = 0
    if beta.any(): # delete if statement ?? OLD_absorption.py does not have it
        try:
            vminindex_for_range = np.min(np.where(beta >= VELOCITY_LIMIT.start)) #index value of the ending point (on the very right) -- index value of VELOCITY_LIMIT.start
        except:
            vminindex_for_range = np.where(beta == np.max(beta)) 

    velocity_range_index = np.arange(vmaxindex_for_range, vminindex_for_range)
    velocity_range_index  = np.array(velocity_range_index[::-1])   # From right to left (reversed list)

    # looping through the velocity ranges
    for current_velocity_index in velocity_range_index:
        C = 0 # C will be 0 or 1 and is the C used in the integral for the calculation of BI

        # ([1 - f(v)/0.9] = bracket) > 0 when there is an absorption feature 
        # bracket is the things inside the bracket from the BI integral calculation 
        bracket = (1. - (normalized_flux[current_velocity_index] / 0.9))
        
        # Handle 3-point spike limit #####################################################
        if bracket > 0:
            non_trough_count = 0
        else:
            non_trough_count += 1
            bracket = 0

        if((bracket > 0) and (non_trough_count <= 3)):
            delta_v = beta[current_velocity_index] - beta[current_velocity_index - 1]
            sum_of_deltas += delta_v
            brac_all.append(bracket)
            delta_v_all.append(delta_v)

            EW = bracket * delta_v
            EW = np.round(EW, 5)
            EW_ind.append(EW)      
        ################################################################################# 
            #----------------------------------------------------------------------------------------------
            ################################# BI calculation ################################
            if sum_of_deltas >= BALNICITY_INDEX_LIMIT: # passing the BALNICITY_INDEX_LIMIT (in this case 2,000 km/s) threshold
                C = 1  #set to 1 only if square bracket is continuously positive over a velocity interval            
                BI = (bracket * C) * (delta_v) #Calculate BAL for this delta_v
                BI_mid.append(np.round(BI, 5)) #Append to intermediate results
                BI_ind.append(np.round(BI, 5)) 

                # plotting the black line
                if non_trough_count == 0: 
                   #plt.plot((beta[current_velocity_index + 1], beta[current_velocity_index]), (1.5,1.5),'k-')
                #----------------------------------------------------------------------------------------------
                ############################# V MIN CALCULATION ######################################################
                if count_v == 0 and non_trough_count == 0:  
                    vmins_index = np.min(np.where(beta >= (beta[current_velocity_index] + BALNICITY_INDEX_LIMIT))) # vmins occurs current beta plus BALNICITY_INDEX_LIMIT
                    vmins.append(np.round(beta[vmins_index], 5))                    
                    #######################################################

                    # Calculate where CIV, CII and OI would be for each pair of VMIN *if* the EHVO absorption found were 
                    # instead not EHVO and due to SiIV: 
                    z_absSiIV = (wavelength[current_velocity_index] / avr_SiIV_doublet) - 1
                     
                    obs_wavelength_C = (z_absSiIV + 1) * (avr_CIV_doublet)
                    obs_wavelength_C_index = np.min(np.where(wavelength > obs_wavelength_C))
                    obs_wavelength_C_vel = beta[obs_wavelength_C_index] + BALNICITY_INDEX_LIMIT
                    #plt.plot((obs_wavelength_C_vel, obs_wavelength_C_vel),(-1,10),'k-')

                    obs_wavelength_CII = (z_absSiIV + 1) * (CII_emitted)
                    obs_wavelength_CII_index = np.min(np.where(wavelength > obs_wavelength_CII))                  
                    obs_wavelength_CII_vel = beta[obs_wavelength_CII_index] + BALNICITY_INDEX_LIMIT
                    #plt.plot((obs_wavelength_CII_vel, obs_wavelength_CII_vel),(-1,10),'b-')

                    obs_wavelength_OI = (z_absSiIV + 1) * (OI_emitted)
                    obs_wavelength_OI_index = np.min(np.where(wavelength > obs_wavelength_OI))                  
                    obs_wavelength_OI_vel = beta[obs_wavelength_OI_index] + BALNICITY_INDEX_LIMIT
                    #plt.plot((obs_wavelength_OI_vel, obs_wavelength_OI_vel),(-1,10),'y-')
                    ############################################################################################

                    count_v = 1
                
                bracket_2 = (1. - (normalized_flux[current_velocity_index - 1] / 0.9))
                bracket_3 = (1. - (normalized_flux[current_velocity_index - 2] / 0.9))
                bracket_4 = (1. - (normalized_flux[current_velocity_index - 3] / 0.9))
                bracket_5 = (1. - (normalized_flux[current_velocity_index - 4] / 0.9))
                
                #----------------------------------------------------------------------------------------------
                ############################# V MAX CALCULATION ######################################################
                if (((bracket > 0 and bracket_2 < 0 and bracket_3 < 0 and bracket_4 < 0 and bracket_5 < 0 and count_v == 1)) or (current_velocity_index == vminindex_for_range)):  

                    vmaxs_index = np.min(np.where (beta >= beta[current_velocity_index]))
                    vmaxs.append(np.round(beta[current_velocity_index], 4))
                                            ###########################################################                   
                    #plt.axvspan(beta[vmins_index], beta[vmaxs_index], alpha = 0.2, color = 'red')

                    z_absSiIV_final = (wavelength[vmaxs_index] / avr_SiIV_doublet) - 1.
                    # Calculate where CIV, CII and OI would be for each pair of VMAX *if* the EHVO absorption found were 
                    # instead not EHVO and due to SiIV: 
                    # if the absorption is SiIV, this finds and plots where CIV, CII and OI would be ###########
                    obs_wavelength_Cfinal = (z_absSiIV_final + 1.) * (avr_CIV_doublet)
                    obs_wavelength_Cfinal_index = np.min(np.where(wavelength > obs_wavelength_Cfinal))
                    obs_wavelength_C_final_vel = beta[obs_wavelength_Cfinal_index]
                    #plt.axvspan(obs_wavelength_C_vel, obs_wavelength_C_final_vel, alpha = 0.2, color = 'grey')

                    obs_wavelength_CIIfinal = (z_absSiIV_final + 1.) * (CII_emitted)
                    obs_wavelength_CIIfinal_index = np.min (np.where (wavelength > obs_wavelength_CIIfinal))
                    obs_wavelength_CII_final_vel = beta[obs_wavelength_CIIfinal_index]
                    #plt.axvspan(obs_wavelength_CII_vel,obs_wavelength_CII_final_vel, alpha = 0.2, color = 'blue')

                    obs_wavelength_OIfinal = (z_absSiIV_final + 1.) * (OI_emitted)
                    obs_wavelength_OIfinal_index = np.min(np.where (wavelength > obs_wavelength_OIfinal))
                    obs_wavelength_OI_final_vel = beta[obs_wavelength_OIfinal_index]
                    #plt.axvspan(obs_wavelength_OI_vel,obs_wavelength_OI_final_vel, alpha = 0.2, color = 'yellow')
                    ############################################################################################

                    ###################### BALnicity Index ###########################
                    BI_ind_sum = round(sum(BI_ind), 2)
                    BI_individual.append(BI_ind_sum) # this array contains one single BI value of each absorption feature in a single spectrum
                    BI_ind = []
                    ######################### Equivalent Width #######################
                    EW_ind_sum = round(sum(EW_ind), 2)
                    EW_individual.append(EW_ind_sum)
                    EW_ind = []
                    ############################ Depth Function #######################  
                    final_depth = round((1. - np.min(normalized_flux[vmaxs_index:vmins_index])), 2)
                    final_depth_individual.append(final_depth)
                    print('depth', final_depth_individual)
                    
"""Per spectrum, we can just take the individual BI for every trough and sum them for total BI;
I'm not sure what differentiates the BI_total and BI_all so I made BI_sum """
#This might be replacing the above work^
BI_sum = np.sum(BI_individual)



######################### Newer Depth function without arrary slicing ###################
"""
def depth(vmins_index,vmaxs_index): #this is most likely slower than the originally given depth function
    #so we can stick with that
"""
"""Returns the depth of each BAL from the given spectra. Must be evaluated on each spectrum.

Parameters
----------
vmins_index : float
    The starting point for the BAL trough
vmaxs_index : float
    The ending border of the BAL trough

Returns
-------
float
    "depth" will return the lowest flux for the region bordered by vmin and vmax. Be aware that broken pixels in the region may
    lead to erroneous depths.

"""
"""
    #so the indices across each BAL need to be looped through:
    depths=[]
    for i in range(len(vmins_index)):
        BALregion=np.linspace(vmins_index[i],vmaxs_index[i],vmaxs_index[i]-vmins_index[i]+1)
        #print(BALregion)
        minimum = np.min(norm_flux_used(BALregion))
        #print(minimum)
        depths.append(minimum)
    return depths
"""
#################################################### output needs to be a textfile ##################
if (len(vmaxs) != 0):
    text = [f"{current_spectrum_file_name}",
        f"BI({VELOCITY_LIMIT.end} > v > {VELOCITY_LIMIT.start}): {BI_total}",
        f"vmins: {vmins}"
        f"vmaxs: {vmaxs}"
        f"BI_individual: {BI_individual}"
        f"EW_individual: {EW_individual}"
        f"Depth: {final_depth_individual}"]
    vlast.extend(['\n'.join(text), '\n'])
else:
    print('there was an error somewhere')
######################################################