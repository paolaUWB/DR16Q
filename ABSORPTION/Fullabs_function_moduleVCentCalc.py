"""
abs_function_module
===================

@author Wendy Garcia Naranjo, Mikel Charles, Nathnael Kahassai, Michael Parker 
based on code prepared by Abdul Khatri and Paola Rodriguez Hidalgo

Notes
-----
    Contains various functions needed for ``absorption.py``. One to convert wavelength to velocity (wavelength_to_velocity), another 
    with the smooth function (``smooth``), and one that calulates the BI, vmin and vmax, equivalent width, and depth for a list of spectra 
    AND plots where CIV, CII, and OI would be *if* the EHVO absorption found was due to SiIV (``abs_parameters_plot_optional``). 
"""

###############################################################################################################################
########################################## IMPORTS ############################################################################
import numpy as np 
import scipy.constants as sc
from scipy import signal
from matplotlib import pyplot as plt
from numpy.lib.function_base import append
from data_types import Range
from abs_plot_module import vmin_plot_IF, vmax_plot_span_IF, vmin_line, span_vmin_vmax
import warnings
###############################################################################################################################
######################################### Functions ###########################################################################

def wavelength_to_velocity(redshift, wavelength):
    """Reads in a list of wavelength values to be converted to velocity.

    Parameters
    ----------
    redshift: list
        The list of redshift values needed for the equation.
    wavelength: list
        The list of wavelength values that will be converted to velocity.
        
    Returns
    -------
    beta: array
        The values of velocity that were converted from the wavelength provided.
    """
    # CIV doublet data from verner table
    avr_CIV_doublet = 1549.0524

    # Transform the wavelength array to velocity (called "beta") based on the CIV doublet: 
    c_in_km = sc.speed_of_light * (10**-3) # speed_of_light is in m/s
    z_absC = (wavelength / avr_CIV_doublet) - 1.
    RC = (1. + redshift) / (1. + z_absC)
    betaC = ((RC**2.) - 1.) / ((RC**2.) + 1.) # betaC is in units of c (speed of light)
    betakm = -betaC * c_in_km # betakm is in km/s
    beta = []

    for velocity in betakm:
        betas = round(velocity, 5)
        beta.append(betas)
    beta = np.array(beta)

    return beta

#############################################################################################################################################
#############################################################################################################################################

def smooth(smooth_this, box_size):   
    """Smooths out values entered.

    Parameters
    ----------
    smooth_this: list
        The value(s) of whatever is going to get smoothed.
    box_size: int
        This is the number of points that are smoothed into one. Always be sure to use an odd 
        number, because we need the same amount of points on each side of the data point to be
        smoothed.

    Returns
    -------
    now_smooth: list
        The values that are now smoothed.
    """   

    now_smooth = signal.savgol_filter(smooth_this, box_size, 2)
    return now_smooth

#############################################################################################################################################
#############################################################################################################################################
def wavelist_calc(beta, wavelength, plot_index, BALNICITY_INDEX_LIMIT):
    """Finds the minimum velocity index and value of where CIV, CII and OI would be *if* the EHVO absorption 
    found were instead not EHVO and due to SiIV, and plots a vertical line for that value.
    
    Parameters
    ----------
    beta: array or list
        The velocity value(s) to be used.
    wavelength: array or list or int
        The wavelength values to be used.
    current_velocity_index: list or int
        When looping through a velocity, the specific velocity index that you are on.
    BALNICITY_INDEX_LIMIT: int
        The minimum value of broad absorption width that we are looking for.

    Returns
    -------
    obs_wavelength_C_vel: array or list
        The minimum value of the velocity for where CIV would be for each pair of VMIN 
        *if* the EHVO absorption found were instead not EHVO and due to SiIV.
    obs_wavelength_CII_vel: array or list
        The minimum value of the velocity for where CII would be for each pair of VMIN 
        *if* the EHVO absorption found were instead not EHVO and due to SiIV.
    obs_wavelength_OI_vel: array or list
        The minimum value of the velocity for where OI would be for each pair of VMIN 
        *if* the EHVO absorption found were instead not EHVO and due to SiIV.
    """
    WAVELENGTH_CIV_EMIT_LIMIT = Range(1548.1950, 1550.7700)                                    #never used?
    AVERAGE_CIV_DOUBLET = 1549.0524 #weighted average
    AVERAGE_SiIV_DOUBLET = 1396.747 # weighted average; individuals: 1402.770, 1393.755
    AVERAGE_NV_DOUBLET = 1240.15 # weighted average; individuals: 1242.80, 1238.82             # not used
    AVERAGE_OVI_DOUBLET=1033.8160 # weighted average; individuals: 1037.6167, 1031. 9261        # also not used
    CII_EMITTED = 1335.313 # (weighted average); individuals:
    OI_EMITTED = 1303.4951 # weighted average; individuals pag 20 in Verner Table
    
    z_absSiIV = (wavelength[plot_index] / AVERAGE_SiIV_DOUBLET) - 1
        
    obs_wavelength_C = (z_absSiIV + 1) * (AVERAGE_CIV_DOUBLET)
    obs_wavelength_C_index = np.min(np.where(wavelength > obs_wavelength_C))
    obs_wavelength_C_vel = beta[obs_wavelength_C_index] + BALNICITY_INDEX_LIMIT
    plt.plot((obs_wavelength_C_vel, obs_wavelength_C_vel), (-1,10), 'k-')

    obs_wavelength_CII = (z_absSiIV + 1) * (CII_EMITTED)
    obs_wavelength_CII_index = np.min(np.where(wavelength > obs_wavelength_CII))                  
    obs_wavelength_CII_vel = beta[obs_wavelength_CII_index] + BALNICITY_INDEX_LIMIT
    plt.plot((obs_wavelength_CII_vel, obs_wavelength_CII_vel), (-1,10), 'b-')

    obs_wavelength_OI = (z_absSiIV + 1) * (OI_EMITTED)
    obs_wavelength_OI_index = np.min(np.where(wavelength > obs_wavelength_OI))                  
    obs_wavelength_OI_vel = beta[obs_wavelength_OI_index] + BALNICITY_INDEX_LIMIT
    plt.plot((obs_wavelength_OI_vel, obs_wavelength_OI_vel), (-1,10), 'y-')

    return [obs_wavelength_C_vel, obs_wavelength_CII_vel, obs_wavelength_OI_vel]
    

#############################################################################################################################################
#############################################################################################################################################


def abs_parameters_plot_optional(z, wavelength, normalized_flux, BALNICITY_INDEX_LIMIT, velocity_limits, trough_detection_percent, limits_percent, plots = 'yes', EW_percent = 1.0, BI_percent = 0.9, absorption_cutoff = 4):
    """Based off and does what find_absorption_parameters does, but also includes plotting.

    Reads in a list of redshift, wavelength, velocity limit (your integral bounds), broad absorption width, and percentage value 
    you want to go below the continuum to calculate BI. Returns BI for each indivdual trough, the total BI from all troughs, 
    vmin/vmaxs of the trough found if there are more than one, the equivalent width, the depth of the trough and the beta values
    that were converted from wavelength. Plots where CIV, CII, and OI would be *if* the EHVO absorption found was due to SiIV and 
    plots the spectra as normalzied flux vs velocity, and error vs velocity.
    
    Parameters
    ----------
    z: list
        The redshift values.
    wavelength: list
        The wavelength values.
    normalized_flux: list
        The normalized flux values.
    BALNICITY_INDEX_LIMIT: int
        The minimum value of broad absorption width that we are looking for. 
    velocity_limits: namedtuple
        The velocity limits that we will be searching for absorption, aka the integral limits of BI calculation.
    trough_detection_percent: float
        The normalized flux threshold required for detection of a trough. Standard is 0.9, aka the trough must be below 0.9 for the program to
        recognize it as a trough (as well as meet the requirements for the BALNICITY_INDEX_LIMIT).
    limits_percent: float
        The normalized flux value that the program identifies the trough limits at (Vmin and Vmax). Standard is 0.9, thus the first point where the trough
        goes below 0.9 is the Vmin and the final point is Vmax. 
    plots: string, default = 'yes'
        Whether you want to plot the values or just want the values, the default is to plot. 
    EW_percent: float, default = '1.0'
        Determines where the Equivalent Width calculations starts, 1.0 is at the continuum line, 0.9 would be 0.1 below the continuum, etc.
    BI_percent: float, default = '0.9'
        The normalized flux value where the BI calculation starts. Balnicity Index is traditionally defined as below 0.9, but you can change this for
        debugging purposes.
    absorption_cutoff: int, default = '4'
        Determines how many of the next points need to be above our limits_percent in order for the trough Vmin and Vmax to be calculated and the trough to end.
        It is ill advised to change from 4, too many can mean troughs never end, too few means an exponential number of trough detections.
        

    Returns
    -------
    BI_total: list
        The total BI for the entire figure added up together. 
    BI_individual: list
        The individual BI value for 1 trough, if there is only 1 trough then ``BI_individual = BI_total``.
    BI_all: list
        Contains ALL of the total BI_total values for all the spectra ran.
    vmins: array
        The value of where the minimum velocity occurs for a trough.
    vmaxs: array
        The value of where the maximum velocity occurs for a trough.
    EW_individual: list
        The value of the equivalent width for each individual trough.
    final_depth_individual: list
        The value of the depth of each individual trough.
    final_depth_all_individual: list
        ALL of the final_depth_individual values for all of the spectra that were ran.
    beta: array
        The velocity values that were converted from wavelength.
    vminindex_for_range: int
        The index of where the minimum velocity is located. This is created with the intent to help scale the y-axis of graphing.
    vmaxindex_for_range: int
        The index of where the maximum velocity is located. This is created with the intent to help scale the y-axis of graphing.
    v_cent: array
        The weighted central velocity for each trough
    avg_depth: array
        The weighted average depth for each trough
    Note
    -----
    ``velocity_limits`` is a namedtuple, in the main code when you call the function make sure you either create your own or manually
    change the velocity limits into ``int``s. An example for how to create a namedtuple (for this specific code) is in ``data_types.py``.
    """
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
            
            # A weighted average is the sum(number * weights)/ sum(weights), in this case the number is our velocity
            V_calc_delta = ((1-normalized_flux[current_velocity_index])*beta[current_velocity_index])
            if V_calc_delta > 0: # Prevents spikes above 1.0 from effecting the location of Vcent
                V_calc_delta = 0
            V_calc += V_calc_delta # Summing (velocity * weights)
            V_weight_delta = (1-normalized_flux[current_velocity_index])
            if V_weight_delta < 0: # Prevents spikes above 1.0 from effecting the location of Vcent
                V_weight_delta = 0
            V_weight += V_weight_delta # Summing weights
            
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

    return BI_total, BI_individual, BI_all, vmins, vmaxs, EW_individual, final_depth_individual, final_depth_all_individual, beta, vminindex_for_range, vmaxindex_for_range, v_cent, avg_depth