"""
==============================
abs_function_module.py
==============================

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
from abs_plot_module import vmin_plot_IF, vmax_plot_span_IF, vmin_line, span_vmin_vmax, black_line
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

    Parameters:
    -----------
    smooth_this: list
        The value(s) of whatever is going to get smoothed.

    box_size: int
        This is the number of points that are smoothed into one. Always be sure to use an odd 
        number, because we need the same amount of points on each side of the data point to be
        smoothed.

    Returns:
    --------
    now_smooth: list
        The values that are now smoothed.
    """   

    now_smooth = signal.savgol_filter(smooth_this, box_size, 2)
    return now_smooth

#############################################################################################################################################
#############################################################################################################################################

def abs_parameters_plot_optional(z, wavelength, normalized_flux, BALNICITY_INDEX_LIMIT, velocity_limits, percent, plots = 'yes'):
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

    percent: float
        The percentage value you want to go below the continuum.

    plots: string, default = 'yes'
        Whether you want to plot the values or just want the values, the default is to plot. 

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
        The index of where the minimum velocity is located. This is created with the intent to help scale the y-axis of grahing.

    vmaxindex_for_range: int
        The index of where the maximum velocity is located. This is created with the intent to help scale the y-axis of grahing.

    Note
    -----
    ``velocity_limits`` is a namedtuple, in the main code when you call the function make sure you either create your own or manually
    change the velocity limits into ``int``s. An example for how to create a namedtuple (for this specific code) is in ``data_types.py``.
    """
    # variables #########################################################################################################
    brac_all = []
    vmins, vmaxs, vmins_all, vmaxs_all, delta_v_all, vmins_all_index, vmaxs_all_index = [], [], [], [], [], [], [] # v = velocity
    final_depth_individual, final_depth_all_individual = [], []
    BI_all, BI_total, BI_ind_sum, BI_individual, BI_all_individual, BI_ind, BI_mid = [], [], [], [], [], [], []
    EW_individual, EW_ind, EW_all_individual = [], [], [] #EW = equivalent width
    non_trough_count = 999 # arbitrary large number that we will never reach
    delta_v = 0 #change in velocity
    sum_of_deltas = 0        
    count_v = 0 # variable initialization to get into vmin/vmax loop

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

    # looping through the velocity ranges ##############################################################################
    for current_velocity_index in velocity_range_index:
        C = 0 # C will be 0 or 1 and is the C used in the integral for the calculation of BI
        # ([1 - f(v)/0.9] = bracket) > 0 when there is an absorption feature 
        # bracket is the things inside the bracket from the BI integral calculation 
        bracket = (1. - (normalized_flux[current_velocity_index] / percent))

        # handle 3-point spike limit ###################################################################################
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
            EW = round(EW, 5)
            EW_ind.append(EW)   
         
            # BI calculation ###########################################################################################
            if sum_of_deltas >= BALNICITY_INDEX_LIMIT: # passing the BALNICITY_INDEX_LIMIT (in this case 2,000 km/s) threshold
                C = 1  #set to 1 only if square bracket is continuously positive over a velocity interval            
                BI = (bracket * C) * (delta_v) #Calculate BAL for this delta_v
                BI_mid.append(round(BI, 5)) #Append to intermediate results
                BI_ind.append(round(BI, 5)) 
                
                if plots == 'yes':
                    # plotting the black line inside the absorption found
                    if non_trough_count == 0: 
                        black_line(beta, current_velocity_index)
                else: 
                    pass
                
                # vMIN calculation + plotting ###############################################################################
                if count_v == 0 and non_trough_count == 0:  
                    vmins_index = np.min(np.where(beta >= (beta[current_velocity_index] + BALNICITY_INDEX_LIMIT))) 
                    vmins.append(round(beta[vmins_index], 5))     
                    
                    if plots == 'yes':
                        # plotting notable vertical line of v min occurance in absorption found
                        vmin_line(beta, vmins_index)

                        # plotting notable vertical line of v min occurance of where CIV would be *if* the EHVO 
                        # absorption found was due to SiIV
                        wavelist = vmin_plot_IF(beta, wavelength, current_velocity_index, BALNICITY_INDEX_LIMIT)
                        carbon_iv = wavelist[0] # the vmin value of where CIV would be *if* the EHVO absorption found was due to SiIV
                        carbon_ii = wavelist[1] # the vmin value of where CII should be *if* the EHVO absorption found was due to SiIV
                        oxygen_i = wavelist[2] # the vmin value of where OI should be *if* the EHVO absorption found was due to SiIV               
                    else: 
                        pass
                    
                    count_v = 1
                
                bracket_1 = (1. - (normalized_flux[current_velocity_index - 1] / 0.9))
                bracket_2 = (1. - (normalized_flux[current_velocity_index - 2] / 0.9))
                bracket_3 = (1. - (normalized_flux[current_velocity_index - 3] / 0.9))
                bracket_4 = (1. - (normalized_flux[current_velocity_index - 4] / 0.9))

                # vMAX calculation + plotting #############################################################################
                if (((bracket > 0 and bracket_1 < 0 and bracket_2 < 0 and bracket_3 < 0 and bracket_4 < 0 and count_v == 1)) or 
                    (current_velocity_index == vmaxindex_for_range)):  

                    vmaxs_index = np.min(np.where(beta >= beta[current_velocity_index]))
                    vmaxs.append(round(beta[current_velocity_index], 5))

                    if plots == 'yes':
                        # plotting the red bar between vmins and vmaxs index
                        span_vmin_vmax(beta, vmins_index, vmaxs_index)
                        # plotting different colored bars based on where CIV, CII, and OI would be *if* the 
                        # EHVO absorption found was due to SiIV
                        vmax_plot_span_IF(beta, wavelength, vmaxs_index, carbon_iv, carbon_ii, oxygen_i)

                        #presentation(beta, normalized_flux, vmins_index, vmaxs_index) made just to create figure for slides
                    else: 
                        pass

                    BI_ind_sum = round(sum(BI_ind), 2)
                    BI_individual.append(BI_ind_sum) # this array contains one single BI value of each absorption feature in a single spectrum
                    BI_ind = []
                    
                    EW_ind_sum = round(sum(EW_ind), 2)
                    EW_individual.append(EW_ind_sum)
                    EW_ind = []
                    
                    # depth calculation ##################################################################################
                    final_depth = round((1. - np.min(normalized_flux[vmaxs_index:vmins_index])), 2)
                    final_depth_individual.append(final_depth)
                    
                    count_v = 0 
        #if the bracket value is not more than zero (so if we don't have absorption feature)
        else: 
            sum_of_deltas = 0 # to reset counting the width of the absorption feature if it is not wider than the BI_index_limit
            count_v = 0 # reset counte in case there is another absorption feature that is wider than 2,000km/s
            EW_ind = []
        
        if current_velocity_index == vmaxindex_for_range:
            BI_total = round(sum(BI_mid), 2)         
            BI_all.append(BI_total)    
            BI_all_individual.append(BI_individual)
            EW_all_individual.append(EW_individual)

    final_depth_all_individual.append(final_depth_individual)   
    BI_all= np.array(BI_all)
    vmins = np.array(vmins)
    vmaxs = np.array(vmaxs)
    '''
    vmins_index = np.array(vmins_index)
    vmaxs_index = np.array(vmaxs_index)

    vmins_all_index.append(vmins_index)
    vmaxs_all_index.append(vmaxs_index)
    '''
    return BI_total, BI_individual, BI_all, vmins, vmaxs, EW_individual, final_depth_individual, final_depth_all_individual, beta, vminindex_for_range, vmaxindex_for_range