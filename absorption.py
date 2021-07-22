"""
=============
absorption.py
=============

@author Wendy Garcia Naranjo, Mikel Charles, Nathnael Kahassai, Michael Parker
based on code prepared by Abdul Khatri and Paola Rodriguez Hidalgo

Creates figure to visually inspect the absorption. Calculates absorption parameters 
(BALnicity Index BI, vmin and vmax) for a list of spectra.

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
from matplotlib import pyplot as plt
from scipy import signal
from numpy.lib.function_base import append
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
from utility_functions import print_to_file, clear_file, read_list_spectra, read_spectra, wavelength_to_velocity
from data_types import Range, RangesData, FigureData, FigureDataOriginal, FlaggedSNRData, DataNormalized 
from abs_plot import draw_abs_figure 
#import basic_absorption_parameters

###############################################################################################################################
############################## CHANGEABLE VARIABLES ###########################################################################

# input which data release you are working with [input the number as a string i.e. '9']
DR = '16'

# defining the config file
CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/OUTPUT_FILES/NORMALIZATION/good_normalization.csv" 

# sets the directory to find the normalized data files
SPEC_DIREC = os.getcwd() + "/DATA/NORM_DR" + DR + "Q/" 

# creates directory for output files
OUT_DIREC = os.getcwd() + "/OUTPUT_FILES/ABSORPTION/"

# do you want to use smoothed norm flux/error
# boxcar_size must always be an odd integer
want_to_smooth = 'no' 
boxcar_size = 101  

# plot all cases or only those with absorption 
plot_all = 'yes'

# lower limit of absorption width to be flagged 
BALNICITY_INDEX_LIMIT = 2000 

# limits on velocity     min,   max
VELOCITY_LIMIT = Range(-30000, -60000.)

# range of spectra you are working with from the good_normalization.csv file
STARTS_FROM, ENDS_AT = 39, 47

# wavelength restframe range
WAVELENGTH_RESTFRAME = Range(1200., 1800.)

###############################################################################################################################
######################################## OUTPUT FILES #########################################################################

# set name of output .txt file with absorption values
ABSORPTION_VALUES = OUT_DIREC + "/" + "absorption_measurements_test.txt"

# set name of output pdf with plots 
ABSORPTION_OUTPUT_PLOT_PDF = PdfPages('absorption_BI' + str(BALNICITY_INDEX_LIMIT) + '_test.pdf') 

###############################################################################################################################
####################################### DO NOT CHANGE #########################################################################

# verner table data
wavelength_CIV_emit1=  1550.7700
wavelength_CIV_emit2 = 1548.1950
avr_CIV_doublet = 1549.0524 #weighted average
avr_SiIV_doublet = 1396.747 # weighted average; individuals: 1402.770, 1393.755
CII_emitted = 1335.313 # (weighted average); individuals:
OI_emitted = 1303.4951 # weighted average; individuals pag 20 in Verner Table
avr_NV_doublet = 1240.15 # weighted average; individuals: 1242.80, 1238.82
avr_OVI_doublet = 1033.8160 # weighted average; individuals: 1037.6167, 1031. 9261

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

brac_all, delta_v_all = [], []
absspeccount = 0
count = 0
BI = 0
vmins, vmaxs, vmins_all, vmaxs_all = [], [], [], [] # v = velocity
final_depth_individual, final_depth_all_individual = [], []
BI_total, BI_ind_sum, BI_individual, BI_all_individual, BI_mid = [], [], [], [], []
EW_individual, EW_all_individual, vlast = [], [], [] #EW = equivalent width

###############################################################################################################################
######################################### MAIN CODE ###########################################################################

# clear files
if __name__ == "__main__":
    clear_file(ABSORPTION_VALUES)
    #clear_file(ABSORPTION_OUTPUT_PLOT) # possibly don't need to to clear pdf, check when runs

# read list of normalized spectra, zem, and calculated snr
norm_spectra_list, redshift_list, calc_snr_list = read_list_spectra(CONFIG_FILE, ["NORM SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR"])

# loops over each spectra
for spectra_index in range(STARTS_FROM, ENDS_AT + 1):
    
    # rounding the numbers of the redshift, calculated snr and setting the norm file name to the current file name
    z = round(redshift_list[spectra_index - 1], 5)
    calc_snr = round(calc_snr_list[spectra_index - 1], 5)
    current_spectrum_file_name = norm_spectra_list[spectra_index - 1]
    
    print(str(spectra_index), "current spectra file name: ", current_spectrum_file_name)
    current_spectra_data = np.loadtxt(SPEC_DIREC + current_spectrum_file_name)

    # reading wavelength, normalized flux, and normalized error from NORM_DRXQ of that paticular spectra
    wavelength, normalized_flux, normalized_error = read_spectra(current_spectra_data)

    # initializing observed wavelength values
    wavelength_observed_from = (z + 1) * WAVELENGTH_RESTFRAME.start
    wavelength_observed_to = (z + 1) * WAVELENGTH_RESTFRAME.end

    # include if statement for smoothing and smooth spectrum.
    if want_to_smooth == 'yes':
        sm_flux = smooth(normalized_flux, boxcar_size)
        sm_error = smooth(normalized_error, boxcar_size) / np.sqrt(boxcar_size)   
        non_sm_flux = normalized_flux
        non_sm_error = normalized_error
        normalized_flux = sm_flux
        normalized_error = sm_error

    # transform the wavelength array to velocity (called "beta") based on the CIV doublet: 
    beta = wavelength_to_velocity(z, wavelength)

    ################################# INITIALIZING  VARIABLES IN LOOP ##############################################################
    #index_depth_final, flux_depth, final_depth_individual = []
    non_trough_count = 999 # arbitrary large number that we will never reach

    delta_v = 0 #change in velocity
    sum_of_deltas = 0
    bb = -1 #???
                        
    countvmins = 0 # variable initialization to get into vmin/vmax loop
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
        
        # Handle 3-point spike limit
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
            EW_individual.append(EW)       

            # BI calculation
            if sum_of_deltas >= BALNICITY_INDEX_LIMIT: # passing the BALNICITY_INDEX_LIMIT (in this case 2,000 km/s) threshold
                C = 1  #set to 1 only if square bracket is continuously positive over a velocity interval            
                BI = (bracket * C) * (delta_v) #Calculate BAL for this delta_v
                BI_mid.append(np.round(BI, 5)) #Append to intermediate results
                BI_individual.append(np.round(BI, 5)) 

                #if non_trough_count == 0: # plotting the black line
                   #plt.plot((beta[current_velocity_index + 1], beta[current_velocity_index]), (1.5,1.5),'k-')

                # vmin calculation               
                if countvmins == 0 and non_trough_count == 0:  
                    vmins_index = np.min(np.where(beta >= (beta[current_velocity_index] + BALNICITY_INDEX_LIMIT)))  # vmins occurs current beta plus countBI
                    vmins.append(np.round(beta[vmins_index], 5))                    
                
                    # Calculate where CIV, CII and OI would be for each pair of vmin and vmax *if* the EHVO absorption found were 
                    # instead not EHVO and due to SiIV: 

                    # If the absorption is SiIV, this finds and plots where CIV, CII and OI would be
                    z_absSiIV = (wavelength[current_velocity_index] / avr_SiIV_doublet) - 1
                     
                    obs_wavelength_C = (z_absSiIV + 1) * (avr_CIV_doublet)
                    obs_wavelength_C_index = np.min(np.where(wavelength > obs_wavelength_C))
                    obs_wavelength_C_vel = beta[obs_wavelength_C_index] + BALNICITY_INDEX_LIMIT

                    obs_wavelength_CII = (z_absSiIV + 1) * (CII_emitted)
                    obs_wavelength_CII_index = np.min(np.where(wavelength > obs_wavelength_CII))                  
                    obs_wavelength_CII_vel = beta[obs_wavelength_CII_index] + BALNICITY_INDEX_LIMIT

                    obs_wavelength_OI = (z_absSiIV + 1) * (OI_emitted)
                    obs_wavelength_OI_index = np.min(np.where(wavelength > obs_wavelength_OI))                  
                    obs_wavelength_OI_vel = beta[obs_wavelength_OI_index] + BALNICITY_INDEX_LIMIT

                    countvmins = 1

    draw_abs_figure(beta, normalized_flux, normalized_error, ABSORPTION_OUTPUT_PLOT_PDF, current_spectrum_file_name, z, calc_snr, obs_wavelength_C_vel, obs_wavelength_CII_vel, obs_wavelength_OI_vel, vmins)
'''      
# ****************************************** NEXT UP ******************************************  
    # Plot figure as if the absorption was SiIV, CII or OI (we will visually inspect this file). 

    if (len(vmaxs) != 0) or (plotall == 'yes'):

        plt.title('Normalized flux vs velocity fhfhkfs')
        plt.xlabel('Velocity (km/s)')
        plt.ylabel('Normalized flux')
        plt.plot((np.min(beta),np.max(beta)),(1,1))
        plt.plot((np.min(beta),np.max(beta)),(0.9,0.9),'r--')
        plt.plot(beta, normalized_flux, 'k-')
    #        plot(beta, sm_flux, 'r-')
    #        plot(beta, norm_flux, 'k-')
        plt.plot(beta, normalized_error,'k--')
        #plot (beta, norm_error,'k--')
        #plot (beta, sm_error,'k--')
        plt.ylim(0,3)

        plt.xlim(-70000,0)
        plt.text(-60000, 2, str(i)+',     z='+str(z)+' snr='+ str(snr), rotation=0, fontsize=9) # <-- Different location?

    # the axvspan using beta need to be defined earlier to plot here f.e., axvspan(beta[current_velocity_index],beta[current_velocity_index-1], alpha=0.05, color='red')

    # Save absorption info in file. 

    if (len(vmaxs) != 0) or (plotall == 'yes'): # < -- I am confused about the "or" there
        absspeccount=absspeccount+1  # <-- I think this is a variable to show how many cases have absorption? 
        yes=(str(count)+';'+str(absspeccount)+'  name: ' + str(i) + '\n' + 'BI (-30000 > v > -60,000): ' + str(BI_total) + '\n' +  'vmins: ' + str(vmins) + '\n' + 'vmaxs: '+str(vmaxs) + '\n' + 'BI_individual: '+ str(BI_individual) + '\n' + 'EW_individual: '+ str(EW_individual) + '\n' + 'Depth: '+ str(final_depth_individual) +'\n'+'\n')
        vlast.append(yes)
    
        pp.savefig()

        s = i.split('-')
        plateid = s[1]
        mjd = s[2]
        s = s[3].split('.')
        fiber = s[0]
        plt.savefig(output_spec + plateid + "-" + mjd + "-" + fiber + ".png") # <-- I don't think we really want a million png... only if doing a single one this makes sense. 

    close(count)

    if (len(vmaxs) != 0) or (plotall == 'yes'):  #<-- Again, I am confused about what the if includes plotall. 
        vmins_all.append(vmins)
        vmaxs_all.append(vmaxs)
        
# Clean below this. Do we need it? Where is it saving it?

BI_total= array(BI_total)

vmins = array(vmins)
vmaxs = array(vmaxs)
pp.close()
vmaxs_final=[]
vmins_final=[]

for loop in range (0, len (vmaxs_all)):
    vmaxs_final.append (str(vmaxs_all[loop])+ ',' )

for loop2 in range (0, len(vmins_all)):
    vmins_final.append (str(vmins_all[loop2])+ ',' )
                    

vmaxs_final = array(vmaxs_final)
vmins_final = array(vmins_final)    
savetxt(ffile,vlast,fmt='%s')
'''
ABSORPTION_OUTPUT_PLOT_PDF.close()