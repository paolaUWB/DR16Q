#####################################################################################################################################
#   Normalization of the quasar spectra
#
# This code normalizes the DRQ spectra with a new algorithm.

# Authors: Paola Rodriguez Hidalgo, Mikel Charles, Wendy Garcia Naranjo, Cora DeFrancesco, Daria K, Can Tosun, David Nguyen, Sean Haas, Abdul Khatri
#
######################################################################################################################################

#############################################################################################
########################################## IMPORTS ##########################################

import os
import sys
import numpy as np 
import pandas as pd
from matplotlib import pyplot as plt
from numpy.lib.function_base import append
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
from utility_functions import print_to_file, clear_file, append_row_to_csv, read_file, read_spectra
from data_types import Range, RangesData, FigureData, FigureDataOriginal, FlaggedSNRData, ColumnIndexes
from useful_wavelength_flux_error_modules import wavelength_flux_error_for_points, wavelength_flux_error_in_range, calculate_snr
from draw_figures import powerlaw, draw_dynamic, draw_dynamic_points, draw_original_figure, draw_normalized_figure
from scipy import signal
import time 
start_time = time.time()
########################################## SPHINX ###########################################
"""
normalization
=============
Normalization module for this project
"""
#############################################################################################

######################################### VARIABLES ######################################### 

DR = '16' ## INPUT WHICH DATA RELEASE YOU ARE WORKING WITH [INPUT NUMBER ONLY i.e. '9']

NORM_FILE_EXTENSION = "norm.dr" + DR

## DEFINES THE CONFIG FILE
CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else "DR" + DR + "_sorted_norm_manual_fit.csv"

## SETS THE DIRECTORY TO FIND THE DATA FILES (DR9, DR16)
SPEC_DIREC = os.getcwd() + "/DATA/DR" + DR + "Q_SNR10/" 

## CREATES DIRECTORY FOR OUTPUT FILES
OUT_DIREC = os.getcwd() + "/OUTPUT_FILES/NORMALIZATION/"

## SETS THE DIRECTORY TO STORE NORMALIZED FILES
NORM_DIREC = os.getcwd() + '/../' + "NORM_DR16Q/"

## RANGE OF SPECTRA YOU ARE WORKING WITH FROM THE DRX_sorted_norm.csv FILE. 
STARTS_FROM, ENDS_AT = 1, 10  ## [1-10, 899-1527 for dr9] [1-21823 [21851 for dr16 (21852-21859 are high redshift cases - must set dynamic = yes to run)] 

## CUTOFF FOR SNR VALUES TO BE FLAGGED; FLAGS VALUES SMALLER THAN THIS
SNR_CUTOFF = 10. 

save_new_output_file = 'yes' ## DO YOU WANT TO SAVE TO THE OUTPUT FILES? 'yes'/'no'
save_new_norm_file = 'yes' ## DO YOU WANT TO CREATE NEW NORM.DRX FILES? 'yes'/'no'
save_figures = 'yes' ## DO YOU WANT TO SAVE PDF FILES OF GRAPHS? 'yes'/'no'

sm = 'no' ## DO YOU WANT TO SMOOTH? 'yes'/'no'

dynamic = 'yes' ## DO YOU WANT TO CHOOSE ANCHOR POINTS? 'yes'/'no'

flag_spectra = 'no' ## DO YOU WANT TO FLAG SPECTRA? 'yes'/'no' [CHANGE TO NO WHEN DYNAMIC PLOTTING]

## VALUE USED IN TEST 1
val1 = 0.15

## VALUE USED IN TEST 2
val2 = 0.05 

BOXCAR_SIZE = 11 ## MUST BE ODD

## INITIAL PARAMETERS OF POWERLAW
b = 1250 
c = -0.5

#############################################################################################
####################################### DO NOT CHANGE #######################################

## RANGES OF WAVELENGTHS IN THE SPECTRA
WAVELENGTH_RESTFRAME = Range(1200., 1800.)
WAVELENGTH_FOR_SNR = Range(1250., 1400.)
WAVELENGTH_RESTFRAME_FOR_LEFT_POINT = Range(1280., 1290.)
WAVELENGTH_RESTFRAME_FOR_MIDDLE_POINT = Range(1420., 1430.)
WAVELENGTH_RESTFRAME_FOR_RIGHT_POINT = Range(1690., 1710.)
WAVELENGTH_RESTFRAME_TEST_1 = Range(1315., 1325.)
WAVELENGTH_RESTFRAME_TEST_2 = Range(1350., 1360.)
WAVELENGTH_RESTFRAME_TESTS = Range(1650., 1700.) ## check these values???? 

#############################################################################################

#############################################################################################
######################################## OUTPUT FILES #######################################

## TEXT OUTPUT FILES
LOG_FILE = OUT_DIREC + "/" + "log.txt"
LOG_NO_LOW_SNR_FILE = OUT_DIREC + "/" + "log_no_low_snr.txt"

FLAGGED_SNR_FILE = OUT_DIREC + "/" + "flagged_snr_in_ehvo.csv"

FLAGGED_ABSORPTION_FILE = OUT_DIREC + "/" + "flagged_absorption.csv"
FLAGGED_BAD_FIT_FILE = OUT_DIREC + "/" + "flagged_bad_fit.csv"
GOOD_FIT_FILE = OUT_DIREC + "/" + "good_fit_manual.csv"
ORIGINAL_FILE = OUT_DIREC + "/" + "original.csv"
UNFLAGGED_FILE = OUT_DIREC + "/" + "unflagged.csv"

## CREATES PDF FOR GRAPHS
FLAGGED_ABSORPTION_PDF = PdfPages('flagged_absorption_graphs.pdf')
FLAGGED_BAD_FIT_PDF = PdfPages('flagged_bad_fit_graphs.pdf')
GOOD_FIT_PDF = PdfPages('good_fit_graphs.pdf')
ORIGINAL_PDF = PdfPages('original_graphs.pdf') 
UNFLAGGED_PDF = PdfPages('unflagged_graphs.pdf')

NORMALIZED_PDF = PdfPages('normalized_graphs.pdf') 

#############################################################################################
######################################### FUNCTIONS #########################################

def smooth(norm_flux, box_size):
    """ Function that smoothes the spectra.

    Parameters:
    -----------
    norm_flux: any
        Normalized flux
    box_size: any
        Boxcar size - how smoothed you want the spectra (higher # = more smooth). This value must be odd.
    
    Returns:
    --------
    y_smooth: any
        Returns the smoothed spectra - as the boxcar size increases, the plot will resemble a smooth line showing the trends of the spectra without extra "noise"
    """

    y_smooth = signal.savgol_filter(norm_flux,box_size,2)
    return y_smooth

def define_three_anchor_points(z: float, spectra_data):
    """ Defines the three anchor points used in the normalization graph.

    Parameters:
    -----------
    z: float
        Values from the data base of the redshift, DR9Q (for now..).
    spectra_data: list
        Current spectra data from files, DR9Q (for now...).
               
    Returns:
    --------
    tuple
        left_point, middle_point, right_point for wavelength_flux_error_for_points.
    """
    ## check these
    WAVELENGTH_OBSERVED_FOR_LEFT_POINT = Range(1280. * (1 + z), 1290. * (1 + z))
    WAVELENGTH_OBSERVED_FOR_MIDDLE_POINT = Range(1420. * (1 + z), 1430. * (1 + z))
    WAVELENGTH_OBSERVED_FOR_RIGHT_POINT = Range(1690. * (1 + z), 1710. * (1 + z))

    left_point = wavelength_flux_error_for_points(
        WAVELENGTH_OBSERVED_FOR_LEFT_POINT.start,
        WAVELENGTH_OBSERVED_FOR_LEFT_POINT.end,
        z,
        spectra_data)

    middle_point = wavelength_flux_error_for_points(
        WAVELENGTH_OBSERVED_FOR_MIDDLE_POINT.start,
        WAVELENGTH_OBSERVED_FOR_MIDDLE_POINT.end,
        z,
        spectra_data)
    
    try: 
        right_point = wavelength_flux_error_for_points(
            WAVELENGTH_OBSERVED_FOR_RIGHT_POINT.start,
            WAVELENGTH_OBSERVED_FOR_RIGHT_POINT.end,
            z,
            spectra_data)
    except:
        right_point = wavelength_flux_error_for_points(
            WAVELENGTH_OBSERVED_FOR_RIGHT_POINT_HIGH_REDSHIFT.start, 
            WAVELENGTH_OBSERVED_FOR_RIGHT_POINT_HIGH_REDSHIFT.end,
            z,
            spectra_data)
    
    return [left_point, middle_point, right_point]

def dynamic_find_anchor_points(spectra_data, number_of_anchor_points):
    """
    Function based on 'define_three_anchor_points'. Defines a user-specified 
        number of anchor points. This function makes use of the function
        'wavelength_flux_error_for_points' to find the closest wavelength bin
        to the user-requested anchor point values.

    Parameters
    ----------
    spectra_data : tuple
        (wavelength, flux, error).
    number_of_anchor_points: array
        [1,2,3,...,n] where n is the number of anchor points defined by the user
    z : float
        redshift.
    user_anchors : list
        User-defined desired anchor points.
    user_delta : float
        Wavelength range to search for anchor points.

    Returns
    -------
    anchor_pts : arr
        List of PointData objects.

    """
    anchor_pts = []
    for i in number_of_anchor_points:
        spec_point = wavelength_flux_error_for_points(wavelength_range[i-1][0], wavelength_range[i-1][1], z, spectra_data)
        anchor_pts.append(spec_point)
    print("User Requested Point: ", anchor_pts)
    return anchor_pts

### Masking points with large errors: 

###    for n in range(1, len(flux_normalized) - 5):  ### !!! Is this too big? 
###     if the change in flux is greater than 0.5 and the error of [n+1] is greater than 0.25   
###     reverts all of the [n+1] back to the n values     
###       if abs(flux_normalized[n + 1] - flux_normalized[n]) > 0.5:
###            if error_normalized[n + 1] > 0.25:
###                error_normalized[n + 1] = error_normalized[n]
###                flux_normalized[n + 1] = flux_normalized[n]
###                error[n + 1] = error[n]
###                flux[n + 1] = flux[n]
###        #if the error is larger than 0.5 then it reverts n back to the [n-1]
###        if error_normalized[n] > 0.5:
###            error_normalized[n] = error_normalized[n - 1]
###            flux_normalized[n] = flux_normalized[n - 1]
###            error[n] = error[n - 1]
###            flux[n] = flux[n - 1]
###        # if the change in flux is greater than 5 then [n+1] reverts back to n  
###        if abs(flux_normalized[n + 1] - flux_normalized[n]) > 5:
###            error_normalized[n + 1] = error_normalized[n]
###            flux_normalized[n + 1] = flux_normalized[n]
###            error[n + 1] = error[n]
###            flux[n + 1] = flux[n]

### --> Create a masking array for flux and error prior to plotting. 

#############################################################################################
######################################### MAIN CODE #########################################

if (__name__ == "__main__"):

    field = ["SPECTRA INDEX", "SPECTRA FILE NAME", "CHI SQUARED"]
    fields_snr = ["SPECTRA INDEX", "SPECTRA FILE NAME", "SDSS SNR", "CALCULATED SNR"]
    fields=["SPECTRA INDEX", "SPECTRA FILE NAME", "NORM SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR", "SDSS SNR", "BF", "CF"]
    
    if save_new_output_file == 'yes':
        clear_file(FLAGGED_ABSORPTION_FILE)
        clear_file(FLAGGED_BAD_FIT_FILE)
        clear_file(FLAGGED_SNR_FILE)
        clear_file(GOOD_FIT_FILE)
        clear_file(ORIGINAL_FILE)
        clear_file(UNFLAGGED_FILE)
        
        append_row_to_csv(FLAGGED_ABSORPTION_FILE, fields)
        append_row_to_csv(FLAGGED_BAD_FIT_FILE, fields)
        append_row_to_csv(FLAGGED_SNR_FILE, fields_snr)
        append_row_to_csv(GOOD_FIT_FILE, fields)
        append_row_to_csv(ORIGINAL_FILE, field)
        append_row_to_csv(UNFLAGGED_FILE, fields)
    
    
    clear_file(LOG_FILE)
    clear_file(LOG_NO_LOW_SNR_FILE)
    
redshift_value_list, snr_value_list, spectra_list = read_file(CONFIG_FILE)

indices, spectra_indices, processed_spectra_file_names, powerlaw_final_b_values, powerlaw_final_c_values = [], [], [], [], []
flagged_indices, flagged_spectra_indices, flagged_spectra_file_names = [], [], []
flagged_snr_indices, flagged_snr_spectra_indices, flagged_snr_spectra_file_names, flagged_snr_in_ehvo_values = [], [], [], []

for spectra_index in range(STARTS_FROM, ENDS_AT + 1):
    z = round(redshift_value_list[spectra_index - 1], 5)
    snr = round(snr_value_list[spectra_index - 1], 5)
    current_spectrum_file_name = spectra_list[spectra_index - 1]
    
    print(str(spectra_index) + ": " + current_spectrum_file_name)
    print_to_file(str(spectra_index) + ": " + current_spectrum_file_name, LOG_FILE)

    current_spectra_data = np.loadtxt(SPEC_DIREC + current_spectrum_file_name)
    
    ## DEFINING WAVELENGTH, FLUX, AND ERROR FOR WHOLE SPECTRA
    wavelength, flux, error = read_spectra(current_spectra_data)
    
    #######
    WAVELENGTH_OBSERVED_FOR_RIGHT_POINT_HIGH_REDSHIFT = Range(np.max(current_spectra_data[:, 0]) - 20., np.max(current_spectra_data[:, 0]))
    #######

    original_ranges = RangesData(wavelength, flux, error)

    ## DYNAMIC PLOTTING - USER INPUT PROVIDES NUMBER OF ANCHOR POINTS TO USE, THEIR LOCATION, AND HOW MUCH OF A RANGE THE ANCHOR POINT CAN BE PLACED IN
    if dynamic == 'yes':
        ### is there a better way to define these? using redshift maybe? 
        wavelength_observed_from = (z + 1) * WAVELENGTH_RESTFRAME.start #np.min(wavelength) #3000 #(1 + z) * 1800
        wavelength_observed_to = (z + 1) * WAVELENGTH_RESTFRAME.end #np.max(wavelength) #6000 #(1 + z) * 3000 

        test1 = wavelength_flux_error_in_range(WAVELENGTH_RESTFRAME_TEST_1.start, WAVELENGTH_RESTFRAME_TEST_1.end, z, current_spectra_data)
        test2 = wavelength_flux_error_in_range(WAVELENGTH_RESTFRAME_TEST_2.start, WAVELENGTH_RESTFRAME_TEST_2.end, z, current_spectra_data)
        max_peak = np.max(flux)
        
        draw_dynamic(wavelength, wavelength_observed_from, wavelength_observed_to, flux, test1, test2, max_peak)
        anchor_point = define_three_anchor_points(z, current_spectra_data)
        print("Anchor Points [OLD]: " + str(anchor_point[0][0]) + ", " + str(anchor_point[1][0]) + ", " + str(anchor_point[2][0]))
        number_of_anchor_points = int(input("How many anchor points would you like to use?: "))
        number_of_anchor_points = [x for x in range(1, number_of_anchor_points + 1)]
        try_again = 'no'
        
        while try_again == 'no':
            user_input_wavelength = []
            wavelength_range = []
            anchor_pts = []
            powerlaw_wavelength = [] 
            powerlaw_flux = []
    
            range_value = int(input("Specify a range of wavelengths you would like used to find an anchor point? (plus or minus this value from your wavelength): "))
            while range_value > 100:
                range_value = int(input("Please choose a range less than 100: "))

            for i in number_of_anchor_points:
                user_guess = int(input("Where would you like anchor point #" + str(i) + " to be?: "))
                user_input_wavelength.append(user_guess)
                range_of_wavelength = [user_input_wavelength[i - 1] - range_value, user_input_wavelength[i - 1] + range_value]
                wavelength_range.append(range_of_wavelength)

            anchor_point = dynamic_find_anchor_points(current_spectra_data, number_of_anchor_points)
            
            for i in number_of_anchor_points:
                anchor_pt = [anchor_point[i-1][0], anchor_point[i-1][1]]
                anchor_pts.append(anchor_pt)
                powerlaw_wavelength.append(anchor_pt[0])
                powerlaw_flux.append(anchor_pt[1])

            try:
                pars, covar = curve_fit(powerlaw, powerlaw_wavelength, powerlaw_flux, p0=[b, c], maxfev=10000)
            except:
                print("Error - curve_fit failed-1st powerlaw " + current_spectrum_file_name)
                print_to_file("Error - curve_fit failed-1st powerlaw " + current_spectrum_file_name, LOG_FILE)
                print_to_file("Error - curve fit failed-1st powerlaw " + current_spectrum_file_name, LOG_NO_LOW_SNR_FILE)
            
            bf, cf = pars[0], pars[1]

            flux_normalized = flux/powerlaw(wavelength, bf, cf)
            error_normalized = error/powerlaw(wavelength, bf, cf)
            
            snr_mean_in_ehvo = calculate_snr(wavelength, z, WAVELENGTH_FOR_SNR, error_normalized)

            norm_w_f_e = (wavelength, flux_normalized, error_normalized) 
            norm_w_f_e = (np.transpose(norm_w_f_e))  
            if save_new_norm_file == 'yes': np.savetxt(NORM_DIREC + current_spectrum_file_name[0:len(current_spectrum_file_name) - 11] + NORM_FILE_EXTENSION, norm_w_f_e)
            norm_spectrum_file_name = current_spectrum_file_name[0: len(current_spectrum_file_name) - 11] + NORM_FILE_EXTENSION

            draw_dynamic_points(spectra_index, wavelength, wavelength_observed_from, wavelength_observed_to, flux, test1, test2, number_of_anchor_points, anchor_pts, max_peak, bf, cf, z, snr, snr_mean_in_ehvo, current_spectrum_file_name, 'null') #ORIGINAL_PDF) 

            try_again = str(input("Are you happy with the fit? 'yes'/'no': "))

            fields = [spectra_index, current_spectrum_file_name, norm_spectrum_file_name, z, snr_mean_in_ehvo, snr, bf, cf]
            
        if try_again == 'yes':
            draw_dynamic_points(spectra_index, wavelength, wavelength_observed_from, wavelength_observed_to, flux, test1, test2, number_of_anchor_points, anchor_pts, max_peak, bf, cf, z, snr, snr_mean_in_ehvo, current_spectrum_file_name, 'null')#GOOD_FIT_PDF)
            if save_new_output_file == 'yes':    
                append_row_to_csv(GOOD_FIT_FILE, fields)

        power_law_data_x = powerlaw_wavelength
        power_law_data_y = powerlaw_flux

    else:
        anchor_point = define_three_anchor_points(z, current_spectra_data)

        ## THE THREE POINTS THAT THE POWER LAW WILL USE (POINTS C, B, AND A)
        power_law_data_x = (anchor_point[0].wavelength, anchor_point[1].wavelength, anchor_point[2].wavelength)
        power_law_data_y = (anchor_point[0].flux, anchor_point[1].flux, anchor_point[2].flux)

        wavelength_observed_from = (z + 1) * WAVELENGTH_RESTFRAME.start
        wavelength_observed_to = (z + 1) * WAVELENGTH_RESTFRAME.end
        
        try:
            pars, covar = curve_fit(powerlaw, power_law_data_x, power_law_data_y, p0=[b, c], maxfev=10000)
        except:
            print("Error - curve_fit failed-1st powerlaw " + current_spectrum_file_name)
            print_to_file("Error - curve_fit failed-1st powerlaw " + current_spectrum_file_name, LOG_FILE)
            print_to_file("Error - curve_fit failed-1st powerlaw " + current_spectrum_file_name, LOG_NO_LOW_SNR_FILE)

        bf, cf = pars[0], pars[1]
        flux_normalized = flux/powerlaw(wavelength, bf, cf)
        error_normalized = error/powerlaw(wavelength, bf, cf)
        snr_mean_in_ehvo = calculate_snr(wavelength, z, WAVELENGTH_FOR_SNR, error_normalized)

        norm_w_f_e = (wavelength, flux_normalized, error_normalized) 
        norm_w_f_e = (np.transpose(norm_w_f_e))  

    fields_snr = [spectra_index, current_spectrum_file_name, snr, snr_mean_in_ehvo]

    ## FLAGGING LOW SNR
    flagged_snr_mean_in_ehvo = False
    if (snr_mean_in_ehvo < SNR_CUTOFF):  
        flagged_snr_mean_in_ehvo = True
        print_to_file('     Flagged - low SNR', LOG_FILE)
        if save_new_output_file == 'yes': 
            append_row_to_csv(FLAGGED_SNR_FILE, fields_snr)
    
    else:
        print_to_file(str(spectra_index) + ": " + current_spectrum_file_name, LOG_NO_LOW_SNR_FILE)
    
    ## SMOOTHING ORIGINAL FIGURES
    if sm == 'yes':
        sm_flux = smooth(flux, BOXCAR_SIZE)
        sm_error = smooth(error, BOXCAR_SIZE) / np.sqrt(BOXCAR_SIZE)   
        non_sm_flux = flux
        non_sm_error = error
        flux = sm_flux
        error = sm_error

    ## SMOOTHING NORMALIZED FIGURES 
    if sm == 'yes':
        sm_flux_norm = smooth(flux_normalized, BOXCAR_SIZE)
        sm_error_norm = smooth(error_normalized, BOXCAR_SIZE) / np.sqrt(BOXCAR_SIZE)
        non_sm_flux_norm = flux_normalized
        non_sm_error_norm = error_normalized 
        flux_normalized = sm_flux_norm
        error_normalized = sm_error_norm

    #############################################################################################
    #################################### TESTING TWO REGIONS ####################################

    flagged = False # BAD FIT
    flagged_anchor_point_bad_fit = False # POWERLAW NOT CLOSE ENOUGH TO ANCHOR POINT
    
    # UNFLAGGING
    flagged_fit_too_high_green = False 
    flagged_fit_too_high_pink = False
    unflagged = False 
    
    # POSSIBLE ABSORPTION
    flagged_completely_below_green = False
    flagged_completely_below_pink = False
    absorption = False

    ##### CHECKING FIT OF CURVE FOR NORMALIZATION 
    ## TEST 1
    if flag_spectra == 'yes':
        index_wavelength_from = np.max(np.where(wavelength <= (z + 1) * WAVELENGTH_RESTFRAME_TESTS.start))
        index_wavelength_to = np.min(np.where(wavelength >= (z + 1) * WAVELENGTH_RESTFRAME_TESTS.end))
        median_normalized_flux = np.median(flux_normalized[index_wavelength_from : index_wavelength_to])
        val = val1 * median_normalized_flux ### why pick this value?

    if not flagged_snr_mean_in_ehvo and flag_spectra == 'yes':
        print_to_file("     Value: " + str(val), LOG_NO_LOW_SNR_FILE)
        
    if flag_spectra == 'yes':
        norm_anchor_pts = define_three_anchor_points(z, norm_w_f_e)

    fields=[spectra_index, current_spectrum_file_name, current_spectrum_file_name[0:len(current_spectrum_file_name) - 11] + NORM_FILE_EXTENSION, z, snr_mean_in_ehvo, snr, bf, cf]

    ## NORMALIZED FLUX AT ANCHOR POINTS
    if flag_spectra == 'yes':
        point_A_flux_norm = norm_anchor_pts[2][1]
        point_B_flux_norm = norm_anchor_pts[1][1]
        point_C_flux_norm = norm_anchor_pts[0][1]
    
    ### Shouldnt point_A_flux_norm be one value? why are we taking the median of it? 
    if flag_spectra == 'yes':    
        flagged_A = abs(1 - np.median(point_A_flux_norm)) >= val
        flagged_B = abs(1 - np.median(point_B_flux_norm)) >= val
        flagged_C = abs(1 - np.median(point_C_flux_norm)) >= val

        if flagged_A and flagged_B and flagged_C:
            flagged_anchor_point_bad_fit = True 

    ## GREEN TEST REGION
    test1 = wavelength_flux_error_in_range(WAVELENGTH_RESTFRAME_TEST_1.start, WAVELENGTH_RESTFRAME_TEST_1.end, z, current_spectra_data)
    normalized_flux_test_1 = test1.flux/powerlaw(test1.wavelength, bf, cf)
    
    ## PINK TEST REGION
    test2 = wavelength_flux_error_in_range(WAVELENGTH_RESTFRAME_TEST_2.start, WAVELENGTH_RESTFRAME_TEST_2.end, z, current_spectra_data)
    normalized_flux_test_2 = test2.flux/powerlaw(test2.wavelength, bf, cf)

    if not flagged_anchor_point_bad_fit and not flagged_snr_mean_in_ehvo and flag_spectra == 'yes':  
        ## TEST 2 
        flagged_by_green_region = abs(np.median(normalized_flux_test_1) - 1) >= val2
        if flagged_by_green_region:
            print("     flagged_by_green_region: ", flagged_by_green_region)
            print_to_file("     flagged_by_green_region: " + str(flagged_by_green_region), LOG_FILE)
            print_to_file("     flagged_by_green_region: " + str(flagged_by_green_region), LOG_NO_LOW_SNR_FILE)

        flagged_by_pink_region = abs(np.median(normalized_flux_test_2) - 1) >= val2
        if flagged_by_pink_region:
            print("     flagged_by_pink_region: ", flagged_by_pink_region)
            print_to_file("     flagged_by_pink_region: " + str(flagged_by_pink_region), LOG_FILE)
            print_to_file("     flagged_by_pink_region: " + str(flagged_by_pink_region), LOG_NO_LOW_SNR_FILE)

        if flagged_by_green_region and flagged_by_pink_region:
            flagged = True
            error_message = "       Flagging figure #" + str(spectra_index) + ", file name: " + current_spectrum_file_name
            print(error_message)
            print_to_file(error_message, LOG_FILE)
            print_to_file(error_message, LOG_NO_LOW_SNR_FILE)

        ## TEST 3
        wavelength_index_from_green_region = np.max(np.where(wavelength <= (z + 1) * WAVELENGTH_RESTFRAME_TEST_1.start))
        wavelength_index_to_green_region = np.min(np.where(wavelength >= (z + 1) * WAVELENGTH_RESTFRAME_TEST_1.end))
        wavelength_index_from_pink_region = np.max(np.where(wavelength <= (z + 1) * WAVELENGTH_RESTFRAME_TEST_2.start))
        wavelength_index_to_pink_region = np.min(np.where(wavelength >= (z + 1) * WAVELENGTH_RESTFRAME_TEST_2.end))

        median_flux_green_region = np.median(flux_normalized[wavelength_index_from_green_region : wavelength_index_to_green_region])
        median_flux_pink_region = np.median(flux_normalized[wavelength_index_from_pink_region : wavelength_index_to_pink_region])

        minimum_flux_green_region = np.min(flux_normalized[wavelength_index_from_green_region : wavelength_index_to_green_region])
        maximum_flux_green_region = np.max(flux_normalized[wavelength_index_from_green_region : wavelength_index_to_green_region])
        minimum_flux_pink_region = np.min(flux_normalized[wavelength_index_from_pink_region : wavelength_index_to_pink_region])
        maximum_flux_pink_region = np.max(flux_normalized[wavelength_index_from_pink_region : wavelength_index_to_pink_region])

        diff_flux_below_green_region = abs(median_flux_green_region - minimum_flux_green_region)
        diff_flux_below_pink_region = abs(median_flux_pink_region - minimum_flux_pink_region)
        diff_flux_above_green_region = abs(maximum_flux_green_region - median_flux_green_region)
        diff_flux_above_pink_region = abs(maximum_flux_pink_region - median_flux_pink_region)

        min_flux_green_region = diff_flux_below_green_region * 0.75
        min_flux_pink_region = diff_flux_below_pink_region * 0.75
        max_flux_green_region = median_flux_green_region + (diff_flux_above_green_region * 0.45)
        max_flux_pink_region = median_flux_pink_region + (diff_flux_above_pink_region * 0.45)

        print_to_file('max flux test - green = ' + str(max_flux_green_region) + ' max flux test - pink = ' + str(max_flux_pink_region), LOG_NO_LOW_SNR_FILE)
        print_to_file('min flux test - green = ' + str(min_flux_green_region) + ' min flux test - pink = ' + str(min_flux_pink_region), LOG_NO_LOW_SNR_FILE)

        print_to_file('avg flux test - green = ' + str(median_flux_green_region) + ' avg flux test - pink = ' + str(median_flux_pink_region), LOG_NO_LOW_SNR_FILE)

        print_to_file('percent green test below = ' + str(np.round(diff_flux_below_green_region/median_flux_green_region, 3)), LOG_NO_LOW_SNR_FILE)
        print_to_file('percent pink test below = ' + str(np.round(diff_flux_below_pink_region/median_flux_pink_region, 3)), LOG_NO_LOW_SNR_FILE)
        print_to_file('percent green above = ' + str(np.round(diff_flux_above_green_region/median_flux_green_region, 3)), LOG_NO_LOW_SNR_FILE)
        print_to_file('percent pink above = ' + str(np.round(diff_flux_above_pink_region/median_flux_pink_region, 3)), LOG_NO_LOW_SNR_FILE)

        if flagged and (max_flux_green_region >= 1 >= min_flux_green_region):
            flagged_fit_too_high_green = True
            print_to_file('     Failed Test 3 [green], percentage: ' + str(np.round(diff_flux_below_green_region/median_flux_green_region, 3)), LOG_FILE)
            print_to_file('     Failed Test 3 [green], percentage: ' + str(np.round(diff_flux_below_green_region/median_flux_green_region, 3)), LOG_NO_LOW_SNR_FILE)

        if flagged and (max_flux_pink_region >= 1 >= min_flux_pink_region):
            flagged_fit_too_high_pink = True
            print_to_file('     Failed Test 3 [pink], percentage: ' + str(np.round(diff_flux_below_pink_region/median_flux_pink_region, 3)), LOG_FILE)
            print_to_file('     Failed Test 3 [pink], percentage: ' + str(np.round(diff_flux_below_pink_region/median_flux_pink_region, 3)), LOG_NO_LOW_SNR_FILE)

        if flagged_fit_too_high_green and flagged_fit_too_high_pink:
            flagged = False
            unflagged = True

        ## TEST 4
        if 1 > maximum_flux_green_region:
            flagged_completely_below_green = True
            print_to_file('     Failed Test 4 [green]', LOG_FILE)
            print_to_file('     Failed Test 4 [green]', LOG_NO_LOW_SNR_FILE)

        if 1 > maximum_flux_pink_region:
            flagged_completely_below_pink = True
            print_to_file('     Failed Test 4 [pink]', LOG_FILE)
            print_to_file('     Failed Test 4 [pink]', LOG_NO_LOW_SNR_FILE)

        if (not flagged_snr_mean_in_ehvo) and (flagged_completely_below_green or flagged_completely_below_pink):
            flagged = False
            absorption = True

    elif not flagged_snr_mean_in_ehvo and flag_spectra == 'yes':
        flagged = True
        print_to_file('     Failed Test 1 [anchor points]', LOG_FILE)
        print_to_file('     Failed Test 1 [anchor points]', LOG_NO_LOW_SNR_FILE)
    
    #### CHECK THIS CHECK THIS CHECK THIS
    if flagged and absorption: 
        flags = ' - ABSORPTION / FLAGGED'
    elif flagged:
        flags = ' - FLAGGED BAD FIT'
    elif unflagged:
        flags = ' - UNFLAGGED'
    elif absorption:
        flags = ' - ABSORPTION / GOOD FIT'
    else:
        flags = ' - GOOD FIT'
        
    ## CHI SQUARED
    residuals_test1 = test1.flux - powerlaw(test1.wavelength, bf, cf)
    residuals_test2 = test2.flux - powerlaw(test2.wavelength, bf, cf)    
    residuals_test1_and_2 = np.concatenate([residuals_test1,residuals_test2])
    wavelength_tests_1_and_2 = np.concatenate([test1.wavelength, test2.wavelength])
    
    chi_sq = sum((residuals_test1_and_2**2)/powerlaw(wavelength_tests_1_and_2, bf, cf))
    
    norm_spectrum_file_name = current_spectrum_file_name[0: len(current_spectrum_file_name) - 11] + NORM_FILE_EXTENSION
    
    field = [spectra_index, current_spectrum_file_name, chi_sq]
    fields=[spectra_index, current_spectrum_file_name, norm_spectrum_file_name, z, snr_mean_in_ehvo, snr, bf, cf]

    if not flagged_snr_mean_in_ehvo:
        append_row_to_csv(ORIGINAL_FILE, field)

    ## SCALING GRAPHS
    if dynamic == 'yes':
        max_peak = np.max(flux)
        max_peak_norm = np.max(flux_normalized)
    else:
        middle_point_from = (z + 1) * WAVELENGTH_RESTFRAME_FOR_MIDDLE_POINT.start
        right_point_to = (z + 1) * WAVELENGTH_RESTFRAME_FOR_RIGHT_POINT.end

        min_wavelength = np.min(np.where(wavelength > middle_point_from))
        max_wavelength = np.max(np.where(wavelength < right_point_to))

        max_peak = np.max(flux[min_wavelength + 1 : max_wavelength + 1])
        max_peak_norm = np.max(flux_normalized[min_wavelength + 1 : max_wavelength + 1])
        
    figure_data = FigureData(current_spectrum_file_name, wavelength_observed_from, wavelength_observed_to, z, snr, snr_mean_in_ehvo)

    ## DRAWING FIGURES
    if flagged_snr_mean_in_ehvo:
        flaggedSNRdata = FlaggedSNRData(figure_data, bf, cf, power_law_data_x, power_law_data_y)
    elif dynamic == 'no':
        original_figure_data = FigureDataOriginal(figure_data, bf, cf, power_law_data_x, power_law_data_y)
        draw_original_figure(spectra_index, original_ranges, original_figure_data, test1, test2, wavelength_observed_from, wavelength_observed_to, max_peak, ORIGINAL_PDF, flags)
        if flagged:
            if save_figures == 'yes':
                draw_original_figure(spectra_index, original_ranges, original_figure_data, test1, test2, wavelength_observed_from, wavelength_observed_to, max_peak, FLAGGED_BAD_FIT_PDF, flags)
            if save_new_output_file == 'yes':
                append_row_to_csv(FLAGGED_BAD_FIT_FILE, fields)
        elif unflagged:
            if save_figures == 'yes':
                draw_original_figure(spectra_index, original_ranges, original_figure_data, test1, test2, wavelength_observed_from, wavelength_observed_to, max_peak, UNFLAGGED_PDF, flags)
                draw_original_figure(spectra_index, original_ranges, original_figure_data, test1, test2, wavelength_observed_from, wavelength_observed_to, max_peak, GOOD_FIT_PDF, flags)
                draw_normalized_figure(spectra_index, original_ranges, figure_data, flux_normalized, error_normalized, test1, test2, normalized_flux_test_1, normalized_flux_test_2, wavelength_observed_from, wavelength_observed_to, max_peak_norm, NORMALIZED_PDF)
            if save_new_output_file == 'yes':
                append_row_to_csv(UNFLAGGED_FILE, fields)
                append_row_to_csv(GOOD_FIT_FILE, fields)
        elif absorption:
            if save_figures == 'yes':
                draw_original_figure(spectra_index, original_ranges, original_figure_data, test1, test2, wavelength_observed_from, wavelength_observed_to, max_peak, FLAGGED_ABSORPTION_PDF, flags)
                draw_original_figure(spectra_index, original_ranges, original_figure_data, test1, test2, wavelength_observed_from, wavelength_observed_to, max_peak, GOOD_FIT_PDF, flags)
                draw_normalized_figure(spectra_index, original_ranges, figure_data, flux_normalized, error_normalized, test1, test2, normalized_flux_test_1, normalized_flux_test_2, wavelength_observed_from, wavelength_observed_to, max_peak_norm, NORMALIZED_PDF)
            if save_new_output_file == 'yes':
                append_row_to_csv(GOOD_FIT_FILE, fields)
                append_row_to_csv(FLAGGED_ABSORPTION_FILE, fields)
        else:
            if save_figures == 'yes':
                draw_normalized_figure(spectra_index, original_ranges, figure_data, flux_normalized, error_normalized, test1, test2, normalized_flux_test_1, normalized_flux_test_2, wavelength_observed_from, wavelength_observed_to, max_peak_norm, NORMALIZED_PDF)
                draw_original_figure(spectra_index, original_ranges, original_figure_data, test1, test2, wavelength_observed_from, wavelength_observed_to, max_peak, GOOD_FIT_PDF, flags)
            if save_new_output_file == 'yes':
                append_row_to_csv(GOOD_FIT_FILE, fields)

    if dynamic == 'yes':
        original_figure_data = FigureDataOriginal(figure_data, bf, cf, power_law_data_x, power_law_data_y)
        draw_original_figure(spectra_index, original_ranges, original_figure_data, test1, test2, wavelength_observed_from, wavelength_observed_to, max_peak, GOOD_FIT_PDF, flags)
        draw_normalized_figure(spectra_index, original_ranges, figure_data, flux_normalized, error_normalized, test1, test2, normalized_flux_test_1, normalized_flux_test_2, wavelength_observed_from, wavelength_observed_to, max_peak_norm, NORMALIZED_PDF)

    if save_new_norm_file == 'yes' and not flagged_snr_mean_in_ehvo and not flagged: np.savetxt(NORM_DIREC + current_spectrum_file_name[0:len(current_spectrum_file_name) - 11] + NORM_FILE_EXTENSION, norm_w_f_e)

    ## OLD END OF PROCESS... 

    # add condition here?
    powerlaw_final_b_values.append(bf)
    powerlaw_final_c_values.append(cf)
    processed_spectra_file_names.append(current_spectrum_file_name)
    indices.append(spectra_index - STARTS_FROM + 1)
    spectra_indices.append(spectra_index)
    if flagged: 
        flagged_spectra_file_names.append(current_spectrum_file_name)
        flagged_indices.append(spectra_index - STARTS_FROM + 1)
        flagged_spectra_indices.append(spectra_index)

    if flagged_snr_mean_in_ehvo:
        flagged_snr_spectra_file_names.append(current_spectrum_file_name)
        flagged_snr_indices.append(spectra_index - STARTS_FROM + 1)
        flagged_snr_spectra_indices.append(spectra_index)
        flagged_snr_in_ehvo_values.append(snr_mean_in_ehvo)

flagged_snr_in_ehvo_graphs = [flagged_snr_spectra_indices, flagged_snr_spectra_file_names, flagged_snr_in_ehvo_values]
flagged_snr_in_ehvo_graphs = (np.transpose(flagged_snr_in_ehvo_graphs))

ORIGINAL_PDF.close()
NORMALIZED_PDF.close()
FLAGGED_BAD_FIT_PDF.close()
UNFLAGGED_PDF.close()
GOOD_FIT_PDF.close()
FLAGGED_ABSORPTION_PDF.close()

print("--- %s seconds" %(time.time()-start_time))
   
