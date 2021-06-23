#####################################################################################################################################
#   Normalization of the quasar spectra
#
# This code normalizes the DRQ spectra with a new algorithm.

# Authors: Paola Rodriguez Hidalgo, Mikel Charles, Wendy Garcia Naranjo, Daria K, Can Tosun, David Nguyen, Sean Haas, Abdul Khatri
#
######################################################################################################################################

#############################################################################################
########################################## IMPORTS ##########################################

import os
import sys
import numpy as np 
from matplotlib import pyplot as plt
from numpy.lib.function_base import append
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
from utility_functions import print_to_file, clear_file, append_row_to_csv
from data_types import Range, RangesData, FigureData, FigureDataOriginal, FlaggedSNRData  ###, DataNormalized
from useful_wavelength_flux_error_modules import wavelength_flux_error_for_points, wavelength_flux_error_in_range, calculate_snr
from file_reader import read_file
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

## READS FILES WITH QUASAR NAMES
CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else "DR" + DR + "_sorted_norm.csv"

## SETS THE DIRECTORY TO FIND THE DATA FILES (DR9, DR16)
SPEC_DIREC = os.getcwd() + "/DATA/DR" + DR + "Q_SNR10/" 

## SETS THE DIRECTORY TO STORE NORMALIZED FILES
NORM_DIREC = os.getcwd() + "/DATA/NORM_DR" + DR + "Q/"

STARTS_FROM, ENDS_AT = 1, 1000 ## [899-1527 for dr9] [1-21823? for dr16] RANGE OF SPECTRA YOU ARE WORKING WITH FROM THE DRX_sorted_norm.csv FILE. 

SNR_CUTOFF = 10. ## CUTOFF FOR SNR VALUES TO BE FLAGGED; FLAGS VALUES SMALLER THAN THIS

sm = 'no' ## DO YOU WANT TO SMOOTH? 'yes'/'no'

BOXCAR_SIZE = 11 ## MUST BE ODD

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

#############################################################################################

#############################################################################################
######################################## OUTPUT FILES #######################################

## CREATES DIRECTORY FOR OUTPUT FILES
OUT_DIREC = os.getcwd() + "/OUTPUT_FILES/"

LOG_FILE = OUT_DIREC + "/" + "log.txt"
FLAGGED_GRAPHS_FILE = OUT_DIREC + "/" + "flagged_bad_fit.csv"
FLAGGED_SNR_GRAPHS_FILE = OUT_DIREC + "/" + "flagged_snr_in_ehvo_graphs.txt"
FLAGGED_ABSORPTION = OUT_DIREC + "/" + "flagged_absorption.csv"
GOOD_NORMALIZATION_FLAGGED_FILE = OUT_DIREC + "/" + "good_normalization.csv"
GOODNESS_OF_FIT_FILE = OUT_DIREC + "/" + "chi_sq_values.csv"

## CREATES PDF FOR GRAPHS
ORIGINAL_PDF = PdfPages('original_graphs.pdf') 
NORMALIZED_PDF = PdfPages('normalized_graphs.pdf') 
FLAGGED_PDF = PdfPages('flagged_spectra.pdf') 
POWERLAW_TEST_PDF = PdfPages('powerlaw_test_graphs.pdf')


#############################################################################################
######################################### FUNCTIONS #########################################

b = 1250 # initial parameter of powerlaw
c = -0.5 # initial parameter of powerlaw


def powerlaw(wavelength, b, c):
    """ Calculates the power law. 

    Parameters:
    -----------
    wavelength: array
        Comes from RangesData().    
    b: int
        Initial parameter of powerlaw. 
    c: float
        Initial parameter of powerlaw.

    Returns:
    --------
    array
        Power law value in the form of an array.
    """
    return b * (np.power(wavelength, c))

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
    left_point = wavelength_flux_error_for_points(
        WAVELENGTH_RESTFRAME_FOR_LEFT_POINT.start,
        WAVELENGTH_RESTFRAME_FOR_LEFT_POINT.end,
        z,
        spectra_data)

    middle_point = wavelength_flux_error_for_points(
        WAVELENGTH_RESTFRAME_FOR_MIDDLE_POINT.start,
        WAVELENGTH_RESTFRAME_FOR_MIDDLE_POINT.end,
        z,
        spectra_data)
   
    right_point = wavelength_flux_error_for_points(
        WAVELENGTH_RESTFRAME_FOR_RIGHT_POINT.start,
        WAVELENGTH_RESTFRAME_FOR_RIGHT_POINT.end,
        z,
        spectra_data)
    
    return (left_point, middle_point, right_point)

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

def draw_original_figure(figure_index: int, original_ranges: RangesData, data: FigureDataOriginal, test1: RangesData, test2: RangesData, max_peak):
    """ Draws the original spectra graph.

    Parameters:
    -----------
    figure_index: int
        Makes a separate graph for each spectra. 
    original_ranges: RangesData
        Ranges of values for the original data.
    data: FigureDataOriginal
        Data from DR9Q (for now...).
    test1: RangesData
        Green highlighted area on graph. 
    test2: RangesData
        Pink highlighted area on graph.
    max_peak: any
        Max peak value of data per spectra.

    Returns:
    --------
    None.

    Note:
    -----
    Returns nothing, but draws the original spectra of the graph.
    """

    main_color = "xkcd:ultramarine"
    test_1_color, test_2_color = "xkcd:green apple", "xkcd:bubblegum"
    subtitle_text = f"z={data.FigureData.z} snr={data.FigureData.snr} snr_mean_in_ehvo={data.FigureData.snr_mean_in_ehvo}"
    plt.figure(figure_index)
    plt.title(data.FigureData.spectrum_file_name)
    plt.xlabel("Wavelength[A]")
    plt.ylabel("Flux[10^[-17]]cgs")
    plt.text(((data.FigureData.wavelength_from + data.FigureData.wavelength_to)/2.3), np.max(original_ranges.flux), subtitle_text)
    plt.plot(original_ranges.wavelength, original_ranges.flux, color = main_color, linestyle = "-")
    plt.plot(data.power_law_data_x, data.power_law_data_y, 'ro')
    plt.plot(original_ranges.wavelength, original_ranges.error, color = "black", linestyle = "-")
    plt.plot(test1.wavelength, test1.flux, color = test_1_color, linestyle = "-")
    plt.plot(test2.wavelength, test2.flux, color = test_2_color, linestyle = "-")
    plt.plot(original_ranges.wavelength, powerlaw(original_ranges.wavelength, data.bf, data.cf), color = "red", linestyle = "--")
    plt.ylim(-2, max_peak + 2)
    ORIGINAL_PDF.savefig()
    plt.close(figure_index)

def draw_normalized_figure(figure_index: int, original_ranges: RangesData, figure_data: FigureData, flux_normalized, error_normalized, #1 param more since I removed the tuple
                            test1: RangesData, test2: RangesData, normalized_flux_test_1, normalized_flux_test_2):
    """ Draws the normalized spectra graph.

    Parameters:
    -----------
    figure_index: int
        Makes a separate graph for each spectra. 
    original_ranges: RangesData
        Ranges of values for the original data.
    figure_data: FigureData
        Data from DR9Q (for now...).
    flux_normalized: array
    error_normalized: array
    test1: RangesData
        Green highlighted area on graph. 
    test2: RangesData
        Pink highlighted area on graph.
    normalized_flux_test_1: any
    normalized_flux_test_2: any

    Returns:
    --------
    None.
    
    Notes:
    ------
    Creates a graph of the spectra and saves to the original_graphs.pdf
    """

    main_color = "xkcd:ultramarine"
    test_1_color, test_2_color = "xkcd:green apple", "xkcd:bubblegum"
    subtitle_text = f"z={figure_data.z} snr={figure_data.snr} snr_mean_in_ehvo={figure_data.snr_mean_in_ehvo}"
    plt.figure(figure_index) 
    plt.text(((figure_data.wavelength_from + figure_data.wavelength_to)/2.3), np.max(flux_normalized)/1.07, figure_data.spectrum_file_name)
    plt.text(((figure_data.wavelength_from + figure_data.wavelength_to)/2.3), np.max(flux_normalized), subtitle_text)
    plt.title(figure_data.spectrum_file_name)
    plt.plot(original_ranges.wavelength, flux_normalized, color = main_color, linestyle = "-")
    plt.plot(original_ranges.wavelength, error_normalized, color = "black", linestyle = "-")
    plt.title("Normalized Data vs. Normalized Error")
    plt.xlabel("Wavelength [A]")
    plt.ylabel("Normalized Flux[10^[-17]]cgs")
    plt.plot(test1.wavelength, normalized_flux_test_1, color = test_1_color, linestyle = "-")
    plt.plot(test2.wavelength, normalized_flux_test_2, color = test_2_color, linestyle = "-")
    plt.plot((original_ranges.wavelength[0], original_ranges.wavelength[-1]), (1, 1), color = "red", linestyle = "-")
    plt.ylim(-2, np.max(flux_normalized) + 1.5)
    NORMALIZED_PDF.savefig()
    plt.close(figure_index)

def draw_flagged_figure(figure_index: int, original_ranges: RangesData, data: FigureDataOriginal, test1: RangesData, test2: RangesData, max_peak):
    """ Draws the spectra graphs for spectra flagged by `test1` and `test2`.

    Parameters:
    -----------
    figure_index: int
        Makes a separate graph for each spectra. 
    original_ranges: RangesData
        Ranges of values for the original data.
    data: FigureDataOriginal
        Data from DR9Q (for now...).
    test1: RangesData
        Green highlighted area on graph. 
    test2: RangesData
        Pink highlighted area on graph.
    max_peak: any
        Max peak value of data per spectra.

    Returns:
    --------
    None.

    Note:
    -----
    Creates a graph of the spectra and saves to the flagged_spectra.pdf
    """

    main_color = "xkcd:ultramarine"
    test_1_color, test_2_color = "xkcd:green apple", "xkcd:bubblegum"
    subtitle_text = f"z={data.FigureData.z} snr={data.FigureData.snr} snr_mean_in_ehvo={data.FigureData.snr_mean_in_ehvo}"
    plt.figure(figure_index)
    plt.title(data.FigureData.spectrum_file_name)
    plt.xlabel("Wavelength[A]")
    plt.ylabel("Flux[10^[-17]]cgs")
    plt.text(((data.FigureData.wavelength_from + data.FigureData.wavelength_to)/2.3), np.max(original_ranges.flux), subtitle_text)
    plt.plot(original_ranges.wavelength, original_ranges.flux, color = main_color, linestyle = "-")
    plt.plot(data.power_law_data_x, data.power_law_data_y, 'ro')
    plt.plot(original_ranges.wavelength, original_ranges.error, color = "black", linestyle = "-")
    plt.plot(test1.wavelength, test1.flux, color = test_1_color, linestyle = "-")
    plt.plot(test2.wavelength, test2.flux, color = test_2_color, linestyle = "-")
    plt.plot(original_ranges.wavelength, powerlaw(original_ranges.wavelength, data.bf, data.cf), color = "red", linestyle = "--")
    plt.ylim(-2, max_peak + 2)
    FLAGGED_PDF.savefig()
    plt.close(figure_index)

def draw_powerlaw_test_figure(figure_index: int, original_ranges: RangesData, data: FigureDataOriginal, test1: RangesData, test2: RangesData, max_peak):
    """ Draws the spectra graphs for spectra flagged by `test1` and `test2` that have a fit line that goes through the anchor points.

    Parameters:
    -----------
    figure_index: int
        Makes a separate graph for each spectra. 
    original_ranges: RangesData
        Ranges of values for the original data.
    data: FigureDataOriginal
        Data from DR9Q (for now...).
    test1: RangesData
        Green highlighted area on graph. 
    test2: RangesData
        Pink highlighted area on graph.
    max_peak: any
        Max peak value of data per spectra.

    Returns:
    --------
    None.

    Note:
    -----
    Creates a graph of the spectra and saves to the powerlaw_test_graphs.pdf
    """

    main_color = "xkcd:ultramarine"
    test_1_color, test_2_color = "xkcd:green apple", "xkcd:bubblegum"
    subtitle_text = f"z={data.FigureData.z} snr={data.FigureData.snr} snr_mean_in_ehvo={data.FigureData.snr_mean_in_ehvo}"
    plt.figure(figure_index)
    plt.title(data.FigureData.spectrum_file_name)
    plt.xlabel("Wavelength[A]")
    plt.ylabel("Flux[10^[-17]]cgs")
    plt.text(((data.FigureData.wavelength_from + data.FigureData.wavelength_to)/2.3), np.max(original_ranges.flux), subtitle_text)
    plt.plot(original_ranges.wavelength, original_ranges.flux, color = main_color, linestyle = "-")
    plt.plot(data.power_law_data_x, data.power_law_data_y, 'ro')
    plt.plot(original_ranges.wavelength, original_ranges.error, color = "black", linestyle = "-")
    plt.plot(test1.wavelength, test1.flux, color = test_1_color, linestyle = "-")
    plt.plot(test2.wavelength, test2.flux, color = test_2_color, linestyle = "-")
    plt.plot(original_ranges.wavelength, powerlaw(original_ranges.wavelength, data.bf, data.cf), color = "red", linestyle = "--")
    plt.ylim(-2, max_peak + 2)
    POWERLAW_TEST_PDF.savefig()
    plt.close(figure_index)

#############################################################################################
######################################### MAIN CODE #########################################

if __name__ == "__main__":
    clear_file(LOG_FILE)
    clear_file(FLAGGED_GRAPHS_FILE)
    clear_file(FLAGGED_SNR_GRAPHS_FILE)
    clear_file(FLAGGED_ABSORPTION)
    clear_file(GOOD_NORMALIZATION_FLAGGED_FILE)
    clear_file(GOODNESS_OF_FIT_FILE)

    field = ["spectra index", "spectra file name", "chi_sq"]
    fields=["spectra index", "spectra file name", "bf", "cf"]
    append_row_to_csv(GOODNESS_OF_FIT_FILE, field)
    append_row_to_csv(FLAGGED_GRAPHS_FILE, fields)
    append_row_to_csv(GOOD_NORMALIZATION_FLAGGED_FILE, fields)

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

    ## DEFINING WAVELENGTH, FLUX, AND ERROR (CHOOSING THEIR RANGE)
    wavelength, flux, error = wavelength_flux_error_in_range(WAVELENGTH_RESTFRAME.start, WAVELENGTH_RESTFRAME.end, z, current_spectra_data)

    def smooth(norm_flux, box_size):
        y_smooth = signal.savgol_filter(norm_flux,box_size,2)
        return y_smooth
        
    ## SMOOTHING ORIGINAL FIGURES
    if sm == 'yes':
        sm_flux = smooth(flux, BOXCAR_SIZE)
        sm_error = smooth(error, BOXCAR_SIZE) / np.sqrt(BOXCAR_SIZE)   
        non_sm_flux = flux
        non_sm_error = error
        flux = sm_flux
        error = sm_error

    wavelength_observed_from = (z + 1) * WAVELENGTH_RESTFRAME.start
    wavelength_observed_to = (z + 1) * WAVELENGTH_RESTFRAME.end

    left_point_from = (z + 1) * WAVELENGTH_RESTFRAME_FOR_LEFT_POINT.start
    middle_point_from = (z + 1) * WAVELENGTH_RESTFRAME_FOR_MIDDLE_POINT.start
    right_point_to = (z + 1) * WAVELENGTH_RESTFRAME_FOR_RIGHT_POINT.end

    point_C, point_B, point_A = define_three_anchor_points(z, current_spectra_data)

    ## THE THREE POINTS THAT THE POWER LAW WILL USE (POINTS C, B, AND A)
    power_law_data_x = (point_C.wavelength, point_B.wavelength, point_A.wavelength)
    power_law_data_y = (point_C.flux, point_B.flux, point_A.flux)

    try:
        pars, covar = curve_fit(powerlaw, power_law_data_x, power_law_data_y, p0=[b, c], maxfev=10000)
    except:
        print("Error - curve_fit failed-1st powerlaw " + current_spectrum_file_name)
        print_to_file("Error - curve_fit failed-1st powerlaw " + current_spectrum_file_name, LOG_FILE)

    bf, cf = pars[0], pars[1]

    flux_normalized = flux/powerlaw(wavelength, bf, cf)
    error_normalized = error/powerlaw(wavelength, bf, cf)

    ## SMOOTHING NORMALIZED FIGURES 
    if sm == 'yes':
        sm_flux_norm = smooth(flux_normalized, BOXCAR_SIZE)
        sm_error_norm = smooth(error_normalized, BOXCAR_SIZE) / np.sqrt(BOXCAR_SIZE)
        non_sm_flux_norm = flux_normalized
        non_sm_error_norm = error_normalized 
        flux_normalized = sm_flux_norm
        error_normalized = sm_error_norm

    original_ranges = RangesData(wavelength, flux, error)

    flagged_snr_mean_in_ehvo = False
    snr_mean_in_ehvo = calculate_snr(wavelength, z, WAVELENGTH_FOR_SNR, error_normalized)

    if snr_mean_in_ehvo < SNR_CUTOFF:  
        flagged_snr_mean_in_ehvo = True

    #############################################################################################
    #################################### TESTING TWO REGIONS ####################################

    ## CHECKING FIT OF CURVE FOR NORMALIZATION 
    flagged = False

    ## GREEN REGION
    test1 = wavelength_flux_error_in_range(WAVELENGTH_RESTFRAME_TEST_1.start, WAVELENGTH_RESTFRAME_TEST_1.end, z, current_spectra_data)
    normalized_flux_test_1 = test1.flux/powerlaw(test1.wavelength, bf, cf)
    
    ## PINK REGION
    test2 = wavelength_flux_error_in_range(WAVELENGTH_RESTFRAME_TEST_2.start, WAVELENGTH_RESTFRAME_TEST_2.end, z, current_spectra_data)
    normalized_flux_test_2 = test2.flux/powerlaw(test2.wavelength, bf, cf)

    ## AVERAGE VALUE OF POWERLAW IN TEST REGIONS
    powerlaw_test1 = powerlaw(test1.wavelength, bf, cf)
    powerlaw_test2 = powerlaw(test2.wavelength, bf, cf)

    ## AVERAGE FLUX VALUE IN TEST REGIONS
    avg_flux_test1 = np.average(test1.flux) + 0.05
    avg_flux_test2 = np.average(test2.flux) + 0.05

    ## MAX FLUX OF TEST REGIONS
    max_test1 = np.max(test1.flux)
    max_test2 = np.max(test2.flux)

    flagged_fit_1 = False
    flagged_fit_2 = False
    
    flagged_t1 = False
    flagged_t2 = False

    if np.average(powerlaw_test1) < avg_flux_test1:
        flagged_fit_1 = True

    if np.average(powerlaw_test2) < avg_flux_test2: 
        flagged_fit_2 = True

    if np.min(powerlaw_test1) > max_test1:
        flagged_t1 = True

    if np.min(powerlaw_test2) > max_test2:
        flagged_t2 = True

    if not flagged_snr_mean_in_ehvo and (flagged_t1 or flagged_t2):
        append_row_to_csv(FLAGGED_ABSORPTION, fields)

    flagged_by_test1 = abs(np.median(normalized_flux_test_1) - 1) >= 0.05
    if flagged_by_test1:
        print("     flagged_by_test1: ", flagged_by_test1)
        print_to_file("     flagged_by_test1: " + str(flagged_by_test1), LOG_FILE)

    flagged_by_test2 = abs(np.median(normalized_flux_test_2) - 1) >= 0.05
    if flagged_by_test2:
        print("     flagged_by_test2: ", flagged_by_test2)
        print_to_file("     flagged_by_test2: " + str(flagged_by_test2), LOG_FILE)


    if flagged_by_test1 and flagged_by_test2:
        flagged = True
        error_message = "       Flagging figure #" + str(spectra_index) + ", file name: " + current_spectrum_file_name
        print(error_message)
        print_to_file(error_message, LOG_FILE)
        point_A_powerlaw = powerlaw(point_A[0], bf, cf)
        point_B_powerlaw = powerlaw(point_B[0], bf, cf)
        point_C_powerlaw = powerlaw(point_C[0], bf, cf)
        point_powerlaw = str(spectra_index) + ": " + str(current_spectrum_file_name) + ", POINT A: " + str(point_A[1]) + ", POINT A PL:" + str(point_A_powerlaw) + ", POINT B: " + str(point_B[1]) + ", POINT B PL:" + str(point_B_powerlaw) + ", POINT C: " + str(point_C[1]) + ", POINT C PL:" + str(point_C_powerlaw)

    residuals_test1 = test1.flux - powerlaw(test1.wavelength, bf, cf)
    residuals_test2 = test2.flux - powerlaw(test2.wavelength, bf, cf)    
    residuals_test1_and_2 = np.concatenate([residuals_test1,residuals_test2])
    wavelength_tests_1_and_2 = np.concatenate([test1.wavelength, test2.wavelength])
    
    chi_sq = sum((residuals_test1_and_2**2)/powerlaw(wavelength_tests_1_and_2, bf, cf))

    field = [spectra_index, current_spectrum_file_name, chi_sq]
    fields=[spectra_index, current_spectrum_file_name, bf, cf]
    
    if not flagged_snr_mean_in_ehvo:
        append_row_to_csv(GOODNESS_OF_FIT_FILE, field)

    if (flagged_by_test1 and flagged_by_test2) and not flagged_snr_mean_in_ehvo:
        append_row_to_csv(FLAGGED_GRAPHS_FILE, fields)
    elif not flagged_snr_mean_in_ehvo:
        append_row_to_csv(GOOD_NORMALIZATION_FLAGGED_FILE, fields)

    ## SCALING GRAPHS
    wavelength_data = current_spectra_data[:,0]
    flux_data = current_spectra_data[:,1]

    min_wavelength = np.min(np.where(wavelength_data > middle_point_from))
    max_wavelength = np.max(np.where(wavelength_data < right_point_to))

    max_peak = np.max(flux_data[min_wavelength + 1 : max_wavelength + 1])

    figure_data = FigureData(current_spectrum_file_name, wavelength_observed_from, wavelength_observed_to, z, snr, snr_mean_in_ehvo)
    
    if flagged_snr_mean_in_ehvo:
        flaggedSNRdata = FlaggedSNRData(figure_data, bf, cf, power_law_data_x, power_law_data_y)
    else:
        original_figure_data = FigureDataOriginal(figure_data, bf, cf, power_law_data_x, power_law_data_y)
        draw_original_figure(spectra_index, original_ranges, original_figure_data, test1, test2, max_peak)
        if flagged:
            draw_flagged_figure(spectra_index, original_ranges, original_figure_data, test1, test2, max_peak)
            val = 0.5
            flagged_A = abs(point_A_powerlaw - point_A[1]) <= val
            flagged_C = abs(point_C_powerlaw - point_C[1]) <= val
            flagged_B = abs(point_B_powerlaw - point_B[1]) <= val
            if (flagged_A and flagged_B and flagged_C) and ((flagged_fit_1 and flagged_fit_2) or (flagged_t1 or flagged_t2)): 
                draw_powerlaw_test_figure(spectra_index, original_ranges, original_figure_data, test1, test2, max_peak)
                append_row_to_csv(GOOD_NORMALIZATION_FLAGGED_FILE, fields)
        else:
            draw_normalized_figure(spectra_index, original_ranges, figure_data, flux_normalized, error_normalized, test1, test2, normalized_flux_test_1, normalized_flux_test_2)


    norm_w_f_e = (wavelength, flux_normalized, error_normalized)
    norm_w_f_e = (np.transpose(norm_w_f_e))  
    np.savetxt(NORM_DIREC + current_spectrum_file_name[0:20] + NORM_FILE_EXTENSION, norm_w_f_e)

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

final_initial_parameters = [spectra_indices, processed_spectra_file_names, powerlaw_final_b_values, powerlaw_final_c_values]
final_initial_parameters = (np.transpose(final_initial_parameters))

flagged_graphs = [flagged_spectra_indices, flagged_spectra_file_names]
flagged_graphs = (np.transpose(flagged_graphs))

flagged_snr_in_ehvo_graphs = [flagged_snr_spectra_indices, flagged_snr_spectra_file_names, flagged_snr_in_ehvo_values]
flagged_snr_in_ehvo_graphs = (np.transpose(flagged_snr_in_ehvo_graphs))
    
ORIGINAL_PDF.close()
NORMALIZED_PDF.close()
FLAGGED_PDF.close()
POWERLAW_TEST_PDF.close()

np.savetxt(FLAGGED_SNR_GRAPHS_FILE, flagged_snr_in_ehvo_graphs, fmt='%s')

## what are we wanting printed to the file?
#else:
#    print(current_spectrum_file_name) ## later will print this to a file

print("--- %s seconds" %(time.time()-start_time))
   
