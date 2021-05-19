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
from scipy.optimize import curve_fit, leastsq
from matplotlib.backends.backend_pdf import PdfPages
from utility_functions import print_to_file, clear_file, append_row_to_csv
from data_types import Range, ColumnIndexes, PointData, RangesData, FigureData, FigureDataOriginal, DataNormalized, FlaggedSNRData
from useful_wavelength_flux_error_modules import wavelength_flux_error_for_points, wavelength_flux_error_in_range, calculate_snr
from file_reader import read_file
from scipy import signal
import time 
start_time = time.time()

#############################################################################################
######################################### VARIABLES ######################################### 

DR = 'dr9' ## Which data release #############################

NORM_FILE_EXTENSION = "norm." + DR

# Reads the file with the quasar names
CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else "sorted_norm.csv"

# How many columns the file with the quasar names has? # XXX Or is it the .dr9? PRH (sorted )
column_index = ColumnIndexes(0, 1, 2)

# Sets the directory to find the data files (dr9, dr16)
SPEC_DIREC = os.getcwd() + "/DATA/" ############################
#SPEC_DIREC = os.getcwd() + "/DATA/" # Set location of input and output spectrum files XXX Set a different one for input & output US LATER

STARTS_FROM, ENDS_AT = 1, 10 # Range of spectra you are working with from the quasar names file. 

SNR_CUTOFF = 10. # Cutoff for SNR values to be flagged; flags values smaller than this

sm = 'no' # Do you want to smooth? yes/no

BOXCAR_SIZE = 71 # Must be odd

# Ranges of wavelengths in the spectra for different tasks
WAVELENGTH_RESTFRAME = Range(1200., 1800.)
WAVELENGTH_FOR_SNR = Range(1250., 1400.)
WAVELENGTH_RESTFRAME_FOR_LEFT_POINT = Range(1280., 1290.)
WAVELENGTH_RESTFRAME_FOR_MIDDLE_POINT = Range(1420., 1430.)
WAVELENGTH_RESTFRAME_FOR_RIGHT_POINT = Range(1690., 1710.)
WAVELENGTH_RESTFRAME_TEST_1 = Range(1315., 1325.)
WAVELENGTH_RESTFRAME_TEST_2 = Range(1350., 1360.)


#############################################################################################
######################################## OUTPUT FILES #######################################

LOG_FILE = "log.txt"
FLAGGED_GRAPHS = SPEC_DIREC + "/" + "flagged_graphs.txt"
FINAL_INIT_PARAMS_FILE = SPEC_DIREC + "/" + "final_initial_parameters.txt"
PROCESSED_SPECTRA_FILE = SPEC_DIREC + "/" + "processed_spectra_filenames.txt"
FLAGGED_GRAPHS_FILE = SPEC_DIREC + "/" + "flagged_for_absorption_or_bad_normalization.txt"
FLAGGED_SNR_GRAPHS_FILE = SPEC_DIREC + "/" + "flagged_snr_in_ehvo_graphs.txt"
GOODNESS_OF_FIT_FILE = SPEC_DIREC + "/" + "chi_sq_values_all.csv"
BAD_NORMALIZATION_FLAGGED_FILE = SPEC_DIREC + "/" + "bad_normalization.csv"
GOOD_NORMALIZATION_FLAGGED_FILE = SPEC_DIREC + "/" + "good_normalization.csv"

ORIGINAL_PDF = PdfPages('original_graphs.pdf') # create pdf
NORMALIZED_PDF = PdfPages('normalized_graphs.pdf') # create pdf
FLAGGED_PDF = PdfPages('flagged_spectra.pdf') # create pdf


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
    plt.ylim(-2, max_peak + 3)
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
    NORMALIZED_PDF.savefig()
    plt.close(figure_index)

def draw_flagged_figure(figure_index: int, original_ranges: RangesData, data: FigureDataOriginal, test1: RangesData, test2: RangesData, max_peak):
    """ Draws the spectra graphs for spectra flagged by test1 and test2.

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
    plt.ylim(-2, max_peak + 3)
    FLAGGED_PDF.savefig()
    plt.close(figure_index)

#############################################################################################
######################################### MAIN CODE #########################################

if __name__ == "__main__":
    clear_file(LOG_FILE)
    clear_file(GOODNESS_OF_FIT_FILE)
    clear_file(BAD_NORMALIZATION_FLAGGED_FILE)
    clear_file(GOOD_NORMALIZATION_FLAGGED_FILE)
    clear_file(FLAGGED_GRAPHS)
    clear_file(FLAGGED_SNR_GRAPHS_FILE)
    print("Hi!")
    
    fields=["index", "spectra index", "chi_sq"] #index and spectra index - will they ever be different? causing a repeat of indexing
    append_row_to_csv(GOODNESS_OF_FIT_FILE, fields)
    append_row_to_csv(BAD_NORMALIZATION_FLAGGED_FILE, fields)
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

    # DEFINING WAVELENGTH, FLUX, AND ERROR (CHOOSING THEIR RANGE)
    wavelength, flux, error = wavelength_flux_error_in_range(WAVELENGTH_RESTFRAME.start, WAVELENGTH_RESTFRAME.end, z, current_spectra_data)
    #original_ranges = RangesData(wavelength, flux, error)

    def smooth(norm_flux, box_size):
        y_smooth = signal.savgol_filter(norm_flux,box_size,2)
        return y_smooth
        
    ### Smoothing original figures
    sm_flux = smooth(flux, BOXCAR_SIZE)
    sm_error = smooth(error, BOXCAR_SIZE) / np.sqrt(BOXCAR_SIZE)

    if sm == 'yes':
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

    # THE THREE POINTS THAT THE POWER LAW WILL USE (Points C, B, and A)
    power_law_data_x = (point_C.wavelength, point_B.wavelength, point_A.wavelength)
    power_law_data_y = (point_C.flux, point_B.flux, point_A.flux)

    try:
        pars, covar = curve_fit(powerlaw, power_law_data_x, power_law_data_y, p0=[b, c], maxfev=10000)
    except:
        print("Error - curve_fit failed-1st powerlaw " + current_spectrum_file_name)
        print_to_file("Error - curve_fit failed-1st powerlaw " + current_spectrum_file_name, LOG_FILE)

    #bf, cf are storing the variance for power_law_data_x and power_law_data_y respectively
    bf, cf = pars[0], pars[1]

    # flux_normalized & error_normalized are used to draw the figure
    flux_normalized = flux/powerlaw(wavelength, bf, cf)
    error_normalized = error/powerlaw(wavelength, bf, cf)

    ### Smoothing normalized figures 
    sm_flux_norm = smooth(flux_normalized, BOXCAR_SIZE)
    sm_error_norm = smooth(error_normalized, BOXCAR_SIZE) / np.sqrt(BOXCAR_SIZE)
    if sm == 'yes':
        non_sm_flux_norm = flux_normalized #to save original normalized flux in case we need it
        non_sm_error_norm = error_normalized 
        flux_normalized = sm_flux_norm #keeps variables consistent
        error_normalized = sm_error_norm

    original_ranges = RangesData(wavelength, flux, error)

    ## flagging spectra with low snr values, we want the high ones
    flagged_snr_mean_in_ehvo = False
    snr_mean_in_ehvo = calculate_snr(wavelength, z, WAVELENGTH_FOR_SNR, error_normalized)

    if snr_mean_in_ehvo < SNR_CUTOFF:  
        flagged_snr_mean_in_ehvo = True

    #############################################################################################
    #################################### TESTING TWO REGIONS ####################################
    ### Checking the fit of the curve for the normalization 
    ## Green and Pink regions in original_graphs.pdf
    flagged = False
    # Green Region
    test1 = wavelength_flux_error_in_range(WAVELENGTH_RESTFRAME_TEST_1.start, WAVELENGTH_RESTFRAME_TEST_1.end, z, current_spectra_data)
    normalized_flux_test_1 = test1.flux/powerlaw(test1.wavelength, bf, cf)
    
    #Pink Region
    test2 = wavelength_flux_error_in_range(WAVELENGTH_RESTFRAME_TEST_2.start, WAVELENGTH_RESTFRAME_TEST_2.end, z, current_spectra_data)
    normalized_flux_test_2 = test2.flux/powerlaw(test2.wavelength, bf, cf)


    flagged_by_test1 = abs(np.median(normalized_flux_test_1) - 1) >= 0.05  ## We tested several values
    if flagged_by_test1:
        print("flagged_by_test1: ", flagged_by_test1)
        print_to_file("flagged_by_test1: " + str(flagged_by_test1), LOG_FILE)

    
    flagged_by_test2 = abs(np.median(normalized_flux_test_2) - 1) >= 0.05
    if flagged_by_test2:
        print("flagged_by_test2: ", flagged_by_test2)
        print_to_file("flagged_by_test2: " + str(flagged_by_test2), LOG_FILE)

    if flagged_by_test1 and flagged_by_test2:
        flagged = True
        error_message = "Flagging figure #" + str(spectra_index) + ", file name: " + current_spectrum_file_name
        print(error_message)
        print_to_file(error_message, LOG_FILE)
        print_to_file(error_message, FLAGGED_GRAPHS)

    ##### residuals: the quantity remaining after other values have been subtracted from it
    ## flux observed, wavelength is the expected value
    # how off the flux is from the values of the power law at that location
    residuals_test1 = test1.flux - powerlaw(test1.wavelength, bf, cf)
    residuals_test2 = test2.flux - powerlaw(test2.wavelength, bf, cf)    
    residuals_test1_and_2 = np.concatenate([residuals_test1,residuals_test2])
    wavelength_tests_1_and_2 = np.concatenate([test1.wavelength, test2.wavelength])
    
    ### chi squared is comparing flux and wavelength
    chi_sq = sum((residuals_test1_and_2**2)/powerlaw(wavelength_tests_1_and_2, bf, cf))
    r_squared = 1 - chi_sq
    print("Chi Squared = ", chi_sq)
    print("R Squared = ", r_squared)

    fields=[spectra_index - STARTS_FROM + 1, spectra_index, chi_sq]
    append_row_to_csv(GOODNESS_OF_FIT_FILE, fields)

    # if chi squared is greater than 8 and meets both flagged tests it is added to bad normalization file
    ########### IS THIS DOING ANYTHING???? ###################
    if chi_sq > 8 and flagged_by_test1 and flagged_by_test2:
        append_row_to_csv(BAD_NORMALIZATION_FLAGGED_FILE, fields)
    else:
        append_row_to_csv(GOOD_NORMALIZATION_FLAGGED_FILE, fields)
    

    ### Finding the highest peak between the middle and right anchor points to scale graphs (in y-direction)
    wavelength_data = current_spectra_data[:,0]
    flux_data = current_spectra_data[:,1]

    min_wavelength = np.min(np.where(wavelength_data > middle_point_from))
    max_wavelength = np.max(np.where(wavelength_data < right_point_to))

    max_peak = np.max(flux_data[min_wavelength + 1 : max_wavelength + 1])

    figure_data = FigureData(current_spectrum_file_name, wavelength_observed_from, wavelength_observed_to, z, snr, snr_mean_in_ehvo)
    
    # Removes low SNR from graphs
    if flagged_snr_mean_in_ehvo:
        flaggedSNRdata = FlaggedSNRData(figure_data, bf, cf, power_law_data_x, power_law_data_y)
    else:
        original_figure_data = FigureDataOriginal(figure_data, bf, cf, power_law_data_x, power_law_data_y)
        draw_original_figure(spectra_index, original_ranges, original_figure_data, test1, test2, max_peak)
        draw_normalized_figure(spectra_index, original_ranges, figure_data, flux_normalized, error_normalized, test1, test2, normalized_flux_test_1, normalized_flux_test_2)
        if flagged:
            draw_flagged_figure(spectra_index, original_ranges, original_figure_data, test1, test2, max_peak)

    norm_w_f_e = (wavelength, flux_normalized, error_normalized)
    norm_w_f_e = (np.transpose(norm_w_f_e))  
    np.savetxt(SPEC_DIREC + current_spectrum_file_name[0:20] + NORM_FILE_EXTENSION, norm_w_f_e)

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

final_initial_parameters = [indices, spectra_indices, processed_spectra_file_names, powerlaw_final_b_values, powerlaw_final_c_values]
final_initial_parameters = (np.transpose(final_initial_parameters))

flagged_graphs = [flagged_indices, flagged_spectra_indices, flagged_spectra_file_names]
flagged_graphs = (np.transpose(flagged_graphs))

flagged_snr_in_ehvo_graphs = [flagged_snr_indices, flagged_snr_spectra_indices, flagged_snr_spectra_file_names, flagged_snr_in_ehvo_values]
flagged_snr_in_ehvo_graphs = (np.transpose(flagged_snr_in_ehvo_graphs))
flagged_snr_in_ehvo_graphs = flagged_snr_in_ehvo_graphs[flagged_snr_in_ehvo_graphs[:,3].argsort()] # sort by snr_mean_in_ehvo column
    
ORIGINAL_PDF.close()
NORMALIZED_PDF.close()
FLAGGED_PDF.close()

np.savetxt(FINAL_INIT_PARAMS_FILE, final_initial_parameters, fmt="%s")
np.savetxt(PROCESSED_SPECTRA_FILE, processed_spectra_file_names, fmt='%s')
np.savetxt(FLAGGED_GRAPHS_FILE, flagged_graphs, fmt='%s')
np.savetxt(FLAGGED_SNR_GRAPHS_FILE, flagged_snr_in_ehvo_graphs, fmt='%s')

## what are we wanting printed to the file?
#else:
#    print(current_spectrum_file_name) ## later will print this to a file



print("--- %s seconds" %(time.time()-start_time))
   
