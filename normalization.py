# Normalization of the quasar spectra
# Please Check README file before changes anything!!!!
import os
import sys
import numpy as np 
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
from utility_functions import print_to_file, clear_file, append_row_to_csv
from data_types import Range, ColumnIndexes, PointData, RangesData, FigureData, FigureDataOriginal, DataNormalized

CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else "sorted_norm.csv"

SPEC_DIREC = os.getcwd() + "/files/" # Set location of input and output spectrum files

LOG_FILE = "log.txt"
FINAL_INIT_PARAMS_FILE = SPEC_DIREC + "/" + "final_initial_parameters.txt"
PROCESSED_SPECTRA_FILE = SPEC_DIREC + "/" + "processed_spectra_filenames.txt"
FLAGGED_GRAPHS_FILE = SPEC_DIREC + "/" + "flagged_for_absorption_or_bad_normalization.txt"
FLAGGED_SNR_GRAPHS_FILE = SPEC_DIREC + "/" + "flagged_snr_in_ehvo_graphs.txt"
GOODNESS_OF_FIT_FILE = SPEC_DIREC + "/" + "chi_sq_values_all.csv"
BAD_NORMALIZATION_FLAGGED_FILE = SPEC_DIREC + "/" + "bad_normalization.csv"
GOOD_NORMALIZATION_FLAGGED_FILE = SPEC_DIREC + "/" + "good_normalization.csv"
NORM_FILE_EXTENSION = "norm.dr9"

ORIGINAL_PDF = PdfPages('original_graphs.pdf') # create pdf
NORMALIZED_PDF = PdfPages('normalized_graphs.pdf') # create pdf

STARTS_FROM, ENDS_AT = 1, 9

WAVELENGTH_RESTFRAME = Range(1200., 1800.)
WAVELENGTH_FOR_SNR = Range(1250., 1400.)
WAVELENGTH_RESTFRAME_FOR_LEFT_POINT = Range(1280., 1290.)
WAVELENGTH_RESTFRAME_FOR_MIDDLE_POINT = Range(1420., 1430.)
WAVELENGTH_RESTFRAME_FOR_RIGHT_POINT = Range(1690., 1710.)
WAVELENGTH_RESTFRAME_TEST_1 = Range(1315., 1325.)
WAVELENGTH_RESTFRAME_TEST_2 = Range(1350., 1360.)

column_index = ColumnIndexes(0, 1, 2)

b = 1250 # powerlaw
c = -0.5 # powerlaw

def powerlaw(wavelength, b, c) -> float:
    return b * (np.power(wavelength, c))

def wavelength_flux_error_for_points(starting_point: float, ending_point: float, z: float, spectra_data) -> RangesData:
    wavelength_column = spectra_data[:, column_index.wavelength]

    wavelength_observed_start = (z + 1) * starting_point
    wavelength_observed_end = (z + 1) * ending_point

    point_from = np.max(np.where(wavelength_column < wavelength_observed_start))
    point_to = np.min(np.where(wavelength_column > wavelength_observed_end))

    wavelength = spectra_data[point_from:point_to, column_index.wavelength]
    flux = spectra_data[point_from:point_to, column_index.flux] 
    error = spectra_data[point_from:point_to, column_index.error] 
  
    return RangesData(wavelength, flux, error)

def define_three_anchor_points(z: float, spectra_data):
    left_point_ranges = wavelength_flux_error_for_points(
        WAVELENGTH_RESTFRAME_FOR_LEFT_POINT.start,
        WAVELENGTH_RESTFRAME_FOR_LEFT_POINT.end,
        z,
        spectra_data)
    left_point = PointData(
        np.median(left_point_ranges.wavelength),
        np.median(left_point_ranges.flux),
        np.median(left_point_ranges.error))

    print(left_point_ranges.flux)
    print_to_file(left_point_ranges.flux, LOG_FILE)

    middle_point_ranges = wavelength_flux_error_for_points(
        WAVELENGTH_RESTFRAME_FOR_MIDDLE_POINT.start,
        WAVELENGTH_RESTFRAME_FOR_MIDDLE_POINT.end,
        z,
        spectra_data)
    middle_point = PointData(
        np.median(middle_point_ranges.wavelength),
        np.median(middle_point_ranges.flux),
        np.median(middle_point_ranges.error))

    right_point_ranges = wavelength_flux_error_for_points(
        WAVELENGTH_RESTFRAME_FOR_RIGHT_POINT.start,
        WAVELENGTH_RESTFRAME_FOR_RIGHT_POINT.end,
        z,
        spectra_data)
    right_point = PointData(
        np.median(right_point_ranges.wavelength),
        np.median(right_point_ranges.flux),
        np.median(right_point_ranges.error))

    return (left_point, middle_point, right_point)

def wavelength_flux_error_in_range(starting_point: float, ending_point: float, z: float, spectra_data) -> RangesData:
    wavelength_column = spectra_data[:, column_index.wavelength]

    wavelength_observed_from = (z + 1) * starting_point
    wavelength_observed_to = (z + 1) * ending_point

    wavelength_lower_limit = np.where(wavelength_column > wavelength_observed_from)
    wavelength_upper_limit = np.where(wavelength_column < wavelength_observed_to)
    
    wavelength = spectra_data[np.min(wavelength_lower_limit[column_index.wavelength]):np.max(wavelength_upper_limit[column_index.wavelength]), column_index.wavelength]
    flux = spectra_data[np.min(wavelength_lower_limit[column_index.wavelength]): np.max(wavelength_upper_limit[column_index.wavelength]), column_index.flux]
    error = spectra_data[np.min(wavelength_lower_limit[column_index.wavelength]): np.max(wavelength_upper_limit[column_index.wavelength]), column_index.error]
    
    return RangesData(wavelength, flux, error)

def draw_original_figure(figure_index: int, original_ranges: RangesData, data: FigureDataOriginal, test1: RangesData, test2: RangesData):
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
    ORIGINAL_PDF.savefig()
    plt.close(figure_index)

def draw_normalized_figure(figure_index: int, original_ranges: RangesData, figure_data: FigureData, normalized_data: DataNormalized,
                            test1: RangesData, test2: RangesData, normalized_flux_test_1, normalized_flux_test_2):
    main_color = "xkcd:ultramarine"
    test_1_color, test_2_color = "xkcd:green apple", "xkcd:bubblegum"
    subtitle_text = f"z={figure_data.z} snr={figure_data.snr} snr_mean_in_ehvo={figure_data.snr_mean_in_ehvo}"
    plt.figure(figure_index) 
    plt.text(((figure_data.wavelength_from + figure_data.wavelength_to)/2.3), np.max(normalized_data.flux_normalized)/1.07, figure_data.spectrum_file_name)
    plt.text(((figure_data.wavelength_from + figure_data.wavelength_to)/2.3), np.max(normalized_data.flux_normalized), subtitle_text)
    plt.title(figure_data.spectrum_file_name)
    plt.plot(original_ranges.wavelength, normalized_data.flux_normalized, color = main_color, linestyle = "-")
    plt.plot(original_ranges.wavelength, normalized_data.error_normalized, color = "black", linestyle = "-")
    plt.title("Normalized Data vs. Normalized Error")
    plt.xlabel("Wavelength [A]")
    plt.ylabel("Normalized Flux[10^[-17]]cgs")
    plt.plot(test1.wavelength, normalized_flux_test_1, color = test_1_color, linestyle = "-")
    plt.plot(test2.wavelength, normalized_flux_test_2, color = test_2_color, linestyle = "-")
    plt.plot((original_ranges.wavelength[0], original_ranges.wavelength[-1]), (1, 1), color = "red", linestyle = "-")
    NORMALIZED_PDF.savefig()
    plt.close(figure_index)

def process_spectra_and_draw_figures(index: int, z, snr, spectrum_file_name):

    print(str(index) + ": " + spectrum_file_name)
    print_to_file(str(index) + ": " + spectrum_file_name, LOG_FILE)

    current_spectra_data = np.loadtxt(SPEC_DIREC + spectrum_file_name)

    wavelength_observed_from = (z + 1) * WAVELENGTH_RESTFRAME.start
    wavelength_observed_to = (z + 1) * WAVELENGTH_RESTFRAME.end

    point_C, point_B, point_A = define_three_anchor_points(z, current_spectra_data)

    # THE THREE POINTS THAT THE POWER LAW WILL USE (Points C, B, and A)
    power_law_data_x = (point_C.wavelength, point_B.wavelength, point_A.wavelength)
    power_law_data_y = (point_C.flux, point_B.flux, point_A.flux)

    # DEFINING WAVELENGTH, FLUX, AND ERROR (CHOOSING THEIR RANGE)
    wavelength, flux, error = wavelength_flux_error_in_range(WAVELENGTH_RESTFRAME.start, WAVELENGTH_RESTFRAME.end, z, current_spectra_data)
    original_ranges = RangesData(wavelength, flux, error)

    print(power_law_data_x)
    print_to_file(power_law_data_x, LOG_FILE)
    print(power_law_data_y)
    print_to_file(power_law_data_y, LOG_FILE)

    # CURVE FIT FOR FIRST POWERLAW
    try:
        pars, covar = curve_fit(powerlaw, power_law_data_x, power_law_data_y, p0=[b, c], maxfev=10000)
    except:
        print("Error - curve_fit failed-1st powerlaw " + spectrum_file_name)
        print_to_file("Error - curve_fit failed-1st powerlaw " + spectrum_file_name, LOG_FILE)

    bf, cf = pars[0], pars[1]

    flux_normalized = flux/powerlaw(wavelength, bf, cf)
    error_normalized = error/powerlaw(wavelength, bf, cf)
    normalized_data = DataNormalized(flux_normalized, error_normalized)

    for n in range(1, len(flux_normalized) - 5):          
        if abs(flux_normalized[n + 1] - flux_normalized[n]) > 0.5:
            if error_normalized[n + 1] > 0.25:
                error_normalized[n + 1] = error_normalized[n]
                flux_normalized[n + 1] = flux_normalized[n]
                error[n + 1] = error[n]
                flux[n + 1] = flux[n]
        if error_normalized[n] > 0.5:
            error_normalized[n] = error_normalized[n-1]
            flux_normalized[n] = flux_normalized[n - 1]
            error[n] = error[n - 1]
            flux[n] = flux[n - 1]
        if abs(flux_normalized[n + 1] - flux_normalized[n]) > 5:
            error_normalized[n + 1] = error_normalized[n]
            flux_normalized[n + 1] = flux_normalized[n]
            error[n + 1] = error[n]
            flux[n + 1] = flux[n]

    ############# TESTING TWO REGIONS ##########################
    flagged = False
    test1 = wavelength_flux_error_in_range(WAVELENGTH_RESTFRAME_TEST_1.start, WAVELENGTH_RESTFRAME_TEST_1.end, z, current_spectra_data)
    normalized_flux_test_1 = test1.flux/powerlaw(test1.wavelength, bf, cf)
    flagged_by_test1 = abs(np.median(normalized_flux_test_1) - 1) >= 0.05
    if flagged_by_test1:
        print("flagged_by_test1: ", flagged_by_test1)
        print_to_file("flagged_by_test1: " + str(flagged_by_test1), LOG_FILE)

    test2 = wavelength_flux_error_in_range(WAVELENGTH_RESTFRAME_TEST_2.start, WAVELENGTH_RESTFRAME_TEST_2.end, z, current_spectra_data)
    normalized_flux_test_2 = test2.flux/powerlaw(test2.wavelength, bf, cf)
    flagged_by_test2 = abs(np.median(normalized_flux_test_2) - 1) >= 0.05
    if flagged_by_test2:
        print("flagged_by_test2: ", flagged_by_test2)
        print_to_file("flagged_by_test2: " + str(flagged_by_test2), LOG_FILE)

    if flagged_by_test1 and flagged_by_test2:
        flagged = True
        error_message = "Flagging figure #" + str(index) + ", file name: " + spectrum_file_name
        print(error_message)
        print_to_file(error_message, LOG_FILE)

    residuals_test1 = test1.flux - powerlaw(test1.wavelength, bf, cf)
    residuals_test2 = test2.flux - powerlaw(test2.wavelength, bf, cf)    
    residuals_test1_and_2 = np.concatenate([residuals_test1,residuals_test2])
    wavelength_tests_1_and_2 = np.concatenate([test1.wavelength, test2.wavelength])
    chi_sq = sum((residuals_test1_and_2**2)/powerlaw(wavelength_tests_1_and_2, bf, cf))
      
    fields=[index - STARTS_FROM + 1, index, chi_sq]
    append_row_to_csv(GOODNESS_OF_FIT_FILE, fields)
    if chi_sq > 8 and flagged_by_test1 and flagged_by_test2:
        append_row_to_csv(BAD_NORMALIZATION_FLAGGED_FILE, fields)
    else:
        append_row_to_csv(GOOD_NORMALIZATION_FLAGGED_FILE, fields)
        
    ##########################################################

    ############# SNR Calculations ##########################
    wavelengths_for_snr_lower = np.where (wavelength/(z + 1.) < WAVELENGTH_FOR_SNR.start)
    wavelengths_for_snr_upper = np.where (wavelength/(z + 1.) > WAVELENGTH_FOR_SNR.end)
    snr_mean_in_ehvo = round(np.mean(1./error_normalized[np.max(wavelengths_for_snr_lower[0]):np.min(wavelengths_for_snr_upper)]), 5)
    ##########################################################

    flagged_snr_mean_in_ehvo = False
    if snr_mean_in_ehvo < 10.:
        flagged_snr_mean_in_ehvo = True

    figure_data = FigureData(spectrum_file_name, wavelength_observed_from, wavelength_observed_to, z, snr, snr_mean_in_ehvo)
    original_figure_data = FigureDataOriginal(figure_data, bf, cf, power_law_data_x, power_law_data_y)
    
    draw_original_figure(index, original_ranges, original_figure_data, test1, test2)
    draw_normalized_figure(index, original_ranges, figure_data, normalized_data, test1, test2, normalized_flux_test_1, normalized_flux_test_2)

    norm_w_f_e = (wavelength, flux_normalized, error_normalized)
    norm_w_f_e = (np.transpose(norm_w_f_e))  
    np.savetxt(SPEC_DIREC + spectrum_file_name[0:20] + NORM_FILE_EXTENSION, norm_w_f_e)
    return bf, cf, flagged, flagged_snr_mean_in_ehvo, snr_mean_in_ehvo


def main(starting_index: int, ending_index: int):
    spectra_list, redshift_value_list, snr_value_list = [], [], []

    # Reading the file and assigning to the specific lists
    with open(CONFIG_FILE) as f:
        for line in f:
            each_row_in_file = line.split(",")
            spectra_list.append(each_row_in_file[0])
            redshift_value_list.append(np.float(each_row_in_file[1]))
            snr_value_list.append(np.float(each_row_in_file[2]))

    indices, spectra_indices, processed_spectra_file_names, powerlaw_final_b_values, powerlaw_final_c_values = [], [], [], [], []
    flagged_indices, flagged_spectra_indices, flagged_spectra_file_names = [], [], []
    flagged_snr_indices, flagged_snr_spectra_indices, flagged_snr_spectra_file_names, flagged_snr_in_ehvo_values = [], [], [], []

    for spectra_index in range(starting_index, ending_index + 1):
        z = round(redshift_value_list[spectra_index - 1], 5)
        snr = round(snr_value_list[spectra_index - 1], 5)
        current_spectrum_file_name = spectra_list[spectra_index - 1]
        b_final, c_final, failed_test, flagged_snr_mean_in_ehvo, snr_mean_in_ehvo = process_spectra_and_draw_figures(spectra_index, z, snr, current_spectrum_file_name)
       
        # add condition here?
        powerlaw_final_b_values.append(b_final)
        powerlaw_final_c_values.append(c_final)
        processed_spectra_file_names.append(current_spectrum_file_name)
        indices.append(spectra_index - starting_index + 1)
        spectra_indices.append(spectra_index)
        if failed_test:
            flagged_spectra_file_names.append(current_spectrum_file_name)
            flagged_indices.append(spectra_index - starting_index + 1)
            flagged_spectra_indices.append(spectra_index)
        
        if flagged_snr_mean_in_ehvo:
            flagged_snr_spectra_file_names.append(current_spectrum_file_name)
            flagged_snr_indices.append(spectra_index - starting_index + 1)
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

    np.savetxt(FINAL_INIT_PARAMS_FILE, final_initial_parameters, fmt="%s")
    np.savetxt(PROCESSED_SPECTRA_FILE, processed_spectra_file_names, fmt='%s')
    np.savetxt(FLAGGED_GRAPHS_FILE, flagged_graphs, fmt='%s')
    np.savetxt(FLAGGED_SNR_GRAPHS_FILE, flagged_snr_in_ehvo_graphs, fmt='%s')

if __name__ == "__main__":
    clear_file(LOG_FILE)
    clear_file(GOODNESS_OF_FIT_FILE)
    clear_file(BAD_NORMALIZATION_FLAGGED_FILE)
    clear_file(GOOD_NORMALIZATION_FLAGGED_FILE)
    
    fields=["index", "spectra index", "chi_sq"]
    append_row_to_csv(GOODNESS_OF_FIT_FILE, fields)
    append_row_to_csv(BAD_NORMALIZATION_FLAGGED_FILE, fields)
    append_row_to_csv(GOOD_NORMALIZATION_FLAGGED_FILE, fields)

    main(STARTS_FROM, ENDS_AT)