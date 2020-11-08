# Normalization of the quasar spectra
# Please Check README file before changes anything!!!!
import os
import sys
import numpy as np 
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
import utility_functions


b = 1250 # powerlaw
c = -0.5 # powerlaw
STARTS_FROM, ENDS_AT = 0, 9
WAVELENGTH_RANGE_RESTFRAME = (1200., 1800.)
WAVELENGTH_RANGE_FOR_SNR = (1250., 1400.)
WAVELENGTH_RESTFRAME_RANGE_POINT_A = (1690., 1710.)
WAVELENGTH_RESTFRAME_RANGE_POINT_B = (1420., 1430.)
WAVELENGTH_RESTFRAME_RANGE_POINT_C = (1280., 1290.)

config_file = sys.argv[1] if len(sys.argv) > 1 else "sorted_norm.csv"

# Set location of spectrum files and create pdfs
specdirec = os.getcwd() + "/files/"
log_file = "log.txt"
original_pdf = PdfPages('original_all_graph_Sean.pdf') 
normalized_pdf = PdfPages('normalized_all_graph_Sean.pdf')

def powerlaw(wavelength, b, c) -> float:
    return b * (np.power(wavelength, c))

def wavelength_flux_error_for_points(starting_point: float, ending_point: float, z: float, spectra_data):
    wavelength_column = spectra_data[:, 0]

    wavelength_observed_start = (z + 1) * starting_point
    wavelength_observed_end = (z + 1) * ending_point

    point_from = np.max(np.where(wavelength_column < wavelength_observed_start))
    point_to = np.min(np.where(wavelength_column > wavelength_observed_end))

    wavelength = spectra_data[point_from:point_to, 0]
    flux = spectra_data[point_from:point_to, 1] 
    error = spectra_data[point_from:point_to, 2] 
  
    return wavelength, flux, error

def wavelength_flux_error_in_range(starting_point: float, ending_point: float, z: float, spectra_data):
    wavelength_column = spectra_data[:, 0]

    wavelength_observed_from = (z + 1) * starting_point
    wavelength_observed_to = (z + 1) * ending_point

    wavelength_lower_limit = np.where(wavelength_column > wavelength_observed_from)
    wavelength_upper_limit = np.where(wavelength_column < wavelength_observed_to)
    
    wavelength = spectra_data[np.min(wavelength_lower_limit[0]):np.max(wavelength_upper_limit[0]), 0]
    flux = spectra_data[np.min(wavelength_lower_limit[0]): np.max(wavelength_upper_limit[0]), 1]
    error = spectra_data[np.min(wavelength_lower_limit[0]): np.max(wavelength_upper_limit[0]), 2]
    
    return wavelength, flux, error

def process_sprectra_and_draw_figure(index: int, z, snr, spectrum_file_name):

    print(str(index + 1) + ": " + spectrum_file_name)
    utility_functions.print_to_file(str(index + 1) + ": " + spectrum_file_name, log_file)

    current_spectra_data = np.loadtxt(specdirec + spectrum_file_name)

    wavelength_observed_from = (z + 1) * WAVELENGTH_RANGE_RESTFRAME[0]
    wavelength_observed_to = (z + 1) * WAVELENGTH_RANGE_RESTFRAME[1]

    # POINT C (LEFTMOST POINT)
    wavelength_point_C, flux_point_C, error_point_C = wavelength_flux_error_for_points(WAVELENGTH_RESTFRAME_RANGE_POINT_C[0], WAVELENGTH_RESTFRAME_RANGE_POINT_C[1], z, current_spectra_data)
    median_wavelength_point_C, median_flux_point_C = np.median(wavelength_point_C), np.median(flux_point_C)

    print(flux_point_C)
    utility_functions.print_to_file(flux_point_C, log_file)
        
    # POINT B (MIDDLE POINT)
    wavelength_point_B, flux_point_B, error_point_B = wavelength_flux_error_for_points(WAVELENGTH_RESTFRAME_RANGE_POINT_B[0], WAVELENGTH_RESTFRAME_RANGE_POINT_B[1], z, current_spectra_data)
    median_wavelength_point_B, median_flux_point_B, median_error_point_B = np.median(wavelength_point_B), np.median(flux_point_B), np.median(error_point_B)
    st_dev_of_flux = np.std(flux_point_B)
        
    # POINT A (RIGHTMOST POINT)
    wavelength_point_A, flux_point_A, error_point_A = wavelength_flux_error_for_points(WAVELENGTH_RESTFRAME_RANGE_POINT_A[0], WAVELENGTH_RESTFRAME_RANGE_POINT_A[1], z, current_spectra_data)
    median_wavelength_point_A, median_flux_point_A = np.median(wavelength_point_A), np.median(flux_point_A)

    # THE THREE POINTS THAT THE POWER LAW WILL USE (Points C, B, and A)
    power_law_datax = (median_wavelength_point_C, median_wavelength_point_B, median_wavelength_point_A)
    power_law_datay = (median_flux_point_C, median_flux_point_B, median_flux_point_A)

    # DEFINING WAVELENGTH, FLUX, AND ERROR (CHOOSING THEIR RANGE)
    wavelength, flux, error = wavelength_flux_error_in_range(WAVELENGTH_RANGE_RESTFRAME[0], WAVELENGTH_RANGE_RESTFRAME[1], z, current_spectra_data)

    print(power_law_datax)
    utility_functions.print_to_file(power_law_datax, log_file)
    print(power_law_datay)
    utility_functions.print_to_file(power_law_datay, log_file)

    # CURVE FIT FOR FIRST POWERLAW
    try:
        pars, covar = curve_fit(powerlaw, power_law_datax, power_law_datay, p0=[b, c], maxfev=10000)
    except:
        print("Error - curve_fit failed-1st powerlaw " + spectrum_file_name)
        utility_functions.print_to_file("Error - curve_fit failed-1st powerlaw " + spectrum_file_name, log_file)

    bf, cf = pars[0], pars[1]

    flux_normalized = flux/powerlaw(wavelength, bf, cf)
    error_normalized = error/powerlaw(wavelength, bf, cf)

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
    wavelength_test_1, flux_test_1, error_test_1 = wavelength_flux_error_in_range(1350., 1360., z, current_spectra_data)
    normalized_flux_test_1 = flux_test_1/powerlaw(wavelength_test_1, bf, cf)
    failed_test_1 = abs(np.median(normalized_flux_test_1) - 1) >= 0.05
    if failed_test_1:
        print("failed_test_1: ", failed_test_1)
        utility_functions.print_to_file("failed_test_1: " + str(failed_test_1), log_file)

    wavelength_test_2, flux_test_2, error_test_2 = wavelength_flux_error_in_range(1315., 1325., z, current_spectra_data)
    normalized_flux_test_2 = flux_test_2/powerlaw(wavelength_test_2, bf, cf)
    failed_test_2 = abs(np.median(normalized_flux_test_2) - 1) >= 0.05
    if failed_test_2:
        print("failed_test_2: ", failed_test_2)
        utility_functions.print_to_file("failed_test_2: " + str(failed_test_2), log_file)

    if failed_test_1 and failed_test_2:
        error_message = "Flagging figure #" + str(index + 1) + ", file name: " + spectrum_file_name
        print(error_message)
        utility_functions.print_to_file(error_message, log_file)
    ##########################################################

    # SNR Calculations:
    wavelengths_for_snr_lower = np.where (wavelength/(z + 1.) < WAVELENGTH_RANGE_FOR_SNR[0])
    wavelengths_for_snr_upper = np.where (wavelength/(z + 1.) > WAVELENGTH_RANGE_FOR_SNR[1])
    snr_mean_in_ehvo = round(np.mean(1./error_normalized[np.max(wavelengths_for_snr_lower[0]):np.min(wavelengths_for_snr_upper)]), 5)

    # Start of Figure 1
    plt.figure(index + 1)

    # PLOT FIGURE
    plt.plot(wavelength, powerlaw(wavelength, bf, cf), color = "red", linestyle = "--")
    plt.plot(median_wavelength_point_B, median_flux_point_B - median_error_point_B, 'yo')
    plt.plot(median_wavelength_point_B, median_flux_point_B - 3 * (median_error_point_B), 'yo')
    plt.plot(median_wavelength_point_B, median_flux_point_B - st_dev_of_flux, color = "green", marker = "o")
    plt.title(spectrum_file_name)
    plt.xlabel("Wavelength[A]")
    plt.ylabel("Flux[10^[-17]]cgs")
    plt.text(((wavelength_observed_from + wavelength_observed_to)/2.17), np.max(flux), f"z= {z} snr={snr} snr_1325={snr_mean_in_ehvo}")
    plt.plot(median_wavelength_point_B, median_flux_point_B, 'yo')
    plt.plot(wavelength, flux, color = "blue", linestyle = "-")
    plt.plot(power_law_datax, power_law_datay, 'ro')
    plt.plot(wavelength, error, color = "black", linestyle = "-")
    plt.plot(wavelength_test_1, flux_test_1, color = "magenta", linestyle = "-")
    plt.plot(wavelength_test_2, flux_test_2, color = "yellow", linestyle = "-")

    original_pdf.savefig()
    plt.close(index + 1)
    # End of Figure 1
        
        
    # Start of Figure 2
    plt.figure(index + 1)
    0.2, "z=" + str(z) + " snr=" + str(snr)

    plt.text(wavelength_observed_from + 1000, np.max(flux_normalized) - 0.2, spectrum_file_name)
    plt.title(spectrum_file_name)
    plt.plot(wavelength, flux_normalized, color = "blue", linestyle = "-")
    plt.plot((wavelength[0], wavelength[-1]), (1, 1), color = "red", linestyle = "-")
    plt.plot(wavelength, error_normalized, color = "black", linestyle = "-")
    plt.title("normalized data vs. normalized error")
    plt.xlabel("Wavelength [A]")
    plt.ylabel("Normalized Flux[10^[-17]]cgs")
    plt.plot(wavelength_test_1, normalized_flux_test_1, color = "magenta", linestyle = "-")
    plt.plot(wavelength_test_2, normalized_flux_test_2, color = "yellow", linestyle = "-")
        
    normalized_pdf.savefig()
    plt.close(index + 1)

    # End of Figure 2

    www = (wavelength, flux_normalized, error_normalized)
    www = (np.transpose(www))
    oo=spectrum_file_name[0:20]
        
    np.savetxt(specdirec + oo +'norm.dr9', www)  # ,fmt='%s')
    return bf, cf


def normalize_spectra(starting_index: int, ending_index: int):
    utility_functions.clear_file(log_file)

    spectra_list, redshift_value_list, snr_value_list = [], [], []

    # Reading the file and assigning to the specific lists
    with open(config_file) as f:
        for line in f:
            each_row_in_file = line.split(",")
            spectra_list.append(each_row_in_file[0])
            redshift_value_list.append(np.float(each_row_in_file[1]))
            snr_value_list.append(np.float(each_row_in_file[2]))

    processed_spectra_file_names, powerlaw_final_b_values, powerlaw_final_c_values = [], [], []

    for spectra_index in range(starting_index, ending_index):
        z = round(redshift_value_list[spectra_index], 5)
        snr = round(snr_value_list[spectra_index], 5)
        current_spectrum_file_name = spectra_list[spectra_index]
        b_final, c_final = process_sprectra_and_draw_figure(spectra_index, z, snr, current_spectrum_file_name)
       
        # add condition here?
        powerlaw_final_b_values.append(b_final)
        powerlaw_final_c_values.append(c_final)
        processed_spectra_file_names.append(current_spectrum_file_name)

    final_initial_parameters = [processed_spectra_file_names, powerlaw_final_b_values, powerlaw_final_c_values]
    final_initial_parameters = (np.transpose(final_initial_parameters))
        
    original_pdf.close()
    normalized_pdf.close()

    np.savetxt(specdirec + "/Final_Initial_Parameters.txt", final_initial_parameters, fmt="%s")
    np.savetxt(specdirec + "/good_spectra.txt", processed_spectra_file_names, fmt='%s')

if __name__ == "__main__":
    normalize_spectra(STARTS_FROM, ENDS_AT)