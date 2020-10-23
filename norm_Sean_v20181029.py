#==============Normalization of the quasar spectra=================
# THIS FILE:
# Imports the neccessary libraries
# Specify the range of spectra list(if needed)
# Reads data from .csv files and assign to the three different lists
# Uses the powerlaw and calculates the normalization of spectra with given SDSS data
# Plots the necessary values to the graph (as figures)
# Saves these figures to the pdf spectrum files. (Given location)
# At the end saves


# ====================================================================



#============Import Files and Libraries========================
import os
import numpy as np 
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
import utility_functions
#================================================================

def normalize_spectra():
    #===========Specifying The Range==============================
    # this is a range that how many spectra you want to check
    # and/or between spectras in the file
    starts_from = 0
    ends_at = 9
    good_categorized_spectra_range = range(starts_from, ends_at)

    #For IDE run, locate the file
    config_file = "sorted_norm.csv" 

    #For command line, activate this config_file
    #config_file = sys.argv[1] 
    #================================================================

    #============Set location of spectrum files========================
    specdirec = os.getcwd() + "/files/"
    original_pdf = PdfPages('original_all_graph_Sean.pdf') 
    normalized_pdf = PdfPages('normalized_all_graph_Sean.pdf')  
    #================================================================

    log_file = "log.txt"
    utility_functions.clear_file(log_file)

    # number_of_spectra = 6760
    number_of_spectra = utility_functions.file_length(config_file)


    spectra_list = list()
    redshift_value_list = list()
    snr_value_list = list()


    #========Reading the file and assignin to the specific lists============
    with open(config_file) as f:
        for line in f:
            each_row_in_file = line.split(",")
            spectra_list.append(each_row_in_file[0])
            redshift_value_list.append(np.float(each_row_in_file[1]))
            snr_value_list.append(np.float(each_row_in_file[2]))
    #============End of Reading and Assigning=============================== 


    b = 1250 #powerlaw
    c = -0.5 #powerlaw

    ee = []
    processed_spectra_file_names = []
    eee = []

    init_pars = [b, c]

    powerlaw_not_made = []

    #==================POWERLAW FUNCTION==============================
    # Powerlaw function takes 3 parameters. Parameters b and c are defined above.
    # The parameter x is calculation of wavelength. The powerlaw calculates a value
    # as in formula and return a floating value. 
    def powerlaw(x, b, c):
        return b * (np.power(x, c))
    #==================END OF POWERLAW FUNCTION=========================
        

    ################################################################################
    ################################################################################
    ###########THIS HUGE FOR LOOP STARTS HERE AND FINISHES AROUND LINE 450s ########

    for index in good_categorized_spectra_range:
    
        current_spectrum_file_name = spectra_list[index]
        
        z = round(redshift_value_list[index], 5)
        snr = round(snr_value_list[index], 5)

        print(str(index + 1) + ": " + current_spectrum_file_name)
        utility_functions.print_to_file(str(index + 1) + ": " + current_spectrum_file_name, log_file)

        
        current_spectra_data = np.loadtxt(specdirec + current_spectrum_file_name)
        wavelength_column = current_spectra_data[:, 0]

        WAVELENGTH_RANGE_RESTFRAME = (1200., 1800.)
        wavelength_observed_from = (z + 1) * WAVELENGTH_RANGE_RESTFRAME[0]
        wavelength_observed_to = (z + 1) * WAVELENGTH_RANGE_RESTFRAME[1]

        # FIRST POINT (Point C)
        WAVELENGTH_RESTFRAME_RANGE_POINT_C = (1280., 1290.) # 1st point for powerlaw, Point C

        wavelength_observed_starting_point_C = (z + 1) * (WAVELENGTH_RESTFRAME_RANGE_POINT_C[0])
        wavelength_observed_ending_point_C = (z + 1) * (WAVELENGTH_RESTFRAME_RANGE_POINT_C[1])
          
        p7 = np.max(np.where(wavelength_column < wavelength_observed_starting_point_C))  
        p9 = np.min(np.where(wavelength_column > wavelength_observed_ending_point_C))
    
        median_flux_point_C = np.median(current_spectra_data[p7:p9, 1])
        median_wavelength_point_C = np.median(current_spectra_data[p7:p9, 0])

        print(current_spectra_data[p7:p9, 1])
        utility_functions.print_to_file(current_spectra_data[p7:p9, 1], log_file)

        ######################### MIDDLE POINT (Between Points A and B)
        wavelength_new_emit2 = 1690.  # Point B (use only this point later(instead of both A and B))
        wavelength_new_emit3 = 1725.669  # Point A
        
        wavelength_new_midpoint = (wavelength_new_emit2 + wavelength_new_emit3)/2

        wavelength_new_obs2 = (z + 1) * wavelength_new_emit2
        wavelength_new_obs3 = (z + 1) * wavelength_new_emit3
        wavelength_new_midpoint_obs = (z + 1) * wavelength_new_midpoint

        p1 = np.max(np.where(wavelength_column < wavelength_new_obs2))
        p2 = np.min(np.where(wavelength_column > wavelength_new_midpoint_obs))
    
        median_flux1 = np.median(current_spectra_data[p1:p2, 1])
        median_wavelength1 = np.median(current_spectra_data[p1:p2, 0])

        # new code
        # WAVELENGTH_RESTFRAME_RANGE_POINT_A = (1690., 1710.)

        # wavelength_observed_starting_point_A = (z + 1) * (WAVELENGTH_RESTFRAME_RANGE_POINT_A[0])
        # wavelength_observed_ending_point_A = (z + 1) * (WAVELENGTH_RESTFRAME_RANGE_POINT_A[1])

        # p1 = np.max(np.where(wavelength_column < wavelength_observed_starting_point_A))  
        # p2 = np.min(np.where(wavelength_column > wavelength_observed_ending_point_A))
    
        # median_flux_point_A = np.median(current_spectra_data[p1:p2, 1])
        # median_wavelength_point_A = np.median(current_spectra_data[p1:p2, 0])

        ########################### END OF MIDDLE POINT ##############################
        
        
        ########################### LAST POINT################################### is Point A
        
        p3 = np.max(np.where(wavelength_column < wavelength_new_midpoint_obs))
        p4 = np.min(np.where(wavelength_column > wavelength_new_obs3))

        flux2 = current_spectra_data[p3:p4, 1]

        median_flux2 = np.median(flux2) 
        wavelength_flux2 = current_spectra_data[p3:p4, 0] 

        median_wavelength2 = np.median(wavelength_flux2)
        ########################### END OF LAST POINT ##############################
        
        
    
        ######################## D POINT AND THREE POINTS #######################
        # range taken in rest frame: 1415-1430
        dpoint_starting_point_restframe = 1415 # Point D
        dpoint_ending_point_restframe = 1430 # Point D


        dpoint_starting_point = (z + 1) * (dpoint_starting_point_restframe)
        dpoint_ending_point = (z + 1) * (dpoint_ending_point_restframe)


        pp7 = np.max(np.where(wavelength_column < dpoint_starting_point))
        pp9 = np.min(np.where(wavelength_column > dpoint_ending_point))

        flux33 = current_spectra_data[pp7:pp9, 1]  
    
        median_flux33 = np.median(flux33) 
        median_wavelength33 = np.median(current_spectra_data[pp7:pp9, 0])

        # AVERAGE ERROR FOR D POINT
        median_flux_error33 = np.median(current_spectra_data[pp7:pp9, 2])  
        st_dev_of_flux = np.std(flux33[:]) #Standart deviation of flux

        # THE THREE POINTS (THE THREE POINTS THAT THE original power law WILL USE), Points C, Point A, Point B
        power_law_datax = (median_wavelength_point_C, median_wavelength1, median_wavelength2)
        power_law_datay = (median_flux_point_C, median_flux1, median_flux2)

        # THE THREE POINTS (THE THREE POINTS THAT THE second power law WILL USE), Points D, Point A, Point B
        power_law_datax2 = (median_wavelength33, median_wavelength1, median_wavelength2)
        power_law_datay2 = (median_flux33, median_flux1, median_flux2)

        ############# END OF D POINT AND THREE POINTS #################################
        

        # BASICALLY, DEFINING MY WAVELENGTH, FLUX, AND ERROR (OR CHOOSING THEIR RANGE)
        wavelength_lower_limit = np.where(wavelength_column > wavelength_observed_from)
        wavelength_upper_limit = np.where(wavelength_column < wavelength_observed_to)
        

        wavelength = current_spectra_data[np.min(wavelength_lower_limit[0])
                                : np.max(wavelength_upper_limit[0]), 0]
        
        flux = current_spectra_data[np.min(wavelength_lower_limit[0]): np.max(
            wavelength_upper_limit[0]), 1] 
        
        
        error = current_spectra_data[np.min(wavelength_lower_limit[0]): np.max(
            wavelength_upper_limit[0]), 2]

        fev = 10000

        print(power_law_datax)
        utility_functions.print_to_file(power_law_datax, log_file)
        print(power_law_datay)
        utility_functions.print_to_file(power_law_datay, log_file)

        # CURVE FIT FOR FIRST POWERLAW
        try:
            pars, covar = curve_fit(powerlaw, power_law_datax,
                                    power_law_datay, p0=init_pars, maxfev=fev)
        except:
            print("Error - curve_fit failed-1st powerlaw " + current_spectrum_file_name)
            utility_functions.print_to_file("Error - curve_fit failed-1st powerlaw " + current_spectrum_file_name, log_file)
            powerlaw_not_made.append(current_spectrum_file_name)

        normalizing = flux/powerlaw(wavelength, *pars)
        error_normalized = error/powerlaw(wavelength, *pars)
    
            
        for n in range(1, len(normalizing)-5):
            
            if abs(normalizing[n + 1] - normalizing[n]) > 0.5:

                if error_normalized[n+1] > 0.25:
                    error_normalized[n+1] = error_normalized[n]
                    normalizing[n+1] = normalizing[n]  # normalized graph
                    error[n+1] = error[n]  # original error
                    flux[n+1] = flux[n]  # original graph

            if error_normalized[n] > 0.5:
                error_normalized[n] = error_normalized[n-1]
                normalizing[n] = normalizing[n-1]  # normalized graph
                error[n] = error[n-1]  # original error
                flux[n] = flux[n-1]  # original graph


            if abs(normalizing[n + 1] - normalizing[n]) > 5:
                
                error_normalized[n+1] = error_normalized[n]
                normalizing[n+1] = normalizing[n]  # normalized graph
                error[n+1] = error[n]  # original error
                flux[n+1] = flux[n]  # original graph

   
        bf = pars[0]
        cf = pars[1]

        # REQUIREMENT FOR USING #POWERLAW.
        if (bf)*(np.power(median_wavelength33, cf)) > (median_flux33) - (3)*(median_flux_error33):

            ee.append(pars[0])
            eee.append(pars[1])
            processed_spectra_file_names.append(current_spectrum_file_name)

        # SNR Calculations:
        WAVELENGTH_RANGE_FOR_SNR = (1250., 1400.)
        wavelengths_for_snr_lower = np.where (wavelength/(z+1.) < WAVELENGTH_RANGE_FOR_SNR[0])
        wavelengths_for_snr_upper = np.where (wavelength/(z+1.) > WAVELENGTH_RANGE_FOR_SNR[1])
        snr_mean_in_ehvo = round(np.mean(1./error_normalized[np.max(wavelengths_for_snr_lower[0]):np.min(wavelengths_for_snr_upper)]), 5)

        # Start of Figure 1
        plt.figure(index + 1)

        # IF STATEMENT IS THERE TO AVOID ANY RANDOM, WEIRD PIXEL PROBLEM WITH THE ERRORS. FOR EXAMPLE: TO AVOID CASES WHERE THE ERROR IS 100
        
        if (bf)*(np.power(median_wavelength33, cf)) < (median_flux33) - (3)*(median_flux_error33):
            
            
            plt.plot(wavelength, powerlaw(wavelength, *pars), color = "red", linestyle = "--")
            plt.plot(median_wavelength33, median_flux33 - median_flux_error33, 'yo')
            plt.plot(median_wavelength33, median_flux33 -
                3*(median_flux_error33), 'yo')
            plt.plot(median_wavelength_point_C, median_flux_point_C, 'yo')
            plt.plot(median_wavelength33, median_flux33 - st_dev_of_flux, color = "green", marker = "o")
            plt.title(current_spectrum_file_name)
            plt.xlabel("Wavelength[A]")
            plt.ylabel("Flux[10^[-17]]cgs")
            plt.text(((wavelength_observed_from+wavelength_observed_to)/2.17), np.max(flux), f"z= {z} snr={snr} snr_1326={snr_mean_in_ehvo}" , style = 'italic')
            plt.plot(median_wavelength33, median_flux33, 'yo')
            plt.plot(wavelength, flux, color = "blue", linestyle = "-")
            plt.plot(power_law_datax2, power_law_datay2, 'ro')
            plt.plot(wavelength, error, color = "black", linestyle = "-")
            
        else:

            plt.plot(wavelength, powerlaw(wavelength, *pars), color = "red", linestyle = "--")
            plt.plot(median_wavelength33, median_flux33 - median_flux_error33, 'yo')
            plt.plot(median_wavelength33, median_flux33 -
                3*(median_flux_error33), 'yo')
            plt.plot(median_wavelength33, median_flux33 - st_dev_of_flux, color = "green", marker = "o")
            plt.title(current_spectrum_file_name)
            plt.xlabel("Wavelength[A]")
            plt.ylabel("Flux[10^[-17]]cgs")
            plt.text(((wavelength_observed_from+wavelength_observed_to)/2.17),np.max(flux), f"z= {z} snr={snr} snr_1325={snr_mean_in_ehvo}")
            #plt.text(wavelength_observed_from + 1000, np.max(flux)-10, current_spectrum_file_name)
            plt.plot(median_wavelength33, median_flux33, 'yo')
            plt.plot(wavelength, flux, color = "blue", linestyle = "-")
            plt.plot(power_law_datax, power_law_datay, 'ro')
            plt.plot(wavelength, error, color = "black", linestyle = "-")


        original_pdf.savefig()
        plt.close(index + 1)
        # End of Figure 1
        
        
        # Start of Figure 2
        plt.figure(index + 1)
        0.2, "z=" + str(z) + " snr=" + str(snr)

        plt.text(wavelength_observed_from + 1000, np.max(normalizing)-0.2, current_spectrum_file_name)
        plt.title(current_spectrum_file_name)
        
        plt.plot(wavelength, normalizing, color = "blue", linestyle = "-")#I CHANGED THIS JUST NOW
        plt.plot((wavelength[0], wavelength[-1]),(1, 1), color = "red", linestyle = "-")
        #plot((wavelength[0], wavelength[-1]),(1, 1))
        plt.plot(wavelength, error_normalized, color = "black", linestyle = "-") #I CHANGED THIS JUST NOW
        plt.title("normalized data vs. normalized error")
        plt.xlabel("Normalized Wavelength [A]")
        plt.ylabel("Flux[10^[-17]]cgs")
        normalized_pdf.savefig()
        plt.close(index + 1)

        # End of Figure 2
        
        www = (wavelength, normalizing, error_normalized)
        www = (np.transpose(www))
        oo=current_spectrum_file_name[0:20]
        
        #########################################################################################
        # This is where normalized files will be saved. Remember to change to your own directory!
        
        np.savetxt(specdirec + oo +'norm.dr9', www)  # ,fmt='%s')
        #########################################################################################


    ############### THE FOR LOOP STARTED AT LINE ~95 FINISHED HERE #####################
    #####################################################################################
    #####################################################################################

    # EVERYTHING BELOW IS NOT IN THE FOR LOOP

    final_initial_parameters = [processed_spectra_file_names, ee, eee]
    final_initial_parameters = (np.transpose(final_initial_parameters))

    # processed_spectra_file_names is the list of all good spectra.

    for l in powerlaw_not_made:
        for lj in processed_spectra_file_names:
            if lj in powerlaw_not_made:
                processed_spectra_file_names.remove(lj)
        
            
    original_pdf.close()
    normalized_pdf.close()


    np.savetxt(specdirec + "/Final_Initial_Parameters.txt", final_initial_parameters, fmt="%s")
    np.savetxt(specdirec + "/good_spectra.txt", processed_spectra_file_names, fmt='%s')
    np.savetxt(specdirec + "/Powerlaw1_did_not_work.txt", powerlaw_not_made, fmt='%s')


def main():
    normalize_spectra()

if __name__ == "__main__":
    main()