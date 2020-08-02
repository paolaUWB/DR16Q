#=====Please Check Readme file before start========
#=====Working on this code=======================
#================================================


#============Import Files and Libraries========================

#import sys #unused
#from pylab import *
#import scipy as sp  #unused
#import csv #unused
import numpy as np 
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

#================================================================


#=====Read, Assign, and  Print the file==========================

#This is how many specturim in the file that we need to process
specnum = 6760

# Set location of spectrum files.
specdirec = 'G:/School/_UWB Classes/ZCapstone/Capstone/Normalized-Spectre/files/'
pp1 = PdfPages('original_all_graph_Sean.pdf')  # Set output PDF file
pp2 = PdfPages('normalized_all_graph_Sean.pdf')  # Set normalized output PDF file

# specname_file = "confs/DR9Q_selection_specnames.lst" #|<-----The normal configs
# zem_file = "confs/DR9Q_selection_zem.lst"            #|

# specname_file = "confs/bad_specnames.lst" #|<------These are for "bad spectra"
# zem_file = "confs/bad_zem.lst"            #|

# Set cfg file path (csv with spec_name,z,snr...)
config_file = "sorted_norm.csv" 

#USE THIS LINE ONLY FOR WHEN COMMAND LINE ACTIVE
#config_file = sys.argv[1] 


#create a list for each column
spectra_action = list() #[row[0] for row in config]
redshifts_action = list() #[row[1] for row in config]
snr_action = list() #[row[2] for row in config]


#assign each column to a different list
for line in open(config_file, 'r'):
	d = line.split(",")
	spectra_action.append(d[0])
	redshifts_action.append(np.float(d[1]))
	snr_action.append(np.float(d[2]))
    


################################################################
##this is a range that how many spectra you want to check
##an/or between spectras in the file
starts_from = 0
ends_at = 5
good_categorized_spectra_action = range(starts_from, ends_at)
################################################################


# COUNTERS FOR THE TWO FIGURES, THEY ARE SET APART (0 AND 20) AS TO NOT OVERLAP BECAUSE OTHERWISE,
# THERE WILL BE TWO GRAPHS IN ONE FIGURE
# (ON TOP OF EACH OTHER)
count_fig1 = 0

# count_fig2 will always start counting from twice the number of total spectra
count_fig2 = 2*specnum


#Initilizations and empty lists
original_graph_number = 0
normalized_graph_number = 0
rows_of_power_law_txt_file = len(spectra_action)


## l = [] ##unused?
ee = []
ll = []
eee = []

# a=25
# b,c are initial parameters for first powerlaw, and bb,cc are intial parameteres for the second powerlaw.
b = 1250
c = -0.5

#==================POWERLAW FUNCTION==============================
#Description: Powerlaw calculation
#PreCondition: Takes 3 parameter as a value
#PostCOndition: Returns a value after calculation
def powerlaw(x, b, c):
    return b*(np.power(x, c))
#==================END OF POWERLAW FUNCTION=========================

#WE HAVE DEFINED THE INITIAL PARAMETERS AS B,C SO ANYTIME WE TELL THE CODE TO 
#PRINT INITIAL PARAMETERS, IT WILL DO SO IN THAT ORDER: B,C
init_pars = [b, c]  # FOR 1ST POWERLAW
init_pars2 = [b, c]  # FOR 2ND POWERLAW


powerlaw1_not_made = []
powerlaw2_not_made = []
i_all = []
counter = 0

# Get spectra names and redshifts

################################################################################
################################################################################
###########THIS HUGE FOR LOOP STARTS HERE AND FINISHES AROUND LINE 580s ########



for index in good_categorized_spectra_action:
    counter += 1 #increment
    # i= THE CURRENT SPECTRA THAT THE LOOP IS GOING THROUGH
    # z= THE CURRENT REDSHIFT THAT THE LOOP IS GOING THROUGH
    # append all the rows and columns in the file and round to 5 decimal point
   
    i = spectra_action[index] #get the spectra name
    i_all.append(i) #append all the values into the list
    z = round(redshifts_action[index], 5)
    snr = round(snr_action[index], 5)

    print(str(counter) + ": " + i)

    
    data = np.loadtxt(specdirec + i)  # Load in spectrum file
    number_rows = len(data[:, :])  # Get the number of points in the spectrum

    # Only calculating spectrum from 1200 - 1800 rest frame wavelength


    # THE POINTS USED FOR POWERLAW (START)
    # THIS IS THE RANGE OF THE ELECTROMAGNETIC SPECTRUM THAT WE WANT OUR SPECTRA TO BE.
    # FROM 1200 to 1800
    wavelength_emit1_initial = 1200
    wavelength_emit2_initial = 1800
    
    
    # Shift start wavelength into frame 
    wavelength_observe1 = (z+1)*wavelength_emit1_initial
    # Shift end wavelength into frame   
    wavelength_observe2 = (z+1)*wavelength_emit2_initial
    
    wavelength_NV_emit = 1242.8040
    # Shift Nitrogen V (NV) line into frame
    wavelength_NV_obs = (z+1)*wavelength_NV_emit

    # REST-FRAME WAVELENGTH RANGE FOR FIRST POINT
    wavelength_restframe_starting_point = 1280.206
    wavelength_restframe_ending_point = 1284.333
    
    # Shift start wavelength of first point into frame |<--This makes the range for the first
    wavelength_observed_starting_point = (z+1)*(wavelength_restframe_starting_point)
    # Shift end wavelength of first point into frame   |     point in the spectrum
    wavelength_observed_ending_point = (z+1)*(wavelength_restframe_ending_point)

    # a good range of rest frame wavelength is = 1686 - 1773
    # their midpoint is 1729, so im gonna take 1686-1729, average them up, and find the midpoint,
    # then find midpoint of 1729-1773
    wavelength_new_emit1 = 1282.398  # 1st point for powerlaw
	
    # Shift first power law point (a) into rest frame
    wavelength_new_obs1 = (z+1)*wavelength_new_emit1
	
    # FIRST POINT
    # Get all points from data with wavelengths less than our starting wavelength for our first point
    q6 = np.where(data[:, 0] < wavelength_observed_starting_point)

    p7 = 0
    try:  # Why is this needed?
        # Get the largest wavelength from q6. Is the index of the closest wavelength in the spectrum to the starting wavelength for our first point
        p7 = np.max(q6)
    except:
        pass

   # Get all points from data with wavelength greater than our ending wavelengt for our first point
    q8 = np.where(data[:, 0] > wavelength_observed_ending_point)
    
    # Get the smallest wavelength from q8. Is the index of the closest wavelength in the sprctrum to the ending wavelength for our first point
    p9 = np.min(q8)

    flux3 = data[p7:p9, 1]  # Get all flux values for our starting point
   
    # Get the median flux in the range.
    median_flux3 = np.median(flux3)
    
    # Get the corresponding wavelengths for this ragion
    wavelength_flux3 = data[p7:p9, 0]
    
    # Get the median wavelength of this region.
    median_wavelength3 = np.median(wavelength_flux3)
    
    # Make a point from the median wavelength and median flux (a)
    first_point = (median_wavelength3, median_flux3)

    print(flux3)

    wavelength_new_emit2 = 1677.938  # 2nd point for powerlaw
    wavelength_new_emit3 = 1725.669  # 3rd point for powerlaw
    
    # Average the points to get a midpoint
    wavelength_new_midpoint = (wavelength_new_emit2 + wavelength_new_emit3)/2

    wavelength_new_obs2 = (z+1)*wavelength_new_emit2  # Shift point 2
    wavelength_new_obs3 = (z+1)*wavelength_new_emit3  # Shift point 3
    wavelength_new_midpoint_obs = (z+1)*wavelength_new_midpoint  # Shift midpoint

    ######################### MIDDLE POINT ##################################
    q1 = np.where(data[:, 0] < wavelength_new_obs2)
    
    # Find the index of the wavelength in data closest to point 2
    p1 = np.max(q1)

    q2 = np.where(data[:, 0] > wavelength_new_midpoint_obs)
    # Find the index of the wavelength in data closest to the midpoint between point 2 and 3
    p2 = np.min(q2)

    # Get all flux values in the region between point 2 and our midpoint
    flux1 = data[p1:p2, 1]
   

    median_flux1 = np.median(flux1)  # Get the median flux in this region
    wavelength_flux1 = data[p1:p2, 0]  # Get the wavelengths in this region
    
    # Get the median wavelength in this region
    median_wavelength1 = np.median(wavelength_flux1)
    # Compose median wavelength and median flux into a point (b)
    middle_point = (median_flux1, median_wavelength1)
    
    ########################### END OF MIDDLE POINT ##############################
    
    ########################### LAST POINT###################################
    
    q3 = np.where(data[:, 0] < wavelength_new_midpoint_obs)
    # Find the index of the wavelength in data closest to the midpoint between point 2 and 3
    p3 = np.max(q3)

    q4 = np.where(data[:, 0] > wavelength_new_obs3)
    # Find the index of the wavelength in data closest to point 3
    p4 = np.min(q4)

    # Get all flux values in the region between our midpoint and point 3
    flux2 = data[p3:p4, 1]

    median_flux2 = np.median(flux2)  # Get the median flux in this region
    wavelength_flux2 = data[p3:p4, 0]  # Get all wavelengths in this region

    # Get the median wavelength in this region
    median_wavelength2 = np.median(wavelength_flux2)
    
    # Compose median wavelength and median flux into a point (c)
    last_point = (median_flux2, median_wavelength2)
    
    ########################### END OF LAST POINT ##############################
    
	
    ######################## D POINT AND THREE POINTS #######################
    # range taken in rest frame: 1415-1430
    dpoint_starting_point_restframe = 1415
    dpoint_ending_point_restframe = 1430

    # Shift start into frame
    dpoint_starting_point = (z+1)*(dpoint_starting_point_restframe)
    # Shift end into frame
    dpoint_ending_point = (z+1)*(dpoint_ending_point_restframe)

    # Get the indexes that correspond to wavelengths in data closest to our start and end points for this region
    qq6 = np.where(data[:, 0] < dpoint_starting_point)
    try:
        pp7 = np.max(qq6)
    except:
        pass

    qq8 = np.where(data[:, 0] > dpoint_ending_point)
    pp9 = np.min(qq8)

    flux33 = data[pp7:pp9, 1]  # Get all flux values in this region
   
    median_flux33 = np.median(flux33)  # Get the median flux for the region
    
    wavelength_flux33 = data[pp7:pp9, 0]  # Get all wavelengths in this region
    
    median_wavelength33 = np.median(wavelength_flux33) # Get the median wavelength for this region
    
    # Compose the median wavelength and median flux into a point (d)
    dpoint = (median_wavelength33, median_flux33)

    # AVERAGE ERROR FOR D POINT (WE HAVE TO DO THIS FOR THE FIRST IF STATEMENT BELOW)
    error33 = data[pp7:pp9, 2]
    
    median_flux_error33 = np.median(error33) #find the median of flux error
    
    st_dev_of_flux = np.std(flux33[:]) #Standart deviation of flux
   
  	# Compose error of dpoint average wavelength and average flux into a point (c)
    dpoint_error = (median_wavelength33, median_flux_error33)

    # THE THREE POINTS (THE THREE POINTS THAT THE original power law WILL USE)
    power_law_datax = (median_wavelength3, median_wavelength1, median_wavelength2)
    power_law_datay = (median_flux3, median_flux1, median_flux2)

    # THE THREE POINTS (THE THREE POINTS THAT THE second power law WILL USE)
    power_law_datax2 = (median_wavelength33, median_wavelength1, median_wavelength2)
    power_law_datay2 = (median_flux33, median_flux1, median_flux2)

    ############# END OF D POINT AND THREE POINTS #################################

    # BASICALLY, DEFINING MY WAVELENGTH, FLUX, AND ERROR (OR CHOOSING THEIR RANGE)
    wavelength_lower_limit = np.where(data[:, 0] > wavelength_observe1)
    wavelength_upper_limit = np.where(data[:, 0] < wavelength_observe2)
    
    # Get wavelengths in out data set that fall into our region of study
    wavelength = data[np.min(wavelength_lower_limit[0])
                             : np.max(wavelength_upper_limit[0]), 0]
    actual_wavelength = wavelength
    
    flux = data[np.min(wavelength_lower_limit[0]): np.max(
        wavelength_upper_limit[0]), 1]  # Get flux values in our region
    
    
    error = data[np.min(wavelength_lower_limit[0]): np.max(
        wavelength_upper_limit[0]), 2]  # Get error values in our region
    
    messed_up_error = np.where(data[np.min(wavelength_lower_limit[0]): np.max(
        wavelength_upper_limit[0]), 2] > 3)  # Get inexes of points with error > 3

    # Unshift(?) the wavelength, back to a rest frame
    wavelength_emit = wavelength/(z+1)

    # SOMETIMES, THERE ARE PIXEL PROBLEMS, AND WE MIGHT GET AN ERROR OF 30 IN FLUX. TO AVOID THAT, WE HAVE DONE THIS. MESSED UP ERROR IS
    # DEFINED ABOVE.
    # if len (messed_up_error[0]) > 0:######################################################original
    #error[messed_up_error[0]]=0####################################################
    # flux[messed_up_error[0]]=0

    fev = 10000

    print(power_law_datax)
    print(power_law_datay)

    # CURVE FIT FOR FIRST POWERLAW
    try:
        pars, covar = curve_fit(powerlaw, power_law_datax,
                                power_law_datay, p0=init_pars, maxfev=fev)
    except:
        print("Error - curve_fit failed-1st powerlaw " + i)
        powerlaw1_not_made.append(i)

    # CURVE FIT FOR SECOND POWERLAW
    try:
        pars2, covar = curve_fit(
            powerlaw, power_law_datax2, power_law_datay2, p0=init_pars2, maxfev=fev)

    # try:
        #popt,pcov = scipy.optimize.curve_fit(f, xdata, ydata, p0=None, sigma=None)

    except:

        print("Error - curve_fit failed-2nd powerlaw " + i)
        powerlaw2_not_made.append(i) 

    # (which are flux values) are on the fitted line, so they are not the 3 points that i chose
    m = []
    normalizing = flux/powerlaw(wavelength, *pars)
    actual_normalizing = normalizing
    error_normalized = error/powerlaw(wavelength, *pars)
    plerror_normalized = error_normalized
    
    for n in range(1, len(normalizing)-5):
        
        if abs(normalizing[n + 1] - normalizing[n]) > 0.5:  # new line

            if plerror_normalized[n+1] > 0.25:
                # normalized error
                plerror_normalized[n+1] = plerror_normalized[n]
                normalizing[n+1] = normalizing[n]  # normalized graph
                error[n+1] = error[n]  # original error
                flux[n+1] = flux[n]  # original graph

        if plerror_normalized[n] > 0.5:
            # normalized error
            plerror_normalized[n] = plerror_normalized[n-1]
            normalizing[n] = normalizing[n-1]  # normalized graph
            error[n] = error[n-1]  # original error
            flux[n] = flux[n-1]  # original graph

        if abs(normalizing[n+1] - normalizing[n]) > 5:
            # normalized error
            plerror_normalized[n+1] = plerror_normalized[n]
            normalizing[n+1] = normalizing[n]  # normalized graph
            error[n+1] = error[n]  # original error
            flux[n+1] = flux[n]  # original graph

   

    # APPEND IS SO TO MAKE AN ARRAY WITH INCREASING ELEMENT AS THE WHOLE CODE RUNS EACH TIME FOR EACH SPECTRA

    bf = pars[0]
    cf = pars[1]

    # REQUIREMENT FOR USING #POWERLAW.
    if (bf)*(np.power(median_wavelength33, cf)) < (median_flux33) - (3)*(median_flux_error33):

        ee.append(pars2[0])
        eee.append(pars2[1])
        ll.append(i)

    else:
        ee.append(pars[0])
        eee.append(pars[1])
        ll.append(i)

    # SNR Calculations:

    less_than_SNR = np.where (wavelength/(z+1.) < 1250.)
    greater_than_SNR = np.where (wavelength/(z+1.) > 1400.)
    
    snr_12001600mean = []
    snr_12001600median = []
    snr_1700 = []
    
   
    if len(less_than_SNR[0]) > 0:
		
        mean_error = np.mean(1./error_normalized[np.max(less_than_SNR[0]):np.min(greater_than_SNR)])
        snr_12001600mean.append(round(mean_error,5))
        
        
        median_error = np.median(1./error_normalized[np.max(less_than_SNR[0]):np.min(greater_than_SNR)])
        snr_12001600median.append(round(median_error,5))
        
        snr_1700.append(snr)


    # Start of Figure 1


    plt.figure(count_fig1)

    # IF STATEMENT IS: IF THE DIFFERENCE BTWN D POINT (RIGHT OF SiIV EMISSION) AND THE FIRST
    # POWER LAW AT THAT POINT IS MORE THAN 3 SIGMA,    #THEN USE THE NEW POWERLAW

    bf = pars[0]
    cf = pars[1]

	  # IF STATEMENT IS THERE TO AVOID ANY RANDOM, WEIRD PIXEL PROBLEM WITH THE ERRORS. FOR EXAMPLE: TO AVOID CASES WHERE THE ERROR IS 100
    
    if (bf)*(np.power(median_wavelength33, cf)) < (median_flux33) - (3)*(median_flux_error33):
        
        plt.plot(wavelength, powerlaw(wavelength, *pars2), 'r--')
        plt.plot(wavelength, powerlaw(wavelength, *pars), 'r--')
        plt.plot(median_wavelength33, median_flux33 - median_flux_error33, 'yo')
        plt.plot(median_wavelength33, median_flux33 -
             3*(median_flux_error33), 'yo')
        plt.plot(median_wavelength3, median_flux3, 'yo')
        plt.plot(median_wavelength33, median_flux33 - st_dev_of_flux, 'go')
        #plt.title(i)
        plt.xlabel("Wavelength[A]")
        plt.ylabel("Flux[10^[-17]]cgs")
        plt.text(wavelength_observe1-50, np.max(flux) - 5, "z = " + str(z) + " snr=" + str(snr)+ " snr_1326=" +str(snr_12001600mean[0]))
        plt.text(wavelength_observe1 + 1000, np.max(flux)-10, i)
        plt.plot(median_wavelength33, median_flux33, 'yo')
        plt.plot(wavelength, flux, 'b-')
        plt.plot(power_law_datax2, power_law_datay2, 'ro')
        plt.plot(wavelength, error, 'k-')
        
    else:

        plt.plot(wavelength, powerlaw(wavelength, *pars), 'r--')
        plt.plot(median_wavelength33, median_flux33 - median_flux_error33, 'yo')
        plt.plot(median_wavelength33, median_flux33 -
             3*(median_flux_error33), 'yo')
        plt.plot(median_wavelength33, median_flux33 - st_dev_of_flux, 'go')
        #plt.title(i)
        plt.xlabel("Wavelength[A]")
        plt.ylabel("Flux[10^[-17]]cgs")
        plt.text(wavelength_observe1 - 50, np.max(flux) - 5, "z = " + str(z) + " snr=" + str(snr)+ " snr_1325=" + str(snr_12001600mean[0]))
        plt.text(wavelength_observe1 + 1000, np.max(flux)-10, i)
        plt.plot(median_wavelength33, median_flux33, 'yo')
        plt.plot(wavelength, flux, 'b-')
        plt.plot(power_law_datax, power_law_datay, 'ro')
        plt.plot(wavelength, error, 'k-')


    pp1.savefig()
    plt.close(count_fig1)

    
    # End of Figure 1
    
    # Start of Figure 2
    """
    0.2, "z=" + str(z) + " snr=" + str(snr)

    plt.text(wavelength_observe1 + 1000, np.max(normalizing)-0.2, i)
    n1 = np.where(normalizing < 1)
    n2 = wavelength[n1] #515
    n3 = normalizing[n1] #515
    pp2.savefig()
    plt.close(count_fig2)
    plt.show()
    """
    
    # End of Figure 2
    
    www = (wavelength, normalizing, error_normalized)
    www = (np.transpose(www))
    oo=i[0:20]
    
    #########################################################################################
    # This is where normalized files will be saved. Remember to change to your own directory!
    
    np.savetxt('G:/School/_UWB Classes/ZCapstone/Capstone/Normalized-Spectre/files/' +
               oo+'norm.dr9', www)  # ,fmt='%s')
		#########################################################################################
    count_fig1 = count_fig1+1
    count_fig2 = count_fig2+1


############### THE FOR LOOP STARTED AT LINE ~115 FINISHED HERE #####################
#####################################################################################
#####################################################################################

# EVERYTHING BELOW IS NOT IN THE FOR LOOP

tt = [ll, ee, eee]
tt = (np.transpose(tt))

# i_all is the list of all good spectra.

for i in powerlaw2_not_made:
    for j in i_all:
        if j in powerlaw2_not_made:
            i_all.remove(j)

for l in powerlaw1_not_made:
    for lj in i_all:
        if lj in powerlaw1_not_made:
            i_all.remove(lj)
      
       
# Closes the PDF that saves all original files        
pp1.close()
# Closes the PDF that saves all normalized files
pp2.close()

# Merges Multiple PDF files



# np.savetxt("initial parameters 1.txt", t,fmt="%s")#INITIAL PARAMETERS FOR POWER LAW 1, THE fmt="%s" MAKES IT SO THAT THE FORMAT(fmt) is a string (=%s), rather than the default setting of float


file_name1 = 'G:/School/_UWB Classes/ZCapstone/Capstone/Normalized-Spectre/files/Final_Initial_Parameters.txt'


# np.savetxt("Final Initial Parameters .txt",tt,fmt="%s")#INITIAL PARAMETERS FOR POWER LAW 2

np.savetxt(file_name1, tt, fmt="%s")
file_name2 = 'G:/School/_UWB Classes/ZCapstone/Capstone/Normalized-Spectre/files/Powerlaw2_did_not_work.txt'
np.savetxt(file_name2, powerlaw2_not_made, fmt='%s')


file_name3 = 'G:/School/_UWB Classes/ZCapstone/Capstone/Normalized-Spectre/files/good_spectra.txt'
np.savetxt(file_name3, i_all, fmt='%s')
file_name4 = 'G:/School/_UWB Classes/ZCapstone/Capstone/Normalized-Spectre/files/Powerlaw1_did_not_work.txt'
np.savetxt(file_name4, powerlaw1_not_made, fmt='%s')

