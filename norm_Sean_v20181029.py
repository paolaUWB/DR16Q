import sys

#from pylab import *
from numpy import *
import scipy as sp  # I'm not using this

import csv

from matplotlib.pylab import*

from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
# import numpy as np# I'm not using this. The np you see in the code is a built in function
#from PyPDF2 import PdfFileMerger, PdfFileReader


# specnum=6760 UNCOMMENT THIS WHEN YOU ARE READY - THIS IS THE TRUE VALUE OF SPECNUM
specnum = 6760
# total number of quasars = 6760
# Set location of spectrum files
specdirec = '/Users/Paola/SDSSIII/DR9Q/FILES_FOR_ABDUL/'
pp1 = PdfPages('original_all_graph_Sean.pdf')  # Set output PDF file
pp2 = PdfPages('normalized_all_graph_Sean.pdf')  # Set normalized output PDF file

# specname_file = "confs/DR9Q_selection_specnames.lst" #|<-----The normal configs
# zem_file = "confs/DR9Q_selection_zem.lst"            #|

# specname_file = "confs/bad_specnames.lst" #|<------These are for "bad spectra"
# zem_file = "confs/bad_zem.lst"            #|


# USE AN ACTUAL CONFIG FILE!!!!!
config_file = sys.argv[1]  # Set cfg file path (csv with spec_name,z,snr...)
#config = loadtxt(config_file,
#                 delimiter=",",
#                 dtype={
#                     'names': ('spec', 'z', 'snr'),
#                     'formats': ('str', np.float, np.float)
#                 })

spectra_action = list() #[row[0] for row in config]
redshifts_action = list() #[row[1] for row in config]
snr_action = list() #[row[2] for row in config]

for line in open(config_file, 'r'):
	d = line.split(",")
	spectra_action.append(d[0])
	redshifts_action.append(np.float(d[1]))
	snr_action.append(np.float(d[2]))
	
# spectra = loadtxt(specname_file, dtype='str') #Load a list of spectra names
# spectra_action = spectra[:] #Make it into an array

# redshifts = loadtxt(zem_file) #Load a list of spectra redshifts
# redshifts_action = redshifts[:] #Make it into an array

# good_categorized_spectra = loadtxt('confs/good_categorized_spectra_through_code_subtracted_one.txt') #Limits spectra to normalize
#good_categorized_spectra_action = good_categorized_spectra[:]

good_categorized_spectra_action = range(0, len(spectra_action))


# COUNTERS FOR THE TWO FIGURES, THEY ARE SET APART (0 AND 20) AS TO NOT OVERLAP BECAUSE OTHERWISE, THERE WILL BE TWO GRAPHS IN ONE FIGURE
# (ON TOP OF EACH OTHER)
count_fig1 = 0
# count_fig2 will always start counting from twice the number of total spectra
count_fig2 = 2*specnum
original_graph_number = 0
normalized_graph_number = 0
rows_of_power_law_txt_file = len(spectra_action)
# these are empty arrays that i will add numbers to by appending them. Some of these are not used.
dr9 = []
e = []
l = []
ee = []
ll = []
eee = []
# a=25
# b,c are initial parameters for first powerlaw, and bb,cc are intial parameteres for the second powerlaw.
b = 1250
c = -0.5
bb = 1250
cc = -0.5
# DEFINING POWERLAW FUNCTIONS AND THEIR INITIAL PARAMETERS (START)
# FIRST POWERLAW FUNCTION


def powerlaw(x, b, c):
    return b*(np.power(x, c))

# SECOND POWERLAW FUNCTION


def powerlaw2(x, bb, cc):
    return bb*(np.power(x, cc))


# DEFINING POWERLAW FUNCTIONS AND THEIR INITIAL PARAMETERS (END)


# WE HAVE DEFINED THE INITIAL PARAMETERS AS B,C OR BB,CC SO ANYTIME WE TELL THE CODE TO PRINT INITIAL
# PARAMETERS, IT WILL DO SO IN THAT ORDER: B,C OR BB,CC
init_pars = [b, c]  # FOR 1ST POWERLAW
init_pars2 = [bb, cc]  # FOR 2ND POWERLAW


powerlaw1_not_made = []
powerlaw2_not_made = []
i_all = []
flux_remove = []
best_counter = 0


#ask = raw_input("Would you like to see potential absorption features? (y/n)  ")


powerlaw1_not_made = []
powerlaw2_not_made = []
i_all = []
flux_remove = []
best_counter = 0


#ask = raw_input("Would you like to see potential absorption features? (y/n)  ")
# for i in spectra_action_again:
# for i,j in zip(spectra_action, redshifts_action):

for qqqqq in good_categorized_spectra_action:
    best_counter = best_counter+1

    # Get spectra names and redshifts
    i = spectra_action[qqqqq]
    i_all.append(i)
    j = redshifts_action[qqqqq]
    snr = round(snr_action[qqqqq], 5)

    print(str(best_counter) + ": " + i)

    j = round(j, 5)  # TO ROUND OFF J(REDSHIFT) TO 5 CHARACTERS
    z = j

    #print "DR9 file: " + `i`,   "-------->  z=" + `j`
    # i= THE CURRENT SPECTRA THAT THE LOOP IS GOING THROUGH
    # j= THE CURRENT REDSHIFT THAT THE LOOP IS GOING THROUGH

    #data = loadtxt('spec-4774-55659-0140.DR9')
    data = np.loadtxt(specdirec + i)  # Load in spectrum file
    number_rows = len(data[:, :])  # Get the number of points in the spectrum
    # i only care about spectrum from 1200 - 1800 rest frame wavelength
    # z=2.102
    j = round(j, 5)  # Redundant?
    z = j

    # THE POINTS USED FOR POWERLAW (START)
    # THIS IS THE RANGE OF THE ELECTROMAGNETIC SPECTRUM THAT WE WANT OUR SPECTRA TO BE. wE DONT CARE ABOUT OTHER REGIONS.
    wavelength_emit1_initial = 1200
    wavelength_emit2_initial = 1800
    # Shift start wavelength into frame |<--This makes our wavelength range for
    wavelength_observe1 = (z+1)*wavelength_emit1_initial
    # Shift end wavelength into frame   |     the region we want to look at
    wavelength_observe2 = (z+1)*wavelength_emit2_initial

    wavelength_NV_emit = 1242.8040
    # Shift Nitrogen V (NV) line into frame
    wavelength_NV_obs = (z+1)*wavelength_NV_emit

    # REST-FRAME WAVELENGTH RANGE FOR FIRST POINT
    wavelength_restframe_starting_point = 1280.206
    wavelength_restframe_ending_point = 1284.333
    # Shift start wavelength of first point into frame |<--This makes the range for the first
    wavelength_observed_starting_point = (
        z+1)*(wavelength_restframe_starting_point)
    # Shift ebd wavelength of first point into frame   |     point in the spectrum
    wavelength_observed_ending_point = (
        z+1)*(wavelength_restframe_ending_point)

    # a good range of rest frame wavelength is = 1686 - 1773
    # their midpoint is 1729, so im gonna take 1686-1729, average them up, and find the midpoint, then find midpoint of 1729-1773
    wavelength_new_emit1 = 1282.398  # 1st point for powerlaw
    # Shift first power law point (a) into rest frame
    wavelength_new_obs1 = (z+1)*wavelength_new_emit1
    # FIRST POINT
    # Get all points from data with wavelengths less than our starting wavelength for our first point
    q6 = np.where(data[:, 0] < wavelength_observed_starting_point)

    p7 = 0
    try:  # Why is this needed?
        # Get the largest wavelength from q6. Is the index of the closest wavelength in the sprctrum to the starting wavelength for our first point
        p7 = np.max(q6)
    except:
        pass

    # Get all points from data with wavelength greater than our ending wavelengt for our first point
    q8 = np.where(data[:, 0] > wavelength_observed_ending_point)
    
    # Get the smallest wavelength from q8. Is the index of the closest wavelength in the sprctrum to the ending wavelength for our first point
    p9 = np.min(q8)

    flux3 = data[p7:p9, 1]  # Get all flux values for our starting point
    length_flux3 = len(flux3)  # Unused, should be commented out or removed.
    #average_flux3 = sum (flux3) / length_flux3
    # Get the median flux in the range. Call it an average?
    average_flux3 = np.median(flux3)
    # Get the corresponding wavelengths for this ragion
    wavelength_flux3 = data[p7:p9, 0]
    #average_wavelength3 = sum(wavelength_flux3) / length_flux3
    # Get the median wavelength of this region. Call is an average?
    average_wavelength3 = np.median(wavelength_flux3)
    # Make a point from the median wavelength and median flux (a)
    first_point = (average_wavelength3, average_flux3)

    print(flux3)

    wavelength_new_emit2 = 1677.938  # 2nd point for powerlaw
    wavelength_new_emit3 = 1725.669  # 3rd point for powerlaw
    # Average the points to get a midpoint
    wavelength_new_midpoint = (wavelength_new_emit2 + wavelength_new_emit3)/2

    wavelength_new_obs2 = (z+1)*wavelength_new_emit2  # Shift point 2
    wavelength_new_obs3 = (z+1)*wavelength_new_emit3  # Shift point 3
    wavelength_new_midpoint_obs = (
        z+1)*wavelength_new_midpoint  # Shift midpoint

    # MIDDLE POINT
    q1 = np.where(data[:, 0] < wavelength_new_obs2)
    # Find the index of the wavelength in data closest to point 2
    p1 = np.max(q1)

    q2 = np.where(data[:, 0] > wavelength_new_midpoint_obs)
    # Find the index of the wavelength in data closest to the midpoint between point 2 and 3
    p2 = np.min(q2)

    # Get all flux values in the region between point 2 and our midpoint
    flux1 = data[p1:p2, 1]
    length_flux1 = len(flux1)  # Unuesd, again?
    #average_flux1 = sum (flux1) / length_flux1
    average_flux1 = np.median(flux1)  # Get the median flux in this region
    wavelength_flux1 = data[p1:p2, 0]  # Get the wavelengths in this region
    #average_wavelength1 = sum(wavelength_flux1) / length_flux1
    # Get the median wavelength in this region
    average_wavelength1 = np.median(wavelength_flux1)
    # Compose median wavelength and median flux into a point (b)
    middle_point = (average_flux1, average_wavelength1)

    # LAST POINT
    q3 = np.where(data[:, 0] < wavelength_new_midpoint_obs)
    # Find the index of the wavelength in data closest to the midpoint between point 2 and 3
    p3 = np.max(q3)

    q4 = np.where(data[:, 0] > wavelength_new_obs3)
    # Find the index of the wavelength in data closest to point 3
    p4 = np.min(q4)

    # Get all flux values in the region between our midpoint and point 3
    flux2 = data[p3:p4, 1]
    length_flux2 = len(flux2)  # Again, unused.
    #average_flux2 = sum (flux2) / length_flux2
    average_flux2 = np.median(flux2)  # Get the median flux in this region
    wavelength_flux2 = data[p3:p4, 0]  # Get all wavelengths in this region

    #average_wavelength2 = sum(wavelength_flux2) / length_flux2
    # Get the median wavelength in this region
    average_wavelength2 = np.median(wavelength_flux2)
    # Compose median wavelength and median flux into a point (c)
    last_point = (average_flux2, average_wavelength2)

    # D POINT
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
    # length_flux33 = len(flux33) #Boom! I knew this was unused.
    #average_flux33 = sum (flux33) / length_flux33
    average_flux33 = np.median(flux33)  # Get the median flux for the region
    wavelength_flux33 = data[pp7:pp9, 0]  # Get all wavelengths in this region
    #average_wavelength33 = sum(wavelength_flux33) / length_flux33
    # Get the median wavelength for this region
    average_wavelength33 = np.median(wavelength_flux33)
    # Compose the median wavelength and median flux into a point (d)
    dpoint = (average_wavelength33, average_flux33)

    # AVERAGE ERROR FOR D POINT (WE HAVE TO DO THIS FOR THE FIRST IF STATEMENT BELOW)
    error33 = data[pp7:pp9, 2]
    #length_error33 = len(error33)
    #average_flux_error33 = sum(error33)/length_error33
    average_flux_error33 = np.median(error33)
    rms = np.std(flux33[:])
    #average_wavelength_error33 =data[pp7:pp9,0]
    #average_flux_error33 = sum (average_wavelength_error33)/length_error33
    dpoint_error = (average_wavelength33, average_flux_error33)

    # THE THREE POINTS (THE THREE POINTS THAT THE original power law WILL USE)
    power_law_datax = (average_wavelength3,
                       average_wavelength1, average_wavelength2)
    power_law_datay = (average_flux3, average_flux1, average_flux2)

    # THE THREE POINTS (THE THREE POINTS THAT THE second power law WILL USE)
    power_law_datax2 = (average_wavelength33,
                        average_wavelength1, average_wavelength2)
    power_law_datay2 = (average_flux33, average_flux1, average_flux2)

    # THE POINTS USED FOR POWERLAW (END)

    # BASICALLY, DEFINING MY WAVELENGTH, FLUX, AND ERROR (OR CHOOSING THEIR RANGE)
    wavelength_lower_limit = np.where(data[:, 0] > wavelength_observe1)
    wavelength_upper_limit = np.where(data[:, 0] < wavelength_observe2)
    # Get wavelengths in out data set that fall into our region of study
    wavelength = data[np.min(wavelength_lower_limit[0])
                             : np.max(wavelength_upper_limit[0]), 0]
    actual_wavelength = wavelength
    flux = data[np.min(wavelength_lower_limit[0]): np.max(
        wavelength_upper_limit[0]), 1]  # Get flux values in our region
    actual_flux = flux
    error = data[np.min(wavelength_lower_limit[0]): np.max(
        wavelength_upper_limit[0]), 2]  # Get error values in our region
    messed_up_error = np.where(data[np.min(wavelength_lower_limit[0]): np.max(
        wavelength_upper_limit[0]), 2] > 3)  # Get inexes of points with error > 3

    plerror = error

    # Unshift(?) the wavelength, back to a rest frame
    wavelength_emit = wavelength/(z+1)

    # SOMETIMES, THERE ARE PIXEL PROBLEMS, AND WE MIGHT GET AN ERROR OF 30 IN FLUX. TO AVOID THAT, WE HAVE DONE THIS. MESSED UP ERROR IS
    # DEFINED ABOVE.
    # if len (messed_up_error[0]) > 0:######################################################original
    #plerror[messed_up_error[0]]=0####################################################
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
            powerlaw2, power_law_datax2, power_law_datay2, p0=init_pars2, maxfev=fev)

    # try:
        #popt,pcov = scipy.optimize.curve_fit(f, xdata, ydata, p0=None, sigma=None)

    except:

        print("Error - curve_fit failed-2nd powerlaw " + i)
        powerlaw2_not_made.append(i) 

     # else:
        # i_all.append(i)

    # the powerlaw(power_law_datax,*pars) (which are flux values) are on the fitted line, so they are not the 3 points that i chose
    mmmoiz = []
    normalizing = flux/powerlaw(wavelength, *pars)
    actual_normalizing = normalizing
    error_normalized = error/powerlaw(wavelength, *pars)
    plerror_normalized = error_normalized
    for mmoiz in range(1, len(normalizing)-5):
        """
        if abs (normalizing [mmoiz + 1]- normalizing[mmoiz] )>0.5:
                if abs (normalizing [mmoiz +2] - normalizing [mmoiz+1]) >0.5 :
                        normalizing [mmoiz +1] = normalizing [mmoiz +7]
                        plerror_normalized [mmoiz +1] = 0
                        mmmoiz.append (mmoiz)
        """
        if abs(normalizing[mmoiz + 1] - normalizing[mmoiz]) > 0.5:  # new line

            if plerror_normalized[mmoiz+1] > 0.25:
                # normalized error
                plerror_normalized[mmoiz+1] = plerror_normalized[mmoiz]
                normalizing[mmoiz+1] = normalizing[mmoiz]  # normalized graph
                plerror[mmoiz+1] = plerror[mmoiz]  # original error
                flux[mmoiz+1] = flux[mmoiz]  # original graph

        if plerror_normalized[mmoiz] > 0.5:
            # normalized error
            plerror_normalized[mmoiz] = plerror_normalized[mmoiz-1]
            normalizing[mmoiz] = normalizing[mmoiz-1]  # normalized graph
            plerror[mmoiz] = plerror[mmoiz-1]  # original error
            flux[mmoiz] = flux[mmoiz-1]  # original graph

        if abs(normalizing[mmoiz+1] - normalizing[mmoiz]) > 5:
            # normalized error
            plerror_normalized[mmoiz+1] = plerror_normalized[mmoiz]
            normalizing[mmoiz+1] = normalizing[mmoiz]  # normalized graph
            plerror[mmoiz+1] = plerror[mmoiz]  # original error
            flux[mmoiz+1] = flux[mmoiz]  # original graph

     # for secondloop in range (0,len (normalizing)):

        """
	if normalizing [secondloop] < plerror_normalized [secondloop]:
		
		plerror_normalized [secondloop]=0   #plerror [secondloop]=0
        """

    #print 'before pixel fix, normalizing [moiz] = ' + `normalizing [moiz]`
    # BAD PIXELS
    """
    if len (messed_up_error [0])> 0:

    	
	for moiz in messed_up_error[0]:
		print 'before pixel fix, normalizing [moiz] = ' + `normalizing [moiz]`
			
		normalizing[moiz]=normalizing[moiz+10]
		print 'normalizing [moiz] = ' + `normalizing [moiz+10]`
    """
    # BAD PIXELS

    # APPEND IS SO TO MAKE AN ARRAY WITH INCREASING ELEMENT AS THE WHOLE CODE RUNS EACH TIME FOR EACH SPECTRA
    #print 'look here ' + `flux[moiz]`

    #print "flux =" +`average_flux33`
    #print "error =" +`average_flux_error33`
    bf = pars[0]
    cf = pars[1]

    # REQUIREMENT FOR USING THE SECOND 		#POWERLAW.
    if (bf)*(np.power(average_wavelength33, cf)) < (average_flux33) - (3)*(average_flux_error33):

        ee.append(pars2[0])
        eee.append(pars2[1])
        ll.append(i)

    else:
        ee.append(pars[0])
        eee.append(pars[1])
        ll.append(i)

    # SNR Calculations:

    oo1 = np.where (wavelength/(z+1.) < 1250.)
    oo2 = np.where (wavelength/(z+1.) > 1400.)
    
   
    if len(oo1[0]) > 0:

        jiji = np.mean(1./error_normalized[np.max(oo1[0]):np.min(oo2)])
        snr_12001600mean.append(round(jiji,5))
        jeje = np.median(1./error_normalized[np.max(oo1[0]):np.min(oo2)])
        snr_12001600median.append(round(jeje,5))
        
        snr_1700.append(snr)


    # FIGURE 1 (START)

    #print ('now graphing')

    plt.figure(count_fig1)
    #original_graph_number = original_graph_number+1

    # IF STATEMENT IS: IF THE DIFFERENCE BTWN D POINT (RIGHT OF SiIV EMISSION) AND THE FIRST
    # POWER LAW AT THAT POINT IS MORE THAN 3 SIGMA,    #THEN USE THE NEW POWERLAW
    # if ((average_flux33  - (b)*(average_wavelength33)**(c)))>=3* (average_flux_error33):#NOT EXACTLY SURE Y WE R DOING THIS
    # if (average_flux33  - 3*(average_flux_error33) > (b)*(average_wavelength33)**(c)):

    bf = pars[0]
    cf = pars[1]
    #print "value powerlaw at D" + `(bf)*(np.power (average_wavelength33,cf))`
    #print "average_flux" + `(average_flux33)  - (3)*(average_flux_error33)`

    if (bf)*(np.power(average_wavelength33, cf)) < (average_flux33) - (3)*(average_flux_error33):
        plt.plot(wavelength, powerlaw2(wavelength, *pars2), 'r--')
        plt.plot(wavelength, powerlaw(wavelength, *pars), 'r--')
        plt.plot(average_wavelength33, average_flux33 - average_flux_error33, 'yo')
        plt.plot(average_wavelength33, average_flux33 -
             3*(average_flux_error33), 'yo')
        plt.plot(average_wavelength3, average_flux3, 'yo')
        plt.plot(average_wavelength33, average_flux33 - rms, 'go')

        #plot(wavelength, powerlaw(wavelength,*pars),'r--')
        plt.title("original data vs error")
        plt.xlabel("Wavelength[A]")
        plt.ylabel("Flux[10^[-17]]cgs")
        #normalizing = flux/ powerlaw(power_law_datax,*pars)
        #normalize = (plot(wavelength, flux))/(plot((power_law_datax, powerlaw(power_law_datax,*pars))))
        #text(wavelength_observe1 +510,np.max(flux)-5,"z = " +`j`)
        plt.text(wavelength_observe1-50, np.max(flux) -
             5, "z = " + str(j) + " snr=" + str(snr)+ " snr_1326=" +str(snr_12001600mean[qqqqq]))
        plt.text(wavelength_observe1 + 1000, np.max(flux)-10, i)
        plt.plot(average_wavelength33, average_flux33, 'yo')
        #print dpoint
        plt.plot(wavelength, flux, 'b-')
        #print 'this is '+`flux[moiz]`
        plt.plot(power_law_datax2, power_law_datay2, 'ro')
        plt.plot(wavelength, plerror, 'k-')
    else:

        plt.plot(wavelength, powerlaw(wavelength, *pars), 'r--')

        #plot(wavelength, powerlaw(wavelength,*pars),'r--')
        #print 'this is '+`flux[moiz]`
        plt.plot(average_wavelength33, average_flux33 - average_flux_error33, 'yo')
        plt.plot(average_wavelength33, average_flux33 -
             3*(average_flux_error33), 'yo')
        plt.plot(average_wavelength33, average_flux33 - rms, 'go')

        plt.title("original data vs error")
        plt.xlabel("Wavelength[A]")
        plt.ylabel("Flux[10^[-17]]cgs")
        #normalizing = flux/ powerlaw(power_law_datax,*pars)
        #normalize = (plot(wavelength, flux))/(plot((power_law_datax, powerlaw(power_law_datax,*pars))))
        plt.text(wavelength_observe1 - 50, np.max(flux) -
             5, "z = " + str(j) + " snr=" + str(snr)+ " snr_1325=" + str(snr_12001600mean[qqqqq]))
        #text(wavelength_observe1 +350,np.max(flux)-5,"z = " +`j`)
        plt.text(wavelength_observe1 + 1000, np.max(flux)-10, i)
        plt.plot(average_wavelength33, average_flux33, 'yo')
        plt.plot(wavelength, flux, 'b-')
        plt.plot(power_law_datax, power_law_datay, 'ro')
        #print 'this is '+`flux[moiz]`
        # IF STATEMENT IS THERE TO AVOID ANY RANDOM, WEIRD PIXEL PROBLEM WITH THE ERRORS. FOR EXAMPLE: TO AVOID CASES WHERE THE ERROR IS 100

        plt.plot(wavelength, plerror, 'k-')

    pp1.savefig()
    plt.close(count_fig1)
    # 'pp'+`original_hellos`.savefig()
    """
    if original_graph_number <=200:
	pp_original.savefig()
    if original_graph_number >200:
	pp_original.close()
	original_hellos=original_hellos+1
	pp_original=PdfPages('original_sixth_many_graphs'+`original_hellos`+'.pdf')  
	original_graph_number=0
    """

    # FIGURE 1 (END)
         0.2, "z=" + str(j) + " snr=" + str(snr))

    text(wavelength_observe1 + 1000, np.max(normalizing)-0.2, i)
    n1 = where(normalizing < 1)

    n2 = wavelength[n1]
    n3 = normalizing[n1]
    pp2.savefig()
    close(count_fig2)
    """
    if normalized_graph_number <=200:
	pp_normalized.savefig()
    if normalized_graph_number >200:
	pp_normalized.close()
	normalized_hellos=normalized_hellos+1
	pp_normalized=PdfPages('normalized_sixth_many_graphs'+`normalized_hellos`+'.pdf')  
	normalized_graph_number=0
    """

    #print ('graphing ended')

    # FIGURE 2 (END)
    www = (wavelength, normalizing, error_normalized)
    www = (np.transpose(www))
    #i=i.replace('\'', '')
    oo=i[0:20]
    np.savetxt('/Users/Paola/QUASAR/Work_EHVO/ROUTINES/norm/' +
               oo+'norm.dr9', www)  # ,fmt='%s')

    count_fig1 = count_fig1+1
    count_fig2 = count_fig2+1
    # show()
# EVERYTHING BELOW IS NOT IN THE FOR LOOP


k = 0
t = []
g = []
# tt=[]
ty = 0
kk = 0


# ll=array(ll)
# ee=array(ee)
tt = [ll, ee, eee]
tt = (np.transpose(tt))


# tt=array(tt)


#whereto=np.where (i_all=powerlaw2_not_made)
for i in powerlaw2_not_made:
    for j in i_all:
        if j in powerlaw2_not_made:
            i_all.remove(j)

for l in powerlaw1_not_made:
    for lj in i_all:
        if lj in powerlaw1_not_made:
            i_all.remove(lj)

# i_all is the list of all good spectra.
pp1.close()
pp2.close()
# pp+pdfs.close()
# trying to merge multiple pdf files
"""  
for earth in range (0,pdfs):
	earth1=('normalized_sixth_no_original_graph' + `earth` +'.pdf')
	merger.append(PdfFileReader(file(earth1, 'rb')))

merger.write("document-output.pdf")
"""

# trying to merge multiple pdf files

# pp1.close()#CLOSES THE PDF THAT SAVES ALL ORIGINAL FILES
# pp2.close()#CLOSES THE PDF THAT SAVES ALL NORMALIZED FILES
# np.savetxt("initial parameters 1.txt", t,fmt="%s")#INITIAL PARAMETERS FOR POWER LAW 1, THE fmt="%s" MAKES IT SO THAT THE FORMAT(fmt) is a string (=%s), rather than the default setting of float

# specdirec='/home/khatri/SDSSIII/SPECTRA/'
file_name1 = '/Users/Paola/QUASAR/Work_EHVO/ROUTINES/Final_Initial_Parameters.txt'


#file_name1 = '/home/khatri/Documents/Python_code/Final_Initial_Parameters.txt'

# np.savetxt("Final Initial Parameters .txt",tt,fmt="%s")#INITIAL PARAMETERS FOR POWER LAW 2
np.savetxt(file_name1, tt, fmt="%s")
file_name2 = '/Users/Paola/QUASAR/Work_EHVO/ROUTINES/Powerlaw2_did_not_work.txt'
np.savetxt(file_name2, powerlaw2_not_made, fmt='%s')

file_name3 = '/Users/Paola/QUASAR/Work_EHVO/ROUTINES/good_spectra.txt'
np.savetxt(file_name3, i_all, fmt='%s')
file_name4 = '/Users/Paola/QUASAR/Work_EHVO/ROUTINES/Powerlaw1_did_not_work.txt'
np.savetxt(file_name4, powerlaw1_not_made, fmt='%s')
# np.savetxt("normalized_spectra.DR9",www,fmt="%s")
# text_file.close()
# np.savetxt("normalized_dr9.dr9",
