# DR16Q 

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Welcome to the DR16Q EHVO Github!

Steps to contribute to this project:

1. Look through the "Welcome! First steps for PHYS students" document on in the Google Drive. This document contains information about the research project [how to register for credits, general resources about quasars & quasar outflows, some of the posters/presentations we have made, etc.].

2. Look through the GitHub Welcome Packet in the Google Drive. This document contains information about how GitHub works [how to get the repository on your own computer, how to push changes, etc.]. 

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### GENERAL NOTES

-In this repository, all code is written in Python. To run the code, you can use either Spyder IDE or Visual Studio Code (VS Code). Regardless of which you choose, be sure to launch the app through Anaconda.

-Any major change should be added to the README file.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### FILES IN DIRECTORY
-abs_function_module.py **WFGN**
-abs_plot_module.py **WFGN**
-absorption.py **WFGN**
-data_types.py
    This file contains definitions of tuples that are used in the normalization.py code. 
-depth_graph_test.py **WFGN**
-depth_testing.py **WFGN**
-DR9_redshift.py
    This code creates the redshift histogram plots for the DR9 data. 
-DR9_sorted_norm.csv
    File that contains all necessary information/data to run normalization.py. This is NOT the parent sample - all spectra we have downloaded from the SDSS are included in this file. 
    **Columns:** Spectra Name, Redshift, SDSS SNR
-DR9Q_selection_minus17.dat
    File that contains DR9 data used to make redshift histogram plots in the DR9_redshift.py code.
-DR16_BAL_parent_sample.csv
    File that contains the DR16 BALQSOs --> these are the BALQSOs that meet our criteria AND are in our parent sample
-DR16_EHVO_sorted_norm.csv
    File that contains the final DR16 EHVO list
-DR16_parent_sample.csv
    File that contains the final DR16 parent sample
-DR16_redshift.py
    This code creates the redshift histogram plots for the DR16 data.
-DR16_sorted_norm.csv
    File that contains all necessary information/data to run normalization.py. This is NOT the parent sample - all spectra we have downloaded from the SDSS are included in this file. 
    **Columns:** Spectra Name, Redshift, SDSS SNR
-draw_figures.py
    This code contains the functions used by the normalization.py code to plot spectra. 
-draw_histogram.py
    This code contains the functions used by the DRX_redshift.py codes to plot the histograms to analyze redshifts of our samples.
-EHVO_DR9.dat
    File that contains the DR9 EHVO data for the DR9_redshift.py code
-normalization.py
    This code normalizes SDSS quasar spectra: fits each spectrum with a power law, performs tests to check the fit of the power law, and normalizes the spectrum if the fit is good. 
-plotindividualspectrum_Figof20.py
    Pretty plots for posters/presentations/etc. --> plots CIV absorption on normalized plot (includes normalized line)
-plotindividualspectrum.py
    Pretty plots for posters/presentations/etc. --> plots all absorption in the plot (must be done manually)
-useful_wavelength_flux_error_modules.py
    Contains the functions used to calculate the location of the anchor points, calculate SNR, etc.
-utility_functions.py
    Contains the functions to read/open/print to files and other general functions. 

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### NORMALIZATION FILE [normalization.py]

-Files needed to run normalization.py: 
    -data_types.py
    -DRX_sorted_norm.csv [X being the current data release; i.e. DR9_sorted_norm.csv]
    -draw_figures.py
        This file contains the functions to plot the spectra when you run the normalization code.
    -useful_wavelength_flux_error_modules.py 
        This file contains functions that calculate SNR and calculate anchor point locations.
    -utility_functions.py
        This file contains basic functions for opening files, printing to files, etc. 


-At the top of the file there is a section containing ranges of wavelengths under the "DO NOT CHANGE" heading. These constant variables are defined by PRH - ask before changing these values. 

-The first variable defined is where you state which Data Release you would like to work with. There are currently files for the DR9 and DR16. 
    To run the code, make sure the DR variable is set to the proper Data Release number (i.e. '16' for DR16). 

-The range you choose for STARTS_FROM, ENDS_AT will depend on the Data Release you are working with. We currently have 6760 files for DR9 and 21859 files for DR16. 
  
  -For DR9 the current range that works is either 1, 10 or 899, 1527 <-- this is because of the data we currently have in the repository. If you would like to run the DR9 data, you will need to get the full data from PRH.
  
  -For DR16 the current range that works is 1, 21823. To run through the range 21824, 21859 set the variable: dynamic = 'yes'. This puts the code into a manual plotting mode where you specify the locations for the anchor points. 

# To download and save the NORM files:

1. In the STUDENT WORK DRQ16 EHVO PROJECT drive open Normalization DR16Q

2. Open Norm Files

3. Right click on NORM_DR16Q.zip and download to your computer

4. The DR16Q repository should be saved in a folder called “GitHub” on your computer. Save the NORM_DR16Q.zip file in the GitHub folder (NOT in the DR16Q repository)
    --> Saving the NORM files in the repository causes issues when pushing to GitHub since there are so many files. 

6. Confirm that the version of the code you are working with has the proper path for the NORM_DIREC:  NORM_DIREC = os.getcwd() + '/../' + "NORM_DR16Q/"

# Variables to change

-save_new_output_file should be set to 'yes' for the csv output files to be created 
    these save in OUTPUT_FILES --> NORMALIZATION

-save_new_norm_file should be set to 'yes' for new norm.drX files to be created

-save_figures should be set to 'yes' for the pdf plots to be created
    these are the plots of the spectra - the code runs faster when it does not plot them, so if you don't need the plots this is a way to speed things up.

-sm should be set to 'yes' for the plots to be smoothed. NOTE: for smoothing, the variable BOXCAR_SIZE must be set to an odd values

-dynamic should be set to 'yes' for manual fitting - this allows for manual choosing of anchor points for the plotting (typically needs to be set to 'yes' for redshift >= 5.11 or the code may throw errors)
    when using manual/dynamic fitting, make sure flag_spectra is set to 'no'

-flag_spectra should be set to 'yes' for the plots to be flagged (good fit/bad fit)

# OUTPUT FILES (csv):

-original.csv
    -spectra index, spectra file name, chi squared
    -the chi squared values for all spectra with SNR>10 are added to this file
    
-good_fit.csv
    -spectra index, spectra file name, norm spectra file name, redshift, calculated SNR, SDSS SNR, bf, cf
    -any spectra deemed to be a good fit originally are added to this file
    -any spectra that has been unflagged by TEST #3 are added to this file
    -the spectra in this file ARE NORMALIZED
    
-unflagged.csv
    -spectra index, spectra file name, norm spectra file name, redshift, calculated SNR, SDSS SNR, bf, cf
    -spectra that were flagged by TEST #2 but have been deemed ok fits when the upper and lower limits of 'closeness' are increased slightly
    -the spectra in this file ARE NORMALIZED [since they are added to the good_fit file]

-flagged_bad_fit.csv
    -spectra index, spectra file name, norm spectra file name, redshift, calculated SNR, SDSS SNR, bf, cf
    -if the powerlaw does not go close enough through the anchor points, the spectra is added to this file
    -if the powerlaw is not close enough through the center of the green and pink test regions, the spectra is added to this file
    -the spectra in this file are NOT NORMALIZED
    
-flagged_absorption.csv
    -spectra index, spectra file name, norm spectra file name, redshift, calculated SNR, SDSS SNR, bf, cf
    -if the green OR pink test region is completely below the powerlaw, the spectra is added to this file
    - TEST #4 determines this
    -NOTE: these spectra are sorted into either good_fit.csv or flagged_bad_fit.csv. These are not included in final count to check.

-flagged_snr_in_ehvo.csv
    -spectra index, spectra file name, SDSS SNR, calculated SNR
    -any spectra with SNR<10 in the region we care about are added to this file
    -spectra are essentially ignored once they are flagged for low SNR

-log.txt
    -a log of the printed outputs for ALL spectra in the current run (including those with SNR<10)
    
-log_no_low_snr.txt
    -same outputs as log.txt but only for spectra with SNR>10
    
original.csv + flagged_snr_in_ehvo.csv = total # spectra in run 
good_fit.csv + flagged_bad_fit.csv = original.csv
    
# PLOTS (pdf files):

-After you run this code, all your graphs will be added to a pdf file in your directory. None of these files contain graphs for spectra with SNR less than 10 
    -original_graphs.pdf
        -all spectra in the range provided that have SNR>10 are plotted and added to this file
    -normalized_graphs.pdf
        -all spefctra that have been deemed good fits are normalized and the normalized spectra are plotted and added to this file [contains good_fit and unflagged]
    -good_fit_graphs.pdf
        -all spectra that have SNR>10 that have been deemed to be a good fit originally are plotted and added to this file
        -all spectra that have been unflagged have been plotted and added to this file
    -unflagged_graphs.pdf
        -any spectra that were previously flagged that have been deemed a good fit are unflagged, plotted and added to this file [TEST #3] 
    -flagged_bad_fit_graphs.pdf
        -all spectra that have been deemed a bad fit by the tests are flagged, plotted, and added to this file [TEST #1 & TEST #2]   
    -flagged_absorption_graphs.pdf
        -any spectra that is flagged for TEST #4 are added to this file due to possible absorption in the test regions [TEST #4]

original_graphs.pdf ~ original.csv
good_fit_graphs.pdf ~ flagged_bad_fit_graphs.pdf ~ original_graphs.pdf
good_fit_graphs.pdf ~ normalized_graphs.pdf

# Normalization Tests:

-TEST #1
    -checks whether the fit of the powerlaw is going closely through the anchor points
        -if the powerlaw goes through all three anchor points well, it proceeds to the other three tests
        -if the powerlaw is too far from the anchor points, two other locations for the middle anchor point are attempted. if the powerlaw is still too far, the spectra is flagged and moved to the flagged_bad_fit file. 

-TEST #2
    -checks whether the powerlaw goes through close to the center of the green and pink test regions
        -if the fit is good, it proceeds to TEST #4 (and skips TEST #3)
        -if the fit is bad, the spectra is flagged but continues on to TEST #3
    
-TEST #3
    -checks if the powerlaw is generally in the center of the pink and green test regions
        -if the fit is good enough in BOTH regions, this test unflags the spectra and adds it to the unflagged file as well as the good_fit file
        -if the fit is still not good, the spectra remains flagged and is added to the flagged_bad_fit file
    
-TEST #4
    -checks if the green or pink test region is completely below the powerlaw (a check for possible absorption in these regions)
        -if either the green or pink region is completely below the powerlaw, the spectra is added to the flagged absorption file
        -if neither region are completely below the powerlaw, nothing is done

-There is a test file in the directory. If you change something wrong, the test file will catch the different results. If you change constant variables, this will cause different output. 
  If changes are correct then update the test file with new results for future tests.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### ABSORPTION FILE

-Follow the GitHub steps above and open the file "absorption.py" using Spyder IDE or Visual Studio Code.

-File currently has DR9 "Data Release 9" extension. In the future, it will be changed to the recent Data Release from SDSS Database
 
-The Absorption file reads the Normalization files. After you generate the normalization file for each spectra, the Absorption file will read that and generate necessary graphs and information. 
  Currently the format of the processed Normalization file is : "spectra name + norm + .dr9"

-Based on "BALNICITY_INDEX_LIMIT", the file will import different results. You can see the results in the same directory under the "Absortion_cleaning" document. 
  Also the same "BALNICITY_INDEX_LIMIT" constant will generate different graph output because it is defining the narrowness of the output.

- At the top of the file, the constant variables are defined by Astrophysicist. Any change of these constants SHOULD discuss with the client.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

BAL_PARENT_SAMPLE.txt:
- Columns: SPECTRA NAME, Redshift, BI Value, BI ERROR, SNR
### .gitignore

-Any files that you would like to be ignored when pushing should be recorded in this file.
** ADD MORE DETAILS **