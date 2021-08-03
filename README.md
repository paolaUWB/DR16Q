# Normalized-Spectra

Welcome to Github!

Steps to contribute to this project:

1. Create a Github account. Once the account is created, the owner of the repository will need to give you access.

2. Once you have access to the repository, Press on the "Pull requests" on the top left tabs of the repository, and press the green button "New pull request on the top right, this allows any changes that you make to be on a separate branch that will not affect the master file until it has been approved. 

3. To make changes, you can choose to use other programs like Git Desktop for an easy interface or use Github itself by clicking on the file you want to make changes to in your branch, click the pencil at the top right, and paste the changed contents into the textbox.
  - Another way of doing this is pressing the "Push" button and it will allow you to paste the updated files into the branch.


## General Notes

-We run this code under Spyder IDE, because Spyder is a powerful scientific environment written in Python. It's designed by and for scientists, engineers, and data analysts. https://www.spyder-ide.org/

-Running on Visual Studio Code: Download Anaconda Navigator and launch the VS Code via Anaconda. This will you give a conda base interpreter.

-Any major change should be added to the README file.


### Normalization File

-Follow the GitHub steps above and open the file "normalization.py" using Spyder IDE or Visual Studio Code.

-The file currently named "DRX_sorted_norm.csv" contains all the necessary information (data) to read [X being the current data release; i.e. 9].

-The "data_types," "draw_figures," "useful_wavelenght_flux_error_modules," and "utility_functions" files must be in your directory for import

-At the top of the file there is a section containing ranges of wavelengths under the "DO NOT CHANGE" heading. These constant variables are defined by Astrophysicist. Any change of these constants SHOULD be discussed with the client.

-The first variable defined is where you state which Data Release you would like to work with. There are currently files for the DR9 and DR16. To run the code, make sure the DR variable is set to the proper Data Release number (i.e. '16' for DR16). 

-The range you choose for STARTS_FROM, ENDS_AT will depend on the Data Release you are working with. We currently have 6760 files for DR9 and 21859 files for DR16. 
  -For DR9 the current range that works is either 1, 10 or 899, 1527
  -For DR16 the current range that works is 1, 21000 (high redshift cases are currently throwing errors in the code)

-After you run this code, all your graphs will be added to a pdf file in your directory. None of these files contain graphs for spectra with SNR less than 10 
  -There are currently 6 pdf files containing graphs:
    -flagged_absorption_graphs.pdf contains all of the graphs that have been flagged for the green or pink test region being completely under the powerlaw (likely absorption in those regions). [TEST #4]
    -flagged_bad_fit_graphs.pdf contains all of the graphs that have been flagged by the code for a poor fit through the pink and green test regions. [TEST #1 & TEST #2]
    -good_fit_graphs.pdf contains all of the graphs that have been fit well. Includes unflagged cases. ****CHECK THIS****
    -normalized_graphs.pdf contains graphs of all of the normalized spectra [contains good_fit and unflagged]
    -original_graphs.pdf contains graphs of all of the spectra run in the range defined
    -unflagged_graphs.pdf contains graphs of all of the spectra that were previously flagged but have been deemed good fits based on a secondary test. [TEST #3]

-There is a test file in the directory. If you change something wrong, the test file will catch the different results. If you change constant variables, this will cause different output. If changes are correct then update the test file with new results for future tests.

OUTPUT FILES (csv):

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
    
OUTPUT FILES (plots):

-original_graphs.pdf
    -all spectra that have SNR>10 are plotted and added to this file
    
-normalized_graphs.pdf
    -all spefctra that have been deemed good fits are normalized and the normalized spectra are plotted and added to this file
    
-good_fit_graphs.pdf
    -all spectra that have SNR>10 that have been deemed to be a good fit originally are plotted and added to this file
    -all spectra that have been unflagged have been plotted and added to this file

-unflagged_graphs.pdf
    -any spectra that were previously flagged that have been deemed a good fit are unflagged, plotted and added to this file
    
-flagged_bad_fit_graphs.pdf
    -all spectra that have been deemed a bad fit by the tests are flagged, plotted, and added to this file
    
-flagged_absorption_graphs.pdf
    -any spectra that is flagged for test#4 are added to this file due to possible absorption in the test regions

original_graphs.pdf = original.csv
good_fit_graphs.pdf + flagged_bad_fit_graphs.pdf = original_graphs.pdf
good_fit_graphs.pdf = normalized_graphs.pdf

NORMALIZATION TESTS:

-TEST #1
    -checks whether the fit of the powerlaw is going closely through the anchor points
        -if the powerlaw goes through all three anchor points well, it proceeds to the other three tests
        -if the powerlaw is too far off of the anchor points they are flagged and moved to the flagged_bad_fit file

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


### Absorption File

-Follow the GitHub steps above and open the file "absorption.py" using Spyder IDE or Visual Studio Code.

-File currently has DR9 "Data Release 9" extension. In the future, it will be changed to the recent Data Release from SDSS Database
 
-The Absorption file reads the Normalization files. After you generate the normalization file for each spectra, the Absorption file will read that and generate necessary graphs and information. Currently the format of the processed Normalization file is : "spectra name + norm + .dr9"

-Based on "BALNICITY_INDEX_LIMIT", the file will import different results. You can see the results in the same directory under the "Absortion_cleaning" document. Also the same "BALNICITY_INDEX_LIMIT" constant will generate different graph output because it is defining the narrowness of the output.

- At the top of the file, the constant variables are defined by Astrophysicist. Any change of these constants SHOULD discuss with the client.


### .gitignore

-Any files that you would like to be ignored when pushing should be recorded in this file.
** ADD MORE DETAILS **