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
    
-There are 6760 spectra (It can be more later on this project). If you don't generate all 6760 files, find variables "STARTS_FROM" and "ENDS_AT" and define a range.

-The "data_types" and "utility_functions" files must be in your directory for import

-The file currently named "DRX_sorted_norm.csv" contains all the necessary information (data) to read [X being the current data release; i.e. 9].

-At the top of the file there is a section containing ranges of wavelengths under the "DO NOT CHANGE" heading. These constant variables are defined by Astrophysicist. Any change of these constants SHOULD be discussed with the client.

-The first variable defined is where you state which Data Release you would like to work with. There are currently files for the DR9 and DR16. To run the code, make sure the DR variable is set to the proper Data Release number (i.e. '16' for DR16). 

-The range you choose for STARTS_FROM, ENDS_AT will depend on the Data Release you are working with. We currently have 6760 files for DR9 and 21859 files for DR16. 
  -For DR9 the current range that works is either 1, 10 or 899, 1527
  -For DR16 the current range that works is 1, 21000 (high redshift cases are currently throwing errors in the code)

-After you run this code, all your graphs will be added to a pdf file in your directory. None of these files contain graphs for spectra with SNR less than 10 
  -There are currently 4 pdf files containing graphs:
    -flagged_spectra.pdf contains all of the graphs that have been flagged by the code for a poor fit through the pink and green test regions
    -normalized_graphs.pdf contains graphs of all of the normalized spectra
    -original_graphs.pdf contains graphs of all of the spectra run in the range defined
    -powerlaw_test_graphs.pdf contains graphs of all of the spectra that were previously flagged but have been deemed good fits based on a secondary test

-There is a test file in the directory. If you change something wrong, the test file will catch the different results. If you change constant variables, this will cause different output. If changes are correct then update the test file with new results for future tests.

OUTPUT FILES:
-bad_normalization.csv
    -spectra index, chi_sq
    -contains spectra that has been flagged as a bad fit by the green and pink test regions
    -these spectra do not move on to be normalized

-chi_sq_values.csv
    -spectra index, chi_sq
    -contains the chi squared values for all spectra run in the range defined

-final_initial_parameters.txt
    -spectra index, spectra file name, bf, cf
    -contains the initial parameters for each spectra that are used for the powerlaw curve fitting

-flagged_absorption.csv
    -spectra file name
    -contains the spectra that have been flagged as having possible absorption in the green and pink test regions

-flagged_bad_fit.csv
    -spectra index, chi_sq
    -contains spectra that have been flagged as a bad fit by the green and pink test regions

-flagged_snr_in_ehvo_graphs.txt
    -spectra index, spectra file name, SNR
    -contains spectra that have SNR less than 10
    -these spectra are not used in the code
    ** sorted by SNR **

-good_normalization.csv
    -spectra index, chi_sq
    -contains spectra that have a good fit and will be normalized

-powerlaw_test2.txt ** DELETE? **

-processed_spectra_file_names.txt
    -contains the spectra file names for all spectra in the defined range


### Absorption File

-Follow the GitHub steps above and open the file "absorption.py" using Spyder IDE or Visual Studio Code.

-File currently has DR9 "Data Release 9" extension. In the future, it will be changed to the recent Data Release from SDSS Database
 
-The Absorption file reads the Normalization files. After you generate the normalization file for each spectra, the Absorption file will read that and generate necessary graphs and information. Currently the format of the processed Normalization file is : "spectra name + norm + .dr9"

-Based on "BALNICITY_INDEX_LIMIT", the file will import different results. You can see the results in the same directory under the "Absortion_cleaning" document. Also the same "BALNICITY_INDEX_LIMIT" constant will generate different graph output because it is defining the narrowness of the output.

- At the top of the file, the constant variables are defined by Astrophysicist. Any change of these constants SHOULD discuss with the client.