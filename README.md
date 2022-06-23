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
### DIRECTORIES IN REPOSITORY [see the README in each directory for more detailed info]
- ABSORPTION
    -contains all files/codes necessary to run the code that searches each normalized quasar spectrum for absorption

- CROSS_CORRELATION
    -?? 

- DATA
    -contains all .drX files; these are the data for each individual spectrum that is downloaded from the SDSS

- DR9Q_EHVO
    -contains all files relevant to the 40 DR9 EHVO cases

- DR16Q_EHVO
    -contains all files relevant to the 98 DR16 EHVO cases (& possibly the candidates later on??? TBD)

- NORMALIZATION
    -contains all files/codes necessary to run the code that fits & normalizes each quasar spectrum
    -this code also contributes to creating our parent sample (flags spectra for SNR, etc.)

- OLD
    -old files/versions of code that we are scared to delete :)

- PRESENTATION_PLOTS
    -codes to make pretty plots for posters/presentations/papers/fun wallpaper/etc. 

- REDSHIFT
    -codes/files to create histogram plots to compare redshifts of different populations for both the DR9 and DR16

- VARIABILITY
    -??


------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### FILES IN DIRECTORY
-data_types.py
    This file contains definitions of tuples that are used in the normalization.py code.  
-DR9_sorted_norm.csv
    File that contains all necessary information/data to run normalization.py. This is NOT the parent sample - all spectra we have downloaded from the SDSS are included in this file. 
    **Columns:** Spectra Name, Redshift, SDSS SNR
-DR16_BAL_parent_sample.csv
    File that contains the DR16 BALQSOs --> these are the BALQSOs that meet our criteria AND are in our parent sample
    -**Columns:** Spectra File Name, Redshift, BI Value, BI Error, SNR
-DR16_parent_sample.csv
    File that contains the final DR16 parent sample
-DR16_sorted_norm.csv
    File that contains all necessary information/data to run normalization.py. This is NOT the parent sample - all spectra we have downloaded from the SDSS are included in this file. 
    **Columns:** Spectra Name, Redshift, SDSS SNR
-draw_figures.py
    This code contains the functions used by the normalization.py code to plot spectra. 
-useful_wavelength_flux_error_modules.py
    Contains the functions used to calculate the location of the anchor points, calculate SNR, etc.
-utility_functions.py
    Contains the functions to read/open/print to files and other general functions. 

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### .gitignore

-Any files that you would like to be ignored when pushing should be recorded in this file.
** ADD MORE DETAILS **