# DR16Q --> NORMALIZATION

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

The normalization code works by first calculating the SNR of each spectrum between 1250-1400 A to determine if any spectra have a lot of noise in the region we are searching for CIV absorption to identify EHVO quasars. Any cases with SNR<10 are flagged and ignored for the remainder of the project. Cases with SNR>10 are fit with a power law using three anchor points in regions of the spectrum that typically lack emission and absorption to best attempt to fit the spectrum through the continuum. A series of tests are performed (described below) to check that the fit goes through the continuum. Any spectrum that is deemed a "good fit" (fit through continuum) is normalized by dividing the flux by the power law and is saved in a file to be run through the absorption.py code to search for absorption. 

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### FILES IN NORMALIZATION DIRECTORY
- normalization.py
    -This is the code that normalizes quasar spectra to analyze and helps to define our parent sample (flags low SNR cases). This is the first step in the project and creates the files that the absorption code uses. 

- OUTPUT_FILES --> pdfFILES
    After you run this code, all your graphs will be added to a pdf file in your directory. None of these files contain graphs for spectra with SNR less than 10 (other than flagged_snr_graphs.pdf)

    - anchor_pt_graphs.pdf
        -saved only when dynamic='yes'
        -these plots contain the full wavelength plot of the spectrum as well as three plots that are zoomed in on the anchor point locations to assist in plotting points
        -the user input anchor points are on these plots which can be nice to reference if you need to redo or recreate the normalization

    - flagged_absorption_graphs.pdf
        -any spectra that is flagged for TEST #4 are added to this file due to possible absorption in the test regions [TEST #4]

    - flagged_bad_fit_graphs.pdf
        -all spectra that have been deemed a bad fit by the tests are flagged, plotted, and added to this file [TEST #1 & TEST #2] 

    - flagged_snr_graphs.pdf
        -all spectra that have been flagged for SNR<10 have been added to this file --> this file has little purpose, it was used for symposium presentations to find examples of noisy vs. high SNR spectra

    - good_fit_graphs.pdf
        -all spectra that have SNR>10 that have been deemed to be a good fit originally are plotted and added to this file
        -all spectra that have been unflagged have been plotted and added to this file

    - normalized_graphs.pdf
        -all spefctra that have been deemed good fits are normalized and the normalized spectra are plotted and added to this file [contains good_fit and unflagged]

    - original_graphs.pdf
        -all spectra in the range provided that have SNR>10 are plotted and added to this file

    - unflagged_graphs.pdf
        -any spectra that were previously flagged that have been deemed a good fit are unflagged, plotted and added to this file [TEST #3] 
    
    --> original_graphs.pdf ~ original.csv
    --> good_fit_graphs.pdf ~ flagged_bad_fit_graphs.pdf ~ original_graphs.pdf
    --> good_fit_graphs.pdf ~ normalized_graphs.pdf

- OUTPUT_FILES --> textFILES
    - anchor_pts.csv
        -saves file with user input anchor points when dynamic='yes'

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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