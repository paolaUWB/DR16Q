# DR16Q --> REDSHIFT

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

This directory contains files that plot histograms to analyze the redshifts of the samples. Four histograms in total can be made from this code, one being a cumulative histogram. These histograms show a comparison of the redshifts in the parent sample, EHVO sample, and BAL sample for each DR. 

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### RUNNING REDSHIFT FILE [DRX_redshift.py]

NOTE: these were separated to account for the different types of files being read (to keep cleaner) - the DR9 provided .dat files and the DR16 was .csv files. The DR9 file has the BAL definition within it, the DR16 requires one of the input files to be a list of the BALs.

For both the DR9 and DR16 Redshift programs, the number of spectra in the parent sample, EHVO sample, and BAL sample needs to be known and input at the top in the corresponding variables. 

-Files needed to run DR9_redshift.py: 
    -DR9Q_selection_minus17.dat
    -EHVO_DR9.dat

-Files needed to run DR16_redshift.py:
    -DR16_parent_sample.csv
    -DR16_EHVO_sorted_norm.csv
    -DR16_BAL_parent_sample

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### FILES IN REDSHIFT DIRECTORY
- DR9_redshift.py
    -Code to plot histograms for the DR9 quasar redshift

- DR16_redshift.py
    -Code to plot histograms for the DR16 quasar redshift

- draw_histogram.py
    -Function that plots the histograms for the DRX_redshift.py codes

- OUTPUT_FILES --> contains all plots created by the redshift codes
    -zem_DR9_cdf.png
        -This file contains the cumulative histogram plot for the parent, EHVOs & BALs in the DR9
    -zem_DR9_DR16_cdf.png
        -This file contains the cumulative histogram plot for the parent, EHVOs, & BALs in the DR9 & DR16 to compare
    -zem_DR9_zem-class.png
        -This file contains the histogram plots for the parent, EHVOs, & BALs, each divided by the total number of spectra in the corresponding class for the DR9
    -zem_DR9_zem-parent.png
        -This file contains the histogram plots for the parent, EHVOs, & BALs, each divided by the total number of spectra in the parent sample per bin for the DR9
    -zem_DR9_zem.png
        -This file contains the histogram plots for the parent, EHVOs, & BALs for the DR9 (with the option to scale the EHVOs and BALs so they appear on the plot better)
    -zem_DR16_cdf.png
        -This file contains the cumulative histogram plot for the parent, EHVOs & BALs in the DR16
    -zem_DR16_zem-class.png
        -This file contains the histogram plots for the parent, EHVOs, & BALs, each divided by the total number of spectra in the corresponding class for the DR16
    -zem_DR16_zem-parent.png
        -This file contains the histogram plots for the parent, EHVOs, & BALs, each divided by the total number of spectra in the parent sample per bin for the DR16
    -zem_DR16_zem.png
        -This file contains the histogram plots for the parent, EHVOs, & BALs for the DR16 (with the option to scale the EHVOs and BALs so they appear on the plot better)