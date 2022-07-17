# DR16Q --> DR16Q_EHVO

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

All information & files for the DR16Q EHVOs.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### FILES IN DR16Q_EHVO DIRECTORY
- NORM_DR16Q_EHVO
    -Contains the normalized files for the EHVOs in the DR16Q

- pdfFILES_abs
    - EHVO_BI2000_smooth11.pdf
        -This file contains all EHVO spectra after they have been run through the absorption code that searches for absorption troughs that dip at least 10% below the continuum for at least 2000km/s when the spectra is smoothed (boxcar=11).  

    - EHVO_BI2000.pdf
        -This file contains all EHVO spectra after they have been run through the absorption code that searches for absorption troughs that dip at least 10% below the continuum for at least 2000km/s when the spectra is not smoothed.

- pdfFILES_norm
    - anchor_pt_graphs_EHVO.pdf
        -This file contains the original plots of all of the EHVO spectra.
        -These have all been manually normalized and includes the full wavelength range with the power law fit and 3 anchor points. 
        -It also contains three zoomed in plots of where the anchor points were chosen to be and has the user input wavelength for each anchor point. 

    - good_fit_graphs_EHVO.pdf
        -This file contains plots of all of the EHVO spectra after they have been fit with a power law, shows the range between 1200-1800A in rest frame. 

    - normalized_graphs_EHVO.pdf
        -This file contains the plots of all EHVO spectra after they have been normalized. 

- DR16_EHVO_sorted_norm.csv
    -This file contains the list of EHVO spectra that is run through the normalization code. 
    -**Columns:** Spectra File Name, Redshift, SNR, Plate, MJD, Fiber

- EHVO_BI2000_smooth11.txt
    -This file contains the outputs for the EHVO spectra after they have been run through the absorption code (for smoothed spectra, boxcar=11).
    -Contains count of spectra that have been flagged for having absorption (and total spectra that have been run through code), the norm spectra file name, the calculated Balnicity Index (BI), vmin and vmax of the CIV trough, individual BI measurement, Equivalent Width (EW) of the CIV trough, and Depth of the CIV trough. 
        -NOTE: normally we cannot immediately assume that the absorption in the red column is CIV, but these have all been visually inspected to confirm this is true in this case. 

- EHVO_BI2000.txt
    -This file contains the outputs for the EHVO spectra after they have been run through the absorption code (for non smoothed spectra).
    -Contains count of spectra that have been flagged for having absorption (and total spectra that have been run through code), the norm spectra file name, the calculated Balnicity Index (BI), vmin and vmax of the CIV trough, individual BI measurement, Equivalent Width (EW) of the CIV trough, and Depth of the CIV trough. 
        -NOTE: normally we cannot immediately assume that the absorption in the red column is CIV, but these have all been visually inspected to confirm this is true in this case. 

- good_fit_EHVO.csv
    -This file contains the outputs for the EHVO spectra after they have been run through the normalization code. This is the file that is used to run the absorption code. 
    -**Columns:** Spectra Index, Spectra File Name, Norm Spectra File Name, Redshift, Calculated SNR, SDSS SNR, BF, CF
