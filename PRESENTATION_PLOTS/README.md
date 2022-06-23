# DR16Q --> PRESENTATION_PLOTS

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Tools/codes to make pretty plots for posters/presentations/papers/etc.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### RUNNING [plotindividualspectrum.py] FILE

-Files needed to run plotindividualspectrum.py: 
    -norm.drX file

1. Make sure path to normalized files is correct (specdirec)

2. Set the pp2 variable to be the name of the plot that is created & saved

3. Put the name of the file you want to plot as the spectra variable 

4. Fill in the following variables with the correct information for the spectra you would like to plot: 
    -zem=z: redshift of the spectrum

    -topy: ylim of the plot (this is sometimes a trial & error to determine)
    -topemlabel: this places the ion labels, typically 0.03 less than the value for topy works well

    -zem_label_x: x coordinate for the redshift label on plot (in rest frame)
    -zem_label_y: y coordinate for the redshift label on plot (in rest frame)

    -wavelength_emit1_initial: xlim for the left side of the plot (in rest frame) - somewhere between 850-1000 usually works well
    -wavelength_emit2_initial: xlim for the right side of the plot (in rest frame) - 1600 works well to include CIV emission

    -vmin: the minimum velocity of the CIV trough (get from BIXXXX.txt file that outputs from the absorption code - this is just a starting point, adjust as necessary)
    -vmax: the maximum velocity of the CIV trough (get from BIXXXX.txt file that outputs from the absorption code - this is just a starting point, adjust as necessary)

    -NVabs: 'yes' if you want to plot the shaded region for NV absorption
    -OVIabs: 'yes' if you want to plot the shaded region for OVI absorption
    -SiIVabs: 'yes' if you want to plot the shaded region for SiIV absorption
    -OVIem: 'yes' if you want to plot the label for OVI emission
    -Lyaabs: 'yes' if you want to plot the shaded region for Lya absorption

5. Run code and observe your beautiful plot! (and adjust where necessary)

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### FILES IN PRESENTATION_PLOTS DIRECTORY