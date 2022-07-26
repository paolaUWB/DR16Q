# EHVO_Variability
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Work on EHVO Variability for DR16Q EHVO quasars.
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
`make_directories.py`
- This file will attempt to make directories in the "DATA_VARIABILITY" folder from any spectra files that are in the main folder of the repository.

`make_info_files.py`
- This file will make info files for all directories within another given directory.

`overplot.py`
- This file will go through a directory of directories and plot the spectra in each directory.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Data File Organization

In the "DATA_VARIABILITY" folder, there are other folders with object IDs for the names. These folders will contain the spectral data for all of thedifferent epochs for each object, as well as a file called `info.txt`,which will have other information about the object, such as its redshift.

- MJD limits:
    - between SDSSI/II and SDSS III/IV is 54663.
    - BOSS DR9 started MJD 55171 2009 December 5 and ended mid July 2011 55752
    - SDSS IV > 55752

