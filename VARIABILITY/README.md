# EHVO_Variability
Work on EHVO Variability for DR16Q EHVO quasars.
---
make_directories.py
    This file will make directories for EHVO spectra that have duplicate observations. The directories are created for each object, and will eventually contain all duplicate spectra of that object (they are named by the SDSS name).
    To run this code, you will need the DR16Q_v4.fits file and the DRX_EHVO_sorted_norm.csv file.
    The code will output a file "data_request.csv" containing all necessary information for each duplicate observation of the spectra to request the data from Pat. 
    It will also output a file "ehvo_duplicate_file.csv" containing all the EHVO spectra and their corresponding duplicates to assist in sorting the spectra files within the directories. 
    

make_info_files.py
    This file will make info files for all directories within another given directory.

overplot.py
    This file will go through a directory of directories and plot the spectra
    in each directory.

---
Data File Organization
    In the "DATA_VARIABILITY" folder, there are other folders with object IDs 
    for the names. These folders will contain the spectral data for all of the
    different epochs for each object, as well as a file called "info.txt",
    which will have other information about the object, such as its redshift.

MJD limits:
between SDSSI/II and SDSS III/IV is 54663.
BOSS DR9 started MJD 55171 2009 December 5 and ended mid July 2011 55752
SDSS IV > 55752

