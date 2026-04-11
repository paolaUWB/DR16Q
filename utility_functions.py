import csv
import numpy as np 
from data_types import ColumnIndexes, RangesData, PointData
import pandas as pd
import scipy.constants as sc
from scipy import signal
from astropy.io import fits

######################################### sphinx ######################################### 
"""
utility_functions
=================
Utility functions for this project.
"""
#############################################################################################

def read_file(FILE: str):
    spectra_list, redshift_value_list, snr_value_list = [], [], []
    
    with open(FILE) as f:  
        for line in f:
            each_row_in_file = line.split(",")
            spectra_list.append(each_row_in_file[0])
            redshift_value_list.append(float(each_row_in_file[1]))
            snr_value_list.append(float(each_row_in_file[2]))


    return(redshift_value_list, snr_value_list, spectra_list)

def print_to_file(text: str, file_name: str):
    print(text, file = open(file_name, 'a'))

def file_length(file_name):
    with open(file_name) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def clear_file(file_name: str):
    open(file_name, 'w').close()

def open_file(file_name: str):
    return open(file_name, 'r')

def append_row_to_csv(file_name: str, fields: list):
    with open(file_name, 'a') as f:
        writer = csv.writer(f)
        writer.writerow(fields)

def read_list_spectra(file_name: str, column_list: list):
    """Reads in a csv file of spectra to allow you to access certain columns of data from the csv file.

    Parameters
    ----------
    file_name: str
        Enter in the name of your csv file as a string. Can also enter the path of where the file
        is as long as the name of the file is included in the pathway.
    column_list: list
        Enter in the names of the columns you want to access from your csv file in the form as
        a list of strings.

    Returns
    -------
    spectra_list: list
        Whatever information is in the spectra column as a list.
    redshift_list: list
        Whatever information is in the redshift column as a list.
    snr_list: list
        Whatever infromation is in the snr column as a list.

    Example
    -------
    >>> CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/OUTPUT_FILES/NORMALIZATION/good_normalization.csv"
    >>> norm_spectra_list, redshift_list, calc_snr_list = read_list_spectra(CONFIG_FILE, ["NORM SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR"])
    [spec-9140-58039-0081norm.dr16  1.9  14.1
    spec-7671-57360-0092norm.dr16   1.9  14.2]

    Notes
    -----
    ``good_normalization.csv`` in this case is a csv file with several headers but we only wanted to select these three.
    Contents of ``good_normalization.csv`` for purpose of example:

    >>> SPECTRA FILE NAME,NORM SPECTRA FILE NAME,REDSHIFT,CALCULATED SNR,SDSS SNR,BF,CF
    >>> spec-9140-58039-0081-dered.dr16,spec-9140-58039-0081norm.dr16,1.9,14.1,31.,9.1,0.1
    >>> spec-7671-57360-0092-dered.dr16,spec-7671-57360-0092norm.dr16,1.9,14.2,20.,9.1,-2.01

    ...

    See Also
    --------
    Pandas is being utilized to read in the csv file, ``pd.read_csv()`` has many different keyword parameters that can be utilized.
    Check pandas api for more details.
    """

    data = pd.read_csv(file_name)

    variable_lists = []
    for i in range(len(column_list)):
        x = data[column_list[i]]
        variable_lists.append(x)

    return variable_lists

def read_spectra(spectra_data, is_fits=False, hdu_index=1):
    """Reads and returns wavelength, flux, and error arrays.
    

 Parameters
 ----------
 spectra_data:
 Either:
 - array-like data (text-loaded spectra), OR
 - path to a FITS file (if is_fits=True)
 is_fits : bool, optional
 If True, reads spectra_data as a FITS file.
 hdu_index : int, optional
 HDU index to read from FITS file (default is 1).

 Returns
 -------
 wavelength: array
 flux: array
 error: array
 """

    column_index = ColumnIndexes(0, 1, 2)

    if is_fits:
        # Open FITS fits
            
        
        with fits.open(spectra_data) as hdul:
            data = hdul[hdu_index].data


            # Adjust column names depending on your FITS structure
            wavelength = data['WAVE']
            flux = data['FLUX']
            error = data['ERROR']
        
        
    else:
        # Assume spectra_data is already a numpy array
        wavelength = spectra_data[:, column_index.wavelength]
        flux = spectra_data[:, column_index.flux]
        error = spectra_data[:, column_index.error]

    return [wavelength, flux, error]