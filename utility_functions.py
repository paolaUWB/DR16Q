import csv
import numpy as np 
from data_types import ColumnIndexes, RangesData, PointData
import pandas as pd
import scipy.constants as sc

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
            redshift_value_list.append(np.float(each_row_in_file[1]))
            snr_value_list.append(np.float(each_row_in_file[2]))
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
    spectra_list = data[column_list[0]]
    redshift_list = data[column_list[1]]
    snr_list = data[column_list[2]]
    return(spectra_list,redshift_list, snr_list)

def read_spectra(spectra_data):
    """Reads in and returns a lists of lists containing the wavelength, flux, and error for each spectra.

    Defines the variables to be used in the code.

    Parameters
    ----------
    spectra_data:
        The drX (X being 9 or 16) files in the form of a text file conataining the wavelength, flux
        and error of that paticular spectra.

    Returns
    -------
    wavelength: list
        All of the wavelength values in a list.
    flux: list
        All of the flux values in a list.
    error: list
        All of the error values in a list.

    Note
    ----
    As shown in the return, the return value are lists within a list.
    """

    column_index = ColumnIndexes(0, 1, 2)
    wavelength = spectra_data[:, column_index.wavelength]
    flux = spectra_data[:, column_index.flux] 
    error = spectra_data[:, column_index.error] 
    
    return [wavelength, flux, error]

def wavelength_to_velocity(redshift, wavelength):
    """Reads in a list of wavelength values to be converted to velocity.

    Parameters
    ----------
    redshift: list
        The list of redshift values needed for the equation.
    wavelength: list
        The list of wavelength values that will be converted to velocity.
        
    Returns
    -------
    beta: array
        The values of velocity that were converted from wavelength.
    """
    # CIV doublet data from verner table
    avr_CIV_doublet = 1549.0524

    # Transform the wavelength array to velocity (called "beta") based on the CIV doublet: 
    c_in_km = sc.speed_of_light**-3
    z_absC = (wavelength / avr_CIV_doublet) - 1.
    RC = (1. + redshift) / (1. + z_absC)
    betaC = ((RC**2.) - 1.) / ((RC**2.) + 1.) # betaC is in units of c (speed of light)
    betakm = -betaC * c_in_km #betakm is in km/s
    beta = []

    for velocity in betakm:
        betas = round(velocity, 5)
        beta.append(betas)
    beta = np.array(beta)

    return beta