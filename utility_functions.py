import csv
import numpy as np 
from data_types import ColumnIndexes, RangesData, PointData
import pandas as pd

######################################### sphinx ######################################### 
"""
utility_functions.py
====================
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

def read_file_abs(FILE: str):
    spectra_index, spectra_filename, norm_spectra_filename, redshift_value, calc_snr_value, sdss_snr, bf, cf = [], [], [], [], [], [], [], []

    with open(FILE) as f:  
        for line in f:
            each_row_in_file = line.split(",")
            norm_spectra_filename.append(each_row_in_file[2])
            redshift_value.append(np.float(each_row_in_file[3]))
            calc_snr_value.append(np.float(each_row_in_file[4]))
    return(redshift_value, calc_snr_value, norm_spectra_filename)

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

def pandas_test(path: str, file_name: str):
    good_norm = pd.read_csv(path, file_name,
        dtype = {"REDSHIFT": float, "CALCULATED SNR": float},
        usecols = ['NORM SPECTRA FILE NAME', 'REDSHIFT', 'CALCULATED SNR']
    ) [['NORM SPECTRA FILE NAME', 'REDSHIFT', 'CALCULATED SNR']]
    return good_norm

def read_list_spectra(file_name: str, column_list: list):
    print("column_list: ", column_list)
    data = pd.read_csv(file_name)
    #data = pd.read_csv(file_name,
    #    dtype = {"REDSHIFT": float, "CALCULATED SNR": float},
     #   names = column_list, 
      #  engine = 'python')
    spectra_list = data[column_list[0]]
    redshift_list = data[column_list[1]]
    snr_list = data[column_list[2]]
    return(spectra_list,redshift_list, snr_list)

def read_spectra(spectra_data):
    column_index = ColumnIndexes(0, 1, 2)
    wavelength = spectra_data[:, column_index.wavelength]
    flux = spectra_data[:, column_index.flux] 
    error = spectra_data[:, column_index.error] 
    
    return [wavelength, flux, error]