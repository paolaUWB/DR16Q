import csv
import numpy as np 
from data_types import ColumnIndexes, RangesData, PointData

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

