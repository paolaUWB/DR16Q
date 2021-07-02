import numpy as np 
from data_types import ColumnIndexes, RangesData, PointData

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