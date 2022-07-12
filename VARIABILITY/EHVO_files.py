import os
import sys
import shutil
import numpy as np
from data_types import ColumnIndexes
sys.path.insert(0, os.getcwd() + '/../' + 'DR16Q') # changes the directory to the DR16Q --> all paths after this will need to be written as if this was in the top level of the DR16Q
from utility_functions import read_file, read_spectra, append_row_to_csv, print_to_file


DR16_SPEC_DIREC = os.getcwd() + '/DATA/DR16Q_SNR10/'
DR16_EHVO_DIREC = os.getcwd() + '/DR16Q_EHVO/DR16Q_EHVO_FILES/'

DR9_EHVO_DIREC = os.getcwd() + '/DR9Q_EHVO/NORM_DR9Q_EHVO/'
# DR9_EHVO_TXT_FILE = os.getcwd() + '/DR9Q_EHVO/EHVO_DR9.dat'
DR9_EHVO_TXT_FILE = np.loadtxt(os.getcwd() + '/DR9Q_EHVO/EHVO_DR9.dat',dtype=bytes,delimiter="\n").astype(str)
DR9_EHVO_SORTED_NORM = os.getcwd() + '/DR9Q_EHVO/DR9_EHVO_sorted_norm.csv'

EHVO_SAMPLE_DR9 = 40

PLATE_EHVO_DR9 = []
MJD_EHVO_DR9 = []
FIBERID_EHVO_DR9 = []

###### for use in dr9_redshift.py code to create dr9_sorted_norm for ehvos
# from utility_functions import append_row_to_csv
# DR9_EHVO_SORTED_NORM = os.getcwd() + '/DR9Q_EHVO/DR9_EHVO_sorted_norm.csv'

# for r in range(EHVO_SAMPLE_DR9):
#     DR9_spec_name = ('spec-' + str(PLATE_EHVO_DR9[r]).zfill(4) + '-' + str(MJD_EHVO_DR9[r]).zfill(4) + '-' + str(FIBERID_EHVO_DR9[r]).zfill(4) + '.dr9')
#     fields = [DR9_spec_name, zem_EHVO_DR9[r], '-1', str(PLATE_EHVO_DR9[r]).zfill(4), str(MJD_EHVO_DR9[r]).zfill(4), str(FIBERID_EHVO_DR9[r]).zfill(4)]
#     append_row_to_csv(DR9_EHVO_SORTED_NORM, fields)
######

# DR9_spec_name = []

# for m in range(0, EHVO_SAMPLE_DR9):
#     DR9_EHVO_row = DR9_EHVO_TXT_FILE[m]
#     DR9_EHVO_columns = DR9_EHVO_row.split()
    
#     PLATE_EHVO_DR9.append(str(DR9_EHVO_columns[1]))
#     MJD_EHVO_DR9.append(str(DR9_EHVO_columns[2]))
#     FIBERID_EHVO_DR9.append(str(DR9_EHVO_columns[3]))
    
#     DR9_spec_name = ('spec-' + PLATE_EHVO_DR9[m] + '-' + MJD_EHVO_DR9[m] + '-' + FIBERID_EHVO_DR9[m] + '.dr9')
    
#     fields = [DR9_spec_name, PLATE_EHVO_DR9[m], MJD_EHVO_DR9[m], FIBERID_EHVO_DR9[m]]
#     append_row_to_csv(DR9_EHVO_SORTED_NORM, fields)
    # print_to_file(DR9_spec_name, DR9_EHVO_SORTED_NORM)

# DR9_SPEC_DIREC = os.getcwd() + '/DATA/DR9Q_SNR10/'
# DR9_EHVO_DIREC = os.getcwd() + '/DR9Q_EHVO/DR9Q_EHVO_FILES/'

# FITS_DUPLICATES_FILE = os.getcwd() + '/fits_duplicates_no-1.txt'


# info_DR16_EHVO = os.getcwd() + '/DR16Q_EHVO/DR16_EHVO_sorted_norm.csv'

# zem_DR16, snr_DR16, spectra_list_DR16 = read_file(info_DR16_EHVO) # parent

# for i in range(len(spectra_list_DR16)):
#     shutil.copyfile(DR16_SPEC_DIREC + spectra_list_DR16[i], DR16_EHVO_DIREC + spectra_list_DR16[i])



# for i in range(len(spectra_list_DR16)): # i starts at 0
#     EHVO_DR16_minus_dered = spectra_list_DR16[i][:-11]) # -11 for dr16, -10 for dr9
#     if 
#         if EHVO_DR16_minus_dered == 
