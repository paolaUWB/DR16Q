import os
import sys
import shutil
import numpy as np
sys.path.insert(0, os.getcwd() + '/../' + 'DR16Q') # changes the directory to the DR16Q --> all paths after this will need to be written as if this was in the top level of the DR16Q
from utility_functions import read_file

DR16_SPEC_DIREC = os.getcwd() + '/DATA/DR16Q_SNR10/'
DR16_EHVO_DIREC = os.getcwd() + '/DR16Q_EHVO/DR16Q_EHVO_FILES/'

# DR9_SPEC_DIREC = os.getcwd() + '/DATA/DR9Q_SNR10/'
# DR9_EHVO_DIREC = os.getcwd() + '/DR9Q_EHVO/DR9Q_EHVO_FILES/'

info_DR16_EHVO = os.getcwd() + '/DR16Q_EHVO/DR16_EHVO_sorted_norm.csv'

zem_DR16, snr_DR16, spectra_list_DR16 = read_file(info_DR16_EHVO) # parent

for i in range(len(spectra_list_DR16)):
    shutil.copyfile(DR16_SPEC_DIREC + spectra_list_DR16[i], DR16_EHVO_DIREC + spectra_list_DR16[i])


