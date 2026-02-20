#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 09:45:45 2026

@author: lilianaflores
"""
import os
import sys
import numpy as np 
import math
import pandas as pd
#from numpy.lib.function_base import append #Remove: Unused in program and incompatible with updated version of numpy
from matplotlib.backends.backend_pdf import PdfPages
sys.path.insert(0, os.getcwd() + '/../' ) # changes the directory to the DR16Q --> all paths after this will need to be written as if this was in the top level of the DR16Q
from utility_functions import clear_file, read_list_spectra, read_spectra, append_row_to_csv
from data_types import Range
from abs_function_module import smooth, abs_parameters_plot_optional
from abs_plot_module import draw_abs_figure

CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/../DR16Q_EHVO/good_fit_EHVO.csv" #"/OUTPUT_FILES/NORMALIZATION/good_fit_EHVO.csv" #good_fit_EHVO.csv" ##_newSNR_flagged_but_ok.csv #_EHVO.csv" 


cross_match = os.getcwd()+'/EHVO_absorption_depths_corrected.csv'


df_og = pd.read_csv(CONFIG_FILE)

df_match = pd.read_csv(cross_match)
print(df_match)

spec_name_og = df_og['NORM SPECTRA FILE NAME']
spec_name_cross = df_match['NORM SPECTRA FILE NAME']

flag_cross = df_match['NEEDS RECALCULATION']



df_og['DEPTH FLAG'] = np.where(spec_name_og == spec_name_cross, flag_cross, 'no match')

print(df_og)