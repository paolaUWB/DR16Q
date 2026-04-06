# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 23:05:15 2026

@author: elijahf
"""
import os
import sys

csv1 = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/spec_lists/xray_list_SNR10_z1.9.csv"
match2 = os.getcwd() + '/SDSS_fullspecs'

from CSV_NAME_MATCH import match_sources 

match_sources(
    source1 = csv1,
    source2 = match2,
    sections1 = [-2, -1],
    sections2 = [-2, -1],
    column_name1="SPECTRA FILE NAME",
    output_csv="test_match_name_match.csv",
    debug=True
)