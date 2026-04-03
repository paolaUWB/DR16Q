# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 15:57:44 2026

@author: elijahf
"""
'''
loop through x-ray_list_snr10_z1.9
loop through folder containing full sdss spectra files
during loop, chop off the first sections to get the last 2 sets of numbers
match the names and export into csv with the original names of sdss_full spec but in order

Look through looping through files in a folder
csv to dataframe
'''

import os
from pathlib import Path
import pandas as pd
import sys
import numpy as np


csv1 = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()+"/spec_lists/xray_list_SNR10_z1.9.csv"
match2 = Path(os.getcwd() + '/SDSS_fullspecs')

file = pd.read_csv(csv1)


SPEC_NAME_1 = file['SPECTRA FILE NAME']
EDITED_NAME_1 = []

for i in range(np.size(SPEC_NAME_1)):
    EDITED_NAME_1.append(SPEC_NAME_1[i][14:])
    

EDITED_NAME_2 = []

for f in match2.iterdir():
    if f.is_file():
        EDITED_NAME_2.append((f.name[13:]))

print(EDITED_NAME_2)


        
    

