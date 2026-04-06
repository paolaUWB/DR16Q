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

# Create edited names from CSV
EDITED_NAME_1 = [name[14:] for name in SPEC_NAME_1]

#for i in range(np.size(SPEC_NAME_1)):
#    EDITED_NAME_1.append(SPEC_NAME_1[i][14:])
    

EDITED_NAME_2 = []

for f in match2.iterdir():
    if f.is_file():
        EDITED_NAME_2.append((f.name[13:]))

#print(EDITED_NAME_2)
'''

        
import os
from pathlib import Path
import pandas as pd
import sys

csv1 = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/spec_lists/xray_list_SNR10_z1.9.csv"
match2 = Path(os.getcwd() + '/SDSS_fullspecs')

file = pd.read_csv(csv1)
SPEC_NAME_1 = file['SPECTRA FILE NAME']

# --- Extraction functions ---
def extract_sdss_key(filename):
    # spec-015005-59193-4399272889.fits
    parts = filename.replace(".fits", "").split("-")
    return f"{parts[-2]}-{parts[-1]}"

def extract_csv_key(name):
    # spec-allepoch-59192-4382235955.fits
    parts = name.replace(".fits", "").split("-")
    return f"{parts[-2]}-{parts[-1]}"

# --- Build SDSS lookup ---
sdss_dict = {}

for f in match2.iterdir():
    if f.is_file():
        key = extract_sdss_key(f.name)
        sdss_dict[key] = f.name

# --- Match in order ---
ordered_sdss_names = []

for name in SPEC_NAME_1:
    key = extract_csv_key(name)
    
    if key in sdss_dict:
        ordered_sdss_names.append(sdss_dict[key])
    else:
        ordered_sdss_names.append(None)

# --- Export ---
output_df = pd.DataFrame({'ordered_sdss_names': ordered_sdss_names})
output_df.to_csv('ordered_sdss_names.csv', index=False)

print("Done.")

