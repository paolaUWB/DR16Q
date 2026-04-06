# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 15:57:44 2026

@author: elijahf
"""
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
# Splits makes it so that it will create a list with each entry being what is between the given parameter
# We then take the last 2 sections and combine them together
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
'''


def match_sources(
    source1,
    source2,
    sections1,
    sections2,
    column_name1=None,
    column_name2=None,
    delimiter="-",
    remove_ext=True,
    debug=False,
    output_csv=None
):
    """
    Match filenames between two sources (CSV or folder) using selected sections.

    Parameters:
        source1 (str): Path to CSV file or folder (primary order is based on this)
        source2 (str): Path to CSV file or folder (lookup source)
        sections1 (list): Indices of parts to extract from source1 filenames eg: [-2, -1] for last two sections
        sections2 (list): Indices of parts to extract from source2 filenames
        column_name1 (str): Column name for source1 if CSV
        column_name2 (str): Column name for source2 if CSV
        delimiter (str): Character to split filenames (default "-")
        remove_ext (bool): Remove file extension before splitting
        debug (bool): Print missing matches if True
        output_csv (str): If provided, saves results to this CSV filename

    Returns:
        list: Matched filenames from source2 in order of source1
    """

    import pandas as pd
    from pathlib import Path

    # --- extract key ---
    def extract_key(name, sections):
        if not isinstance(name, str):
            return None
        #remove whitespaces
        name = name.strip()
        #removes the file extension
        if remove_ext:
            name = name.rsplit(".", 1)[0]
        
        #splits into each subsequent piece
        parts = name.split(delimiter)
        
        #takes our sections and then grabs the ones indicated
        try:
            selected = [parts[i] for i in sections]
        except IndexError:
            return None
        #joins a final key of each selected sections
        return delimiter.join(selected)

    # --- Load source1 ---
    # checks if csv
    if source1.endswith(".csv"):
        df1 = pd.read_csv(source1)
        if column_name1 is None:
            raise ValueError("column_name1 must be provided for CSV source1")
        list1 = df1[column_name1].dropna().astype(str).tolist()
    # read as folder
    else:
        list1 = [f.name for f in Path(source1).iterdir() if f.is_file()]

    # --- Load source2 ---
    # checks if csv
    if source2.endswith(".csv"):
        df2 = pd.read_csv(source2)
        if column_name2 is None:
            raise ValueError("column_name2 must be provided for CSV source2")
        list2 = df2[column_name2].dropna().astype(str).tolist()
    # read as folder
    else:
        list2 = [f.name for f in Path(source2).iterdir() if f.is_file()]

    # --- Build lookup dictionary ---
    # creates a quick key lookup with the joined selected files and its full name
    # makes it so we don't need to loop
    lookup = {}
    for name in list2:
        key = extract_key(name, sections2)
        if key:
            lookup[key] = name

    # --- Match in order of source1 ---
    results = []
    # takes our dictionary created for lookup and makes a key that matches the named selection from source1
    # tells you where error is otherwise
    for name in list1:
        key = extract_key(name, sections1)

        if key in lookup:
            results.append(lookup[key])
        else:
            if debug:
                print(f"NO MATCH: {key} (from {name})")
            results.append(None)

    # --- Optional export ---
    # create csv if we have a given name
    if output_csv:
        pd.DataFrame({'matched_names': results}).to_csv(output_csv, index=False)
        if debug:
            print(f"Saved results to {output_csv}")

    return results