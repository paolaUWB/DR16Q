# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 15:57:44 2026

@author: elijahf
"""

"""
TO DO:

- Make column header changable parameter 

"""
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


'''
# This is an example of how it can be used. In this case we are matching file names that have different structures but contain the same 
two number sequences. We need a csv with the new file name structure in a folder 'match2' to be in a csv in the order seen in 'csv1'.:

# imports
import os
import sys
from CSV_NAME_MATCH import match_sources 


csv1 = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/spec_lists/xray_list_SNR10_z1.9.csv" # using a csv as reference to match to
match2 = os.getcwd() + '/Downloading_SDSS_specs/SDSS_fullspecs' # folder that contains files that need to be ordered

# implementing function, defining sections of the file names that can be matched up, what the new column name and csv name should be
match_sources(
    source1 = csv1,
    source2 = match2,
    sections1 = [-2, -1],
    sections2 = [-2, -1],
    column_name1="SPECTRA FILE NAME",
    output_csv="test_match_name_match.csv",
    debug=True
)

'''


