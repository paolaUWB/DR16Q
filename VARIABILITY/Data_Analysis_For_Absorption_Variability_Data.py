# -*- coding: utf-8 -*-
"""
Created on Mon May 29 20:38:00 2023

@author: phyfr
"""

import numpy as np
import pandas as pd
data = pd.read_csv('Variability_Results - Sheet1 - Variability_Results - Sheet1.csv',',')
numpy_data = data.to_numpy()

data['MJD'] = data['Spec Files'].astype(str).str[-10:-5]


grouped = data.groupby('SDSS name')

lister = []

quasars = {}

for name,group in grouped:
    quasars[name] = group
    
# def find_smallest_difference(df, column_name): #Find smallest values between the MJD (you can pass in any column though)
#     sorted_column = df[column_name].sort_values() #Sorting allows us to subtract subquential rows, ensuring that for ex: row 1 and row 4 will never be the smallest diff
#     smallest_diff = float('inf')  # Initialize with infinity

#     for i in range(len(sorted_column) - 1):
#         diff = sorted_column.iloc[i + 1] - sorted_column.iloc[i]
#         smallest_diff = min(smallest_diff, diff)

#     return smallest_diff

def find_smallest_difference(df, column_name):
    sorted_df = df.sort_values(column_name)
    smallest_diff = float('inf')  
    indices = None  

    for i in range(len(sorted_df) - 1):
        diff = sorted_df[column_name].iloc[i + 1] - sorted_df[column_name].iloc[i]
        if diff < smallest_diff:
            smallest_diff = diff
            indices = (sorted_df.index[i], sorted_df.index[i + 1])

    return smallest_diff, indices

for name, df in quasars.items():
    df['MJD'] = df['MJD'].astype(int)
    df = df.sort_values('MJD')
    df['Tmin'] = find_smallest_difference(df, 'MJD')[0]/(df['Redshift'] + 1)
    lister.append(find_smallest_difference(df, 'MJD')[1])
    # df['Tmin'] = df['MJD'].nsmallest(2, keep = 'all')
    # df['Tmin'] = df['Tmin'].max() - df['Tmin'].min()

    quasars[name] = df


recombined = pd.concat(quasars.values())

        
