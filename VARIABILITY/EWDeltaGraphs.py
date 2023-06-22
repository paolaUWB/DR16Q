# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 01:33:03 2023

@author: phyfr
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('Variability_Results - Sheet1 - Unsmoothed.csv', ',')

for column in df.columns:
    for i in range(len(df[column])):
        # Check if the cell contains brackets
        if isinstance(df[column][i], str) and '[' in df[column][i]:
            # Remove the brackets using the replace() method
            df[column][i] = df[column][i].replace('[', '').replace(']', '')


# Iterate through each cell in the "EW INDIVIDUAL 500" column
for i in range(len(df['EW INDIVIDUAL 500'])):
    # Check if the cell contains commas
    if isinstance(df['EW INDIVIDUAL 500'][i], str) and ',' in df['EW INDIVIDUAL 500'][i]:
        try:
            # Split the values by commas, strip whitespace, and convert them to integers
            values = [float(x.strip()) for x in df['EW INDIVIDUAL 500'][i].split(',')]
            # Calculate the sum of the values
            total = sum(values)
            # Update the cell with the sum
            df.at[i, 'EW INDIVIDUAL 500'] = total
        except ValueError:
            # Handle the exception when encountering non-numeric values
            df.at[i, 'EW INDIVIDUAL 500'] = 'Error: Non-numeric value'
            
            
df.sort_values('OG EHVO Obs', inplace=True)
df['OG EHVO Obs'] = df['OG EHVO Obs'].astype(float)
df['EW INDIVIDUAL 500'] = df['EW INDIVIDUAL 500'].astype(float)
unique_values = df['Tmin From EHVO'].unique().astype(float)
# tryer = df.sort_values('Redshift')
tryer = df['Redshift'].copy()
print(tryer)
tryer = tryer.unique()

if 'OG EHVO Obs' in df.columns and 'SDSS name' in df.columns:
    # Sort the DataFrame based on 'OG EHVO Obs' column in ascending order
    df.sort_values('OG EHVO Obs', inplace=True)

    # Create a new DataFrame to store the differences
    new_df = pd.DataFrame(columns=['SDSS name', 'Difference'])

    # Iterate through each unique 'SDSS name'
    for name in df['SDSS name'].unique():
        # Get the rows corresponding to the current 'SDSS name'
        rows = df[df['SDSS name'] == name]

        # Get the values from the 'EW INDIVIDUAL 500' column for the current 'SDSS name'
        values = rows['EW INDIVIDUAL 500'].tolist()

        # Calculate the differences based on the order determined by 'OG EHVO Obs'
        differences = [values[i + 1] - values[i] for i in range(len(values) - 1)]

        # Create a new row with 'SDSS name' and 'Difference' values
        new_row = pd.DataFrame({'SDSS name': [name] * len(differences), 'Difference': differences})

        # Append the new row to the new DataFrame
        new_df = new_df.append(new_row, ignore_index=True)
unique_values = df['Tmin From EHVO'].unique().astype(float)
new_df['Tmin From EHVO'] = tryer
new_df['Difference'] = new_df['Difference']*-1

# deltaEWman = [-172.42,164.98,-176.22,1291.03,-212.2,918.75,-422.05,-27.18,304.59,-642.96,
#             712.9,-914.78,-388.64,-51.1,-621.5,-126.52,-527.44,-14.73,189.02,-58.96,-188.88,4.73]
# deltaTminman = [1772.477694,738.8794567,1084.214003,683.2682292,1064.306593,578.7261558,540.8194234,1016.928658,471.6354858,1180.41543,
#               223.550606,180.2232855,858.0343214,204.4041451,564.6053293,442.818,631.0562016,0.713606,1.589825119,196.1797753,158.7713606,1.095090345]

deltaEW = new_df['Difference'].tolist()


# Convert the unique values to a list
Tmin = unique_values.tolist()

# redshift =

print(np.sum(np.abs(deltaEW)))
print(np.sum(Tmin))

plt.scatter(Tmin,deltaEW, label = 'EHVO quasars')
plt.title('Changes in Absorption vs Time')
plt.xlabel('$T_{min}$ from EHVO (Restframe days)')
plt.ylabel('$\Delta$ EW (km*$s^{-1}$)')
# plt.scatter(deltaTminman,deltaEWman, color = 'r')
plt.legend()
