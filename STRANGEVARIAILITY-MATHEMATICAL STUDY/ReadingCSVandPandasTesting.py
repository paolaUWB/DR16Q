# -*- coding: utf-8 -*-
"""
Created on Wed May  1 14:40:39 2024

@author: tayta
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def calculations(x,y):
    return x+y

def plotting(x,y,title,xaxis,yaxis):
    plt.figure(1)
    plt.plot(x,y)
    plt.title(title)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    
#Reads the CSV file given by the strange Variability finder code by Alex
filename = 'StrangeVariability.csv'

data = pd.read_csv(filename)


#Code for Testing out pandas
for index,row in data.iterrows():
    #print('I01: ' + str(row['I01']), 'I02: ' + str(row['I02']))
    print('Calculating: ' + str(row['I01']) + ' + ' + str(row['I02']) + ' = ' + str(calculations(row['I01'],row['I02'])))
    #plotting(row['I01'], row['I02'], 'Testing', 'I01', 'I02')
    
coverage = np.array([1,0.5])
#data = data.assign(Cf= [1,0.5])
#data = data.assign(tau = [0.65,0.89])

'''
The plan for the Strange Variability code is to have it so that it is the most efficient by using 
pandas as a way to refer to the data.

The spec files are more of a reference, less of a need.
We most need to use the values of intensity. 

The last column will be added to the data by calculating
the mintau value. I don't think that there will be a column for the optical depth or coverage
fraction unless we want the minimum value for that or whatever else.

The data will be directly called using a for loop to go through each of the situations.
ie.:
    for loop:
        calls coverage fraction calculator(row['I01'],row['I02'],row['I2'])
        calls optical depth calculator(row['I01'],row['I02'],row['I2'])
        calls graphing function(....)
        
then moves through the dataframe till the end of the file.
'''


