 # -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 22:32:42 2024

@author: Original Code written by: Dr. Paola Rodriguez Hidalgo
        Modified by Taylor Gibbons
"""


#imports
import numpy as np

import pandas as pd

import time

import os

import StrangeVariabilityCalculationFunctions as svc

import StrangeVariabilityPlottingFunctions as svp


##########################################################################################################################################################

#Main Programming
'''
If there is going to be anything added, make sure that it is being added here
(unless it is a new function then whatever works bruh)
'''


graph = 'No'

#yes = saving figures to CovOpt_Plots folder
save_figs = 'no'

start_time = time.time()

#Name of file being read
filename = os.getcwd() + '/CSV Files/out.csv'
#filename = 'SVTesting.csv'
#filename = 'Alpha_Group_Template.csv'

#Name of file being saved to
#out_filename = 'Alpha_Template.csv'
out_filename = os.getcwd() + '/CSV Files/StrangeVariabilityCalculations.csv'

#Makes a dataframe of all the data from csv using PANDAS
og_data = pd.read_csv(filename)


#Removes all of the detections where In > I0n, or I01 = I02 (which goes against the idea of
#   Strange Variability), then reindexes the dataframe
data = og_data[og_data['I1'] < og_data['I01']]
data = data[data['I01'] != data['I02']]
data = data[data['I2'] < data['I02']]
data.reset_index(inplace = True, drop = True)


#Initializes all of the variable lists that will eventually be saved into the dataframe
minOpt = []
minCov1 = []
minCov2 = []
alpha = []
#new_alpha = []
Spec_Cf1 = []
Spec_Cf2 = []
Spec_Cf2Cf1 = []


#How we determine if we have an object that has a repeated detection
prev_name = ''
prev_num = 0
number = 1


#Starts going through the dataframe row by row
for index,row in data.iterrows():    
    
    #print('Working on plotting ' + row['objectName'])
    
    #If there is a repeated name, then it increases the prev_num so that it will
    #not overwrite the graph that is already been ploted
    if prev_name == row['objectName']:
        prev_num += 1
    else:
        prev_num = 0

    
    alp,inv_alpha, mintau, Spec_v_1, Spec_v_2, Spectral_Cf1, Spectral_Cf2, Spectral_Cf2Cf1 = svc.min_value_finder(row['I01'], row['I02'], row['I1'], row['I2'], prev_num)
    prev_name = row['objectName']
    minOpt.append(mintau)
    minCov1.append(Spec_v_1)
    minCov2.append(Spec_v_2)
    alpha.append(alp)
    Spec_Cf1.append(Spectral_Cf1)
    Spec_Cf2.append(Spectral_Cf2)
    Spec_Cf2Cf1.append(Spectral_Cf2Cf1)
    number += 1
     

data['alpha'] = alpha  
#data['new_alpha'] = new_alpha #New alpha is a ratio of change (smaller / larger)
data['minOpt'] = minOpt
data['minCov1'] = minCov1
data['minCov2'] = minCov2
data['Spec_Cf1'] = Spec_Cf1
data['Spec_Cf2'] = Spec_Cf2
data['Spec_Cf2Cf1'] = Spec_Cf2Cf1
data['save_fig'] = save_figs

data_alpha = data.sort_values(by = ['alpha'])
data_alpha.reset_index(inplace = True, drop = True)

data_alpha.to_csv(out_filename,index =False)


if graph == 'Yes':
    svp.Cf_grapher(data)
    


prog_time = time.time() - start_time
print('Program Finished\n Took: ' + str(prog_time) + ' seconds')


    #Start figuring a way to move the legend per case so that we are not 
    #covering the important info
