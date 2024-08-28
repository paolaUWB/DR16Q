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

#import openpyxl

##########################################################################################################################################################

#Main Programming
'''
If there is going to be anything added, make sure that it is being added here
(unless it is a new function then whatever works bruh)
'''

#yes = saving figures to CovOpt_Plots folder
graph = 'no'
alpha_grouping = 'no'
save_figs = 'no'

alpha_template = 'no'
start_time = time.time()

#Name of file being read
filename = os.getcwd() + '/CSV Files/out.xlsx'
filename = os.getcwd() + '/CSV Files/out.csv'
filename = os.getcwd() + '/CSV Files/'
#Name of file being saved to
out_filename = os.getcwd() + '/CSV Files/StrangeVariabilityCalculations.xlsx'
alpha_out_filename = os.getcwd() + '/CSV Files/StrangeVariabilityCalculationsAlphaSort.xlsx'
alpha_pickle_out = os.getcwd() + '/CSV Files/StrangeVariabilityCalculationsAlphaSort.pkl'


#Used for alpha grouping! 

if alpha_template == 'yes':
    #sets the filenames of input and output
    filename = os.getcwd() + '/CSV Files/Alpha_Group_Template.xlsx'
    out_filename = os.getcwd() + '/CSV Files/Alpha_Template.pkl'
    #yes will make a new template file that will be based on the interval wanted: IN PROGRESSSSSSSSSS
    make_temp = 'no'

    #the +- value
    epsilon = 0.05

    #where you want the alpha grouping to start
    alpha = 0.25

    #Random value to start getting a template datafram

    if make_temp == 'yes':
         I01 = 10
         I02 = []
         alpha_list = []
    
         while alpha < 1:
             I02.append(I01 * alpha)
             alpha_list.append(alpha)
             alpha = alpha + ( 2 * epsilon)






#Makes a dataframe of all the data from csv using PANDAS
og_data = pd.read_excel(filename)


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
Spec_Cf1 = []
Spec_Cf2 = []
Spec_Cf2Cf1 = []


#How we determine if we have an object that has a repeated detection
prev_name = ''
prev_num = 0



#Starts going through the dataframe row by row
for index,row in data.iterrows():    
    
    #If there is a repeated name, then it increases the prev_num so that it will
    #not overwrite the graph that is already been ploted
    if prev_name == row['objectName']:
        prev_num += 1
    else:
        prev_num = 0

    #calls the calculation function to get values of Cf
    alp,inv_alpha, mintau, Spec_v_1, Spec_v_2, Spectral_Cf1, Spectral_Cf2, Spectral_Cf2Cf1 = svc.min_value_finder(row['I01'], row['I02'], row['I1'], row['I2'], prev_num)

    #keeps track of the previous name so that if there was more than one detection of the same name, then it wouldn't save over
    prev_name = row['objectName']

    #filling the lists with valuse
    minOpt.append(mintau)
    minCov1.append(Spec_v_1)
    minCov2.append(Spec_v_2)
    alpha.append(alp)
    Spec_Cf1.append(Spectral_Cf1)
    Spec_Cf2.append(Spectral_Cf2)
    Spec_Cf2Cf1.append(Spectral_Cf2Cf1)

    




     
#Fills the dataframe with the newly filled lists
data['alpha'] = alpha  
data['minOpt'] = minOpt
data['minCov1'] = minCov1
data['minCov2'] = minCov2
data['Spec_Cf1'] = Spec_Cf1
data['Spec_Cf2'] = Spec_Cf2
data['Spec_Cf2Cf1'] = Spec_Cf2Cf1
data['save_fig'] = save_figs



print(type(data.Spec_Cf1[0]))

#sorts the data frame by alpha values and reindexes
data_alpha = data.sort_values(by = ['alpha'])
data_alpha.reset_index(inplace = True, drop = True)

#This is to make sure all of the detections do not have alpha = 1
data_alpha = data_alpha[data_alpha.alpha < 1]

#saves dataframes to a csv
#data.to_excel(out_filename, index=False)
data_alpha.to_excel(alpha_out_filename,index =False)

data_alpha.to_pickle(out_filename)






#this is what will call the graphing function and send the data frame
if graph == 'yes':
    svp.Cf_tau_grapher(data_alpha)
    
#this is to plot the data based on a selected alpha interval
if alpha_grouping == 'yes':
    grouped_df = svp.alpha_group_grapher(data_alpha)
    grouped_df.to_pickle(os.getcwd() + '/CSV Files/groupedDfPickle.pkl')
    print(type(grouped_df))
    svp.Cf_tau_grapher(grouped_df,'yes')
    
    
#To see how long the program takes to run/graph
prog_time = time.time() - start_time
print('Program Finished\n Took: ' + str(prog_time) + ' seconds, or ' + str(round(prog_time/60,2)) + ' minutes')

