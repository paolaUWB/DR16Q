# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 17:14:55 2024

@author: tayta
"""
import os
import pandas as pd
import numpy as np

filename = os.getcwd() + '/CSV Files/alphaPickle.pkl'
temp_filename = os.getcwd() + '/CSV Files/templateDF.pkl'
grouped_filename = os.getcwd() + '/CSV Files/groupedDfPickle.pkl'


data = pd.read_pickle(filename)
temp_df = pd.read_pickle(temp_filename)
grouped_df = pd.read_pickle(grouped_filename)

#print(type(data.Spec_Cf1[0]))

if isinstance(temp_df,pd.DataFrame):
    print('DATAFRAME')
    
    if isinstance(temp_df.obs_Cf1, pd.Series):
        print('obs_Cf1 is a Series')
        
        if isinstance(temp_df.obs_Cf1[2],np.ndarray):
            print('obs_Cf1[2] is an array')
            
            
        
        
for index,rows in grouped_df.iterrows(): #BIGGEST LOOP
    #print(rows)
    
    
    
    if isinstance(rows.obs_Cf1,np.ndarray):
        print('This Row has an array for its observational Data')
        
    if isinstance(rows.obs_Cf1, pd.Series):
        print('This row is a series object holding observational Data')
        
        #for k in rows.obs_Cf1:
            #print(k)