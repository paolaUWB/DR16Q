# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 11:08:16 2024

@author: Imadm
"""
import numpy as np

import pandas as pd

import matplotlib.pyplot as plt 

import time

#import UpdatedPlotting as UP





#t = float(time.time())

calcs = 'TestCalculations.csv'
#calcs = 'OpticalCoverageCalculations.csv'
calc_data = pd.read_csv(calcs)

alpha = calc_data['alpha']
tau = np.arange(0.1,5,0.1)
Cf = 0.1
Cf1 = []


Sv1 = calc_data['minCov1']

for m in range(len(calc_data)): 
    if alpha[m] > 1:
        alpha[m] = alpha[m] ** -1
    fig, graph = plt.subplots(nrows=1, ncols=2, figsize = (12,7))
#    plt.subplots_adjust(space=0.5)

#    plt.figure(figsize = (15,10))
    graph[1].set_xlim(-1,5)
    graph[1].set_ylim(0,1)
    graph[1].set_ylabel(r'$C_{f1}$')
    graph[1].set_xlabel(r'$\tau$') 
    graph[1].set_title(r'$C_{f2}$' + ' vs. ' + r'$\tau$')

    graph[0].set_xlim(-1,5)
    graph[0].set_ylim(0,1)
    graph[0].set_ylabel(r'$C_{f2}$')
    graph[0].set_xlabel(r'$\tau$')       #copied the latex lol
    graph[0].set_title(r'$C_{f1}$' + ' vs. ' + r'$\tau$')
    
    
    
    
    mincov12 = np.zeros(13)+ (calc_data.minCov2[m] / calc_data.minCov1[m])
    
    minOpt = calc_data['minOpt']
    Cf = 0.1
    mintau = calc_data.loc[m]['minOpt']
    
    plt.plot(tau, calc_data['Spec_Cf1'][m], color = 'r')
    
#    calc_data['Spec_Cf1'].str.split('').apply(set)
#    calc_data['Spec_Cf1'] = calc_data['Spec_Cf1'].str.split(' ')

    for n in range(10):
        Cf1 = (Cf / alpha[m]) - ((alpha[m] - 1) / (alpha[m] * ((np.exp(-tau) -1) )))
        Cf2 = (Cf / (alpha[m] ** -1)) - (((alpha[m] ** -1) - 1) / ((alpha[m] ** -1) * ((np.exp(-tau) -1) )))
        graph[0].plot(tau,Cf1, color=[0,0.5,1-Cf])
        graph[1].plot(tau,Cf2, color=[0,0.5,1-Cf])
#        print(Cf)
        Cf += 0.1
        
        
#Vertical Lines  
    graph[0].plot(np.zeros(11)+minOpt[m], np.arange(0,1.1,0.1), 'r--')
    graph[1].plot(np.zeros(11)+minOpt[m], np.arange(0,1.1,0.1), 'r--')
    
#Horizontal Lines
    graph[0].plot(np.arange(-1,5.5,0.5), mincov12, 'r')
    graph[1].plot(np.arange(-1,5.5,0.5), mincov12, 'r')
#    graph[0].plot(np.arange(0.1,5,1), Sv1, 'r')
#    graph[1].plot(tau, minSv2[m], 'r')
    


    plt.show()
    
    




#t = time.time() - t
#print(t)
# UP.Cf_grapher(calc_data)
    