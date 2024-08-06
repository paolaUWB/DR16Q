# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 11:07:01 2024

@author: tayta
"""

import numpy as np
import matplotlib.pyplot as plt

def epsilon_finder(I01, I02, mintau, cf1):
    '''

    Parameters
    ----------
    I01 : float, int
        This is the value for the average flux around the trough on the first observation
    I02 : float, int
        This is the value for the average flux around the trough on the second observation
    mintau : float
        This is the minimum value of the optical depth based on the initial fluxes.

    Returns
    -------
    Epsilon: float
        A value that we can consider to be the most pessimistic value of error in order to determine if the 
        observations are considered useful data. 
    '''
    
    alpha = I01/I02
    epsilon = (alpha*cf1 - (alpha - 1))*(I02*(np.exp(-mintau)-1))
    
    return epsilon



#testing Epsilon Finder

I01 = np.arange(1,100,1)
I02 = reversed(np.array(1,100,1))

alpha = I01/I02

Error = []
mintau = 0.687

cf1 = .50

for i in range(len(I01)):
    Error.append(epsilon_finder(float(I01[i]),float(I02[i]),mintau,cf1))
    
plt.figure(1)
plt.title('Error Value vs Alpha')
plt.xlabel('Alpha Ratio')
plt.ylabel('Error Value')
plt.plot(alpha,Error)




