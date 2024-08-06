# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 13:36:07 2024

@author: tayta
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import StrangeVariabilityCoverageOptical.py as STVC



filename = 'StrangeVariability.csv'
data = pd.read_csv(filename)


#CASE 1
a12 = data['I01']/data['I02']
tau = np.arange(0.1,5,0.1)
Cf1 = 0.1
index = 0

plt.figure(figsize=(12,5))
plt.subplot(1,2,1)


 
 




#CASE 2
a21 = data['I02']/data['I01']

tau = np.arange(0.1,5.1,0.1)

    

