# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 15:07:27 2024

@author: tayta
"""

import numpy as np
import pandas as pd
import StrangeVariabilityPlottingFunctions as svp

filename = 'StrangeVariability.csv'
filename = 'OpticalCoverageCalculations.csv'
data = pd.read_csv(filename)


graph = 'yes'




for index, row in data.iterrows():
    
    svp.alpha_Testing(row['I01'],row['I02'],row['I1'],row['I2'], row['minCov1'], 
                      row['minCov2'], row['Spec_Cf1'], row['Spec_Cf2'], graph)