# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 14:24:02 2024

@author: phyfr
"""

import numpy as np
import matplotlib.pyplot as plt

lister = [71,72,72,73,71,73,75,74,72,73,71,74,71,71,71,72,72]

print(np.std(lister))

standards = []
window_size = 5
for i in range(len(lister)-window_size):
    
    current = np.nanstd(lister[i:i+window_size], ddof = 1)
    standards.append(current)



print(np.nanstd(lister, ddof = 1))


np.linspace()
