# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 13:11:38 2024

@author: tayta
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

filename = 'OpticalCoverageCalculations.csv'

data = pd.read_csv(filename)



#Pivot Table Testing and data initialization

table = pd.pivot_table(data, values = ['new_alpha','minOpt','minCov1','minCov2'],index = ['objectName'], aggfunc= 'mean')
xdata = table['new_alpha']
cov1 = table['minCov1']
cov2 = table['minCov2']
opt = table['minOpt']
#Bin Histogram example testing

xbins = np.array([0,0.25,0.5,0.75,1])

style = {'facecolor': 'none', 'edgecolor': 'C0','linewidth':3}

fig,ax = plt.subplots()

ax.hist(xdata,bins = xbins, **style)
ax.set_ylabel('Number of Objects with new alpha values')
ax.set_xlabel('New Alpha')
ax.set_title('Number of Objects vs Value of New Alpha')

#Making labels: WORK IN PROGRESS

rects = ax.patches
labels = ['Detections: %d' % i for i in range(len(rects))]

for rect, label in zip(rects, labels):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width() / 2, height + 0.01, label, ha= 'center', va = 'bottom')

#Hexbin testing
fig,ax2 = plt.subplots()
ax2.set_ylabel('Minimum Coverage Fraction n = 2')
ax2.set_xlabel('Minimum Coverage Fraction n = 1')
ax2.hexbin(cov1,cov2, xdata)
