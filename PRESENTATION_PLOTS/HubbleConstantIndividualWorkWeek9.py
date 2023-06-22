# -*- coding: utf-8 -*-
"""
Created on Sat May 27 16:49:34 2023

@author: phyfr
"""

import numpy as np
import pandas as pd
import math
import seaborn as sns
import scipy.stats as sci
import matplotlib.pyplot as plt
data = pd.read_csv('values.csv')
data = data.sort_values(by =['distance'], axis = 0)
data = data.replace(-1000,math.nan)
data = data.drop(labels = 2, axis = 0)
data['Calculated Redshift']*=3e5
# data.plot(x = 'distance', y = 'Calculated Redshift', xlabel = 'Distance (Mpc)', ylabel = 'Velocity (m/s)', kind = 'scatter' )

graph = sns.regplot(x = data['distance'], y = data['Calculated Redshift'], scatter = True, label = 'Line Of Best Fit', color = 'g')
plt.scatter(data['distance'], data['Calculated Redshift'], label = 'Data')
plt.title('Distance vs Radial Velocity in Nearby Galatic Objects')
plt.xlabel('Distance (Mpc)')
plt.ylabel('Velocity (km/s)')
plt.legend(fontsize = 9, loc = 'lower right')
slope, intercept, r, graph, sterr = sci.linregress(x=graph.get_lines()[0].get_xdata(),
                                                       y=graph.get_lines()[0].get_ydata())
print(slope)
