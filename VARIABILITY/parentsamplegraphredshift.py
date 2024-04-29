# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 14:35:51 2024

@author: phyfr
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

if os.getcwd()[-5:] != 'DR16Q': # Sets directory to the default location for DR16Q
    user = os.getlogin()
    drive = os.getenv('SYSTEMDRIVE')
    path = f"{drive}/Users/{user}"

    os.chdir(path + '/Documents/Github/DR16Q')

# Importing csv data
EHVO = pd.read_csv("DR16_EHVO_sorted_norm.csv", header = None)
EHVO_variable = pd.read_csv("Variability_Results  - 1.0EWUnsmoothedFiltered.csv")
parent = pd.read_csv("DR16_parent_sample.csv", header=None)
BAL_parent = pd.read_csv("DR16_BAL_parent_sample.csv", header = None)
rankine = pd.read_csv("DR16parent_DR14RankineInfo.csv", header = None)

# Assigning column names
EHVO.columns = ["Spectra Name", "Redshift", "SNR", "Plate", "MJD", "Fiber"]
parent.columns = ["Spectra Name", "Redshift", "SNR"]
BAL_parent.columns = ["Spectra Name", "Redshift", "BI Value", "BI Error","Unknown?", "SNR"]
rankine.columns = ["Object Name", "Plate", "MJD", "Fiber", "Redshift", "CIV_blue0", "CIV_EW0", "HeII_EW0", "BAL0" , "nBAL0", "BI_BI0", 'BI_VMAX0', 'BI_VMIN0', 'BI_FMIN0',
'CIVdist0', 'SNR0', 'log_MBH_Ran0', 'log_Lbol_Ran0', 'log_Redd_Ran0', 'log_MBH_Rak140','log_MBH_err_Rak140' ,'Q_MBH_Rak140' ,'log_Lbol_Rak140' ,'Q_Lbol_Rak140' ,
'log_Redd_Rak140' ,'Q_Redd_Rak140' ,'BI_CIV_Rak140', 'BI_CIV_err_Rak140', 'BAL_flag_Rak140']

rankine.insert(4, 'Spectra Name', 0)

EHVO_min = EHVO["Redshift"].min()
EHVO_max = EHVO["Redshift"].max()

parent_min = parent["Redshift"].min()
parent_max = parent["Redshift"].max()

BAL_parent_min = BAL_parent["Redshift"].min()
BAL_parent_max = BAL_parent["Redshift"].max()

interval = np.linspace(1.9, parent_max,25)

unique = EHVO_variable.drop_duplicates(subset=["Redshift"])

rankine['Plate'], rankine['MJD'], rankine['Fiber'] = rankine['Plate'].astype(str) , rankine['MJD'].astype(str), rankine['Fiber'].astype(str)

# Filling in the spectra name column to match the others so we can use them as a key later.
for i in range(len(rankine["Redshift"])):
    
    if len(rankine['Fiber'][i]) == 1:
        rankine['Fiber'][i] = '000' + rankine['Fiber'][i]
    
    if len(rankine['Fiber'][i]) == 2:
        rankine['Fiber'][i] = '00' + rankine['Fiber'][i]
    
    if len(rankine['Fiber'][i]) == 3:
        rankine['Fiber'][i] = '0' + rankine['Fiber'][i]
    
    rankine['Spectra Name'][i] = 'spec-' + rankine['Plate'][i] + '-' + rankine['MJD'][i] + '-' + rankine['Fiber'][i] + '-dered.dr16'



# This merges parent and the 'Spectra Name' and 'Redshift' columns of rankine, it matches them via the 'Spectra Name' in on = 'Spectra Name', how = left, tells it to
# use the keys from the parent dataframe. If a redshift doesnt exist in rankine for parent, then redshift_y = NaN. redshift_x = the parent redshifts 
fixed = pd.merge(parent,rankine[['Spectra Name', 'Redshift']], on = 'Spectra Name', how = 'left')

# This takes the redshifts from redshift_y (rankine redshifts) and puts them in parent. If a redshift_y doesnt exist (NaN) fill that with redshift_x (the original)
# parent redshift. Thus rankine redshifts get put in when available and if no rankine exists, the parent is put back in (redshift_x) via .fillna
parent['Redshift'] = fixed['Redshift_y'].fillna(fixed['Redshift_x'])

fixers = pd.merge(BAL_parent,rankine[['Spectra Name', 'Redshift']], on = 'Spectra Name', how = 'left')
BAL_parent['Redshift'] = fixers['Redshift_y'].fillna(fixers['Redshift_x'])

fixest = pd.merge(EHVO,rankine[['Spectra Name', 'Redshift']], on = 'Spectra Name', how = 'left')
EHVO['Redshift'] = fixest['Redshift_y'].fillna(fixest['Redshift_x'])



parent["Redshift"].plot.hist(alpha = 0.5, bins = 25, color = "#96948c", label = "Parent Sample")
BAL_parent["Redshift"].plot.hist(alpha = 0.5, bins = 25, color = "#363333", label = "BAL Parent Sample")
EHVO["Redshift"].plot.hist(alpha = 0.5, bins = 25, color = 'k', label = "RH+ EHVOs")
unique["Redshift"].plot.hist(alpha = 0.6, bins = 25, label = "This Study", color = "k")



plt.title("")
plt.yscale('log')
plt.gca().yaxis.set_major_formatter(ScalarFormatter())
plt.ylabel("Number of Quasars")
plt.xlabel("Redshift")
plt.legend()
plt.show()

# EHVO.to_csv('EHVO_WithRankineRedshift.csv')
# parent.to_csv('DR16Parent_WithRankineRedshift.csv')
# BAL_parent.to_csv('BALParent_WithRankineRedshift.csv')

# EHVO["Bucket"] = pd.cut(EHVO["Redshift"], bins = interval, include_lowest=True)
# group = EHVO.groupby("Bucket").size()

# parent["Bucket"] = pd.cut(parent["Redshift"], bins = interval, include_lowest=True)
# grouper = parent.groupby("Bucket").size().reset_index()
# grouper.columns = ["bucket", "count"]
# # group.plot(kind = "bar", color = 'g')
# # grouper.plot(kind = "bar", color = 'g')

# print(grouper["count"])
# plt.plot(grouper["bucket"], grouper["count"])
# plt.yscale("log")
# plt.show()