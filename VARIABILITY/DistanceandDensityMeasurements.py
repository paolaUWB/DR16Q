# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 14:42:16 2024

@author: phyfr
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

if os.getcwd()[-11:] != 'VARIABILITY': # Sets directory to the default location for Variability
    user = os.getlogin()
    drive = os.getenv('SYSTEMDRIVE')
    path = f"{drive}/Users/{user}"

    os.chdir(path + '/Documents/Github/DR16Q/VARIABILITY')

# t_recomb = tmin

rankine_info = pd.read_csv('DR16EHVO_DR14RankineInfo_wEHVOspeed_ONLYVARIABLE.csv', header = None)
variability_info = pd.read_csv('FULLCorrected Variability Results - 1.0EWUnsmoothedFiltered.csv')
rankine_info = rankine_info.drop(range(29,37), axis = 1)

rankine_info.columns = ["SDSS name", "Plate", "MJD", "Fiber", "Redshift", "CIV_blue0", "CIV_EW0", "HeII_EW0", "BAL0" , "nBAL0", "BI_BI0", 'BI_VMAX0', 'BI_VMIN0', 'BI_FMIN0',
'CIVdist0', 'SNR0', 'log_MBH_Ran0', 'log_Lbol_Ran0', 'log_Redd_Ran0', 'log_MBH_Rak140','log_MBH_err_Rak140' ,'Q_MBH_Rak140' ,'log_Lbol_Rak140' ,'Q_Lbol_Rak140' ,
'log_Redd_Rak140' ,'Q_Redd_Rak140' ,'BI_CIV_Rak140', 'BI_CIV_err_Rak140', 'BAL_flag_Rak140']

variability_info['SDSS name'] = variability_info['SDSS name'].str[1:]

# lbol is column 17 (for now, check with Paola)

lbol = rankine_info['log_Lbol_Ran0'].tolist()
combined = pd.merge(rankine_info, variability_info, on = 'SDSS name', how = 'left')
tmin = combined['Tmin From EHVO'].unique()
U = 0.0251188643   # Log(U) = -1.6 -> U = 10^-1.6 -> 0.0251188643

# tmin = np.array(tmin)





def ncalc(t, alpha = 2.8e-12):
    n = 1/(alpha*(t*86400))
    return n
n = ncalc(tmin)
r = np.sqrt(lbol/(U*n))


