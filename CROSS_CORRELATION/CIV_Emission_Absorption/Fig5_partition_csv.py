#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 11:52:56 2026

@author: lilianaflores

This program produces a csv with the file names, mjd, plate, fiber, HeII EW, and CIV Distance 
of cases above the partition line in figure 5, 'EHVOs_HeII_CIVDist_above_partition.csv', and cases below the partition line, 'EHVOs_HeII_CIVDist_below_partition.csv'.

Note that the cases will not add up to 99 as we do not have HeII EW values for all cases?
Some are -999999 or 0000000 

"""


import os
import sys
from os.path import exists
import csv
import pandas as pd
import scipy.stats as stat
import numpy as np
from pylab import*
from sympy import sympify
from matplotlib.backends.backend_pdf import PdfPages
from scipy import*
from astropy import*
from astropy.table import Table
from scipy.stats import ks_2samp
from astropy import constants as const
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import pearsonr
from astropy.io import fits
from astropy import stats
import re
import ndtest #pip installed from here: https://github.com/syrte/ndtest/blob/master/README.md
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'


# plotting function imports
from plot_functions import plot_CIVBlue_hexbin_ehvo, figure_2_outline, plot_EHVO_CIV_scatter_hist, figure_5_outline, figure_1

##############################################

#Rankines info file PRH selected the good cases
infoRankineEHVO = os.getcwd() + "/../DR16EHVO_DR14RankineInfo_wEHVOspeed_vmax_reordered.csv" #adding BAL flag corrected version with reordered vmax
infoRankineEHVO_DR9 = os.getcwd() + "/../DR9EHVO_DR14RankineInfo_wEHVOspeed_vmax_reordered.csv" #adding BAL flag corrected version with reordered vmax





#Extracting values from Rankines info files
dfRHV = pd.read_csv(infoRankineEHVO, header=0)
dfRHV_9 = pd.read_csv(infoRankineEHVO_DR9, header=0)


#EHVOs
plate_EHVOR=dfRHV[dfRHV.columns[1]].to_numpy()
mjd_EHVOR=dfRHV[dfRHV.columns[2]].to_numpy()
fiber_EHVOR=dfRHV[dfRHV.columns[3]].to_numpy()

CivBlue_EHVOR=dfRHV[dfRHV.columns[5]].to_numpy()
CivEW_EHVOR=dfRHV[dfRHV.columns[6]].to_numpy()
CivDist_EHVOR=dfRHV[dfRHV.columns[14]].to_numpy()
HeiiEW_EHVOR=dfRHV[dfRHV.columns[7]].to_numpy()
Lbol_EHVOR=dfRHV[dfRHV.columns[17]].to_numpy()
Edd_EHVOR=dfRHV[dfRHV.columns[18]].to_numpy()
MBH_EHVOR=dfRHV[dfRHV.columns[16]].to_numpy()


plate_EHVOR_9=dfRHV_9[dfRHV_9.columns[1]].to_numpy()
mjd_EHVOR_9=dfRHV_9[dfRHV_9.columns[2]].to_numpy()
fiber_EHVOR_9=dfRHV_9[dfRHV_9.columns[3]].to_numpy()

CivBlue_EHVOR_9=dfRHV_9[dfRHV_9.columns[5]].to_numpy()
CivEW_EHVOR_9=dfRHV_9[dfRHV_9.columns[6]].to_numpy()
CivDist_EHVOR_9=dfRHV_9[dfRHV_9.columns[14]].to_numpy()
HeiiEW_EHVOR_9=dfRHV_9[dfRHV_9.columns[7]].to_numpy()
Lbol_EHVOR_9=dfRHV_9[dfRHV_9.columns[17]].to_numpy()
vmin_EHVO=dfRHV[dfRHV.columns[29]].to_numpy()
vmin_EHVO_9=dfRHV_9[dfRHV_9.columns[29]].to_numpy()
vmin_EHVO_rel=dfRHV[dfRHV.columns[31]].to_numpy()
vmin_EHVO_9_rel=dfRHV_9[dfRHV_9.columns[31]].to_numpy()
vmax_EHVO=dfRHV[dfRHV.columns[33]].to_numpy()
vmax_EHVO_9=dfRHV_9[dfRHV_9.columns[32]].to_numpy()
Edd_EHVOR_9=dfRHV_9[dfRHV_9.columns[18]].to_numpy()
MBH_EHVOR_9=dfRHV_9[dfRHV_9.columns[16]].to_numpy()

#combining 9 with 16
plate_EHVO_combined = np.concatenate((plate_EHVOR_9, plate_EHVOR))
mjd_EHVO_combined = np.concatenate((mjd_EHVOR_9, mjd_EHVOR))
fiber_EHVO_combined = np.concatenate((fiber_EHVOR_9, fiber_EHVOR))


CivBlue_EHVO_combined = np.concatenate((CivBlue_EHVOR_9, CivBlue_EHVOR))
CivEW_EHVO_combined = np.concatenate((CivEW_EHVOR_9, CivEW_EHVOR))
CivDist_EHVOR_combined = np.concatenate((CivDist_EHVOR_9, CivDist_EHVOR))
HeiiEW_EHVOR_combined = np.concatenate((HeiiEW_EHVOR_9, HeiiEW_EHVOR))
Lbol_EHVOR_combined = np.concatenate((Lbol_EHVOR_9, Lbol_EHVOR))
vmin_EHVOR_combined = np.concatenate((vmin_EHVO_9, vmin_EHVO))
vmin_EHVOR_rel_combined = np.concatenate((vmin_EHVO_rel, vmin_EHVO_9_rel))
vmax_EHVOR_combined = np.concatenate((vmax_EHVO_9, vmax_EHVO))
Edd_EHVO_combined = np.concatenate((Edd_EHVOR_9, Edd_EHVOR))
MBH_EHVOR_combined = np.concatenate((MBH_EHVOR_9, MBH_EHVOR))

#EHVO non-BAL vs EHVO BAL
EHVO_BAL_9 = np.where(dfRHV_9['BAL0']==True)
EHVO_BAL_16 = np.where(dfRHV['BAL0']==True)

EHVO_nonBAL_9 = np.where(dfRHV_9['BAL0']==False)
EHVO_nonBAL_16 = np.where(dfRHV['BAL0']==False)

CivDist_EHVOR_combined_BAL = np.concatenate((CivDist_EHVOR_9[EHVO_BAL_9], CivDist_EHVOR[EHVO_BAL_16]))
CivDist_EHVOR_combined_nonBAL = np.concatenate((CivDist_EHVOR_9[EHVO_nonBAL_9], CivDist_EHVOR[EHVO_nonBAL_16]))
vmin_EHVOR_combined_BAL = np.concatenate((vmin_EHVO_9[EHVO_BAL_9], vmin_EHVO[EHVO_BAL_16]))
vmin_EHVOR_combined_nonBAL = np.concatenate((vmin_EHVO_9[EHVO_nonBAL_9], vmin_EHVO[EHVO_nonBAL_16]))
vmax_EHVOR_combined_BAL = np.concatenate((vmax_EHVO_9[EHVO_BAL_9], vmax_EHVO[EHVO_BAL_16]))
vmax_EHVOR_combined_nonBAL = np.concatenate((vmax_EHVO_9[EHVO_nonBAL_9], vmax_EHVO[EHVO_nonBAL_16]))


outlier = np.where(np.log10(CivEW_EHVOR)>2.1)
outlier
print(f'The outlier is: {outlier}')

# Remove outlier
outlier = np.where(np.log10(CivEW_EHVOR)>2.1)
copyCivBlue = np.copy(CivBlue_EHVOR)
no_outlier_CivBlue = np.delete(copyCivBlue, [40])
copyCivEW = np.copy(CivEW_EHVOR)
no_outlier_CivEW = np.delete(copyCivEW, [40])
copyCivDist = np.copy(CivDist_EHVOR)
no_outlier_CivDist = np.delete(copyCivDist, [40])
copyHeiiEW = np.copy(HeiiEW_EHVOR)
no_outlier_HeiiEW = np.delete(copyHeiiEW, [40])
copyEdd_EHVO = np.copy(Edd_EHVOR)
no_outlier_Edd_EHVO = np.delete(copyEdd_EHVO, [40])
copyvmin_EHVO = np.copy(vmin_EHVO)
no_outlier_vmin_EHVO = np.delete(copyvmin_EHVO, [40])
copyLbol_EHVO = np.copy(Lbol_EHVOR)
no_outlier_Lbol_EHVO = np.delete(copyLbol_EHVO, [40])
copyMBH_EHVO = np.copy(MBH_EHVOR)
no_outlier_MBH_EHVO = np.delete(copyMBH_EHVO, [40])

copyvmax_EHVO = np.copy(vmax_EHVO)
no_outlier_vmax_EHVO = np.delete(copyvmax_EHVO, [40])

copymjd_EHVO = np.copy(mjd_EHVOR)
no_outlier_mjd_EHVO = np.delete(copymjd_EHVO, [40])
copyplate_EHVO = np.copy(plate_EHVOR)
no_outlier_plate_EHVO = np.delete(copyplate_EHVO, [40])
copyfiber_EHVO = np.copy(fiber_EHVOR)
no_outlier_fiber_EHVO = np.delete(copyfiber_EHVO, [40])

outlier_idx = 40

dfRHV_noOut = dfRHV.drop(index=outlier_idx).copy()
dfRHV_noOut = dfRHV_noOut.reset_index(drop=True)


EHVO_BAL_9 = np.where(dfRHV_9['BAL0']==True)
EHVO_BAL_16 = np.where(dfRHV_noOut['BAL0']==True)

EHVO_nonBAL_9 = np.where(dfRHV_9['BAL0']==False)
EHVO_nonBAL_16 = np.where(dfRHV_noOut['BAL0']==False)

#combining 9 with 16 and remove outlier
CivBlue_EHVO_combined_noOut = np.concatenate((CivBlue_EHVOR_9, no_outlier_CivBlue))
CivEW_EHVO_combined_noOut = np.concatenate((CivEW_EHVOR_9, no_outlier_CivEW))
CivDist_EHVOR_combined_noOut = np.concatenate((CivDist_EHVOR_9, no_outlier_CivDist))
HeiiEW_EHVOR_combined_noOut = np.concatenate((HeiiEW_EHVOR_9, no_outlier_HeiiEW))
Edd_EHVO_combined_noOut = np.concatenate((Edd_EHVOR_9, no_outlier_Edd_EHVO))
Lbol_EHVOR_combined_noOut = np.concatenate((Lbol_EHVOR_9, no_outlier_Lbol_EHVO))
MBH_EHVOR_combined_noOut = np.concatenate((MBH_EHVOR_9, no_outlier_MBH_EHVO))
vmin_EHVOR_combined_noOut = np.concatenate((vmin_EHVO_9, no_outlier_vmin_EHVO))
vmax_EHVOR_combined_noOut = np.concatenate((vmax_EHVO_9, no_outlier_vmax_EHVO))

mjd_EHVOR_combined_noOut = np.concatenate((mjd_EHVOR_9, no_outlier_mjd_EHVO))
plate_EHVOR_combined_noOut = np.concatenate((plate_EHVOR_9, no_outlier_plate_EHVO))
fiber_EHVOR_combined_noOut = np.concatenate((fiber_EHVOR_9, no_outlier_fiber_EHVO))



# line option 2
HeiiEW_bisect_above = np.array([0, 1, 3, 9, 11, 19, 27, 29, 42, 47, 53, 60, 63, 65, 67, 69, 70, 73, 77, 89, 91])
#HeiiEW_bisect_below = np.array([4, 5, 6, 8, 10, 13, 14, 15, 16, 17, 18, 20, 21, 22, 24, 25, 26, 28, 30, 31, 32, 34, 35, 37,38, 39, 40, 43, 44, 45, 46, 48, 49, 50, 51, 52, 55, 56, 57, 58, 59, 61, 62, 64, 66, 68, 71, 72, 74, 75, 76, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 93, 94, 95, 97, 98])
# new below bisect array vvv
HeiiEW_bisect_below = np.array([4, 5, 6, 8, 10, 13, 14, 15, 16, 17, 18, 20, 21, 22, 24, 25, 26, 28, 30, 31, 32, 34, 35, 38,
 39, 40, 44, 45, 46, 48, 49, 50, 51, 52, 55, 56, 57, 58, 59, 61, 62, 64, 66, 68, 71, 72, 74,
 75, 76, 78, 79, 81, 83, 84, 85, 86, 87, 94, 95, 97, 98])
# note there are less than 99 EHVOs on this plot as 11 log10(HeII) values are nan
# there are 6 EHVOs with log10(HeII) = -inf that needed to be changed to nan values with np.log10(np.where(HeiiEW_EHVOR_combined_noOut> 0, HeiiEW_EHVOR_combined_noOut, np.nan))
# The new array adds up to 82 EHVOs plotted in figure 5

spec_file_name = []
for i in range(len(mjd_EHVOR_combined_noOut)):
    fiber_EHVOR_combined_noOut_str = str(fiber_EHVOR_combined_noOut[i]).zfill(4)
    spec_file_name.append('spec-' + str(plate_EHVOR_combined_noOut[i]) + '-' + str(mjd_EHVOR_combined_noOut[i]) + '-' + fiber_EHVOR_combined_noOut_str + 'norm.dr')

spec_file_name = np.array(spec_file_name)
above_dict = {'Spec File Name': spec_file_name[HeiiEW_bisect_above],
              'mjd': mjd_EHVOR_combined_noOut[HeiiEW_bisect_above],
              'plate': plate_EHVOR_combined_noOut[HeiiEW_bisect_above],
              'fiber': fiber_EHVOR_combined_noOut[HeiiEW_bisect_above],
              'HeII EW': HeiiEW_EHVOR_combined_noOut[HeiiEW_bisect_above],
              'CIV Distance': CivDist_EHVOR_combined_noOut[HeiiEW_bisect_above]
              }
  

above_df = pd.DataFrame(above_dict, index=None)
above_df.to_csv('EHVOs_HeII_CIVDist_above_partition.csv', index=False)  
    
below_dict = {'Spec File Name': spec_file_name[HeiiEW_bisect_below],
              'mjd': mjd_EHVOR_combined_noOut[HeiiEW_bisect_below],
              'plate': plate_EHVOR_combined_noOut[HeiiEW_bisect_below],
              'fiber': fiber_EHVOR_combined_noOut[HeiiEW_bisect_below],
              'HeII EW': HeiiEW_EHVOR_combined_noOut[HeiiEW_bisect_below],
              'CIV Distance': CivDist_EHVOR_combined_noOut[HeiiEW_bisect_below]
              }


below_df = pd.DataFrame(below_dict)
below_df.to_csv('EHVOs_HeII_CIVDist_below_partition.csv', index=False)


all_dict = {'Spec File Name': spec_file_name,
              'mjd': mjd_EHVOR_combined_noOut,
              'plate': plate_EHVOR_combined_noOut,
              'fiber': fiber_EHVOR_combined_noOut,
              'HeII EW': HeiiEW_EHVOR_combined_noOut,
              'CIV Distance': CivDist_EHVOR_combined_noOut
              }


all_df = pd.DataFrame(all_dict)
all_df.to_csv('EHVOs_HeII_CIVDist_all.csv', index=False)