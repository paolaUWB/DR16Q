#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 12:03:36 2026

@author: lilianaflores
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

#################################################
############################## Changeable Variables
# [Fig 1]
plot_fig1 = True
plot_fig1 = False

# [Fig 2]
plot_fig2 = True
plot_fig2 = False

# [Fig 3] made and provided by Amy Rankine


# [Fig 4]
plot_fig4 = True
plot_fig4 = False

# [Fig 5]
plot_fig5 = True
plot_fig5 = False

# [Fig 6] HeII colormap hex/scatter plot in physical property parameter spaces
plot_fig6 = True
plot_fig6 = False

HeII_flat_distance_cut = True #plotting with a flat CIV Distance cutoff provided by dist_cut
dist_cut = 0.83

HeII_flat_distance_cut = 'line' #plotting with sample separation via line spearation of cluster in Fig 5 (HeII vs CIV Distance scat/hist)
dist_cut = None

hexbin_function6 = np.median
#hexbin_function6 = np.max

gridsize = 35 # this applies to both fig 6 and 8 as they should have the same gridsize to be easily compared

# [Fig 7] corner plot found in CornerPlotsCIVBlueshift.py
gridsize = 35

# [Fig 8] CIV Blueshift colormap hex/scatter plot in physical property parameter spaces
plot_fig8 = True
#plot_fig8 = False

Blueshift_cut = 'flat'
dist_cut2 = 0.9

#Blueshift_cut = 'line'
#dist_cut2 = None

hexbin_function8 = np.median
hexbin_function8 = np.max
##############################################

#The inputs in this program should be:
rankineAll = os.getcwd() + "/../Rankine20_CIV_BAL.fits"

#Rankines info file PRH selected the good cases
infoRankineparent = os.getcwd() + "/../DR16parent_DR14RankineInfo.csv"
infoRankineEHVO = os.getcwd() + "/../DR16EHVO_DR14RankineInfo_wEHVOspeed_vmax_reordered.csv" #adding BAL flag corrected version with reordered vmax


infoRankineparent_DR9 = os.getcwd() + "/../DR9parent_DR14RankineInfo.csv" #LEF: I added as these were undefined in code
infoRankineEHVO_DR9 = os.getcwd() + "/../DR9EHVO_DR14RankineInfo_wEHVOspeed_vmax_reordered.csv" #adding BAL flag corrected version with reordered vmax


#Reading RankineAll fit table for Emission Parameters
rankine_all = fits.open(rankineAll) # reading in Rankine20_CIV_BAL.fits data
data_rankAll = rankine_all[1].data
CivBlue_rankAll = data_rankAll['CIV_blue']
CivEW_rankAll = data_rankAll['CIV_EW']
HeIIBlue_rankAll = data_rankAll['HeII_blue']
HeIIEW_rankAll = data_rankAll['HeII_EW']
CivDistance_rankAll = data_rankAll['CIV_distance']
SNR_rankAll = data_rankAll['S/N']
bi_bi = data_rankAll['BI_BI']
vmin_rankAll = data_rankAll['BI_VMIN']
vmax_rankAll = data_rankAll['BI_VMAX']
good = data_rankAll['good']
LogLbol_rankAll = data_rankAll['LOG_LBOL']
sdss_name_rankAll = data_rankAll['SDSS_NAME']
z_rankAll = data_rankAll['z']
edd_rankAll = data_rankAll['LOG_REDD']
MBH_rankAll = data_rankAll['LOG_BHM_CIV']

good_rankAll=np.where(good==1)
snr_10 = np.where((SNR_rankAll > 10)&(good==1))
bal_rankAll = np.where((bi_bi > 0)&(SNR_rankAll > 10)&(good == 1))
nonbal_rankAll = np.where((bi_bi==0)&(SNR_rankAll > 10)&(good==1))

CivBlue_rankAll_SNR = CivBlue_rankAll[snr_10]
CivBlue_rankAll_good = CivBlue_rankAll[good_rankAll]
CivEW_rankAll_good = CivEW_rankAll[good_rankAll]
CivEW_rankAll_SNR = CivEW_rankAll[snr_10]
CivBlue_rankAll_bal = CivBlue_rankAll[bal_rankAll]
CivBlue_rankAll_nonbal = CivBlue_rankAll[nonbal_rankAll]
CivEW_rankAll_bal = CivEW_rankAll[bal_rankAll]
CivEW_rankAll_nonbal = CivEW_rankAll[nonbal_rankAll]
HeIIBlue_rankAll_bal = HeIIBlue_rankAll[bal_rankAll]
HeIIBlue_rankAll_nonbal = HeIIBlue_rankAll[nonbal_rankAll]
HeIIEW_rankAll_bal = HeIIEW_rankAll[bal_rankAll]
HeIIEW_rankAll_nonbal = HeIIEW_rankAll[nonbal_rankAll]
CivDistance_rankAll_bal = CivDistance_rankAll[bal_rankAll]
CivDistance_rankAll_nonbal = CivDistance_rankAll[nonbal_rankAll]
LogLbol_rankAll_bal = LogLbol_rankAll[bal_rankAll]
LogLbol_rankAll_nonbal = LogLbol_rankAll[nonbal_rankAll]
vmin_rankAll_bal = vmin_rankAll[bal_rankAll]
sdss_rankAll_bal = sdss_name_rankAll[bal_rankAll]
z_rankAll_bal = z_rankAll[bal_rankAll]
vmax_rankAll_bal = vmax_rankAll[bal_rankAll]
edd_rankAll_bal = edd_rankAll[bal_rankAll]
edd_rankAll_nonbal = edd_rankAll[nonbal_rankAll]
MBH_rankAll_bal = MBH_rankAll[bal_rankAll]
MBH_rankAll_nonbal = MBH_rankAll[nonbal_rankAll]

rankine_all = fits.open(rankineAll)
data_rankAll = rankine_all[1].data
vmin_rankAll = data_rankAll['BI_VMIN']
vmax_rankAll = data_rankAll['BI_VMAX']
bi_bi = data_rankAll['BI_BI']
good_rankAll=np.where(good==1)
bal_rankAll_25000 = np.where((bi_bi > 0)&(SNR_rankAll > 10)&(good == 1)&(vmin_rankAll>25000))
bal_rankAll_25000

#Extracting values from Rankines info files
dfRPA = pd.read_csv(infoRankineparent, header=None)
dfRHV = pd.read_csv(infoRankineEHVO, header=0)
dfRPA_9 = pd.read_csv(infoRankineparent_DR9, header=None)
dfRHV_9 = pd.read_csv(infoRankineEHVO_DR9, header=0)

#Parent sample
CivBlue_parentR=dfRPA[dfRPA.columns[5]].to_numpy()
CivEW_parentR=dfRPA[dfRPA.columns[6]].to_numpy()
CivDist_parentR=dfRPA[dfRPA.columns[14]].to_numpy()
HeiiEW_parentR=dfRPA[dfRPA.columns[7]].to_numpy()
Lbol_parentR=dfRPA[dfRPA.columns[17]].to_numpy()
vmin_parentR=dfRPA[dfRPA.columns[12]].to_numpy()
vmax_parentR=dfRPA[dfRPA.columns[11]].to_numpy()
Edd_parentR=dfRPA[dfRPA.columns[18]].to_numpy()
MBH_parentR=dfRPA[dfRPA.columns[16]].to_numpy()


CivBlue_parentR_9=dfRPA_9[dfRPA_9.columns[5]].to_numpy()
CivEW_parentR_9=dfRPA_9[dfRPA_9.columns[6]].to_numpy()
CivDist_parentR_9=dfRPA_9[dfRPA_9.columns[14]].to_numpy()
HeiiEW_parentR_9=dfRPA_9[dfRPA_9.columns[7]].to_numpy()
Lbol_parentR_9=dfRPA_9[dfRPA_9.columns[17]].to_numpy()
Edd_parentR_9=dfRPA_9[dfRPA_9.columns[18]].to_numpy()
vmin_parentR_9=dfRPA_9[dfRPA_9.columns[12]].to_numpy()
vmax_parentR_9=dfRPA_9[dfRPA_9.columns[11]].to_numpy()
MBH_parentR_9=dfRPA_9[dfRPA_9.columns[16]].to_numpy()

Lbol_parentR_combined = np.concatenate((Lbol_parentR_9, Lbol_parentR))
Edd_parentR_combined = np.concatenate((Edd_parentR_9, Edd_parentR))
vmin_parentR_combined = np.concatenate((vmin_parentR_9, vmin_parentR))
HeiiEW_parentR_combined = np.concatenate((HeiiEW_parentR_9, HeiiEW_parentR))
CivBlue_parentR_combined = np.concatenate((CivBlue_parentR_9, CivBlue_parentR))
CivEW_parentR_combined = np.concatenate((CivEW_parentR_9, CivEW_parentR))
MBH_parentR_combined = np.concatenate((MBH_parentR_9, MBH_parentR))

bi_bi0 = dfRPA[dfRPA.columns[10]].to_numpy()
bi_vmax = dfRPA[dfRPA.columns[11]].to_numpy()
bi_vmin = dfRPA[dfRPA.columns[12]].to_numpy()
SNR = dfRPA[dfRPA.columns[15]].to_numpy()
bi_bi0_9 = dfRPA_9[dfRPA_9.columns[10]].to_numpy()
bi_vmax_9 = dfRPA_9[dfRPA_9.columns[11]].to_numpy()
bi_vmin_9 = dfRPA_9[dfRPA_9.columns[12]].to_numpy()
SNR_9 = dfRPA_9[dfRPA_9.columns[15]].to_numpy()
non_bal = np.where((bi_bi0==0)&(SNR > 10))
pos_bal = np.where((bi_bi0>0) & (SNR > 10))
pos_10k = np.where((bi_vmax<10000)&(bi_bi0>0)&(SNR>10))
pos_20k = np.where((bi_vmax>=20000)&(bi_bi0>0)&(SNR>10))
pos_25k = np.where((bi_vmax>10000)&(bi_vmax<25000)&(bi_bi0>0)&(SNR>10))
non_bal_9 = np.where((bi_bi0_9==0)&(SNR_9 > 10))
pos_bal_9 = np.where((bi_bi0_9>0) & (SNR_9 > 10))
pos_10k_9 = np.where((bi_vmax_9<10000)&(bi_bi0_9>0)&(SNR_9>10))
pos_20k_9 = np.where((bi_vmax_9>10000)&(bi_vmax_9<20000)&(bi_bi0_9>0)&(SNR_9>10))
pos_25k_9 = np.where((bi_vmax_9>10000)&(bi_vmax_9<25000)&(bi_bi0_9>0)&(SNR_9>10))

CivBlue_bal = CivBlue_parentR[pos_bal]
CivBlue_bal_20k = CivBlue_parentR[pos_20k]
CivEW_bal = CivEW_parentR[pos_bal]
CivDist_bal=CivDist_parentR[pos_bal]
HeiiEW_bal = HeiiEW_parentR[pos_bal]
CivBlue_nonbal=CivBlue_parentR[non_bal]
CivEW_nonbal = CivEW_parentR[non_bal]
CivDist_nonbal=CivDist_parentR[non_bal]
HeiiEW_nonbal = HeiiEW_parentR[non_bal]
Lbol_bal = Lbol_parentR[pos_bal]
Lbol_nonbal = Lbol_parentR[non_bal]
vmax_bal = vmax_parentR[pos_bal]
vmin_bal = vmin_parentR[pos_bal]
vmax_nonbal = vmax_parentR[non_bal]
vmin_nonbal = vmin_parentR[non_bal]
Edd_bal = Edd_parentR[pos_bal]
Edd_nonbal = Edd_parentR[non_bal]
MBH_bal = MBH_parentR[pos_bal]
MBH_bal_20k = MBH_parentR[pos_20k]
MBH_nonbal = MBH_parentR[non_bal]
Lbol_bal_20k = Lbol_parentR[pos_20k]


CivBlue_bal_9 = CivBlue_parentR_9[pos_bal_9]
CivBlue_bal_20k_9 = CivBlue_parentR_9[pos_20k_9]
CivEW_bal_9 = CivEW_parentR_9[pos_bal_9]
CivDist_bal_9=CivDist_parentR_9[pos_bal_9]
HeiiEW_bal_9 = HeiiEW_parentR_9[pos_bal_9]
CivBlue_nonbal_9=CivBlue_parentR_9[non_bal_9]
CivEW_nonbal_9 = CivEW_parentR_9[non_bal_9]
CivDist_nonbal_9=CivDist_parentR_9[non_bal_9]
HeiiEW_nonbal_9 = HeiiEW_parentR_9[non_bal_9]
Lbol_bal_9 = Lbol_parentR_9[pos_bal_9]
Lbol_nonbal_9 = Lbol_parentR_9[non_bal_9]
vmax_bal_9 = vmax_parentR_9[pos_bal_9]
vmin_bal_9 = vmin_parentR_9[pos_bal_9]
vmax_nonbal_9 = vmax_parentR_9[non_bal_9]
vmin_nonbal_9 = vmin_parentR_9[non_bal_9]
Edd_bal_9 = Edd_parentR_9[pos_bal_9]
Edd_nonbal_9 = Edd_parentR_9[non_bal_9]
MBH_bal_9 = MBH_parentR_9[pos_bal_9]
MBH_bal_20k_9 = MBH_parentR_9[pos_20k_9]
MBH_nonbal_9 = MBH_parentR_9[non_bal_9]
Lbol_bal_20k_9 = Lbol_parentR_9[pos_20k_9]

vmax_bal_combined = np.concatenate((vmax_bal, vmax_bal_9))
vmin_bal_combined = np.concatenate((vmin_bal, vmin_bal_9))
Lbol_bal_combined = np.concatenate((Lbol_bal, Lbol_bal_9))
CivBlue_bal_combined = np.concatenate((CivBlue_bal, CivBlue_bal_9))
CivBlue_bal_20k_combined = np.concatenate((CivBlue_bal_20k, CivBlue_bal_20k_9))
CivEW_bal_combined = np.concatenate((CivEW_bal, CivEW_bal_9))
CivDist_bal_combined = np.concatenate((CivDist_bal, CivDist_bal_9))
HeiiEW_bal_combined = np.concatenate((HeiiEW_bal, HeiiEW_bal_9))
Lbol_bal_combined = np.concatenate((Lbol_bal, Lbol_bal_9))
Lbol_nonbal_combined = np.concatenate((Lbol_nonbal, Lbol_nonbal_9))
CivBlue_nonbal_combined = np.concatenate((CivBlue_nonbal, CivBlue_nonbal_9))
CivEW_nonbal_combined = np.concatenate((CivEW_nonbal, CivEW_nonbal_9))
CivDist_nonbal_combined = np.concatenate((CivDist_nonbal, CivDist_nonbal_9))
HeiiEW_nonbal_combined = np.concatenate((HeiiEW_nonbal, HeiiEW_nonbal_9))
Edd_bal_combined = np.concatenate((Edd_bal, Edd_bal_9))
Edd_nonbal_combined = np.concatenate((Edd_nonbal, Edd_nonbal_9))
MBH_bal_combined = np.concatenate((MBH_bal, MBH_bal_9))
MBH_nonbal_combined = np.concatenate((MBH_nonbal, MBH_nonbal_9))
MBH_bal_20k_combined = np.concatenate((MBH_bal_20k, MBH_bal_20k_9))
Lbol_bal_20k_combined = np.concatenate((Lbol_bal_20k, Lbol_bal_20k_9))

#EHVOs
CivBlue_EHVOR=dfRHV[dfRHV.columns[5]].to_numpy()
CivEW_EHVOR=dfRHV[dfRHV.columns[6]].to_numpy()
CivDist_EHVOR=dfRHV[dfRHV.columns[14]].to_numpy()
HeiiEW_EHVOR=dfRHV[dfRHV.columns[7]].to_numpy()
Lbol_EHVOR=dfRHV[dfRHV.columns[17]].to_numpy()
Edd_EHVOR=dfRHV[dfRHV.columns[18]].to_numpy()
MBH_EHVOR=dfRHV[dfRHV.columns[16]].to_numpy()

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

#plate & mjds 
plate_EHVOR = dfRHV[dfRHV.columns[1]].to_numpy()
mjd_EHVOR = dfRHV[dfRHV.columns[2]].to_numpy()

plate_EHVOR_9 = dfRHV_9[dfRHV_9.columns[1]].to_numpy()
mjd_EHVOR_9 = dfRHV_9[dfRHV_9.columns[2]].to_numpy()

#combining 9 with 16
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

copyplate_EHVO = np.copy(plate_EHVOR)
no_outlier_plate_EHVO = np.delete(copyplate_EHVO, [40])

copymjd_EHVO = np.copy(mjd_EHVOR)
no_outlier_mjd_EHVO = np.delete(copymjd_EHVO, [40])


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

plate_EHVOR_combined_noOut = np.concatenate((plate_EHVOR_9, no_outlier_plate_EHVO))
mjd_EHVOR_combined_noOut = np.concatenate((mjd_EHVOR_9, no_outlier_mjd_EHVO))


# === BAL subset of the no-outlier combined dataset ===
CivBlue_EHVO_BAL_combined_noOut = np.concatenate((CivBlue_EHVOR_9[EHVO_BAL_9], no_outlier_CivBlue[EHVO_BAL_16]))
CivEW_EHVO_BAL_combined_noOut = np.concatenate((CivEW_EHVOR_9[EHVO_BAL_9], no_outlier_CivEW[EHVO_BAL_16]))
CivDist_EHVOR_BAL_combined_noOut = np.concatenate((CivDist_EHVOR_9[EHVO_BAL_9], no_outlier_CivDist[EHVO_BAL_16]))
HeiiEW_EHVOR_BAL_combined_noOut = np.concatenate((HeiiEW_EHVOR_9[EHVO_BAL_9], no_outlier_HeiiEW[EHVO_BAL_16]))
Edd_EHVO_BAL_combined_noOut = np.concatenate((Edd_EHVOR_9[EHVO_BAL_9], no_outlier_Edd_EHVO[EHVO_BAL_16]))
Lbol_EHVOR_BAL_combined_noOut = np.concatenate((Lbol_EHVOR_9[EHVO_BAL_9], no_outlier_Lbol_EHVO[EHVO_BAL_16]))
MBH_EHVOR_BAL_combined_noOut = np.concatenate((MBH_EHVOR_9[EHVO_BAL_9], no_outlier_MBH_EHVO[EHVO_BAL_16]))
vmin_EHVOR_BAL_combined_noOut = np.concatenate((vmin_EHVO_9[EHVO_BAL_9], no_outlier_vmin_EHVO[EHVO_BAL_16]))
vmax_EHVOR_BAL_combined_noOut = np.concatenate((vmax_EHVO_9[EHVO_BAL_9], no_outlier_vmax_EHVO[EHVO_BAL_16]))


# === non-BAL subset of the no-outlier combined dataset (NEW BLOCK) ===
CivBlue_EHVO_nonBAL_combined_noOut = np.concatenate((CivBlue_EHVOR_9[EHVO_nonBAL_9], no_outlier_CivBlue[EHVO_nonBAL_16]))
CivEW_EHVVO_nonBAL_combined_noOut = np.concatenate((CivEW_EHVOR_9[EHVO_nonBAL_9], no_outlier_CivEW[EHVO_nonBAL_16]))
CivDist_EHVOR_nonBAL_combined_noOut = np.concatenate((CivDist_EHVOR_9[EHVO_nonBAL_9], no_outlier_CivDist[EHVO_nonBAL_16]))
HeiiEW_EHVOR_nonBAL_combined_noOut = np.concatenate((HeiiEW_EHVOR_9[EHVO_nonBAL_9], no_outlier_HeiiEW[EHVO_nonBAL_16]))
Edd_EHVO_nonBAL_combined_noOut = np.concatenate((Edd_EHVOR_9[EHVO_nonBAL_9], no_outlier_Edd_EHVO[EHVO_nonBAL_16]))
Lbol_EHVOR_nonBAL_combined_noOut = np.concatenate((Lbol_EHVOR_9[EHVO_nonBAL_9], no_outlier_Lbol_EHVO[EHVO_nonBAL_16]))
MBH_EHVOR_nonBAL_combined_noOut = np.concatenate((MBH_EHVOR_9[EHVO_nonBAL_9], no_outlier_MBH_EHVO[EHVO_nonBAL_16]))
vmin_EHVOR_nonBAL_combined_noOut = np.concatenate((vmin_EHVO_9[EHVO_nonBAL_9], no_outlier_vmin_EHVO[EHVO_nonBAL_16]))
vmax_EHVOR_nonBAL_combined_noOut = np.concatenate((vmax_EHVO_9[EHVO_nonBAL_9], no_outlier_vmax_EHVO[EHVO_nonBAL_16]))

# [Fig 1] CIV parameter space scatter/hist w/parent+rankine --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if plot_fig1 == True:
    figure_1(CivBlue_rankAll_good,np.log10(CivEW_rankAll_good),CivBlue_rankAll_SNR,np.log10(CivEW_rankAll_SNR),CivBlue_parentR,np.log10(CivEW_parentR),'Rankine+2020 all','Rankine+S/N>10','CS+in prep')
else:
    print('Fig 1 not plotted')
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# [Fig 2] CIV parameter space scatter/hist ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if plot_fig2 == True:
    figure_2_outline(CivBlue_rankAll_nonbal, np.log10(CivEW_rankAll_nonbal), CivBlue_rankAll_bal, np.log10(CivEW_rankAll_bal),
             CivBlue_EHVO_combined, np.log10(CivEW_EHVO_combined), 'Rankine+2020 non-BALs', 'Rankine+2020 BALS',
             'EHVO (RH+2020 & CS+in prep)', CivBlue_EHVOR_9, np.log10(CivEW_EHVOR_9))
else:
    print('Fig 2 not plotted')
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# [Fig 3] plot is made and provided by Amy Rankine

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# [Fig 4] EHVO subsamples separated by vmedian scatter/hist w/velocity colormap -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

print(f'vmin median: {np.median(vmin_EHVOR_combined)}')
print(f'vmax median: {np.median(vmax_EHVOR_combined)}')


large_speed_vmin = np.where(abs(vmin_EHVOR_combined_noOut) >= 36700)#50 quasars   #36515.8231 #NEW median: 36685.20573
small_speed_vmin = np.where(abs(vmin_EHVOR_combined_noOut) < 36700)#49 quasars

large_speed_vmax = np.where(abs(vmax_EHVOR_combined_noOut) >= 43300)#50 quasars     #43378.17873 #NEW median: 43294.439150000006
small_speed_vmax = np.where(abs(vmax_EHVOR_combined_noOut) < 43300)#49 quasars


#2D 2sample K-S Test EHVOs split by median speed comparing CIV distance vs EHVO velocity (Vmin and Vmax)
print('Samples: EHVOs split by median speed (Vmin & Vmax)     Parameters: CIV EW vs CIV Blueshift')
P, D = ndtest.ks2d2s(CivEW_EHVO_combined_noOut[large_speed_vmin], CivBlue_EHVO_combined_noOut[large_speed_vmin], CivEW_EHVO_combined_noOut[small_speed_vmin], CivBlue_EHVO_combined_noOut[small_speed_vmin], extra=True) #Vmin
print(f"{P=:.3g}, {D=:.3g}")
P, D = ndtest.ks2d2s(CivEW_EHVO_combined_noOut[large_speed_vmax], CivBlue_EHVO_combined_noOut[large_speed_vmax], CivEW_EHVO_combined_noOut[small_speed_vmax], CivBlue_EHVO_combined_noOut[small_speed_vmax], extra=True) #Vmin
print(f"{P=:.3g}, {D=:.3g}")
print()
print()

if plot_fig4 == True:
    # vmin
    plot_EHVO_CIV_scatter_hist(
        CivBlue_EHVO_combined_noOut,
        CivEW_EHVO_combined_noOut,
        vmin_EHVOR_combined_noOut,
        small_speed_vmin,
        large_speed_vmin,
        r'$\mathrm{EHVO}\ V_{\min} < 36700$',
        r'$\mathrm{EHVO}\ V_{\min} \geq 36700$',
        v_label=r'$V_{\min} [\mathrm{km\ s^{-1}}]$')
    
    # vmax
    plot_EHVO_CIV_scatter_hist(
        CivBlue_EHVO_combined_noOut,
        CivEW_EHVO_combined_noOut,
        vmax_EHVOR_combined_noOut,
        small_speed_vmax,
        large_speed_vmax,
        r'$\mathrm{EHVO}\ V_{\max} < 43300$',
        r'$\mathrm{EHVO}\ V_{\max} \geq 43300$',
        v_label=r'$V_{\max} [\mathrm{km\ s^{-1}}]$')
else:
    print('Fig 4 not plotted')
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# [Fig 5] HeII vs CIV distance scatter/hist plot  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if plot_fig5 == True:
    figure_5_outline(CivDistance_rankAll_nonbal, np.log10(np.where(HeIIEW_rankAll_nonbal > 0, HeIIEW_rankAll_nonbal, np.nan)), CivDistance_rankAll_bal, np.log10(np.where(HeIIEW_rankAll_bal > 0, HeIIEW_rankAll_bal, np.nan)),
             CivDist_EHVOR_combined_noOut, np.log10(np.where(HeiiEW_EHVOR_combined_noOut> 0, HeiiEW_EHVOR_combined_noOut, np.nan)), 'Rankine+2020 non-BALs', 'Rankine+2020 BALS',
             'EHVO (RH+2020 & CS+in prep)', CivDist_EHVOR_9, np.log10(HeiiEW_EHVOR_9))
else:
    print('Fig 5 not plotted')
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- needed cut offs when plotting in physical property parameter space ------------------------------------------------------------------
#'Good' Values, Coatman et al. 2017 cutoff for Black hole mass (and thus eddington ratio) plots
def bluecut(sample, cut=500):
    return np.where(sample > cut) 

parentR_Bluecut = bluecut(CivBlue_parentR_combined) #cuts points with CIV Blueshift less than 500 km/s



EHVOcomb_Bluecut = bluecut(CivBlue_EHVO_combined_noOut)
EHVO_bluecut_cond = np.zeros(len(CivBlue_EHVO_combined_noOut), dtype=bool)
EHVO_bluecut_cond[EHVOcomb_Bluecut] = True
#all EHVOs have CIV blueshift > 500 km/s so no condition for this needs to be applied to EHVOs when plotting
# ^^ this is because the outlier which is removed prior to here is the only EHVO case with CIV blueshift < 500 km/s
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# cutoff option along cluster separation line in HeII EW vs CIV Distance scatter/hist
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
# LEF TO DO implement np.log10(np.where(HeiiEW_EHVOR_combined_noOut> 0, HeiiEW_EHVOR_combined_noOut, np.nan)) in any plots that have HeII DONE


#print(f'EHVO plate of cases above partition in Fig 5: {plate_EHVOR_combined_noOut[HeiiEW_bisect_above]}')
#print(f'EHVO mjd of cases above partition in Fig 5: {mjd_EHVOR_combined_noOut[HeiiEW_bisect_above]}')

'''
la = np.column_stack((plate_EHVOR_combined_noOut[HeiiEW_bisect_above], mjd_EHVOR_combined_noOut[HeiiEW_bisect_above]))
df = pd.DataFrame(la, columns=["plate", "mjd"])
df.to_csv('EHVOs_HeII_CIVDist_cluster.csv', index=False)
'''

# [Fig 6] HeII colormap hex/scatter plots in physical property parameter spaces -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if plot_fig6 == True:
    
    if HeII_flat_distance_cut == True:
        condition_less = CivDist_EHVOR_combined_noOut > dist_cut
        
        plot_CIVBlue_hexbin_ehvo(
            MBH_parentR_combined[parentR_Bluecut],
            Lbol_parentR_combined[parentR_Bluecut],
            HeiiEW_parentR_combined[parentR_Bluecut],
            MBH_EHVOR_combined_noOut,
            Lbol_EHVOR_combined_noOut,
            HeiiEW_EHVOR_combined_noOut,
            condition_less,
            "log$_{10}$($M_{BH}/M_\odot$)",
            "log$_{10}$($L_{bol}$/erg s$^{-1}$)",
            xlim = (8.4,10.5),
            ylim = (46.1,48.2),
            gridsize=gridsize,
            show_left=True,
            cmap_sample='HeIIEW',
            hex_func=hexbin_function6,
            left_label=fr'EHVOs ($\mathrm{{C\,IV}}\ \mathrm{{ distance}} \leq {dist_cut}$)',
            right_label=fr'EHVOs ($\mathrm{{C\,IV}}\ \mathrm{{ distance}} > {dist_cut}$)')
        
        plot_CIVBlue_hexbin_ehvo(
            MBH_parentR_combined[parentR_Bluecut],
            Edd_parentR_combined[parentR_Bluecut],
            HeiiEW_parentR_combined[parentR_Bluecut],
            MBH_EHVOR_combined_noOut,
            Edd_EHVO_combined_noOut,
            HeiiEW_EHVOR_combined_noOut,
            condition_less,
            "log$_{10}$($M_{BH}/M_\odot$)",
            r'Eddington Ratio $\log_{10}(L_{\mathrm{bol}}/L_{\mathrm{edd}})$',
            xlim = (8.4,10.5),
            ylim = (-1.5, 0.5),
            gridsize=gridsize,
            show_left=True,
            cmap_sample='HeIIEW',
            hex_func=hexbin_function6,
            left_label=fr'EHVOs ($\mathrm{{C\,IV}}\ \mathrm{{ distance}} \leq {dist_cut}$)',
            right_label=fr'EHVOs ($\mathrm{{C\,IV}}\ \mathrm{{ distance}} > {dist_cut}$)')
        
        plot_CIVBlue_hexbin_ehvo(
            Lbol_parentR_combined[parentR_Bluecut],
            Edd_parentR_combined[parentR_Bluecut],
            HeiiEW_parentR_combined[parentR_Bluecut],
            Lbol_EHVOR_combined_noOut,
            Edd_EHVO_combined_noOut,
            HeiiEW_EHVOR_combined_noOut,
            condition_less,
            "log$_{10}$($L_{bol}$/erg s$^{-1}$)",
            r'Eddington Ratio $\log_{10}(L_{\mathrm{bol}}/L_{\mathrm{edd}})$',
            xlim = (46.1,48.2),
            ylim = (-1.5, 0.5),
            gridsize=gridsize,
            show_left=True,
            cmap_sample='HeIIEW',
            hex_func=hexbin_function6,
            left_label=fr'EHVOs ($\mathrm{{C\,IV}}\ \mathrm{{ distance}} \leq {dist_cut}$)',
            right_label=fr'EHVOs ($\mathrm{{C\,IV}}\ \mathrm{{ distance}} > {dist_cut}$)')
    
    elif HeII_flat_distance_cut == 'line':
        
        #np.arange(len(HeiiEW_EHVOR_combined_noOut)) creates the indices of the data array.
        HeII_bisect_condition_more = np.isin(np.arange(len(HeiiEW_EHVOR_combined_noOut)), HeiiEW_bisect_below)
        HeII_bisect_condition_less = np.isin(np.arange(len(HeiiEW_EHVOR_combined_noOut)), HeiiEW_bisect_above)
    
    
        plot_CIVBlue_hexbin_ehvo(
            MBH_parentR_combined[parentR_Bluecut],
            Lbol_parentR_combined[parentR_Bluecut],
            HeiiEW_parentR_combined[parentR_Bluecut],
            MBH_EHVOR_combined_noOut,
            Lbol_EHVOR_combined_noOut,
            HeiiEW_EHVOR_combined_noOut,
            HeII_bisect_condition_more,
            "log$_{10}$($M_{BH}/M_\odot$)",
            "log$_{10}$($L_{bol}$/erg s$^{-1}$)",
            xlim = (8.4,10.5),
            ylim = (46.1,48.2),
            gridsize=gridsize,
            show_left=True,
            cmap_sample='HeIIEW',
            hex_func=hexbin_function6,
            left_label='EHVOs Above Partition Line',
            right_label='EHVOs Below Partition Line')
        
        plot_CIVBlue_hexbin_ehvo(
            MBH_parentR_combined[parentR_Bluecut],
            Edd_parentR_combined[parentR_Bluecut],
            HeiiEW_parentR_combined[parentR_Bluecut],
            MBH_EHVOR_combined_noOut,
            Edd_EHVO_combined_noOut,
            HeiiEW_EHVOR_combined_noOut,
            HeII_bisect_condition_more,
            "log$_{10}$($M_{BH}/M_\odot$)",
            r'Eddington Ratio $\log_{10}(L_{\mathrm{bol}}/L_{\mathrm{edd}})$',
            xlim = (8.4,10.5),
            ylim = (-1.5, 0.5),
            gridsize=gridsize,
            show_left=True,
            cmap_sample='HeIIEW',
            hex_func=hexbin_function6,
            left_label='EHVOs Above Partition Line',
            right_label='EHVOs Below Partition Line')
        
        plot_CIVBlue_hexbin_ehvo(
            Lbol_parentR_combined[parentR_Bluecut],
            Edd_parentR_combined[parentR_Bluecut],
            HeiiEW_parentR_combined[parentR_Bluecut],
            Lbol_EHVOR_combined_noOut,
            Edd_EHVO_combined_noOut,
            HeiiEW_EHVOR_combined_noOut,
            HeII_bisect_condition_more,
            "log$_{10}$($L_{bol}$/erg s$^{-1}$)",
            r'Eddington Ratio $\log_{10}(L_{\mathrm{bol}}/L_{\mathrm{edd}})$',
            xlim = (46.1,48.2),
            ylim = (-1.5, 0.5),
            gridsize=gridsize,
            show_left=True,
            cmap_sample='HeIIEW',
            hex_func=hexbin_function6,
            left_label='EHVOs Above Partition Line',
            right_label='EHVOs Below Partition Line')
else:
    print('Fig 6 not plotted')

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# [Fig 8] CIV blueshift colormap hex/scatter plots in physical property parameter spaces -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if plot_fig8 == True:
    if Blueshift_cut == 'line':
        Blue_bisect_condition_more = np.isin(np.arange(len(CivBlue_EHVO_combined_noOut)), HeiiEW_bisect_below)
        
        
        # ---- double panels ----------------------------
        # MBH vs Lbol
        plot_CIVBlue_hexbin_ehvo(
            MBH_parentR_combined[parentR_Bluecut],
            Lbol_parentR_combined[parentR_Bluecut],
            CivBlue_parentR_combined[parentR_Bluecut],
            MBH_EHVOR_combined_noOut,
            Lbol_EHVOR_combined_noOut,
            CivBlue_EHVO_combined_noOut,
            Blue_bisect_condition_more,
            "log$_{10}$($M_{BH}/M_\odot$)",
            "log$_{10}$($L_{bol}$/erg s$^{-1}$)",
            xlim = (8.4,10.5),
            ylim = (46.1,48.2),
            gridsize=gridsize,
            show_left=True,
            hex_func=hexbin_function8,
            left_label='EHVOs Above Partition Line',
            right_label='EHVOs Below Partition Line')
        
        # MBH vs Edd
        plot_CIVBlue_hexbin_ehvo(
            MBH_parentR_combined[parentR_Bluecut],
            Edd_parentR_combined[parentR_Bluecut],
            CivBlue_parentR_combined[parentR_Bluecut],
            MBH_EHVOR_combined_noOut,
            Edd_EHVO_combined_noOut,
            CivBlue_EHVO_combined_noOut,
            Blue_bisect_condition_more,
            "log$_{10}$($M_{\mathrm{BH}}/M_{\odot})$",
            'Eddington Ratio (log$_{10}$($L_{bol}$/$L_{Edd}$)',
            xlim = (8.4,10.5),
            ylim = (-1.5, 0.5),
            gridsize=gridsize,
            show_left=True,
            hex_func=hexbin_function8,
            left_label='EHVOs Above Partition Line',
            right_label='EHVOs Below Partition Line')
        
        # Lbol vs EDD
        plot_CIVBlue_hexbin_ehvo(
            Lbol_parentR_combined[parentR_Bluecut],
            Edd_parentR_combined[parentR_Bluecut],
            CivBlue_parentR_combined[parentR_Bluecut],
            Lbol_EHVOR_combined_noOut,
            Edd_EHVO_combined_noOut,
            CivBlue_EHVO_combined_noOut,
            Blue_bisect_condition_more,
            "log$_{10}$($L_{bol}$/erg s$^{-1}$)",
            'Eddington Ratio (log$_{10}$($L_{bol}$/$L_{Edd}$)',
            xlim = (46.1,48.2),
            ylim = (-1.5, 0.5),
            gridsize=gridsize,
            show_left=True,
            hex_func=hexbin_function8,
            left_label='EHVOs Above Partition Line',
            right_label='EHVOs Below Partition Line')
    
    elif Blueshift_cut == 'flat':
        condition_less = CivDist_EHVOR_combined_noOut > dist_cut2
    
        # ---- double panels ----------------------------
        # MBH vs Lbol
        plot_CIVBlue_hexbin_ehvo(
            MBH_parentR_combined[parentR_Bluecut],
            Lbol_parentR_combined[parentR_Bluecut],
            CivBlue_parentR_combined[parentR_Bluecut],
            MBH_EHVOR_combined_noOut,
            Lbol_EHVOR_combined_noOut,
            CivBlue_EHVO_combined_noOut,
            condition_less,
            "log$_{10}$($M_{BH}/M_\odot$)",
            "log$_{10}$($L_{bol}$/erg s$^{-1}$)",
            xlim = (8.4,10.5),
            ylim = (46.1,48.2),
            gridsize=gridsize,
            show_left=True,
            hex_func=hexbin_function8,
            left_label=fr'EHVOs ($\mathrm{{C\,IV}}\ \mathrm{{ distance}} \leq {dist_cut2}$)',
            right_label=fr'EHVOs ($\mathrm{{C\,IV}}\ \mathrm{{ distance}} > {dist_cut2}$)')
        
        # MBH vs Edd
        plot_CIVBlue_hexbin_ehvo(
            MBH_parentR_combined[parentR_Bluecut],
            Edd_parentR_combined[parentR_Bluecut],
            CivBlue_parentR_combined[parentR_Bluecut],
            MBH_EHVOR_combined_noOut,
            Edd_EHVO_combined_noOut,
            CivBlue_EHVO_combined_noOut,
            condition_less,
            "log$_{10}$($M_{\mathrm{BH}}/M_{\odot})$",
            'Eddington Ratio (log$_{10}$($L_{bol}$/$L_{Edd}$)',
            xlim = (8.4,10.5),
            ylim = (-1.5, 0.5),
            gridsize=gridsize,
            show_left=True,
            hex_func=hexbin_function8,
            left_label=fr'EHVOs ($\mathrm{{C\,IV}}\ \mathrm{{ distance}} \leq {dist_cut2}$)',
            right_label=fr'EHVOs ($\mathrm{{C\,IV}}\ \mathrm{{ distance}} > {dist_cut2}$)')
        
        # Lbol vs EDD
        plot_CIVBlue_hexbin_ehvo(
            Lbol_parentR_combined[parentR_Bluecut],
            Edd_parentR_combined[parentR_Bluecut],
            CivBlue_parentR_combined[parentR_Bluecut],
            Lbol_EHVOR_combined_noOut,
            Edd_EHVO_combined_noOut,
            CivBlue_EHVO_combined_noOut,
            condition_less,
            "log$_{10}$($L_{bol}$/erg s$^{-1}$)",
            'Eddington Ratio (log$_{10}$($L_{bol}$/$L_{Edd}$)',
            xlim = (46.1,48.2),
            ylim = (-1.5, 0.5),
            gridsize=gridsize,
            show_left=True,
            hex_func=hexbin_function8,
            left_label=fr'EHVOs ($\mathrm{{C\,IV}}\ \mathrm{{distance}} \leq {dist_cut2}$)',
            right_label=fr'EHVOs ($\mathrm{{C\,IV}}\ \mathrm{{ distance}} > {dist_cut2}$)')
else:
    print('Fig 8 not plotted')
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
