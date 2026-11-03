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
from scipy.stats import pearsonr, spearmanr, kendalltau
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
from astropy.io import fits
from astropy import stats
import re
import ndtest #pip installed from here: https://github.com/syrte/ndtest/blob/master/README.md
import matplotlib as mpl
from CIVfunctions import project,CIV_distance
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
#plot_fig5 = False

# [Fig 6] HeII colormap hex/scatter plot in physical property parameter spaces
plot_fig6 = True
plot_fig6 = False

HeII_flat_distance_cut = True #plotting with a flat CIV Distance cutoff provided by dist_cut
dist_cut = 0.83

HeII_flat_distance_cut = 'line' #plotting with sample separation via line spearation of cluster in Fig 5 (HeII vs CIV Distance scat/hist)
HeII_flat_distance_cut = 'line w/parent' #When doing this you MUST set plot_fig5 = True
dist_cut = None



hexbin_function6 = np.median
#hexbin_function6 = np.max

gridsize = 35 # this applies to both fig 6 and 8 as they should have the same gridsize to be easily compared

# [Fig 7] corner plot found in CornerPlotsCIVBlueshift.py
gridsize = 35

# [Fig 8] CIV Blueshift colormap hex/scatter plot in physical property parameter spaces
plot_fig8 = True
plot_fig8 = False

Blueshift_cut = 'flat'
dist_cut2 = 0.9

#Blueshift_cut = 'line'
#dist_cut2 = None

hexbin_function8 = np.median
hexbin_function8 = np.max

#Plot testing plots
plot_tests = True
#plot_tests = False

##############################################

#The inputs in this program should be:
rankineAll = os.getcwd() + "/../Rankine20_CIV_BAL.fits"

#Rankines info file PRH selected the good cases
infoRankineparent = os.getcwd() + "/../DR16parent_DR14RankineInfo.csv"
infoRankineEHVO = os.getcwd() + "/../DR16EHVO_DR14RankineInfo_wEHVOspeed_vmax_reordered.csv" #adding BAL flag corrected version with reordered vmax


infoRankineparent_DR9 = os.getcwd() + "/../DR9parent_DR14RankineInfo.csv" #LEF: I added as these were undefined in code
infoRankineEHVO_DR9 = os.getcwd() + "/../DR9EHVO_DR14RankineInfo_wEHVOspeed_vmax_reordered.csv" #adding BAL flag corrected version with reordered vmax


infoBIEHVO = "https://raw.githubusercontent.com/paolaUWB/DR16Q/Liliana/CROSS_CORRELATION/EHVO_DR16_list_wBIspeed_reordered.csv"
infoBIEHVO_DR9 = "https://raw.githubusercontent.com/paolaUWB/DR16Q/Liliana/CROSS_CORRELATION/DR9Q_EHVO_list_amy_wBIspeed_reordered.csv"


#Extracting values 
dfBI = pd.read_csv(infoBIEHVO, header=0)
dfBI_9 = pd.read_csv(infoBIEHVO_DR9, header=0)

#Reading BI_EHVO Values
EHVOBI_total=dfBI[dfBI.columns[4]].to_numpy()
copyEHVOBI_total = np.copy(EHVOBI_total)
no_outlier_EHVOBI_total = np.delete(copyEHVOBI_total, [40])
EHVOBI_total_9=dfBI_9[dfBI_9.columns[8]].to_numpy()
EHVOBI_total_all = np.concatenate([EHVOBI_total_9, no_outlier_EHVOBI_total])

#individual BI values as there are cases with multiple EHVOs
'''
EHVOBI_indv=dfBI[dfBI.columns[7]].to_numpy()
EHVOBI_indv_9=dfBI_9[dfBI_9.columns[]].to_numpy()
EHVOBI_indv_all = np.concatenate([EHVOBI_indv_9, EHVOBI_indv])
print(EHVOBI_all)
stop
'''

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
BI_rankAll_bal = bi_bi[bal_rankAll]
BI_rankAll_nonbal = bi_bi[nonbal_rankAll]

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
CivDist_parentR_combined = np.concatenate((CivDist_parentR_9, CivDist_parentR))


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
vmin_EHVO_rel=dfRHV[dfRHV.columns[31]].to_numpy() #what is rel and why do we have these?
vmin_EHVO_9_rel=dfRHV_9[dfRHV_9.columns[31]].to_numpy() # ^^^
vmax_EHVO=dfRHV[dfRHV.columns[33]].to_numpy() #double checked and col 33 is correct for the DR16 EHVO sample csv containing speeds 
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
print(f'CIV EW of outlier: {CivEW_EHVOR[40]}')

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


#'Good' Values, Coatman et al. 2017 cutoff for Black hole mass (and thus eddington ratio) plots
def bluecut(sample, cut=500):
    return np.where(sample > cut) 

parentR_Bluecut = bluecut(CivBlue_parentR_combined) #cuts points with CIV Blueshift less than 500 km/s




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


#1D 2sample K-S Test EHVOs split by median speed comparing CIV distance vs EHVO velocity (Vmin and Vmax)
print('1D K-S test CIV Distance')
print('Samples: EHVOs split by median speed (Vmin & Vmax)     Parameters: CIV Distance')
ks_stat, p_value = ks_2samp(CivDist_EHVOR_combined_noOut[small_speed_vmin], CivDist_EHVOR_combined_noOut[large_speed_vmin], alternative="two-sided", method="auto") #Vmin
print(f"Vmin: {ks_stat=:.3g}, {p_value=:.3g}")
ks_stat, p_value = ks_2samp(CivDist_EHVOR_combined_noOut[small_speed_vmax], CivDist_EHVOR_combined_noOut[large_speed_vmax], alternative="two-sided", method="auto") #Vmax
print(f"Vmax: {ks_stat=:.3g}, {p_value=:.3g}")
print()
print()

#----------- Abby original Statistical testing for Table 3 "Bootstrap" test ----------------------------
# copied at the end of this notebook https://colab.research.google.com/drive/1zmikhb-BpbKL4EbkqeTjNU8lmRVgASKo?usp=sharing
'''
# replacement NO
# Select 100 random indices from the range of the BAL sample size
num_samples = 100
total_samples = len(CivBlue_rankAll_bal_10000)
random_indices = np.random.choice(total_samples, num_samples, replace=False)

# Extract the corresponding entries using the random indices
CivBlue_sampled = CivBlue_rankAll_bal_10000[random_indices]
vmin_sampled = vmin_rankAll_bal_10000[random_indices]
vmax_sampled = vmax_rankAll_bal_10000[random_indices]

# Calculate Pearson correlation coefficients and p-values
corr_vmin, p_value_vmin = pearsonr(CivBlue_sampled, vmin_sampled)
corr_vmax, p_value_vmax = pearsonr(CivBlue_sampled, vmax_sampled)

print(f'Correlation between CivBlue and vmin: {corr_vmin:.3f}, p-value: {p_value_vmin:.3e}')
print(f'Correlation between CivBlue and vmax: {corr_vmax:.3f}, p-value: {p_value_vmax:.3e}')

# Plot the correlations
plt.figure(figsize=(12, 5))

# Plot CivBlue vs vmin
plt.subplot(1, 2, 1)
plt.scatter(CivBlue_sampled, vmin_sampled, color='purple')
plt.xlabel('C IV Blueshift (km/s)')
plt.ylabel('vmin (km/s)')
plt.title(f'C IV Blueshift vs vmin\nCorrelation: {corr_vmin:.3f}, p-value: {p_value_vmin:.3e}')

# Plot CivBlue vs vmax
plt.subplot(1, 2, 2)
plt.scatter(CivBlue_sampled, vmax_sampled, color='green')
plt.xlabel('C IV Blueshift (km/s)')
plt.ylabel('vmax (km/s)')
plt.title(f'C IV Blueshift vs vmax\nCorrelation: {corr_vmax:.3f}, p-value: {p_value_vmax:.3e}')

plt.tight_layout()
plt.show()


# replacement YES
# i believe this was the version initially used in the paper and presented as a Bootstrap test
from scipy.stats import pearsonr
# Define the number of bootstrap samples
num_bootstrap_samples = 100

# Initialize arrays to store the correlation coefficients for each bootstrap sample
corr_vmin_list = []
corr_vmax_list = []

# Perform bootstrap resampling
for _ in range(num_bootstrap_samples):
    # Randomly select indices with replacement
    random_indices = np.random.choice(len(CivBlue_rankAll_bal_10000), len(CivBlue_rankAll_bal_10000), replace=True)

    # Extract the corresponding entries using the random indices
    CivBlue_sampled = CivBlue_rankAll_bal_10000[random_indices]
    vmin_sampled = vmin_rankAll_bal_10000[random_indices]
    vmax_sampled = vmax_rankAll_bal_10000[random_indices]

    # Calculate Pearson correlation coefficients for each resample
    corr_vmin, _ = pearsonr(CivBlue_sampled, vmin_sampled)
    corr_vmax, _ = pearsonr(CivBlue_sampled, vmax_sampled)

    # Append correlation coefficients to the lists
    corr_vmin_list.append(corr_vmin)
    corr_vmax_list.append(corr_vmax)

# Convert lists to arrays for easier analysis
corr_vmin_array = np.array(corr_vmin_list)
corr_vmax_array = np.array(corr_vmax_list)

# Calculate summary statistics
mean_corr_vmin = np.mean(corr_vmin_array)
std_corr_vmin = np.std(corr_vmin_array)
ci_vmin = np.percentile(corr_vmin_array, [2.5, 97.5])  # 95% CI

mean_corr_vmax = np.mean(corr_vmax_array)
std_corr_vmax = np.std(corr_vmax_array)
ci_vmax = np.percentile(corr_vmax_array, [2.5, 97.5])  # 95% CI

# Print the results
print(f'vmin vs C IV Blueshift: Mean correlation = {mean_corr_vmin:.3f}, Std = {std_corr_vmin:.3f}, 95% CI = {ci_vmin}')
print(f'vmax vs C IV Blueshift: Mean correlation = {mean_corr_vmax:.3f}, Std = {std_corr_vmax:.3f}, 95% CI = {ci_vmax}')

# Plot histograms of the bootstrap correlation coefficients
plt.figure(figsize=(14, 6))

# Histogram for vmin correlations
plt.subplot(1, 2, 1)
plt.hist(corr_vmin_array, bins=15, color='purple', alpha=0.7, edgecolor='black')
plt.axvline(mean_corr_vmin, color='blue', linestyle='dashed', linewidth=1)
plt.axvline(ci_vmin[0], color='red', linestyle='dashed', linewidth=1)
plt.axvline(ci_vmin[1], color='red', linestyle='dashed', linewidth=1)
plt.xlabel('Correlation Coefficient (vmin vs C IV Blueshift)')
plt.ylabel('Frequency')
plt.title(f'Bootstrap Distribution of Correlations\nMean: {mean_corr_vmin:.3f}, 95% CI: [{ci_vmin[0]:.3f}, {ci_vmin[1]:.3f}]')

# Histogram for vmax correlations
plt.subplot(1, 2, 2)
plt.hist(corr_vmax_array, bins=15, color='green', alpha=0.7, edgecolor='black')
plt.axvline(mean_corr_vmax, color='blue', linestyle='dashed', linewidth=1)
plt.axvline(ci_vmax[0], color='red', linestyle='dashed', linewidth=1)
plt.axvline(ci_vmax[1], color='red', linestyle='dashed', linewidth=1)
plt.xlabel('Correlation Coefficient (vmax vs C IV Blueshift)')
plt.ylabel('Frequency')
plt.title(f'Bootstrap Distribution of Correlations\nMean: {mean_corr_vmax:.3f}, 95% CI: [{ci_vmax[0]:.3f}, {ci_vmax[1]:.3f}]')

plt.tight_layout()
plt.show()
'''

# Table 3 redo testing ------------------------------------------------------------
# - I will run with 99 BALs selected to better compare with the sample of 99 EHVOs
# - I will run both with Pearson and Spearman correlation
# - maybe run with Kendall's Tau [More resarch on this..]
# - I am leaning towards no replacement: This means that one quasar cannot be randomly selected twice 
#   in the same resampling group of 99 quasars
'''
def statistical_correlations(
    CivBlue,
    vmin,
    vmax,
    parameter,
    num_samples=100,
    resample_size=None,
    replace=True,
    correlation_type='Pearson',
    EHVO_ref=None
):
    # Set the resampling size (defaults to the original sample size)
    if resample_size is None:
        resample_size = len(CivBlue)

    # Check that all arrays have the same length
    if not (len(CivBlue) == len(vmin) == len(vmax)):
        raise ValueError("***All input arrays must have the same length***")

    # Initialize arrays to store the correlation coefficients
    corr_vmin_list = []
    corr_vmax_list = []
    med_compare = []
    
    # Perform resampling
    for _ in range(num_samples):

        # Randomly select indices with or without replacement
        random_indices = np.random.choice(
            len(CivBlue),
            size=resample_size,
            replace=replace
        )

        # Extract the corresponding entries using the random indices
        CivBlue_sampled = CivBlue[random_indices]
        vmin_sampled = vmin[random_indices]
        vmax_sampled = vmax[random_indices]
        
        if correlation_type.lower() == 'pearson':
            # Calculate Pearson correlation coefficients
            corr_vmin, _ = pearsonr(CivBlue_sampled, vmin_sampled)
            corr_vmax, _ = pearsonr(CivBlue_sampled, vmax_sampled)
        elif correlation_type.lower() == 'spearman':
            corr_vmin, _ = spearmanr(CivBlue_sampled, vmin_sampled)
            corr_vmax, _ = spearmanr(CivBlue_sampled, vmax_sampled)
        elif correlation_type.lower() == 'kendall':
            corr_vmin, _ = kendalltau(CivBlue_sampled, vmin_sampled)
            corr_vmax, _ = kendalltau(CivBlue_sampled, vmax_sampled)
        elif correlation_type.lower() == 'median':
            med = np.median(CivBlue_sampled)
            if EHVO_ref is not None:
                l = med>np.median(EHVO_ref)
                #print(f'{med} vs {np.median(EHVO_ref)}') #used to check if working properly
                med_compare.append(l)
                corr_vmin = np.nan
                corr_vmax = np.nan
                
            else:
                print('***Must provide reference sample for median comparison***')

        else:
            print('***Please specify a correlation type: "Pearson" or "Spearman"***')
        
        
        # Append correlation coefficients to the lists
        corr_vmin_list.append(corr_vmin)
        corr_vmax_list.append(corr_vmax)
    
    if correlation_type.lower() == 'median':
        print(
            f'Median {parameter} resample vs Reference sample median, n = {num_samples}, resample size = {resample_size}'
            f'\nNumber of times resample {parameter} median > reference sample median = {sum(med_compare)}')
        return 

    # Convert lists to arrays for easier analysis
    corr_vmin_array = np.array(corr_vmin_list)
    corr_vmax_array = np.array(corr_vmax_list)

    # Calculate summary statistics
    mean_corr_vmin = np.mean(corr_vmin_array)
    std_corr_vmin = np.std(corr_vmin_array)
    ci_vmin = np.percentile(corr_vmin_array, [2.5, 97.5])

    mean_corr_vmax = np.mean(corr_vmax_array)
    std_corr_vmax = np.std(corr_vmax_array)
    ci_vmax = np.percentile(corr_vmax_array, [2.5, 97.5])

    # Print the results
    print()
    print(f'BAL {parameter} vs Speed Statistics')
    print(f'Conditions: \n - number of resamplings = {num_samples} \n - resample size = {resample_size} \n - replacement = {replace} \n - correlation type = {correlation_type}')
    print(
        f'vmin vs {parameter}: Mean correlation = '
        f'{mean_corr_vmin:.3f}, Std = {std_corr_vmin:.3f}, '
        f'95% CI = {ci_vmin}'
    )

    print(
        f'vmax vs {parameter}: Mean correlation = '
        f'{mean_corr_vmax:.3f}, Std = {std_corr_vmax:.3f}, '
        f'95% CI = {ci_vmax}'
    )

    # Return results for further analysis
    return {
        "corr_vmin_array": corr_vmin_array,
        "corr_vmax_array": corr_vmax_array,
        "mean_corr_vmin": mean_corr_vmin,
        "std_corr_vmin": std_corr_vmin,
        "ci_vmin": ci_vmin,
        "mean_corr_vmax": mean_corr_vmax,
        "std_corr_vmax": std_corr_vmax,
        "ci_vmax": ci_vmax
    }
'''

def statistical_correlations(
    parameter_values,
    vmin=None,
    vmax=None,
    parameter_name='parameter',
    n_resamples=100,
    sample_size=None,
    with_replacement=True,
    test_type='Pearson',
    reference_sample=None
):
    # Set the resampling size (defaults to the original sample size)
    if sample_size is None:
        sample_size = len(parameter_values)

    # Normalize and validate the test type
    test_type = test_type.lower()

    valid_tests = ['pearson', 'spearman', 'kendall', 'median']
    if test_type not in valid_tests:
        raise ValueError(
            f"Invalid test_type '{test_type}'. "
            f"Choose from {valid_tests}."
        )

    # Median comparisons require a reference sample
    if test_type == 'median':
        if reference_sample is None:
            raise ValueError(
                "***Must provide reference_sample for median comparison***"
            )

    # Correlation tests require vmin and vmax
    else:
        if vmin is None or vmax is None:
            raise ValueError(
                "***vmin and vmax must be provided for correlation tests***"
            )

        # Check that all arrays have the same length
        if not (
            len(parameter_values) == len(vmin) == len(vmax)
        ):
            raise ValueError(
                "***All input arrays must have the same length***"
            )

    # Initialize arrays to store the correlation coefficients
    corr_vmin_list = []
    corr_vmax_list = []
    med_compare = []

    # Perform resampling
    for _ in range(n_resamples):

        # Randomly select indices with or without replacement
        random_indices = np.random.choice(
            len(parameter_values),
            size=sample_size,
            replace=with_replacement
        )

        # Extract the corresponding entries using the random indices
        parameter_sampled = parameter_values[random_indices]

        if test_type == 'median':
            # Compare the resampled median with the reference median
            sampled_median = np.median(parameter_sampled)
            reference_median = np.median(reference_sample)

            med_compare.append(sampled_median > reference_median)

        else:
            vmin_sampled = vmin[random_indices]
            vmax_sampled = vmax[random_indices]

            if test_type == 'pearson':
                corr_vmin, _ = pearsonr(
                    parameter_sampled, vmin_sampled
                )
                corr_vmax, _ = pearsonr(
                    parameter_sampled, vmax_sampled
                )

            elif test_type == 'spearman':
                corr_vmin, _ = spearmanr(
                    parameter_sampled, vmin_sampled
                )
                corr_vmax, _ = spearmanr(
                    parameter_sampled, vmax_sampled
                )

            elif test_type == 'kendall':
                corr_vmin, _ = kendalltau(
                    parameter_sampled, vmin_sampled
                )
                corr_vmax, _ = kendalltau(
                    parameter_sampled, vmax_sampled
                )

            # Append correlation coefficients
            corr_vmin_list.append(corr_vmin)
            corr_vmax_list.append(corr_vmax)

    # Handle median comparison separately
    if test_type == 'median':
        n_greater = sum(med_compare)

        print()
        print(f'Median comparison: {parameter_name}')
        print(
            f'Test: Resampled {parameter_name} median > '
            f'reference sample median'
        )
        print(
            f'Conditions: n_resamples = {n_resamples}, '
            f'sample_size = {sample_size}, '
            f'with_replacement = {with_replacement}'
        )
        print(
            f'Original sample size = {len(parameter_values)}, '
            f'Reference sample size = {len(reference_sample)}'
        )
        print(f'Reference median = {np.median(reference_sample):.3f}')
        print(
            f'Resampled median > reference median: '
            f'{n_greater}/{n_resamples} '
            f'({100 * n_greater / n_resamples:.1f}%)'
        )
        print(
            f'Resampled median <= reference median: '
            f'{n_resamples - n_greater}/{n_resamples}'
        )
        return

    # Convert lists to arrays for easier analysis
    corr_vmin_array = np.array(corr_vmin_list)
    corr_vmax_array = np.array(corr_vmax_list)

    # Calculate summary statistics
    mean_corr_vmin = np.mean(corr_vmin_array)
    std_corr_vmin = np.std(corr_vmin_array)
    ci_vmin = np.percentile(corr_vmin_array, [2.5, 97.5])

    mean_corr_vmax = np.mean(corr_vmax_array)
    std_corr_vmax = np.std(corr_vmax_array)
    ci_vmax = np.percentile(corr_vmax_array, [2.5, 97.5])

    # Print the results
    print()
    print(f'Statistical correlation analysis: {parameter_name}')
    print(f'Test type = {test_type.capitalize()} correlation')
    print(
        f'Conditions: n_resamples = {n_resamples}, '
        f'sample_size = {sample_size}, '
        f'with_replacement = {with_replacement}'
    )
    print(f'Original sample size = {len(parameter_values)}')

    print(
        f'vmin vs {parameter_name}: Mean correlation = '
        f'{mean_corr_vmin:.3f}, Std = {std_corr_vmin:.3f}, '
        f'95% CI = [{ci_vmin[0]:.3f}, {ci_vmin[1]:.3f}]'
    )

    print(
        f'vmax vs {parameter_name}: Mean correlation = '
        f'{mean_corr_vmax:.3f}, Std = {std_corr_vmax:.3f}, '
        f'95% CI = [{ci_vmax[0]:.3f}, {ci_vmax[1]:.3f}]'
    )

    # Return results for further analysis
    return {
        "corr_vmin_array": corr_vmin_array,
        "corr_vmax_array": corr_vmax_array,
        "mean_corr_vmin": mean_corr_vmin,
        "std_corr_vmin": std_corr_vmin,
        "ci_vmin": ci_vmin,
        "mean_corr_vmax": mean_corr_vmax,
        "std_corr_vmax": std_corr_vmax,
        "ci_vmax": ci_vmax
    }



statistical_correlations(
    CivDistance_rankAll_bal,
    vmin_rankAll_bal,
    vmax_rankAll_bal,
    'C IV Distance',
    n_resamples=100,
    sample_size=99,
    with_replacement=False,
    test_type='kendall'
)

statistical_correlations(
    CivDistance_rankAll_bal,
    vmin_rankAll_bal,
    vmax_rankAll_bal,
    'C IV Distance',
    n_resamples=100,
    sample_size=99,
    with_replacement=False,
    test_type='Pearson'
)

statistical_correlations(
    CivDistance_rankAll_bal,
    vmin_rankAll_bal,
    vmax_rankAll_bal,
    'C IV Distance',
    n_resamples=100,
    sample_size=99,
    with_replacement=False,
    test_type='Spearman'
)

statistical_correlations(
    CivDistance_rankAll_bal,
    vmin_rankAll_bal,
    vmax_rankAll_bal,
    'C IV Distance',
    n_resamples=100,
    sample_size=99,
    with_replacement=True,
    test_type='Pearson'
)

statistical_correlations(
    CivDistance_rankAll_bal,
    vmin_rankAll_bal,
    vmax_rankAll_bal,
    'C IV Distance',
    n_resamples=100,
    sample_size=99,
    with_replacement=True,
    test_type='Spearman'
)

print()
print('******* C IV Blueshift *********')
print('********************************')
ok = np.where(vmin_rankAll_bal>10000)[0]
statistical_correlations(
    CivBlue_rankAll_bal[ok],
    vmin_rankAll_bal[ok],
    vmax_rankAll_bal[ok],
    'C IV Blueshift',
    n_resamples=100,
    sample_size=99,
    with_replacement=False,
    test_type='Kendall'
)
stop
statistical_correlations(
    CivBlue_rankAll_bal,
    vmin_rankAll_bal,
    vmax_rankAll_bal,
    'C IV Blueshift',
    n_resamples=100,
    sample_size=99,
    with_replacement=False,
    test_type='Spearman'
)

statistical_correlations(
    CivBlue_rankAll_bal,
    vmin_rankAll_bal,
    vmax_rankAll_bal,
    'C IV Blueshift',
    n_resamples=100,
    sample_size=99,
    with_replacement=True,
    test_type='Pearson'
)

statistical_correlations(
    CivBlue_rankAll_bal,
    vmin_rankAll_bal,
    vmax_rankAll_bal,
    'C IV Blueshift',
    n_resamples=100,
    sample_size=99,
    with_replacement=True,
    test_type='Spearman'
)

print('TESTING MEDIAN BLUESHIFT RESAMPLING')
print('non-BAL')
statistical_correlations(
    CivBlue_rankAll_nonbal,
    parameter_name='C IV Blueshift',
    n_resamples=10000,
    sample_size=99,
    with_replacement=False,
    test_type='median',
    reference_sample=CivBlue_EHVO_combined_noOut
)

print('BAL')
statistical_correlations(
    CivBlue_rankAll_bal,
    parameter_name='C IV Blueshift',
    n_resamples=10000,
    sample_size=99,
    with_replacement=False,
    test_type='median',
    reference_sample=CivBlue_EHVO_combined_noOut
)

la = np.where(vmax_rankAll_bal > 20000)[0]

print('BAL vmax > 20,000')
statistical_correlations(
    CivBlue_rankAll_bal[la],
    parameter_name='C IV Blueshift (BAL vmax > 20,000)',
    n_resamples=10000,
    sample_size=99,
    with_replacement=False,
    test_type='median',
    reference_sample=CivBlue_EHVO_combined_noOut
)

print('*** WITH REPLACEMENT***')
print('TESTING MEDIAN BLUESHIFT RESAMPLING')
print('non-BAL')
statistical_correlations(
    CivBlue_rankAll_nonbal,
    parameter_name='C IV Blueshift',
    n_resamples=10000,
    sample_size=99,
    with_replacement=True,
    test_type='median',
    reference_sample=CivBlue_EHVO_combined_noOut
)

print('BAL')
statistical_correlations(
    CivBlue_rankAll_bal,
    parameter_name='C IV Blueshift',
    n_resamples=10000,
    sample_size=99,
    with_replacement=True,
    test_type='median',
    reference_sample=CivBlue_EHVO_combined_noOut
)

la = np.where(vmax_rankAll_bal > 20000)[0]

print('BAL vmax > 20,000')
statistical_correlations(
    CivBlue_rankAll_bal[la],
    parameter_name='C IV Blueshift (BAL vmax > 20,000)',
    n_resamples=10000,
    sample_size=99,
    with_replacement=True,
    test_type='median',
    reference_sample=CivBlue_EHVO_combined_noOut
)

print('all Rankine+2020')
statistical_correlations(
    CivBlue_rankAll,
    parameter_name='C IV Blueshift all Rankine+2020',
    n_resamples=10000,
    sample_size=99,
    with_replacement=False,
    test_type='median',
    reference_sample=CivBlue_EHVO_combined_noOut
)
#stop

print('###-------CIV DISTANCE------###')

print('TESTING MEDIAN BLUESHIFT RESAMPLING')
print('non-BAL')
statistical_correlations(
    CivDistance_rankAll_nonbal,
    parameter_name='C IV Distance',
    n_resamples=10000,
    sample_size=99,
    with_replacement=False,
    test_type='median',
    reference_sample=CivDist_EHVOR_combined_noOut
)

print('BAL')
statistical_correlations(
    CivDistance_rankAll_bal,
    parameter_name='C IV Distance',
    n_resamples=10000,
    sample_size=99,
    with_replacement=False,
    test_type='median',
    reference_sample=CivDist_EHVOR_combined_noOut
)

la = np.where(vmax_rankAll_bal > 20000)[0]

print('BAL vmax > 20,000')
statistical_correlations(
    CivDistance_rankAll_bal[la],
    parameter_name='C IV Distance (BAL vmax > 20,000)',
    n_resamples=10000,
    sample_size=99,
    with_replacement=False,
    test_type='median',
    reference_sample=CivDist_EHVOR_combined_noOut
)

print('*** WITH REPLACEMENT***')
print('TESTING MEDIAN BLUESHIFT RESAMPLING')
print('non-BAL')
statistical_correlations(
    CivDistance_rankAll_nonbal,
    parameter_name='C IV Distance',
    n_resamples=10000,
    sample_size=99,
    with_replacement=True,
    test_type='median',
    reference_sample=CivDist_EHVOR_combined_noOut
)

print('BAL')
statistical_correlations(
    CivDistance_rankAll_bal,
    parameter_name='C IV Distance',
    n_resamples=10000,
    sample_size=99,
    with_replacement=True,
    test_type='median',
    reference_sample=CivDist_EHVOR_combined_noOut
)

la = np.where(vmax_rankAll_bal > 20000)[0]

print('BAL vmax > 20,000')
statistical_correlations(
    CivDistance_rankAll_bal[la],
    parameter_name='C IV Distance (BAL vmax > 20,000)',
    n_resamples=10000,
    sample_size=99,
    with_replacement=True,
    test_type='median',
    reference_sample=CivDist_EHVOR_combined_noOut
)

print('----------------- FULL BAL sample correlation Stats ------------------------------')
print("\nBAL (5841) C IV Distance vs Vmin/Vmax")
corr_vmin_p, _ = pearsonr(CivDistance_rankAll_bal, abs(vmin_rankAll_bal))
corr_vmax_p, _ = pearsonr(CivDistance_rankAll_bal, abs(vmax_rankAll_bal))
print(f'Pearson correlation coefficient (Vmin) = {corr_vmin_p}')
print(f'Pearson correlation coefficient (Vmax) = {corr_vmax_p}')

corr_vmin_s, _ = spearmanr(CivDistance_rankAll_bal, abs(vmin_rankAll_bal))
corr_vmax_s, _ = spearmanr(CivDistance_rankAll_bal, abs(vmax_rankAll_bal))
print(f'Spearman correlation coefficient (Vmin) = {corr_vmin_s}')
print(f'Spearman correlation coefficient (Vmax) = {corr_vmax_s}')

corr_vmin_k, _ = kendalltau(CivDistance_rankAll_bal, abs(vmin_rankAll_bal))
corr_vmax_k, _ = kendalltau(CivDistance_rankAll_bal, abs(vmax_rankAll_bal))
print(f'Kendall Tau (Vmin) = {corr_vmin_k}')
print(f'Kendall Tau (Vmax) = {corr_vmax_k}')
print()
print()

print('----------------- EHVO Stats ------------------------------')
print("\nEHVO C IV Distance vs Vmin/Vmax")
corr_vmin_p, _ = pearsonr(CivDist_EHVOR_combined_noOut, abs(vmin_EHVOR_combined_noOut))
corr_vmax_p, _ = pearsonr(CivDist_EHVOR_combined_noOut, abs(vmax_EHVOR_combined_noOut))
print(f'Pearson correlation coefficient (Vmin) = {corr_vmin_p}')
print(f'Pearson correlation coefficient (Vmax) = {corr_vmax_p}')

corr_vmin_s, _ = spearmanr(CivDist_EHVOR_combined_noOut, abs(vmin_EHVOR_combined_noOut))
corr_vmax_s, _ = spearmanr(CivDist_EHVOR_combined_noOut, abs(vmax_EHVOR_combined_noOut))
print(f'Spearman correlation coefficient (Vmin) = {corr_vmin_s}')
print(f'Spearman correlation coefficient (Vmax) = {corr_vmax_s}')

corr_vmin_k, _ = kendalltau(CivDist_EHVOR_combined_noOut, abs(vmin_EHVOR_combined_noOut))
corr_vmax_k, _ = kendalltau(CivDist_EHVOR_combined_noOut, abs(vmax_EHVOR_combined_noOut))
print(f'Kendall Tau (Vmin) = {corr_vmin_k}')
print(f'Kendall Tau (Vmax) = {corr_vmax_k}')
print()
print()

print("\nEHVO C IV Blueshift vs Vmin/Vmax")
corr_vmin_p, _ = pearsonr(CivBlue_EHVO_combined_noOut, abs(vmin_EHVOR_combined_noOut))
corr_vmax_p, _ = pearsonr(CivBlue_EHVO_combined_noOut, abs(vmax_EHVOR_combined_noOut))
print(f'Pearson correlation coefficient (Vmin) = {corr_vmin_p}')
print(f'Pearson correlation coefficient (Vmax) = {corr_vmax_p}')

corr_vmin_s, _ = spearmanr(CivBlue_EHVO_combined_noOut, abs(vmin_EHVOR_combined_noOut))
corr_vmax_s, _ = spearmanr(CivBlue_EHVO_combined_noOut, abs(vmax_EHVOR_combined_noOut))
print(f'Spearman correlation coefficient (Vmin) = {corr_vmin_s}')
print(f'Spearman correlation coefficient (Vmax) = {corr_vmax_s}')


corr_vmin_k, _ = kendalltau(CivBlue_EHVO_combined_noOut, abs(vmin_EHVOR_combined_noOut))
corr_vmax_k, _ = kendalltau(CivBlue_EHVO_combined_noOut, abs(vmax_EHVOR_combined_noOut))
print(f'Kendall Tau (Vmin) = {corr_vmin_k}')
print(f'Kendall Tau (Vmax) = {corr_vmax_k}')
print()
print()




########### OLD VVVVVVVVVVVVVVVVVVVVVVV ###########################################################################################

#################### These were initial tests done trying to recreate Abby's tests ################################################
#----------- Bootstrap resampling stats -------------------------------
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bootstrap.html

#Example from ^^^^^
from scipy.stats import bootstrap
# ============================================================
# METHOD 1: Standard Bootstrap Resampling
# ============================================================
# - Uses the FULL BAL sample (N = 5841) for every resample.
# - Each resample contains 5841 quasars, sampled WITH replacement.
# - Measures the uncertainty/stability of the Pearson r estimated
#   from the full BAL population.
# - paired=True preserves the CIV–velocity pairing for each quasar.
# - n_resamples=100 gives 100 bootstrap correlation coefficients.
#
# PRO:
#   - Standard method for estimating uncertainty in r.
#   - Uses all available BAL data.
#   - Provides a relatively precise estimate of the full sample r.
#
# CON:
#   - Not directly matched to the ~100 object EHVO sample.
#   - BAL samples contain 5841 objects, not 100.
#
# Best for: estimating uncertainty in the correlation of the
#           COMPLETE BAL sample.
print('----------------- BAL Bootstrap Test ---------------------------')
CIV  = np.asarray(CivBlue_rankAll_bal)
VMIN = np.asarray(vmin_rankAll_bal)
VMAX = np.asarray(vmax_rankAll_bal)

# Pearson correlation functions
def corr_civ_vel(civblue, v, axis=0):
    return pearsonr(civblue, v, axis=axis).statistic



#CIV vs VMIN
res_vmin = bootstrap(
    data=(CIV, VMIN),
    statistic=corr_civ_vel,
    paired=True,
    n_resamples=100,
    confidence_level=0.95,
    method='percentile',
    random_state=42
)


#CIV vs VMAX
res_vmax = bootstrap(
    data=(CIV, VMAX),
    statistic=corr_civ_vel,
    paired=True,
    n_resamples=100,
    confidence_level=0.95,
    method='percentile',
    random_state=42
)


#Full-sample Pearson r
r_vmin = corr_civ_vel(CIV, VMIN)
r_vmax = corr_civ_vel(CIV, VMAX)

#The 100 bootstrap r values
boot_r_vmin = res_vmin.bootstrap_distribution
boot_r_vmax = res_vmax.bootstrap_distribution

#Mean and standard deviation of bootstrap r
mean_vmin = np.mean(boot_r_vmin)
std_vmin  = np.std(boot_r_vmin, ddof=1)

mean_vmax = np.mean(boot_r_vmax)
std_vmax  = np.std(boot_r_vmax, ddof=1)


print("\nBAL C IV blueshift vs Vmin")
print(f"Original Pearson r:       {r_vmin:.4f}")
print(f"Mean bootstrap r:         {mean_vmin:.4f}")
print(f"Bootstrap standard dev:   {std_vmin:.4f}")
print(f"95% CI:                   "
      f"[{res_vmin.confidence_interval.low:.4f}, "
      f"{res_vmin.confidence_interval.high:.4f}]")


print("\nBAL C IV blueshift vs Vmax")
print(f"Original Pearson r:       {r_vmax:.4f}")
print(f"Mean bootstrap r:         {mean_vmax:.4f}")
print(f"Bootstrap standard dev:   {std_vmax:.4f}")
print(f"95% CI:                   "
      f"[{res_vmax.confidence_interval.low:.4f}, "
      f"{res_vmax.confidence_interval.high:.4f}]")



# ============================================================
# METHOD 2: Repeated Random Subsampling
# ============================================================
# - Randomly selects 100 BAL quasars from the 5841 BALs.
# - Sampling is WITHOUT replacement within each iteration.
# - Pearson r is calculated for each 100 object sample.
# - Repeated 100 times to obtain a distribution of r values.
#
# PRO:
#   - Directly matches the sample size of an ~100 object EHVO sample.
#   - Shows how much the BAL correlation can vary for samples
#     of the same size as the EHVO sample.
#   - Useful for comparing BAL and EHVO correlations at comparable
#     sample sizes.
#
# CON:
#   - Throws away most of the available BAL data in each iteration.
#   - Correlations from n=100 have substantially larger sampling
#     variation than correlations from n=5841.
#   - This is repeated random subsampling, NOT conventional bootstrap.
#
# Best for: comparing BAL and EHVO correlations when the EHVO sample
#           is much smaller (~100 objects).

print('\n-------- Random select 100 BAL 100 times -----------------')

CIV  = np.asarray(CivBlue_rankAll_bal)
VMIN = np.asarray(vmin_rankAll_bal)
VMAX = np.asarray(vmax_rankAll_bal)


sample_size = 100
n_iterations = 100

rng = np.random.default_rng(42)

#Pearson r from each random BAL sample
r_vmin = np.empty(n_iterations)
r_vmax = np.empty(n_iterations)

N = len(CIV)

for i in range(n_iterations):

    #Randomly select 100 DIFFERENT BAL quasars
    idx = rng.choice(N, size=sample_size, replace=False)

    #Correlation for 100 BAL sample
    r_vmin[i] = pearsonr(
        CIV[idx],
        VMIN[idx]
    ).statistic

    r_vmax[i] = pearsonr(
        CIV[idx],
        VMAX[idx]
    ).statistic


#Summary statistics
mean_vmin = np.mean(r_vmin)
std_vmin = np.std(r_vmin, ddof=1)
ci_vmin = np.percentile(r_vmin, [2.5, 97.5])

mean_vmax = np.mean(r_vmax)
std_vmax = np.std(r_vmax, ddof=1)
ci_vmax = np.percentile(r_vmax, [2.5, 97.5])

print("BAL C IV blueshift vs Vmin")
print(f"Mean Pearson r:         {mean_vmin:.4f}")
print(f"Standard deviation:     {std_vmin:.4f}")
print(f"95% interval:           [{ci_vmin[0]:.4f}, {ci_vmin[1]:.4f}]")

print("\nBAL C IV blueshift vs Vmax")
print(f"Mean Pearson r:         {mean_vmax:.4f}")
print(f"Standard deviation:     {std_vmax:.4f}")
print(f"95% interval:           [{ci_vmax[0]:.4f}, {ci_vmax[1]:.4f}]")







