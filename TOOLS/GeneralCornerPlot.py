# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 00:23:24 2024

@author: annam
"""

import numpy as np
from matplotlib import pyplot as plt
from astropy import stats
import pandas as pd



#XQR-30 Properties File
xqrdata = pd.read_csv('https://raw.githubusercontent.com/paolaUWB/DR16Q/Anna/CROSS_CORRELATION/XQR30_Property_Data.csv')

#Rankines info file: PRH selected the good cases
infoRankineparent = 'https://raw.githubusercontent.com/paolaUWB/DR16Q/master/CROSS_CORRELATION/DR16parent_DR14RankineInfo.csv'
infoRankineEHVO = "https://raw.githubusercontent.com/paolaUWB/DR16Q/master/CROSS_CORRELATION/DR16EHVO_DR14RankineInfo.csv"
infoRankineparent_DR9 = "https://raw.githubusercontent.com/paolaUWB/DR16Q/Liliana/CROSS_CORRELATION/DR9parent_DR14RankineInfo.csv"
infoRankineEHVO_DR9 = "https://raw.githubusercontent.com/paolaUWB/DR16Q/Liliana/CROSS_CORRELATION/DR9EHVO_DR14RankineInfo_withEHVOspeeds.csv"
RH_parent_16= 'https://raw.githubusercontent.com/paolaUWB/DR16Q/masterCROSS_CORRELATION/DR16_parent_sample.csv'

#Extracting values from Rankines info file for EHVO Parent
dfRPA = pd.read_csv(infoRankineparent, header=None)
dfRHV = pd.read_csv(infoRankineEHVO, header=None)
dfRPA_9 = pd.read_csv(infoRankineparent_DR9, header=None)
dfRHV_9 = pd.read_csv(infoRankineEHVO_DR9, header=None)

#Extracting Values from XQR-30 Properties Files
xqr_loglbol = np.array(xqrdata["Log_Lbol"])
xqr_logmbh = np.array(xqrdata["Log_bhm_civ"])
xqr_edd = np.array(xqrdata["edd_civ"])
xqr_logedd = np.log10(xqr_edd)
xqr_EHVOstatus = np.array(xqrdata["status"])
xqr_bs = np.array(xqrdata["CIV_Emission_Blueshift"])

#XQR30: Separating Quasars by type using the midZstatus column of the csv
xqr_EHVO = (xqr_EHVOstatus == "EHVO")
xqr_nonBAL = (xqr_EHVOstatus == "nBAL")
xqr_BAL = (xqr_EHVOstatus == 'BAL')
xqr_unknown = (xqr_EHVOstatus == "N/A")

def safe_log10(arr):
  arr = np.array(arr)
  arr[arr <= 0] = np.nan
  return np.log10(arr)


#EHVO Physical Properties-----------------------------------------------------------
xqr_lbol_EHVO = xqr_loglbol[xqr_EHVO]
xqr_mbh_EHVO = xqr_logmbh[xqr_EHVO]
xqr_edd_EHVO = xqr_logedd[xqr_EHVO]
xqr_bs_EHVO = xqr_bs[xqr_EHVO]



#XQR-30 Non-EHVO Physical Properties-----------------------------------------------------------

xqr_nonBALlbol = xqr_loglbol[xqr_nonBAL]
xqr_nonBALmbh = xqr_logmbh[xqr_nonBAL]
xqr_nonBALedd = xqr_logedd[xqr_nonBAL]
xqr_bs_nonBAL = xqr_bs[xqr_nonBAL]
xqr_BALlbol = xqr_loglbol[xqr_BAL]
xqr_BALmbh = xqr_logmbh[xqr_BAL]
xqr_BALedd = xqr_logedd[xqr_BAL]
xqr_bs_BAL = xqr_bs[xqr_BAL]


xqr_lbol_Parent = np.concatenate([xqr_nonBALlbol, xqr_BALlbol])
xqr_mbh_Parent = np.concatenate([xqr_nonBALmbh, xqr_BALmbh])
xqr_edd_Parent = np.concatenate([xqr_nonBALedd, xqr_BALedd])
xqr_bs_Parent = np.concatenate([xqr_bs_nonBAL, xqr_bs_BAL])

xqr_bs_both = np.concatenate([xqr_bs_Parent, xqr_bs_EHVO])

###DR16###
#for parent sample
Parentin16Rank_mbh=dfRPA[dfRPA.columns[16]].to_numpy()
Parentin16Rank_lbol=dfRPA[dfRPA.columns[17]].to_numpy()
Parentin16Rank_redd=dfRPA[dfRPA.columns[18]].to_numpy()
Parentin16Rank_CivBlue=dfRPA[dfRPA.columns[5]].to_numpy()
Parentin16Rank_CivEW=dfRPA[dfRPA.columns[6]].to_numpy()
Parentin16Rank_HeIIEW=dfRPA[dfRPA.columns[7]].to_numpy()


 #for EHVOs
EHVOin16Rank_mbh=dfRHV[dfRHV.columns[16]].to_numpy()
EHVOin16Rank_lbol=dfRHV[dfRHV.columns[17]].to_numpy()
EHVOin16Rank_redd=dfRHV[dfRHV.columns[18]].to_numpy()
EHVOin16Rank_CivBlue=dfRHV[dfRHV.columns[5]].to_numpy()
EHVOin16Rank_CivEW=dfRHV[dfRHV.columns[6]].to_numpy()
EHVOin16Rank_HeIIEW=dfRPA[dfRPA.columns[7]].to_numpy()

###DR9###

#for parent sample
Parentin9Rank_mbh=dfRPA_9[dfRPA_9.columns[16]].to_numpy()
Parentin9Rank_lbol=dfRPA_9[dfRPA_9.columns[17]].to_numpy()
Parentin9Rank_redd=dfRPA_9[dfRPA_9.columns[18]].to_numpy()
Parentin9Rank_CivBlue=dfRPA_9[dfRPA_9.columns[5]].to_numpy()
Parentin9Rank_CivEW=dfRPA_9[dfRPA_9.columns[6]].to_numpy()
Parentin9Rank_HeIIEW=dfRPA[dfRPA.columns[7]].to_numpy()

#for EHVOs
EHVOin9Rank_mbh=dfRHV_9[dfRHV_9.columns[16]].to_numpy()
EHVOin9Rank_lbol=dfRHV_9[dfRHV_9.columns[17]].to_numpy()
EHVOin9Rank_redd=dfRHV_9[dfRHV_9.columns[18]].to_numpy()
EHVOin9Rank_CivBlue=dfRHV_9[dfRHV_9.columns[5]].to_numpy()
EHVOin9Rank_CivEW=dfRHV_9[dfRHV_9.columns[6]].to_numpy()
EHVOin9Rank_HeIIEW=dfRPA[dfRPA.columns[7]].to_numpy()

###Combined DR9 and DR16 Arrays###

#Parent Sample
parentboth_mbh = np.concatenate([Parentin9Rank_mbh, Parentin16Rank_mbh])
parentboth_lbol = np.concatenate([Parentin9Rank_lbol, Parentin16Rank_lbol])
parentboth_redd = np.concatenate([Parentin9Rank_redd, Parentin16Rank_redd])
parentboth_CivBlue = np.concatenate([Parentin9Rank_CivBlue, Parentin16Rank_CivBlue])
parentboth_CivEW = np.concatenate([Parentin9Rank_CivEW, Parentin16Rank_CivEW])
parentboth_HeIIEW = np.concatenate([Parentin9Rank_HeIIEW, Parentin16Rank_HeIIEW])
parentboth_logHeIIEW = np.log10(parentboth_HeIIEW)

#EHVOs
EHVOboth_mbh = np.concatenate([EHVOin9Rank_mbh, EHVOin16Rank_mbh])
EHVOboth_lbol = np.concatenate([EHVOin9Rank_lbol, EHVOin16Rank_lbol])
EHVOboth_redd = np.concatenate([EHVOin9Rank_redd, EHVOin16Rank_redd])
EHVOboth_CivBlue = np.concatenate([EHVOin9Rank_CivBlue, EHVOin16Rank_CivBlue])
EHVOboth_CivEW = np.concatenate([EHVOin9Rank_CivEW, EHVOin16Rank_CivEW])
EHVOboth_HeIIEW = np.concatenate([EHVOin9Rank_HeIIEW, EHVOin16Rank_HeIIEW])
EHVOboth_logHeIIEW = np.log10(EHVOboth_HeIIEW)

#'Good' Values, Coatman et al. 2017 cutoff for Black hole mass (and thus eddington ratio) plots
bb9 = np.where(Parentin9Rank_CivBlue > 500)
bb9_EHVO = np.where(EHVOin9Rank_CivBlue > 500)
bb16 = np.where(Parentin16Rank_CivBlue > 500)
bb16_EHVO = np.where(EHVOin16Rank_CivBlue > 500)
bb_both = np.where(parentboth_CivBlue > 500)
bbboth_EHVO = np.where(EHVOboth_CivBlue > 500)

bbXQR_nonBAL = np.where(xqr_bs_nonBAL > 500)
bbXQR_BAL = np.where(xqr_bs_BAL > 500)
bbXQR_EHVO = np.where(xqr_bs_EHVO > 500)



#highbsEHVO = np.where(EHVOboth_CivBlue > 2500)



#Combined EHVO and Parent Blueshift Arrays for Color Map
c_DR16_CivBlue = np.concatenate([Parentin16Rank_CivBlue[bb16], EHVOin16Rank_CivBlue[bb16_EHVO]])
c_DR9_CivBlue = np.concatenate([Parentin9Rank_CivBlue[bb9], EHVOin9Rank_CivBlue])
c_both_CivBlue = np.concatenate([parentboth_CivBlue[bb_both], EHVOboth_CivBlue[bbboth_EHVO]])

#COMBINED EHVO and Parent CIV Equivalent Width Arrays for Color Map
c_DR16_CivEW = np.concatenate([Parentin16Rank_CivEW[bb16], EHVOin16Rank_CivEW[bb16_EHVO]])
c_DR9_CivEW = np.concatenate([Parentin9Rank_CivEW[bb9], EHVOin9Rank_CivEW])
c_both_CivEW = np.concatenate([parentboth_CivEW[bb_both], EHVOboth_CivEW[bbboth_EHVO]])

#COMBINED EHVO and Parent He II Equivalent Width Arrays for Color Map
c_DR16_HeIIEW = np.concatenate([Parentin16Rank_HeIIEW[bb16], EHVOin16Rank_HeIIEW[bb16_EHVO]])
c_DR9_HeIIEW = np.concatenate([Parentin9Rank_HeIIEW[bb9], EHVOin9Rank_HeIIEW])
c_both_HeIIEW = np.concatenate([parentboth_HeIIEW[bb_both], EHVOboth_HeIIEW[bbboth_EHVO]])
c_both_logHeIIEW = safe_log10(c_both_HeIIEW)

#Color Map
cm = plt.get_cmap('viridis')
cmap = plt.get_cmap('viridis')

#Histogram------------------------------------------------
#MBH EHVO
c_both_mbh = np.concatenate([parentboth_mbh[bb_both], EHVOboth_mbh[bbboth_EHVO]])
knuthbinMBH = (stats.knuth_bin_width(c_both_mbh))
bin_MBH = np.arange(min(c_both_mbh), max(c_both_mbh), knuthbinMBH)

#Lbol EHVO
c_both_lbol = np.concatenate([parentboth_lbol[bb_both], EHVOboth_lbol[bbboth_EHVO]])
knuthbinLbol = (stats.knuth_bin_width(c_both_lbol))
bin_Lbol = np.arange(min(c_both_lbol), max(c_both_lbol), knuthbinLbol)

#Edd EHVO
c_both_edd = np.concatenate([parentboth_redd[bb_both], EHVOboth_redd[bbboth_EHVO]])
knuthbinEdd = (stats.knuth_bin_width(c_both_edd))
bin_Edd = np.arange(min(c_both_edd), max(c_both_edd), knuthbinEdd)

#Histogram Weights

weightsEHVOmbh = [np.ones_like(EHVOboth_mbh[bbboth_EHVO])/float(len(EHVOboth_mbh[bbboth_EHVO]))]
weightsEHVOlbol = [np.ones_like(EHVOboth_lbol[bbboth_EHVO])/float(len(EHVOboth_lbol[bbboth_EHVO]))]
weightsEHVOedd = [np.ones_like(EHVOboth_redd[bbboth_EHVO])/float(len(EHVOboth_redd[bbboth_EHVO]))]
weightsParentmbh = [np.ones_like(parentboth_mbh[bb_both])/float(len(parentboth_mbh[bb_both]))]
weightsParentlbol = [np.ones_like(parentboth_lbol[bb_both])/float(len(parentboth_lbol[bb_both]))]
weightsParentedd = [np.ones_like(parentboth_redd[bb_both])/float(len(parentboth_redd[bb_both]))]

#Setting up cornerplot positioning, works like coordinates starting from the top left

fig = plt.figure(figsize=(16.5, 10.5), dpi=500)

gs = plt.GridSpec(3, 4)

#Middle left: LBOL vs MBH
ax1 = fig.add_subplot(gs[1, 0]) 

#Bottom left: EDD vs MBH
ax2 = fig.add_subplot(gs[2, 0])

#Bottom middle: EDD vs LBOL
ax3 = fig.add_subplot(gs[2, 1])

#Top left: MBH Histogram
ax4 = fig.add_subplot(gs[0, 0])

#Middle middle: LBOL Histogram
ax5 = fig.add_subplot(gs[1, 1])

#Bottom right: EDD Histogram
ax6 = fig.add_subplot(gs[2, 2])

#Top middle: Legend for Histograms
ax7 = fig.add_subplot(gs[0,1])

#--------------LBOL vs MBH--------------------

#Put hexagon data here
cm = ax1.hexbin(parentboth_mbh[bb_both], parentboth_lbol[bb_both],  C = parentboth_CivBlue[bb_both], gridsize = 35, cmap = 'viridis', reduce_C_function = np.median, mincnt = 4, vmin=500, vmax=max(EHVOboth_CivBlue), edgecolor = 'None', zorder= 99)

cbarax = fig.add_axes([-0.03, 0.055, 0.02, 0.93])
cbar = plt.colorbar(cm, cax=cbarax)
cbarax.set_ylabel('CIV Blueshift (km/s)', rotation = 90, fontname='serif')

#Tick positioning
cbarax.yaxis.set_ticks_position('left')
cbarax.yaxis.set_label_position('left')

#Put scatterplot data here
sc = ax1.scatter(EHVOboth_mbh[bbboth_EHVO], EHVOboth_lbol[bbboth_EHVO], c = EHVOboth_CivBlue[bbboth_EHVO], zorder = 100, cmap = 'viridis', edgecolor = 'k', s = 30, vmin=cbar.vmin, vmax=cbar.vmax, label = 'EHVO (RH+ 2020 and RH+ in prep)')


ax1.set_xlim(8.4,10.5)
ax1.set_ylim(46.1,48.1)
ax1.set_xticklabels([])
ax1.set_ylabel('log($L_{bol}$/erg s$^{-1}$)', fontname='serif')
ax1.tick_params(top=True, labeltop=False, bottom=True, right = True, labelright = False)

#-------------------EDD vs MBH------------------

#Put hexagon data here
cm = ax2.hexbin(parentboth_mbh[bb_both], parentboth_redd[bb_both],  C = parentboth_CivBlue[bb_both], gridsize = 35, cmap = 'viridis', reduce_C_function = np.median, mincnt = 4, vmin=500, vmax=max(EHVOboth_CivBlue), edgecolor = 'None', zorder= 99)
#Put scatterplot data here
sc = ax2.scatter(EHVOboth_mbh[bbboth_EHVO], EHVOboth_redd[bbboth_EHVO], c = EHVOboth_CivBlue[bbboth_EHVO], zorder = 100, cmap = 'viridis', edgecolor = 'k', s = 30, vmin=cbar.vmin, vmax=cbar.vmax, label = 'EHVO (RH+ 2020 and RH+ in prep)')

ax2.set_xlabel("log($M_{\mathrm{BH}}/M_{\odot})$", fontname='serif')
ax2.set_xlim(8.4,10.5)
ax2.set_ylim(-1.5,0.5)
ax2.legend(loc='upper right',prop={'size': 9.7})
ax2.set_ylabel('Eddington Ratio (log($L_{bol}$/$L_{Edd}$))',fontname='serif')
ax2.tick_params(top=True, labeltop=False, bottom=True, right = True, labelright = False)

#---------------EDD vs LBOL---------------------

#Put hexagon data here
cm = ax3.hexbin(parentboth_lbol[bb_both], parentboth_redd[bb_both],  C = parentboth_CivBlue[bb_both], gridsize = 35, cmap = 'viridis', reduce_C_function = np.median, mincnt = 4, vmin=500, vmax=max(EHVOboth_CivBlue), edgecolor = 'None', zorder= 99)
#Put scatterplot data here
sc = ax3.scatter(EHVOboth_lbol[bbboth_EHVO], EHVOboth_redd[bbboth_EHVO], c = EHVOboth_CivBlue[bbboth_EHVO], zorder = 100, cmap = 'viridis', edgecolor = 'k', s = 30, vmin=cbar.vmin, vmax=cbar.vmax, label = 'EHVO (RH+ 2020 and RH+ in prep)')
ax3.set_xlabel('log($L_{bol}$/erg s$^{-1}$)', fontname='serif')
#Removing yticks for second plot
ax3.set_yticklabels([])
ax3.set_xlim(46,48.5)
ax3.set_ylim(-1.5,0.5)
#ax3.legend(loc='upper right')
ax3.tick_params(top=True, labeltop=False, bottom=True, right = True, labelright = False)

#---------------MBH Histogram---------------------

#Removing x and y labels
ax4.set_xticklabels([])
ax4.set_yticklabels([])
ax4.set_xlim(8.4,10.5)

#Parent histogram data
ax4.hist(parentboth_mbh[bb_both], bins = bin_MBH, weights= weightsParentmbh, histtype = 'step', color = 'blue', label = 'Parent RH+ 2020 and RH+ in prep)')

#EHVO histogram data
ax4.hist(EHVOboth_mbh[bbboth_EHVO], bins = bin_MBH, weights = weightsEHVOmbh, histtype = 'step', color = 'red', label='EHVO (RH+ 2020 and RH+ in prep)')

ax4.tick_params(top=True, labeltop=False, bottom=True, right = True, labelright = False)


#---------------LBOL Histogram--------------------

#Removing x and y labels
ax5.set_xticklabels([])
ax5.set_yticklabels([])

#Parent histogram data
ax5.hist(parentboth_lbol[bb_both], bins = bin_Lbol, weights = weightsParentlbol, histtype = 'step', color = 'blue', label = 'Parent RH+ 2020 and RH+ in prep)')

#EHVO histogram data
ax5.hist(EHVOboth_lbol[bbboth_EHVO], bins = bin_Lbol, weights = weightsEHVOlbol, histtype = 'step', color = 'red', label='EHVO (RH+ 2020 and RH+ in prep)')

ax5.set_xlim(46,48.5)
ax5.tick_params(top=True, labeltop=False, bottom=True, right = True, labelright = False)

#---------------EDD Histogram----------------------
ax6.set_yticklabels([])
ax6.set_xlabel('Eddington Ratio (log($L_{bol}$/$L_{Edd}$))',fontname='serif')
ax6.hist(parentboth_redd[bb_both], bins = bin_Edd, weights = weightsParentedd, histtype = 'step', color = 'blue', label = 'Parent RH+ 2020 and RH+ in prep)')
ax6.hist(EHVOboth_redd[bbboth_EHVO], bins = bin_Edd, weights = weightsEHVOedd, histtype = 'step', color = 'red', label='EHVO (RH+ 2020 and RH+ in prep)')
ax6.set_xlim(-1.5,0.5)
ax6.tick_params(top=True, labeltop=False, bottom=True, right = True, labelright = False)

#HISTOGRAM LEGEND-----------------------------------------------------
#Legend: Everything in here is in an effort to make everything but the legend invisible, a way to force where the legend goes

#Sets up the colored lines on the legend, but makes the lines not appear on the graph
parent_line = plt.Line2D([0], [0], color='blue', label='Parent RH+ 2020 and RH+ in prep')
EHVO_line = plt.Line2D([0], [0], color='red', label='EHVO (RH+ 2020 and RH+ in prep)')

#no ticks
ax7.set_xticks([])
ax7.set_yticks([])

#Legend location
ax7.legend(handles=[parent_line, EHVO_line], loc='lower left', prop={'size': 9.7})

#Making graph lines invisible
ax7.spines['top'].set_visible(False)
ax7.spines['right'].set_visible(False)
ax7.spines['bottom'].set_visible(False)
ax7.spines['left'].set_visible(False)




#plt.subplots_adjust(wspace=0.16, hspace=0.16)
plt.tight_layout()
#filepath = "C:/Users/annam/Desktop/Research/Emission_PhysicalProperties/Plots/CORNERPLOTSV3.pdf"
#plt.savefig(filepath, dpi=400, bbox_inches='tight')