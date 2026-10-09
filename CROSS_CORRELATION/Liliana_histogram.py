#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 13:16:38 2026

@author: awjaustin
"""

#Display histograms?
DisplayHistogram1 = "y"
#DisplayHistogram1 = "n"

DisplayHistogram2 = "y"
#DisplayHistogram2 = "n"

DisplayHistogram3 = "y"
#DisplayHistogram3 = "n"

#sample settings
sample = "DR9 & DR16"
#sample = "DR9"
#sample = "DR16"

#velocity values
velocity = max
#velocity = min

#values for speed separation
BAL_value = 15_000
EHVO_value = 42_000
# note: BAl min median - 6,000    BAL max median - 15,000
#       EHVO min median - 38,000  EHVO max median - 42,000

#display median lines?
med_lines = 'yes'
#med_lines = 'no'

# colors
bal_lowv_color = 'aqua'
bal_highv_color = 'cornflowerblue'

EHVO_lowv_color = 'plum'
EHVO_highv_color = 'purple'

# Location of saved images (currently not in use as of 10/18 - LEF)
#images_dir = '/content/drive/MyDrive/Images/'


import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from astropy import stats
from astropy.io import fits
import os
from collections import defaultdict
from matplotlib.bezier import math
#from numpy.lib.function_base import average
from numpy.core.fromnumeric import sort
from numpy.ma.extras import median



#If you are trying to add a new csv, or the location of one of these files changed, remember that all you have to do is:
#"copy link address" of the file in GitHub, then paste here in this format: "https://github.com/paolaUWB/DR16Q/blob/branch/folder_name/file_name.csv"
#then just replace "github" with "raw.githubusercontent" and you're good to go.
#highZfile = pd.read_csv('https://raw.githubusercontent.com/paolaUWB/DR16Q/Anna/CROSS_CORRELATION/highRedshiftQuasarCSV.csv')
infoRankineparent = 'https://raw.githubusercontent.com/paolaUWB/DR16Q/master/CROSS_CORRELATION/DR16parent_DR14RankineInfo.csv'
infoRankineEHVO = 'https://raw.githubusercontent.com/paolaUWB/DR16Q/Liliana/CROSS_CORRELATION/DR16EHVO_DR14RankineInfo_wEHVOspeed_vmax_reordered.csv'


#adding in DR9 sample
infoRankineparent_DR9 = "https://raw.githubusercontent.com/paolaUWB/DR16Q/Liliana/CROSS_CORRELATION/DR9parent_DR14RankineInfo.csv"
infoRankineEHVO_DR9 = "https://raw.githubusercontent.com/paolaUWB/DR16Q/Liliana/CROSS_CORRELATION/DR9EHVO_DR14RankineInfo_wEHVOspeed_vmax_reordered.csv"


#Extracting values from files -------------------------
dfRPA = pd.read_csv(infoRankineparent, header=None)
dfRHV = pd.read_csv(infoRankineEHVO, header=0)
dfRPA_9 = pd.read_csv(infoRankineparent_DR9, header=None)
dfRHV_9 = pd.read_csv(infoRankineEHVO_DR9, header=0)

#DR16............................................
bi_bi0_DR16 = dfRPA[dfRPA.columns[10]].to_numpy()
bi_vmax_DR16 = dfRPA[dfRPA.columns[11]].to_numpy()
bi_vmin_DR16 = dfRPA[dfRPA.columns[12]].to_numpy()
SNR_DR16 = dfRPA[dfRPA.columns[15]].to_numpy()

mbh_values_DR16 = dfRPA[dfRPA.columns[16]].to_numpy()
lbol_values_DR16 = dfRPA[dfRPA.columns[17]].to_numpy()
redd_values_DR16 = dfRPA[dfRPA.columns[18]].to_numpy()

z_DR16 = dfRPA[dfRPA.columns[4]].to_numpy()

     #for EHVOs
EHVO_mbh_values_DR16 = dfRHV[dfRHV.columns[16]].to_numpy()
EHVO_lbol_values_DR16 = dfRHV[dfRHV.columns[17]].to_numpy()
EHVO_redd_values_DR16 = dfRHV[dfRHV.columns[18]].to_numpy()

EHVO_z_DR16 = dfRHV[dfRHV.columns[4]].to_numpy()

EHVO_vmin_DR16 = dfRHV[dfRHV.columns[29]].to_numpy()
EHVO_vmax_DR16 = dfRHV[dfRHV.columns[33]].to_numpy()

#DR9.................................................
bi_bi0_DR9 = dfRPA_9[dfRPA_9.columns[10]].to_numpy()
bi_vmax_DR9 = dfRPA_9[dfRPA_9.columns[11]].to_numpy()
bi_vmin_DR9 = dfRPA_9[dfRPA_9.columns[12]].to_numpy()
SNR_DR9 = dfRPA_9[dfRPA_9.columns[15]].to_numpy()

mbh_values_DR9 = dfRPA_9[dfRPA_9.columns[16]].to_numpy()
lbol_values_DR9 = dfRPA_9[dfRPA_9.columns[17]].to_numpy()
redd_values_DR9 = dfRPA_9[dfRPA_9.columns[18]].to_numpy()

z_DR9 = dfRPA_9[dfRPA_9.columns[4]].to_numpy()

      #for EHVOs
EHVO_mbh_values_DR9 = dfRHV_9[dfRHV_9.columns[16]].to_numpy()
EHVO_lbol_values_DR9 = dfRHV_9[dfRHV_9.columns[17]].to_numpy()
EHVO_redd_values_DR9 = dfRHV_9[dfRHV_9.columns[18]].to_numpy()

EHVO_z_DR9 = dfRHV_9[dfRHV_9.columns[4]].to_numpy()

EHVO_vmax_DR9 = dfRHV_9[dfRHV_9.columns[29]].to_numpy()
EHVO_vmin_DR9 = dfRHV_9[dfRHV_9.columns[32]].to_numpy()

#if statements for sample choice ------------------------------
if sample == "DR16":
  bi_bi0 = bi_bi0_DR16
  bi_vmax = bi_vmax_DR16
  bi_vmin = bi_vmin_DR16
  SNR = SNR_DR16

  mbh_values = mbh_values_DR16
  lbol_values = lbol_values_DR16
  redd_values = redd_values_DR16

  EHVO_mbh_values = EHVO_mbh_values_DR16
  EHVO_lbol_values = EHVO_lbol_values_DR16
  EHVO_redd_values = EHVO_redd_values_DR16

  EHVO_vmin = EHVO_vmin_DR16
  EHVO_vmax = EHVO_vmax_DR16

  plot_title = 'DR16'

if sample == "DR9":
  bi_bi0 = bi_bi0_DR9
  bi_vmax = bi_vmax_DR9
  bi_vmin = bi_vmin_DR9
  SNR = SNR_DR9

  mbh_values = mbh_values_DR9
  lbol_values = lbol_values_DR9
  redd_values = redd_values_DR9

  EHVO_mbh_values = EHVO_mbh_values_DR9
  EHVO_lbol_values = EHVO_lbol_values_DR9
  EHVO_redd_values = EHVO_redd_values_DR9

  EHVO_vmin = EHVO_vmin_DR9
  EHVO_vmax = EHVO_vmax_DR9

  plot_title = 'DR9'

if sample == "DR9 & DR16":
  bi_bi0 = np.concatenate((bi_bi0_DR16, bi_bi0_DR9))
  bi_vmax = np.concatenate((bi_vmax_DR16, bi_vmax_DR9))
  bi_vmin = np.concatenate((bi_vmin_DR16, bi_vmin_DR9))
  SNR = np.concatenate((SNR_DR16, SNR_DR9))

  mbh_values = np.concatenate((mbh_values_DR16, mbh_values_DR9))
  lbol_values = np.concatenate((lbol_values_DR16, lbol_values_DR9))
  redd_values = np.concatenate((redd_values_DR16, redd_values_DR9))

  EHVO_mbh_values = np.concatenate((EHVO_mbh_values_DR16, EHVO_mbh_values_DR9))
  EHVO_lbol_values = np.concatenate((EHVO_lbol_values_DR16, EHVO_lbol_values_DR9))
  EHVO_redd_values = np.concatenate((EHVO_redd_values_DR16, EHVO_redd_values_DR9))

  EHVO_vmin = np.concatenate((EHVO_vmin_DR16, EHVO_vmin_DR9))
  EHVO_vmax = np.concatenate((EHVO_vmax_DR16, EHVO_vmax_DR9))

  plot_title = 'DR9 & DR16'

#Defining BAL vs Non-BAL -------------------------------------------------------------
non_bal = np.where((bi_bi0==0)&(SNR > 10))
pos_bal = np.where((bi_bi0>0) & (SNR > 10))

#Z value variables---------------------------------------------------------------------
if sample == "DR16":
  z_EHVO = EHVO_z_DR16
  z_bal = z_DR16[pos_bal]

if sample == "DR9":
  z_EHVO = EHVO_z_DR9
  z_bal = z_DR9[pos_bal]

if sample == "DR9 & DR16":
  z_EHVO = np.concatenate((EHVO_z_DR16, EHVO_z_DR9))
  z_DR16_DR9 = np.concatenate((z_DR16, z_DR9))
  z_bal = z_DR16_DR9[pos_bal]

#Velocity if statement --------------------------------
if velocity == min:
  bi_v = bi_vmin
  EHVO_v = EHVO_vmin
  text_v = 'vmin'
elif velocity == max:
  bi_v = bi_vmax
  EHVO_v = EHVO_vmax
  text_v = 'vmax'
else:
  print('please input "min" or "max"' )

#speed analysis -------------------------------------
#BAL
BAL_low_velo_index =  np.where((bi_v < BAL_value)&(bi_bi0>0)&(SNR>10))
BAL_high_velo_index = np.where((bi_v > BAL_value)&(bi_bi0>0)&(SNR>10))

mbh_bal_lowv, mbh_bal_highv   = mbh_values[BAL_low_velo_index], mbh_values[BAL_high_velo_index]
redd_bal_lowv, redd_bal_highv = redd_values[BAL_low_velo_index], redd_values[BAL_high_velo_index]
lbol_bal_lowv, lbol_bal_highv = lbol_values[BAL_low_velo_index], lbol_values[BAL_high_velo_index]

#EHVO
EHVO_low_velo_index =  np.where(abs(EHVO_v) < EHVO_value)
EHVO_high_velo_index = np.where(abs(EHVO_v) > EHVO_value)

mbh_EHVO_lowv, mbh_EHVO_highv   = EHVO_mbh_values[EHVO_low_velo_index], EHVO_mbh_values[EHVO_high_velo_index]
redd_EHVO_lowv, redd_EHVO_highv = EHVO_redd_values[EHVO_low_velo_index], EHVO_redd_values[EHVO_high_velo_index]
lbol_EHVO_lowv, lbol_EHVO_highv = EHVO_lbol_values[EHVO_low_velo_index], EHVO_lbol_values[EHVO_high_velo_index]

#bins ------------------------------------------------------------------------
lbol_values = lbol_values[pos_bal]
mbh_values = mbh_values[pos_bal]
redd_values = redd_values[pos_bal]

knuthbin1 = 0.25*(stats.knuth_bin_width(lbol_bal_lowv) + stats.knuth_bin_width(lbol_bal_highv) + stats.knuth_bin_width(lbol_EHVO_lowv) + stats.knuth_bin_width(lbol_EHVO_highv))
bin_lbol = np.arange(min(min(lbol_values), min(EHVO_lbol_values)), max(max(lbol_values), max(EHVO_lbol_values)), knuthbin1)

knuthbin2 = 0.25*(stats.knuth_bin_width(mbh_bal_lowv) + stats.knuth_bin_width(mbh_bal_highv) + stats.knuth_bin_width(mbh_EHVO_lowv) + stats.knuth_bin_width(mbh_EHVO_highv))
bin_mbh = np.arange(min(min(mbh_values), min(EHVO_mbh_values)), max(max(mbh_values), max(EHVO_mbh_values)), knuthbin2)

knuthbin3 = 0.25*(stats.knuth_bin_width(redd_bal_lowv) + stats.knuth_bin_width(redd_bal_highv) + stats.knuth_bin_width(redd_EHVO_lowv) + stats.knuth_bin_width(redd_EHVO_highv))
bin_redd = np.arange(min(min(redd_values), min(EHVO_redd_values)), max(max(redd_values), max(EHVO_redd_values)), knuthbin3)

#Histograms--------------------------------------------------------------------------------
#histogram 1 ..................
if DisplayHistogram1 == "y":
  weightslbol_bal_lowv = [np.ones_like(lbol_bal_lowv)/float(len(lbol_bal_lowv))]
  plt.hist(lbol_bal_lowv, bin_lbol, weights = weightslbol_bal_lowv, label = 'BAL '+ str(text_v) + '<'+ str(BAL_value) + ' km/s', color = bal_lowv_color, alpha = 1,)

  weightslbol_bal_highv = [np.ones_like(lbol_bal_highv)/float(len(lbol_bal_highv))]
  plt.hist(lbol_bal_highv, bin_lbol, weights = weightslbol_bal_highv, label = 'BAL '+ str(BAL_value) + ' km/s<'+ str(text_v) +'<25,000 km/s', color = bal_highv_color,alpha = 0.6,)

  weightslbol_EHVO_lowv = [np.ones_like(lbol_EHVO_lowv)/float(len(lbol_EHVO_lowv))]
  plt.hist(lbol_EHVO_lowv, bin_lbol, weights = weightslbol_EHVO_lowv, label = 'EHVO 30,000 km/s<'+ str(text_v) +'<'+ str(EHVO_value) + ' km/s', color = EHVO_lowv_color, alpha=0.8,)

  weightslbol_EHVO_highv = [np.ones_like(lbol_EHVO_highv)/float(len(lbol_EHVO_highv))]
  plt.hist(lbol_EHVO_highv, bin_lbol, weights = weightslbol_EHVO_highv, label = 'EHVO '+ str(text_v) +'>'+ str(EHVO_value) +' km/s', color = EHVO_highv_color, alpha = 0.6,)

  if med_lines == 'yes':
    plt.axvline(x = median(lbol_bal_lowv), color='aqua', linestyle='dashed', linewidth=2)
    plt.axvline(x = median(lbol_bal_highv), color='cornflowerblue', linestyle='dashed', linewidth=2)
    plt.axvline(x = median(lbol_EHVO_lowv), color='plum', linestyle='dashed', linewidth=2)
    plt.axvline(x = median(lbol_EHVO_highv), color='purple', linestyle='dashed', linewidth=2)

  plt.xlabel('log($L_{bol}$/erg s$^{-1}$)')
  plt.ylabel('N')
  plt.legend(title = plot_title, bbox_to_anchor=(0.8, 1), loc='upper left', prop={'size': 7})
  plt.legend(title = plot_title, loc='best',prop={'size': 7})
  plt.show()

#histogram 2 ......................
if DisplayHistogram2 == "y":
  weights_mbh_bal_lowv = [np.ones_like(mbh_bal_lowv)/float(len(mbh_bal_lowv))]
  plt.hist(mbh_bal_lowv, bins = bin_mbh, weights = weights_mbh_bal_lowv, label = 'BAL '+ str(text_v) + '<'+ str(BAL_value) + ' km/s', color = bal_lowv_color, alpha=1 )

  weights_mbh_bal_highv = [np.ones_like(mbh_bal_highv)/float(len(mbh_bal_highv))]
  plt.hist(mbh_bal_highv, bins = bin_mbh, weights = weights_mbh_bal_highv, label = 'BAL '+ str(BAL_value) + ' km/s<'+ str(text_v) +'<25,000 km/s', color = bal_highv_color, alpha=0.7)

  weights_mbh_EHVO_lowv = [np.ones_like(mbh_EHVO_lowv)/float(len(mbh_EHVO_lowv))]
  plt.hist(mbh_EHVO_lowv, bins = bin_mbh, weights = weights_mbh_EHVO_lowv, label = 'EHVO 30,000 km/s<'+ str(text_v) +'<'+ str(EHVO_value) + ' km/s', color = EHVO_lowv_color, alpha=0.9,)

  weights_mbh_EHVO_highv = [np.ones_like(mbh_EHVO_highv)/float(len(mbh_EHVO_highv))]
  plt.hist(mbh_EHVO_highv, bins = bin_mbh, weights = weights_mbh_EHVO_highv, label = 'EHVO '+ str(text_v) +'>'+ str(EHVO_value) +' km/s', color = EHVO_highv_color, alpha = 0.6,)

  if med_lines == 'yes':
    plt.axvline(x = median(mbh_bal_lowv), color='aqua', linestyle='dashed', linewidth=2)
    plt.axvline(x = median(mbh_bal_highv), color='cornflowerblue', linestyle='dashed', linewidth=2)
    plt.axvline(x = median(mbh_EHVO_lowv), color='plum', linestyle='dashed', linewidth=2)
    plt.axvline(x = median(mbh_EHVO_highv), color='purple', linestyle='dashed', linewidth=2)

  plt.xlabel('log($M_{\mathrm{BH}}/M_{\odot})$')
  plt.ylabel('N')
  plt.legend(title = plot_title, loc='best',prop={'size': 7})
  plt.show()

#histogram 3 .........................
if DisplayHistogram3 == "y":
  weights_redd_bal_lowv = [np.ones_like(redd_bal_lowv)/float(len(redd_bal_lowv))]
  plt.hist(redd_bal_lowv, bins = bin_redd, weights = weights_redd_bal_lowv, label = 'BAL '+ str(text_v) + '<'+ str(BAL_value) + ' km/s', color = bal_lowv_color,  alpha = 1 )

  weights_redd_bal_highv = [np.ones_like(redd_bal_highv)/float(len(redd_bal_highv))]
  plt.hist(redd_bal_highv, bins = bin_redd, weights = weights_redd_bal_highv, label = 'BAL '+ str(BAL_value) + ' km/s<'+ str(text_v) +'<25,000 km/s', color = bal_highv_color, alpha=0.7)

  weights_redd_EHVO_lowv = [np.ones_like(redd_EHVO_lowv)/float(len(redd_EHVO_lowv))]
  plt.hist(redd_EHVO_lowv, bins = bin_redd, weights = weights_mbh_EHVO_lowv, label = 'EHVO 30,000 km/s<'+ str(text_v) +'<'+ str(EHVO_value) + ' km/s', color = EHVO_lowv_color, alpha=1)

  weights_redd_EHVO_highv = [np.ones_like(redd_EHVO_highv)/float(len(redd_EHVO_highv))]
  plt.hist(redd_EHVO_highv,bins = bin_redd, weights = weights_redd_EHVO_highv, label = 'EHVO '+ str(text_v) +'>'+ str(EHVO_value) +' km/s', color = EHVO_highv_color, alpha = 0.7)

  if med_lines == 'yes':
    plt.axvline(x = median(redd_bal_lowv), color='aqua', linestyle='dashed', linewidth=2)
    plt.axvline(x = median(redd_bal_highv), color='cornflowerblue', linestyle='dashed', linewidth=2)
    plt.axvline(x = median(redd_EHVO_lowv), color='plum', linestyle='dashed', linewidth=2)
    plt.axvline(x = median(redd_EHVO_highv), color='purple', linestyle='dashed', linewidth=2)

  plt.xlabel('Eddington Ratio (log($L_{bol}$/$L_{Edd}$)')
  plt.ylabel('N')
  plt.legend(title = plot_title, loc='best',prop={'size': 7})
  plt.show()

print(np.sort(mbh_EHVO_highv))
print(np.sort(mbh_EHVO_lowv))

print(median(mbh_bal_lowv))
print(median(mbh_bal_highv))
print(median(redd_EHVO_lowv))
print(median(redd_EHVO_highv))
print(median(lbol_EHVO_lowv))
print(median(lbol_EHVO_highv))
print(median(mbh_EHVO_lowv))
print(median(mbh_EHVO_highv))