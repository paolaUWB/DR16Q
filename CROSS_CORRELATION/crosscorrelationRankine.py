#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 01:31:21 2023

@author: alex
"""

#%%

import sys
import scipy.stats as stat
import numpy as np
from pylab import*
from sympy import sympify
from matplotlib.backends.backend_pdf import PdfPages
from scipy import*
from astropy import*
from scipy.stats import ks_2samp
from astropy import constants as const
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import stats
from utility_functions import read_spectra, append_row_to_csv, clear_file, print_to_file, read_list_spectra
import os
from os.path import exists
import csv
import pandas as pd

#Number of samples 

#specnumDR16=XXX
#specnumDR14=XXX
specnumEHVO=98 
specnumPARENT=int(18165)

wnewfile = 'no'
# wnewfile = 'yes'

#The inputs in this program should be:
infoDR16 = os.getcwd() + "/../../DR16Q_v4.fits" #DR16 fits file
infoDR14 = os.getcwd() + "/../../dr14q_spec_prop.fits"  #DR14Q table fits file from Rakshit+2020

#Rankines info file PRH selected the good cases
infoRankineparent = os.getcwd() + "/DR16parent_DR14RankineInfo.csv"
infoRankineEHVO = os.getcwd() + "/DR16EHVO_DR14RankineInfo.csv"


#%%
#Reading DR14 & 16 fits table files

hdu_16 = fits.open(infoDR16)
data_16 = hdu_16[1].data
#print(data_16.columns)
SDSS_name_16 = data_16['SDSS_NAME']
plate_16 = data_16['PLATE  ']
mjd_16 = data_16['MJD']
fiber_16 = data_16['FIBERID ']
SNR_16 = data_16['SN_MEDIAN_ALL']
redshift_16 = data_16['z']
fits_16__duplicate_PLATE = data_16['PLATE_DUPLICATE']
fits_16__duplicate_MJD = data_16['MJD_DUPLICATE']
fits_16__duplicate_FIBER = data_16['FIBERID_DUPLICATE']
hdu_16.close()




hdu_14 = fits.open(infoDR14)
data_14 = hdu_14[1].data
# print(data_14.columns)
SDSS_name_14 = data_16['SDSS_NAME']
plate_14 = data_14['PLATE  ']
mjd_14 = data_14['MJD     ']
fiber_14 = data_14['FIBERID ']
log_mbh_14 = data_14['LOG_MBH ']
log_mbh_err_14 = data_14['LOG_MBH_ERR']
q_mbh_14 = data_14['QUALITY_MBH']
log_lbol_14 = data_14['LOG_LBOL']
q_lbol_14 = data_14['QUALITY_LBOL']
log_redd_14 = data_14['LOG_REDD']
q_redd_14 = data_14['QUALITY_REDD']
bi_civ_14 = data_14['BI_CIV  ']
bi_civ_err_14 = data_14['ERR_BI_CIV']
bal_flag_14 = data_14['BAL_FLAG']
hdu_14.close()



# BALQSO_index = np.where(bi_civ_14 > 0)
# BALQSO = 'no'


#%%
#Extracting values from Rankines info file for EHVO Parent
dfRPA = pd.read_csv(infoRankineparent, header=None)
dfRHV = pd.read_csv(infoRankineEHVO, header=None)

#for parent sample
Parentin14Rank_mbh=dfRPA[dfRPA.columns[16]].to_numpy()
Parentin14Rank_lbol=dfRPA[dfRPA.columns[17]].to_numpy()
Parentin14Rank_redd=dfRPA[dfRPA.columns[18]].to_numpy()

#for EHVOs
EHVOin14Rank_mbh=dfRHV[dfRHV.columns[16]].to_numpy()
EHVOin14Rank_lbol=dfRHV[dfRHV.columns[17]].to_numpy()
EHVOin14Rank_redd=dfRHV[dfRHV.columns[18]].to_numpy()


MBH_parentR = Parentin14Rank_mbh
Lbol_parentR = Parentin14Rank_lbol
Redd_parentR = Parentin14Rank_redd

MBH_EHVOR = EHVOin14Rank_mbh
Lbol_EHVOR = EHVOin14Rank_lbol
Redd_EHVOR = EHVOin14Rank_redd 

#%%
#Rankines info file for BAL (10k, 25k)
# bals:  16 17 18 
#bal pos:
bi_bi0 = dfRPA[dfRPA.columns[10]].to_numpy()
bi_vmax = dfRPA[dfRPA.columns[11]].to_numpy()
bi_vmin = dfRPA[dfRPA.columns[12]].to_numpy()
SNR = dfRPA[dfRPA.columns[15]].to_numpy()

pos_bal = np.where((bi_bi0>0) & (SNR > 10))
pos_10k = np.where((bi_vmax<10000)&(bi_bi0>0)&(SNR>10))
pos_25k = np.where((bi_vmax>10000)&(bi_vmax<25000)&(bi_bi0>0)&(SNR>10))

mbh_bal = Parentin14Rank_mbh[pos_bal]
lbol_bal = Parentin14Rank_lbol[pos_bal]
redd_bal = Parentin14Rank_redd[pos_bal]

mbh_bal10k,mbh_bal25k   = Parentin14Rank_mbh[pos_10k], Parentin14Rank_mbh[pos_25k]
redd_bal10k, redd_bal25k = Parentin14Rank_redd[pos_10k], Parentin14Rank_redd[pos_25k]
lbol_bal10k, lbol_bal25k = Parentin14Rank_lbol[pos_10k], Parentin14Rank_lbol[pos_25k]

# pos_bal = np.where(bi0>0, SNR>10)
                   # (vmin>10k, vmin<vmax<25k)

                  
#%% 
#Optimal Bin 
# comb_mbhdata = np.concatenate([MBH_EHVOR , MBH_parentR])
# comb_edddata = np.concatenate([Redd_EHVOR , Redd_parentR])

# Parent sample knuth  
op_binmbh_knuthp = (stats.knuth_bin_width(MBH_EHVOR)+stats.knuth_bin_width(MBH_parentR))/2
op_binedd_knuthp = (stats.knuth_bin_width(Redd_EHVOR)+stats.knuth_bin_width(Redd_parentR))/2
op_binlbol_knuthp = (stats.knuth_bin_width(Lbol_EHVOR)+stats.knuth_bin_width(Lbol_parentR))/2

# Bal sample knuth
op_binmbh_knuthbal = (stats.knuth_bin_width(MBH_EHVOR)+stats.knuth_bin_width(mbh_bal))/2
op_binedd_knuthbal = (stats.knuth_bin_width(Redd_EHVOR)+stats.knuth_bin_width(redd_bal))/2
op_binlbol_knuthbal = (stats.knuth_bin_width(Lbol_EHVOR)+stats.knuth_bin_width(lbol_bal))/2

# op_binmbh_knuthbal = (stats.knuth_bin_width(MBH_EHVOR)+stats.knuth_bin_width(mbh_bal10k)+stats.knuth_bin_width(mbh_bal25k))/3
# op_binedd_knuthbal = (stats.knuth_bin_width(Redd_EHVOR)+stats.knuth_bin_width(redd_bal10k)+stats.knuth_bin_width(redd_bal25k))/3
# op_binlbol_knuthbal = (stats.knuth_bin_width(Lbol_EHVOR)+stats.knuth_bin_width(lbol_bal10k)+stats.knuth_bin_width(lbol_bal25k))/3

# BAL sample freedman
op_binmbh_freedmanbal = (stats.freedman_bin_width(MBH_EHVOR)+stats.freedman_bin_width(mbh_bal10k)+stats.freedman_bin_width(mbh_bal10k))/3
op_binedd_freedmanbal = (stats.freedman_bin_width(Redd_EHVOR)+stats.freedman_bin_width(redd_bal10k)+stats.freedman_bin_width(redd_bal25k))/3
op_binlbol_freedmanbal = (stats.freedman_bin_width(Lbol_EHVOR)+stats.freedman_bin_width(lbol_bal10k)+stats.freedman_bin_width(lbol_bal25k))/3

# parent sample freedman
op_binmbh_freedmanp = (stats.freedman_bin_width(MBH_EHVOR)+stats.freedman_bin_width(MBH_parentR))/2
op_binedd_freedmanp = (stats.freedman_bin_width(Redd_EHVOR)+stats.freedman_bin_width(Redd_parentR))/2
op_binlbol_freedmanp = (stats.freedman_bin_width(Lbol_EHVOR)+stats.freedman_bin_width(Lbol_parentR))/2


#%%
#Histogram:


def scatter_hist2(x, y, ax, ax_histx, ax_histy, color, area, mult, factor, ax_set, limx, limy,binwidth_x,binwidth_y):
    # no labels
    ax_histx.tick_params(axis='x', labelbottom=False)
    ax_histy.tick_params(axis='y', labelleft=False)
    # the scatter plot:
    #ax.set_xlim([7.5, 11.25]) #MBH
    #ax.set_ylim([-2.0, 1.0]) #Redd
    ax.set_ylim([ylowlim, yuplim])
    ax.set_xlim([xlowlim, xuplim]) 
    ax.scatter(x, y, s = area, color = color)
    ax.text(10.5,0.5,'')
    # now determine nice limits by hand:
    # binwidth_x = op_binmbh
    # binwidth_y = op_binlbol
    #limx limy redefine lines
    # limx = 11.25
    # limx = 48.3
    # limy = 48.3
    binsx = np.arange(-limx, limx, binwidth_x)
    binsy = np.arange(-limy, limy, binwidth_y)
    
    if mult == 'yes':
        for i in range(0, factor):
            x = np.append(x,x)
            y = np.append(y,y)
        x = np.hstack(x)
        y = np.hstack(y)
    #kwargs = dict(histtype=‘stepfilled’, alpha=0.5, normed=True, bins=binsx)
    weightsx = [np.ones_like(x)/float(len(x))]
    print('Sum =', np.sum(weightsx))
    weightsy = [np.ones_like(y)/float(len(y))]
    print('Sum =', np.sum(weightsy))
    ax_histx.hist(x, bins = binsx, weights = weightsx, color = color, alpha = 0.5)
    ax_histy.hist(y, bins = binsy, weights = weightsy, orientation='horizontal', color=color, alpha = 0.5)
 # definitions for the axes



# PLOTS paramters

left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.005
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]


reddlim = 48.3
masslim = 11.25
lbollim = 48.3

# reddlimbal = (max(Redd_EHVOR) + max(redd_bal10k) + max(redd_bal25k))/3
# masslimbal = (max(MBH_EHVOR) + max(mbh_bal10k) + max(mbh_bal25k))/3
# lbollimbal = (max(Lbol_EHVOR) + max(lbol_bal10k)+ max(lbol_bal25k))/3

# stop
#%%
#Property values:
xpm = MBH_parentR
xem = MBH_EHVOR
x10km = mbh_bal10k
x25km = mbh_bal25k

ypr = Redd_parentR
yer = Redd_EHVOR
y10kr = redd_bal10k
y25kr = redd_bal25k

xpl = Lbol_parentR
xel = Lbol_EHVOR
x10kl = lbol_bal10k
x25kl = lbol_bal25k

#%%
# Plot of  redd vs mass

fig = plt.figure(1)
fig = plt.figure(figsize=(8, 8))


#setting axis limits based on properties
# xlowlim = 46.
# xuplim = 48.3
xlowlim = 8.20
xuplim = 10.25


# ylowlim = 7.5
# yuplim = 11.25
ylowlim = -2.0
yuplim = 1.0

ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

# ax_set_2 = ax.set(ylabel=r'log($M_{\mathrm{BH}}/M_{\odot})$', xlabel='log($L_{bol}$/erg s$^{-1}$)')
ax_set_2 = ax.set(ylabel='Eddington Ratio (log($L_{bol}$/$L_{Edd}$)', xlabel=r'log($M_{\mathrm{BH}}/M_{\odot})$')

ax_set_ylim_Redd = ax.set_ylim([ylowlim, yuplim])
ax_set_xlim_Lbol = ax.set_xlim([xlowlim, xuplim])


# use the previously defined function
scatter_hist2(x10km, y10kr, ax, ax_histx, ax_histy,'aqua', 5, 'yes', 2, ax_set_2, limx=masslim, limy=reddlim, binwidth_x = op_binmbh_knuthbal, binwidth_y = op_binedd_knuthbal )
scatter_hist2(x25km, y25kr, ax, ax_histx, ax_histy ,'cornflowerblue', 5, 'yes', 2, ax_set_2,limx=masslim, limy=reddlim, binwidth_x = op_binmbh_knuthbal, binwidth_y = op_binedd_knuthbal)

# scatter_hist2(xpm, ypr, ax, ax_histx, ax_histy,'cornflowerblue', 5, 'yes', 2, ax_set_2, limx=11.25, limy=48.3,binwidth_x = op_binmbh_knuthp,binwidth_y = op_binedd_knuthp )
scatter_hist2(xem, yer, ax, ax_histx, ax_histy ,'purple', 60, 'yes', 10, ax_set_2,limx=11.25, limy=48.3,binwidth_x = op_binmbh_knuthp, binwidth_y = op_binedd_knuthp)

# ax.add_artist(ax.legend(title='test'))
ax.legend(['Rankine+2020 BALs vmin<10,000 km/s','Rankine+2020 BALs 10,000 km/s<vmin<25,000 km/s','RH+(in prep) EHVO'],loc='lower left')
# ax.legend(['Rankine+2020 Parent sample','RH+(in prep) EHVO'],loc='lower left')
ax.set_title("Title for second plot")
# plt.savefig(PLOT_DIREC + 'Redd_Lbol.png') 
plt.show()
plt.close()

#%%
# Plot of  redd vs Lbol

fig = plt.figure(1)
fig = plt.figure(figsize=(8, 8))

xlowlim = 46.
xuplim = 48.3
# xlowlim = 8.20
# xuplim = 10.25


# ylowlim = 7.5
# yuplim = 11.25
ylowlim = -2.0
yuplim = 1.0

ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

# ax_set_2 = ax.set(ylabel=r'log($M_{\mathrm{BH}}/M_{\odot})$', xlabel='log($L_{bol}$/erg s$^{-1}$)')
ax_set_2 = ax.set(ylabel='Eddington Ratio (log($L_{bol}$/$L_{Edd}$)', xlabel='log($L_{bol}$/erg s$^{-1}$)')

ax_set_ylim_Redd = ax.set_ylim([ylowlim, yuplim])
ax_set_xlim_Lbol = ax.set_xlim([xlowlim, xuplim])


# use the previously defined function
scatter_hist2(x10kl, y10kr, ax, ax_histx, ax_histy,'aqua', 5, 'yes', 2, ax_set_2,limx=lbollim, limy=reddlim,binwidth_x = op_binlbol_knuthbal, binwidth_y = op_binedd_knuthbal)
scatter_hist2(x25kl, y25kr, ax, ax_histx, ax_histy ,'cornflowerblue', 5, 'yes', 2, ax_set_2,limx=lbollim, limy=reddlim,binwidth_x = op_binlbol_knuthbal, binwidth_y = op_binedd_knuthbal)

# scatter_hist2(xpl, ypr, ax, ax_histx, ax_histy,'cornflowerblue', 5, 'yes', 2, ax_set_2, limx=48.3, limy=48.3, binwidth_x = op_binmbh_knuthp, binwidth_y = op_binedd_knuthp)
scatter_hist2(xel, yer, ax, ax_histx, ax_histy ,'purple', 60, 'yes', 10, ax_set_2,limx=lbollim, limy=reddlim, binwidth_x = op_binlbol_knuthp, binwidth_y = op_binedd_knuthp)

# ax.add_artist(ax.legend(title='test'))
ax.legend(['Rankine+2020 BALs vmin<10,000 km/s','Rankine+2020 BALs 10,000 km/s<vmin<25,000 km/s','RH+(in prep) EHVO'],loc='lower left')
# ax.legend(['Rankine+2020 Parent sample','RH+(in prep) EHVO'],loc='lower left')
ax.set_title("Title for second plot")
# plt.savefig(PLOT_DIREC + 'Redd_Lbol.png') 
plt.show()
plt.close()



