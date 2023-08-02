#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 01:31:21 2023

@author: alex
"""

#%%

import ast
import sys
import scipy.stats as stat
import numpy as np
from pylab import*
from sympy import sympify
from matplotlib.backends.backend_pdf import PdfPages
from scipy import*
from astropy import*
from scipy.stats import ks_2samp, anderson_ksamp
from astropy import constants as const
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import stats
# from utility_functions import read_spectra, append_row_to_csv, clear_file, print_to_file, read_list_spectra
import os
from os.path import exists
import csv
import pandas as pd
from tabulate import tabulate
#Number of samples 

#specnumDR16=XXX
#specnumDR14=XXX
# specnumEHVO=98 
# specnumPARENT=int(18165)

# wnewfile = 'no'
# wnewfile = 'yes'

#The inputs in this program should be:
infoDR16 = os.getcwd() + "/../../DR16Q_v4.fits" #DR16 fits file
infoDR14 = os.getcwd() + "/../../dr14q_spec_prop.fits"  #DR14Q table fits file from Rakshit+2020

#Rankines info file PRH selected the good cases
infoRankineparent = os.getcwd() + "/DR16parent_DR14RankineInfo.csv"
infoRankineEHVO = os.getcwd() + "/DR16EHVO_DR14RankineInfo_wEHVOspeed.csv"

infoPRHEHVO = os.getcwd() + "/DR16_EHVO_sorted_norm.csv"
infoPRHparent = os.getcwd() + "/DR16_parent_sample.csv"

filebal = '/Users/alex/Documents/GitHub_copy1/DR16Q/DR16_BAL_parent_sample.csv'
dfbal= pd.read_csv(filebal, header=None)

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
# filebal = '/Users/alex/Documents/GitHub/DR16Q/DR16_BAL_parent_sample.csv'
# dfbal= pd.read_csv(filebal, header=None)


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
EHVOin14Rank_z=dfRHV[dfRHV.columns[4]].to_numpy()

MBH_parentR = Parentin14Rank_mbh
Lbol_parentR = Parentin14Rank_lbol
Redd_parentR = Parentin14Rank_redd


MBH_EHVOR = EHVOin14Rank_mbh
Lbol_EHVOR = EHVOin14Rank_lbol
Redd_EHVOR = EHVOin14Rank_redd 

z_parentR = dfRPA[dfRPA.columns[4]].to_numpy()
z_EHVOR = EHVOin14Rank_z





#%%
#Rankines info file for BAL (10k, 25k)
# bals:  16 17 18 
#bal pos:
bi_bi0 = dfRPA[dfRPA.columns[10]].to_numpy()
bi_vmax = dfRPA[dfRPA.columns[11]].to_numpy()
bi_vmin = dfRPA[dfRPA.columns[12]].to_numpy()
SNR = dfRPA[dfRPA.columns[15]].to_numpy()

pos_bal = np.where((bi_bi0>0) & (SNR > 10))
pos_10k = np.where((bi_vmin<6000)&(bi_bi0>0)&(SNR>10))
pos_25k = np.where((bi_vmin>6000)&(bi_vmin<25000)&(bi_bi0>0)&(SNR>10))

pos_15k = np.where((bi_vmin<13000)&(bi_bi0>0)&(SNR>10))
pos_15kplus = np.where((bi_vmin>13000)&(bi_vmin<24000)&(bi_bi0>0)&(SNR>10))

mbh_bal = Parentin14Rank_mbh[pos_bal]
lbol_bal = Parentin14Rank_lbol[pos_bal]
redd_bal = Parentin14Rank_redd[pos_bal]
z_bal = z_parentR[pos_bal]
# stop

mbh_bal10k,mbh_bal25k   = Parentin14Rank_mbh[pos_10k], Parentin14Rank_mbh[pos_25k]
redd_bal10k, redd_bal25k = Parentin14Rank_redd[pos_10k], Parentin14Rank_redd[pos_25k]
lbol_bal10k, lbol_bal25k = Parentin14Rank_lbol[pos_10k], Parentin14Rank_lbol[pos_25k]

mbh_bal15k,mbh_bal15kplus   = Parentin14Rank_mbh[pos_15k], Parentin14Rank_mbh[pos_15kplus]
redd_bal15k, redd_bal15kplus = Parentin14Rank_redd[pos_15k], Parentin14Rank_redd[pos_15kplus]
lbol_bal15k, lbol_bal15kplus = Parentin14Rank_lbol[pos_15k], Parentin14Rank_lbol[pos_15kplus]

# stop()

# pos_bal = np.where(bi0>0, SNR>10)
                   # (vmin>10k, vmin<vmax<25k)
                   
#%%
# Analysis with the speed of EHVOs 

EHVOin14Rank_maxv=dfRHV[dfRHV.columns[33]].to_numpy()
EHVOin14Rank_minv=dfRHV[dfRHV.columns[29]].to_numpy()
# xx = list(ast.literal_eval(x) for x in EHVOin14Rank_maxv)

# ff = pd.DataFrame(xx)

pos_lbol40k = np.where(EHVOin14Rank_minv > -37000)
# pos_lbol48k = np.where((EHVOin14Rank_maxv < -40000) & (EHVOin14Rank_maxv > -45000))
pos_lbol48kmore = np.where((EHVOin14Rank_minv < -37000))

ehvoin14rank_speedlbol40K = EHVOin14Rank_lbol[pos_lbol40k]
# ehvoin14rank_speedlbol48K = EHVOin14Rank_lbol[pos_lbol48k]
ehvoin14rank_speedlbol48Kmore = EHVOin14Rank_lbol[pos_lbol48kmore]

ehvoin14rank_edd43k = EHVOin14Rank_redd[pos_lbol40k]
ehvoin14rank_edd43kbeyond = EHVOin14Rank_redd[pos_lbol48kmore]

ehvoin14rank_mbh43k = EHVOin14Rank_mbh[pos_lbol40k]
ehvoin14rank_mbh43kbeyond = EHVOin14Rank_mbh[pos_lbol48kmore]

ehvoin14rank_speed40K = EHVOin14Rank_maxv[pos_lbol40k]
# ehvoin14rank_speed48K = EHVOin14Rank_maxv[pos_lbol48k]
ehvoin14rank_speed48Kmore = EHVOin14Rank_maxv[pos_lbol48kmore]



op_bin40k_knuthp = ((stats.knuth_bin_width(ehvoin14rank_speedlbol40K) 
                    + stats.knuth_bin_width(ehvoin14rank_speedlbol48Kmore)) /2)
# op_bin48k_knuthp = (stats.knuth_bin_width(ehvoin14rank_speedlbol48K))
op_bin48kmore_knuthp = (stats.knuth_bin_width(ehvoin14rank_speedlbol48Kmore))



op_binspd40k_knuthp = (stats.knuth_bin_width(ehvoin14rank_speed40K/-1000))
# op_binspd48k_knuthp = (stats.knuth_bin_width(ehvoin14rank_speed48K/-1000))
op_binspd48kmore_knuthp = (stats.knuth_bin_width(ehvoin14rank_speed48Kmore/-1000))

# jj = np.arange(0,70)

scale_fac = 15

ehvorescale_lbol43k = np.repeat(ehvoin14rank_speedlbol40K, scale_fac)
ehvorescale_lbol43kplus = np.repeat(ehvoin14rank_speedlbol48Kmore, scale_fac)

ehvorescale_mbh43k = np.repeat(ehvoin14rank_mbh43k, scale_fac)
ehvorescale_mbh43kplus = np.repeat(ehvoin14rank_mbh43kbeyond, scale_fac)

ehvorescale_edd43k = np.repeat(ehvoin14rank_edd43k, scale_fac)
ehvorescale_edd43kplus = np.repeat(ehvoin14rank_edd43kbeyond, scale_fac)

xx = np.concatenate((lbol_bal15k, lbol_bal15kplus, ehvoin14rank_speedlbol40K, ehvoin14rank_speedlbol48Kmore))

sturges_bin1 = (1+3.22*(math.log(len(lbol_bal15kplus))))
knuthbin1 = 0.5*(stats.knuth_bin_width(lbol_bal15k) + stats.knuth_bin_width(ehvorescale_lbol43kplus))

# bin_spd = np.linspace(min(lbol_bal15k), max(ehvorescale_lbol43kplus), int(sturges_bin1))
bin_spd1 = np.arange(min(lbol_bal15k), max(ehvorescale_lbol43kplus)+knuthbin1, knuthbin1 )

# bin_spd1 = np.arange(min(lbol_bal15k), max(ehvorescale_lbol43kplus) + knuthbin1, knuthbin1)

# sturges_bin2 = (1+3.22*(math.log(len(mbh_bal15kplus))))
knuthbin2 = 0.5*(stats.knuth_bin_width(mbh_bal15k) + stats.knuth_bin_width(ehvorescale_mbh43kplus))
# bin_mbh = np.linspace(min(mbh_bal15k), max(ehvorescale_mbh43kplus), int(sturges_bin2))
bin_mbh1 = np.arange(min(mbh_bal15k), max(ehvorescale_mbh43kplus)+knuthbin2, knuthbin2 )


# sturges_bin3 = (1+3.22*(math.log(len(redd_bal15kplus))))
knuthbin3 = 0.5*(stats.knuth_bin_width(redd_bal15k) + stats.knuth_bin_width(ehvorescale_edd43kplus))
# bin_redd = np.linspace(min(redd_bal15k), max(ehvorescale_edd43kplus), int(sturges_bin3))

bin_redd1 = np.arange(min(redd_bal15k), max(ehvorescale_edd43kplus)+knuthbin3, knuthbin3 )

# plt.hist(lbol_bal15k,bins=bin_spd1, label='BAL <15,000 km/s', color='lightskyblue', alpha = 0.5,
#           )
# plt.hist(lbol_bal15kplus, bins=bin_spd1,label='BAL 15,000 km/s<vmin<25,000 km/s', color='dodgerblue',alpha = 0.6,
#           )
# plt.hist(ehvorescale_lbol43k, bins=bin_spd1, label='EHVO 30,000 km/s<vmin<37,000 km/s', color= 'magenta', alpha=0.6,
#           )
# plt.hist(ehvorescale_lbol43kplus, bins=bin_spd1, label='EHVO vmin>37,000 km/s', color='darkmagenta', alpha = 0.8,
#           )
# plt.xlabel('log($L_{bol}$/erg s$^{-1}$)')
# plt.ylabel('N')
# plt.legend(loc='upper right',prop={'size': 6})    
# plt.savefig('spdlbol.png', dpi=1200)



# plt.hist(mbh_bal15k,bins=bin_mbh1, label='BAL <15,000 km/s', color='lightskyblue', alpha=0.6 )
# plt.hist(mbh_bal15kplus,bins=bin_mbh1, label='BAL 15,000 km/s<vmin<25,000 km/s', color='dodgerblue', alpha=0.6)
# plt.hist(ehvorescale_mbh43k, bins=bin_mbh1, label='EHVO 30,000 km/s<vmin<37,000 km/s', color= 'magenta', alpha=0.6,
#           )
# plt.hist(ehvorescale_mbh43kplus, bins=bin_mbh1, label='EHVO vmin>37,000 km/s', color='darkmagenta', alpha = 0.8,
#           )
# plt.xlabel('log($M_{\mathrm{BH}}/M_{\odot})$')
# plt.ylabel('N')
# plt.legend(loc='upper left',prop={'size': 6})    
# plt.savefig('spdmbh.png', dpi=1200)



plt.hist(redd_bal15k, bins=bin_redd1, label='BAL <15,000 km/s',color='lightskyblue',  alpha = 0.6 )
plt.hist(redd_bal15kplus, bins=bin_redd1,label='BAL 15,000 km/s<vmin<25,000 km/s',color='dodgerblue', alpha=0.6)
plt.hist(ehvorescale_edd43k, bins=bin_redd1, label='EHVO 30,000 km/s<vmin<37,000 km/s', color= 'magenta', alpha=0.6,
          )
plt.hist(ehvorescale_edd43kplus,bins=bin_redd1,  label='EHVO vmin>37,000 km/s', color='darkmagenta', alpha = 0.8,
          )
plt.xlabel('Eddington Ratio (log($L_{bol}$/$L_{Edd}$)')
plt.ylabel('N')
plt.legend(loc='upper left',prop={'size': 6})    
plt.savefig('spdredd.png', dpi=1200)



#%% Median 
import statistics as sta

medehvo37lbol = sta.median(ehvoin14rank_speedlbol40K)
medehvo37pluslbol = sta.median(ehvoin14rank_speedlbol48Kmore)

medehvo37mbh = sta.median(ehvoin14rank_mbh43k)
medehvo37plusmbh = sta.median(ehvoin14rank_mbh43kbeyond)

medehvo37redd = sta.median(ehvoin14rank_edd43k)
medehvo37plusredd = sta.median(ehvoin14rank_edd43kbeyond)

medbal37lbol = sta.median(lbol_bal15k)
medbal37pluslbol = sta.median(lbol_bal15kplus)

medbal37mbh = sta.median(mbh_bal15k)
medbal37plusmbh = sta.median(mbh_bal15kplus)

medbal37redd = sta.median(redd_bal15k)
medbal37plusredd = sta.median(redd_bal15kplus)

datarow0 = [['vmin<15,000 km/s','15,000km/s<vmin<25,000 km/s', '30,000km/s<vmin<37,000 km/s', 'vmin>37,000 km/s'],
            [medbal37lbol, medbal37pluslbol, medehvo37lbol, medehvo37pluslbol],
            [medbal37mbh, medbal37plusmbh, medehvo37mbh, medehvo37plusmbh],
            [medbal37redd, medbal37plusredd, medehvo37redd, medehvo37plusredd]
            ]
    

# dataheaders = ['','Lbol','Mbh','Eddington ratio']
datacol = ['Outflow velocity','Lbol', 'MBH', 'Eddington Ratio']

plt.table(cellText= datarow0, rowLabels = ['BAL15k-BAL15k+', 'EHVO37k-EHVO37k+'])

plt.rcParams["figure.figsize"] = [10, 5]
plt.rcParams["figure.autolayout"] = True
fig, axs = plt.subplots(1, 1)
# data = np.random.random((10, 3))
# columns = ("Column I", "Column II", "Column III")
axs.axis('tight')
axs.axis('off')
the_table = axs.table(cellText=datarow0, 
                      rowLabels = datacol,  loc='center')
# the_table.set_fontsize(24)

# plt.show()
# plt.savefig('mediantable.png', dpi = 400)
stop()

#%%
from draw_histogram import plot_zem_histograms

dfDR16Q = pd.read_csv(infoPRHparent, header=None)
dfEHVO = pd.read_csv(infoPRHEHVO, header=None)

pos_z = np.where((z_parentR == min(z_parentR)))
plate_z, mjd_z, fibre_z = dfRPA[dfRPA.columns[1]].to_numpy()[pos_z], dfRPA[dfRPA.columns[2]].to_numpy()[pos_z], dfRPA[dfRPA.columns[3]].to_numpy()[pos_z]
outlier = [plate_z, mjd_z, fibre_z]



z_dr16q = dfDR16Q[dfDR16Q.columns[1]].to_numpy()
z_EHVO = dfEHVO[dfEHVO.columns[4]].to_numpy()

pos_out = np.where(dfDR16Q[dfDR16Q.columns[0]].to_numpy() =='spec-8160-57071-0650-dered.dr16')
z_out = z_dr16q[pos_out]

z_checkpos = np.where(z_dr16q == z_out)
z_outsdssinfo = dfDR16Q[dfDR16Q.columns[0]].to_numpy()[z_checkpos]

print('z outlier Rankine value:', min(z_parentR), 'outlier sdss info:', outlier)
print('z outlier PRH value:', z_out, 'outlier sdss info:', z_outsdssinfo)

# stop()

#### rankine z valie
OUT_CONFIG = os.getcwd()

zem_parentr = np.hstack(z_parentR)
zem_ehvo = np.hstack(z_EHVOR)
zem_bal = np.hstack(z_bal)


#bin size selection
# sturges_bin = (1+3.22*(math.log(len(zem_parentr))))
# bin_z = np.linspace(min(z_parentR), max(zem_ehvo), int(sturges_bin))

# x_vals = [zem_ehvo, zem_parentr, zem_bal]
# x_label = '$z_{em}$'
# color = ['xkcd:shocking pink','black', 'xkcd:purpleish blue']
# label = ['EHVO', 'parent', 'BAL']

# zem_bal_size = 1
# zem_ehvo_size = 1

# plot_zem_histograms(x_vals, bin_z, color, label, zem_bal_size, zem_ehvo_size, 
#                     x_label, 'zem/class', 'right', OUT_CONFIG + 'zem_DR16_zem_test.png')




#######  PRH z value
# sturges_bin = (1+3.22*(math.log(len(z_dr16q))))
# bin_z = np.linspace(min(z_dr16q), max(z_EHVO), int(sturges_bin))

# x_vals = [z_EHVO, z_dr16q, zem_bal]
# x_label = '$z_{em}$'
# color = ['xkcd:shocking pink','black', 'xkcd:purpleish blue']
# label = ['EHVO', 'parent', 'BAL']

# zem_bal_size = 1
# zem_ehvo_size = 1

# plot_zem_histograms(x_vals, bin_z, color, label, zem_bal_size, zem_ehvo_size, 
#                     x_label, 'zem/class', 'right', OUT_CONFIG + 'zem_DR16_zem_test.png')

# sturges_bin = (1+3.22*(math.log(len(lbol_bal15kplus))))
# bin_spd = np.linspace(min(lbol_bal15k), max(ehvorescale_lbol43kplus), int(sturges_bin))


# plt.plot(z_parentR,Lbol_parentR)
# plt.legend()

# stop()


#%%
# K-S test
pa_ehvo_lbol = stat.ks_2samp(Lbol_parentR, Lbol_EHVOR, mode = 'auto')
bal_ehvo_lobl = stat.ks_2samp(lbol_bal, Lbol_EHVOR, mode = 'auto')
# pa_bal = stat.ks_2samp(lbol_bal, Lbol_parentR, mode = 'auto')
pa_ehvo_edd = stat.ks_2samp(Redd_parentR, Redd_EHVOR, mode = 'auto')
bal_ehvo_edd = stat.ks_2samp(redd_bal, Redd_EHVOR, mode = 'auto')

pa_ehvo_mbh = stat.ks_2samp(MBH_parentR, MBH_EHVOR, mode = 'auto')
bal_ehvo_mbh = stat.ks_2samp(mbh_bal, MBH_EHVOR, mode = 'auto')

## K-S test Parent vs EHVO vs BALQSOs

# tabledata = [['EHVO-Parent', pa_ehvo_lbol[1], pa_ehvo_mbh[1], pa_ehvo_edd[1]],
#              ['EHVO-BALQSOs',bal_ehvo_lobl[1], bal_ehvo_mbh[1], bal_ehvo_edd[1]]]

# table_columns = ['samples', 'Lbol', 'Mbh', 'Lbol/Ledd' ]

# print(tabulate(tabledata, headers =table_columns, tablefmt = 'fancy_grid'))
# datarow1 = [[pa_ehvo_lbol[1], pa_ehvo_mbh[1], pa_ehvo_edd[1]],
#            [bal_ehvo_lobl[1], bal_ehvo_mbh[1], bal_ehvo_edd[1]]]

# dataheaders = ['Lbol','Mbh','Eddington ratio']

# plt.table(cellText= datarow1, rowLabels = ['EHVO-Parent', 'EHVO-BALQSOs'], colLabels= dataheaders)
# plt.box(on=None)
# plt.show

# Anderson Darling Test
ad_pa_ehvo = stat.anderson_ksamp([Lbol_parentR, Lbol_EHVOR])
ad_bal_ehvo = stat.anderson_ksamp([lbol_bal, Lbol_EHVOR])
# ad_pa_bal = stat.anderson_ksamp([lbol_bal, Lbol_parentR])

ad_pa_ehvo_mbh = stat.anderson_ksamp([MBH_parentR, MBH_EHVOR])
ad_bal_ehvo_mbh = stat.anderson_ksamp([mbh_bal, MBH_EHVOR])

ad_pa_ehvo_redd = stat.anderson_ksamp([Redd_parentR, Redd_EHVOR])
ad_bal_ehvo_redd = stat.anderson_ksamp([redd_bal, Redd_EHVOR])

datarow2 = [[ad_pa_ehvo[2], ad_pa_ehvo_mbh[2], ad_pa_ehvo_redd[2]],
           [ad_bal_ehvo[2], ad_pa_ehvo_mbh[2], ad_bal_ehvo_redd[2]]]



## K-S test Outflow Speed Analysis 
bal15_bal15plus_lbol = stat.ks_2samp(lbol_bal15k, lbol_bal15kplus, mode = 'auto')
ehvo37_ehvo37plus_lobl = stat.ks_2samp(ehvoin14rank_speedlbol40K, ehvoin14rank_speedlbol48Kmore, mode = 'auto')
# pa_bal = stat.ks_2samp(lbol_bal, Lbol_parentR, mode = 'auto')
bal15_bal15plus_edd = stat.ks_2samp(redd_bal15k, redd_bal15kplus, mode = 'auto')
ehvo37_ehvo37plus_edd = stat.ks_2samp(ehvoin14rank_edd43k, ehvoin14rank_edd43kbeyond, mode = 'auto')

bal15_bal15plus_mbh = stat.ks_2samp(mbh_bal15k, mbh_bal15kplus, mode = 'auto')
ehvo37_ehvo37plus_mbh = stat.ks_2samp(ehvoin14rank_mbh43k, ehvoin14rank_mbh43kbeyond, mode = 'auto')

datarow1 = [['Lbol','Mbh','Eddington ratio'],
            [round(bal15_bal15plus_lbol[1],11), round(bal15_bal15plus_mbh[1],2), round(bal15_bal15plus_edd[1],11)],
           [round(ehvo37_ehvo37plus_lobl[1],11), round(ehvo37_ehvo37plus_mbh[1],11), round(ehvo37_ehvo37plus_edd[1],11)]]

# dataheaders = ['','Lbol','Mbh','Eddington ratio']
datarow = ['Outflow class','BAL15k-BAL15k+', 'EHVO37k-EHVO37k+']

# plt.table(cellText= datarow1, rowLabels = ['BAL15k-BAL15k+', 'EHVO37k-EHVO37k+'], colLabels= dataheaders)

plt.rcParams["figure.figsize"] = [10.00, 5.0]
plt.rcParams["figure.autolayout"] = True
fig, axs = plt.subplots(1, 1)
# data = np.random.random((10, 3))
# columns = ("Column I", "Column II", "Column III")
axs.axis('tight')
axs.axis('off')
the_table = axs.table(cellText=datarow1, 
                      rowLabels = datarow,  loc='center')
plt.show()
plt.savefig('KSoutflowspd.png')

stop

# print(ad_bal_ehvo)
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

#ehvo 10k 25k bin
# op_bin_knuthehvo10k = (stats.knuth_bin_width(lbol_ehvo_10k))
# op_bin_knuthehvo25k = (stats.knuth_bin_width(lbol_ehvo_25k))


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
    # limy = 
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
spacing = 0.01
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]


reddlim = 48.3
masslim = 11.25
lbollim = 48.3
speedlim = 60

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
xlowlim = 8.0
xuplim = 10.49


# ylowlim = 7.5
# yuplim = 11.25
ylowlim = -2.0
yuplim = 0.99

ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

# ax_set_2 = ax.set(ylabel=r'log($M_{\mathrm{BH}}/M_{\odot})$', xlabel='log($L_{bol}$/erg s$^{-1}$)')
ax_set_2 = ax.set(ylabel='Eddington Ratio (log($L_{bol}$/$L_{Edd}$)', xlabel=r'log($M_{\mathrm{BH}}/M_{\odot})$')

ax_set_ylim_Redd = ax.set_ylim([ylowlim, yuplim])
ax_set_xlim_Lbol = ax.set_xlim([xlowlim, xuplim])


# use the previously defined function
# scatter_hist2(x10km, y10kr, ax, ax_histx, ax_histy,'aqua', 5, 'yes', 2, ax_set_2, limx=masslim, limy=reddlim, binwidth_x = op_binmbh_knuthbal, binwidth_y = op_binedd_knuthbal )
# scatter_hist2(x25km, y25kr, ax, ax_histx, ax_histy ,'cornflowerblue', 5, 'yes', 2, ax_set_2,limx=masslim, limy=reddlim, binwidth_x = op_binmbh_knuthbal, binwidth_y = op_binedd_knuthbal)

# scatter_hist2(xpm, ypr, ax, ax_histx, ax_histy,'cornflowerblue', 5, 'yes', 2, ax_set_2, limx=11.25, limy=48.3,binwidth_x = op_binmbh_knuthp,binwidth_y = op_binedd_knuthp )
# scatter_hist2(xem, yer, ax, ax_histx, ax_histy ,'purple', 60, 'yes', 10, ax_set_2,limx=11.25, limy=48.3,binwidth_x = op_binmbh_knuthp, binwidth_y = op_binedd_knuthp)

# ax.add_artist(ax.legend(title='test'))
# ax.legend(['Rankine+2020 BALs vmin<10,000 km/s','Rankine+2020 BALs 10,000 km/s<vmin<25,000 km/s','RH+(in prep) EHVO'],loc='lower left')
ax.legend(['Rankine+2020 Parent sample','RH+(in prep) EHVO'],loc='lower left')
ax.set_title("Title for second plot")
# plt.savefig(PLOT_DIREC + 'Redd_Lbol.png') 
# plt.show()
# plt.close()

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
yuplim = 0.99

ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

# ax_set_2 = ax.set(ylabel=r'log($M_{\mathrm{BH}}/M_{\odot})$', xlabel='log($L_{bol}$/erg s$^{-1}$)')
ax_set_2 = ax.set(ylabel='Eddington Ratio (log($L_{bol}$/$L_{Edd}$)', xlabel='log($L_{bol}$/erg s$^{-1}$)')

ax_set_ylim_Redd = ax.set_ylim([ylowlim, yuplim])
ax_set_xlim_Lbol = ax.set_xlim([xlowlim, xuplim])


# use the previously defined function
# scatter_hist2(x10kl, y10kr, ax, ax_histx, ax_histy,'aqua', 5, 'yes', 2, ax_set_2,limx=lbollim, limy=reddlim,binwidth_x = op_binlbol_knuthbal, binwidth_y = op_binedd_knuthbal)
# scatter_hist2(x25kl, y25kr, ax, ax_histx, ax_histy ,'cornflowerblue', 5, 'yes', 2, ax_set_2,limx=lbollim, limy=reddlim,binwidth_x = op_binlbol_knuthbal, binwidth_y = op_binedd_knuthbal)

# scatter_hist2(xpl, ypr, ax, ax_histx, ax_histy,'cornflowerblue', 5, 'yes', 2, ax_set_2, limx=48.3, limy=48.3, binwidth_x = op_binmbh_knuthp, binwidth_y = op_binedd_knuthp)
# scatter_hist2(xel, yer, ax, ax_histx, ax_histy ,'purple', 60, 'yes', 10, ax_set_2,limx=lbollim, limy=reddlim, binwidth_x = op_binlbol_knuthp, binwidth_y = op_binedd_knuthp)

# ax.add_artist(ax.legend(title='test'))
# ax.legend(['Rankine+2020 BALs vmin<10,000 km/s','Rankine+2020 BALs 10,000 km/s<vmin<25,000 km/s','RH+(in prep) EHVO'],loc='lower left')
ax.legend(['Rankine+2020 Parent sample','RH+(in prep) EHVO'],loc='lower left')
ax.set_title("Title for second plot")
# plt.savefig(PLOT_DIREC + 'Redd_Lbol.png') 
# plt.show()
# plt.close()


#%%
# ploting ehvo 25k vs 10k lbol
fig = plt.figure(1)
fig = plt.figure(figsize=(8, 8))

xlowlim = 1.6959
xuplim = 3.499

ylowlim = 46.
yuplim = 48.3

ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

# ax_set_2 = ax.set(ylabel=r'log($M_{\mathrm{BH}}/M_{\odot})$', xlabel='log($L_{bol}$/erg s$^{-1}$)')
ax_set_2 = ax.set(xlabel='z', ylabel='log($L_{bol}$/erg s$^{-1}$)')

ax_set_ylim_Redd = ax.set_ylim([ylowlim, yuplim])
ax_set_xlim_Lbol = ax.set_xlim([xlowlim, xuplim])


# use the previously defined function
# scatter_hist2(x10kl, y10kr, ax, ax_histx, ax_histy,'aqua', 5, 'yes', 2, ax_set_2,limx=lbollim, limy=reddlim,binwidth_x = op_binlbol_knuthbal, binwidth_y = op_binedd_knuthbal)
# scatter_hist2(x25kl, y25kr, ax, ax_histx, ax_histy ,'cornflowerblue', 5, 'yes', 2, ax_set_2,limx=lbollim, limy=reddlim,binwidth_x = op_binlbol_knuthbal, binwidth_y = op_binedd_knuthbal)

scatter_hist2(z_parentR, Lbol_parentR, ax, ax_histx, ax_histy,'cornflowerblue', 60, 'yes', 10, ax_set_2, limx=lbollim, limy=speedlim , binwidth_x = 0.21, binwidth_y = op_binspd40k_knuthp)
# scatter_hist2(ehvoin14rank_speedlbol48K, ehvoin14rank_speed48K/-1000, ax, ax_histx, ax_histy ,'purple', 60, 'yes', 10, ax_set_2,limx=lbollim, limy=speedlim, binwidth_x = op_bin48k_knuthp, binwidth_y = op_binspd48k_knuthp)
# scatter_hist2(ehvoin14rank_speedlbol48Kmore, ehvoin14rank_speed48Kmore/-1000, ax, ax_histx, ax_histy ,'magenta', 60, 'yes', 10, ax_set_2,limx=lbollim, limy=speedlim, binwidth_x = 0.21, binwidth_y = op_binspd48kmore_knuthp)

# ax.add_artist(ax.legend(title='test'))
# ax.legend(['Rankine+2020 BALs vmin<10,000 km/s','Rankine+2020 BALs 10,000 km/s<vmin<25,000 km/s','RH+(in prep) EHVO'],loc='lower left')
ax.legend(['EHVO < 43,000 km/s','EHVO > 43,000 km/s', 'EHVO >48k'],loc='lower left')
ax.set_title("Title for second plot")
# plt.savefig(PLOT_DIREC + 'Redd_Lbol.png') 
plt.show()
plt.close()

stop
#%%
# ploting ehvo 25k vs 10k eddington 
# fig = plt.figure(1)
# fig = plt.figure(figsize=(8, 8))

# xlowlim = -2.0
# xuplim = 0.99

# ylowlim = 30
# yuplim = 61.

# ax = fig.add_axes(rect_scatter)
# ax_histx = fig.add_axes(rect_histx, sharex=ax)
# ax_histy = fig.add_axes(rect_histy, sharey=ax)

# # ax_set_2 = ax.set(ylabel=r'log($M_{\mathrm{BH}}/M_{\odot})$', xlabel='log($L_{bol}$/erg s$^{-1}$)')
# ax_set_2 = ax.set(ylabel='Outflow Velocity (km/s)', xlabel='Eddington Ratio (log($L_{bol}$/$L_{Edd}$)')

# ax_set_ylim_Redd = ax.set_ylim([ylowlim, yuplim])
# ax_set_xlim_Lbol = ax.set_xlim([xlowlim, xuplim])


# # use the previously defined function
# # scatter_hist2(x10kl, y10kr, ax, ax_histx, ax_histy,'aqua', 5, 'yes', 2, ax_set_2,limx=lbollim, limy=reddlim,binwidth_x = op_binlbol_knuthbal, binwidth_y = op_binedd_knuthbal)
# # scatter_hist2(x25kl, y25kr, ax, ax_histx, ax_histy ,'cornflowerblue', 5, 'yes', 2, ax_set_2,limx=lbollim, limy=reddlim,binwidth_x = op_binlbol_knuthbal, binwidth_y = op_binedd_knuthbal)

# scatter_hist2(ehvoin14rank_edd43k, ehvoin14rank_speed40K/-1000, ax, ax_histx, ax_histy,'cornflowerblue', 60, 'yes', 10, ax_set_2, limx=lbollim, limy=speedlim , binwidth_x = 0.21, binwidth_y = op_binspd40k_knuthp)
# # scatter_hist2(ehvoin14rank_speedlbol48K, ehvoin14rank_speed48K/-1000, ax, ax_histx, ax_histy ,'purple', 60, 'yes', 10, ax_set_2,limx=lbollim, limy=speedlim, binwidth_x = op_bin48k_knuthp, binwidth_y = op_binspd48k_knuthp)
# scatter_hist2(ehvoin14rank_edd43kbeyond, ehvoin14rank_speed48Kmore/-1000, ax, ax_histx, ax_histy ,'magenta', 60, 'yes', 10, ax_set_2,limx=lbollim, limy=speedlim, binwidth_x = 0.21, binwidth_y = op_binspd48kmore_knuthp)

# # ax.add_artist(ax.legend(title='test'))
# # ax.legend(['Rankine+2020 BALs vmin<10,000 km/s','Rankine+2020 BALs 10,000 km/s<vmin<25,000 km/s','RH+(in prep) EHVO'],loc='lower left')
# ax.legend(['EHVO < 43,000 km/s','EHVO > 43,000 km/s', 'EHVO >48k'],loc='lower left')
# ax.set_title("Title for second plot")
# # plt.savefig(PLOT_DIREC + 'Redd_Lbol.png') 
# plt.show()
# plt.close()

#%%
# ploting ehvo 25k vs 10k mbh
# fig = plt.figure(1)
# fig = plt.figure(figsize=(8, 8))

# xlowlim = 8.20
# xuplim = 10.2

# ylowlim = 30 
# yuplim = 61.

# ax = fig.add_axes(rect_scatter)
# ax_histx = fig.add_axes(rect_histx, sharex=ax)
# ax_histy = fig.add_axes(rect_histy, sharey=ax)

# # ax_set_2 = ax.set(ylabel=r'log($M_{\mathrm{BH}}/M_{\odot})$', xlabel='log($L_{bol}$/erg s$^{-1}$)')
# ax_set_2 = ax.set(ylabel='Outflow Velocity (km/s)', xlabel='log($M_{\mathrm{BH}}/M_{\odot})$')

# ax_set_ylim_Redd = ax.set_ylim([ylowlim, yuplim])
# ax_set_xlim_Lbol = ax.set_xlim([xlowlim, xuplim])


# # use the previously defined function
# # scatter_hist2(x10kl, y10kr, ax, ax_histx, ax_histy,'aqua', 5, 'yes', 2, ax_set_2,limx=lbollim, limy=reddlim,binwidth_x = op_binlbol_knuthbal, binwidth_y = op_binedd_knuthbal)
# # scatter_hist2(x25kl, y25kr, ax, ax_histx, ax_histy ,'cornflowerblue', 5, 'yes', 2, ax_set_2,limx=lbollim, limy=reddlim,binwidth_x = op_binlbol_knuthbal, binwidth_y = op_binedd_knuthbal)

# scatter_hist2(ehvoin14rank_mbh43k, ehvoin14rank_speed40K/-1000, ax, ax_histx, ax_histy,'cornflowerblue', 60, 'yes', 10, ax_set_2, limx=lbollim, limy=speedlim , binwidth_x = 0.22, binwidth_y = op_binspd40k_knuthp)
# # scatter_hist2(ehvoin14rank_speedlbol48K, ehvoin14rank_speed48K/-1000, ax, ax_histx, ax_histy ,'purple', 60, 'yes', 10, ax_set_2,limx=lbollim, limy=speedlim, binwidth_x = op_bin48k_knuthp, binwidth_y = op_binspd48k_knuthp)
# scatter_hist2(ehvoin14rank_mbh43kbeyond, ehvoin14rank_speed48Kmore/-1000, ax, ax_histx, ax_histy ,'magenta', 60, 'yes', 10, ax_set_2,limx=lbollim, limy=speedlim, binwidth_x = 0.22, binwidth_y = op_binspd48kmore_knuthp)

# # ax.add_artist(ax.legend(title='test'))
# # ax.legend(['Rankine+2020 BALs vmin<10,000 km/s','Rankine+2020 BALs 10,000 km/s<vmin<25,000 km/s','RH+(in prep) EHVO'],loc='lower left')
# ax.legend(['EHVO < 43,000 km/s','EHVO > 43,000 km/s', 'EHVO >48k'],loc='lower left')
# ax.set_title("Title for second plot")
# # plt.savefig(PLOT_DIREC + 'Redd_Lbol.png') 
# plt.show()
# plt.close()