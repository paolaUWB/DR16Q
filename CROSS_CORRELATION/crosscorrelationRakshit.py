#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 16:05:09 2022
@author: Tzitzi
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
infoEHVO = os.getcwd() + "/../DR16Q_EHVO/DR16_EHVO_sorted_norm.csv" #the list of EHVOs (?? with more info or just the list)
infoPARENT = os.getcwd() + "/DR16_parent_sample.csv" #the list of DR16 parent sample
#Outputs
OUTFILE = os.getcwd() + "/DR16PA_DR14INFO.csv"
OUTFILE2 = os.getcwd() + "/DR16EHVOs_inDR14INFO.csv"

infoPARENT = os.getcwd() + "/../DR16_parent_sample.csv" #the list of DR16 parent sample

fieldnames = ['SDSS_Names', 'Mbh', 'L_bol', 'Redd']

PLOT_DIREC = os.getcwd() + "/OUTPUT_FILES/pngFILES/"

#for new csv: columns 16 mass 17 lbol  18 redd
#SDSS NAME 

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
#Writing csv file for DR16 parent with DR14 info 

if wnewfile == 'yes':
    print('I am here!')
    #DR16Parent sample read:
    df = pd.read_csv(infoPARENT, header=None)
    parent_first_col = df[df.columns[0]].to_numpy()
    spec_name_PARENT = parent_first_col.tolist()
    #print(spec_name_PARENT)
    
    #ehvo csv file read:
    df = pd.read_csv(infoEHVO)
    EHVO_first_col = df[df.columns[1]].to_numpy()
    spec_name_EHVO = EHVO_first_col.tolist()
    #print(spec_name_EHVO)
    
    OUTFILE = os.getcwd() + "/DR16parent_DR14INFO.csv"
    clear_file(OUTFILE)
    
    MBH_parent=[] ; errMBH_parent=[] ; Redd_parent=[] ; Lbol_parent= [] ; 
    #SDSS_Names_14, Plate_14, Mjd_14, Fiber_14, Mbh_14,Mbh_err_14, Lbol_14, Quality_Lbol_14, Redd_14, Quality_Redd_14 = ['SDSS_Names_14', 'Plate_14', 'Mjd_14', 'Fiber_14', 'Mbh_14','Mbh_err_14', 'Lbol_14', 'Quality_Lbol_14', 'Redd_14', 'Quality_Redd_14']
    Parentheaders = ['SDSS_Names_14', 'Plate_14', 'Mjd_14', 'Fiber_14', 'Mbh_14','Mbh_err_14', 'Lbol_14', 'Quality_Lbol_14', 'Redd_14', 'Quality_Redd_14', "BALQSO?"]
    
    with open(OUTFILE, 'w') as file:
        dw = csv.DictWriter(file, fieldnames = Parentheaders)
        dw.writeheader() 
    
    
        
    # reader = csv.reader(file)
    # writer = csv.writer(csvfile, fmtparams)


    #DR16parent and DR14 fits CC:
    for ii in range(specnumPARENT):
        aa=spec_name_PARENT[ii]
        #print(aa)
        bb=aa.split('-')
        #print(bb)
        cc=bb[3].split('.')
        #print(cc)
        plate_PARENT = int (bb[1])
       # print(plate_PARENT)
        mjd_PARENT = int (bb[2])
        fiber_PARENT = int (cc[0])
        vv, = np.where((plate_14[:] == plate_PARENT) & (mjd_14[:] == mjd_PARENT) & (fiber_14[:] == fiber_PARENT))
        #print(plate_14[:])
        if (len(vv) != 0):
            MBH_parent.append(log_mbh_14[vv,])
            Lbol_parent.append(log_lbol_14[vv,])
            Redd_parent.append(log_redd_14[vv,])
            #print(vv)
            #print(SDSS_name_14[vv,],plate_14[vv,],mjd_14[vv,],fiber_14[vv,],log_mbh_14[vv,],log_mbh_err_14[vv,],log_lbol_14[vv,], q_lbol_14[vv,], log_redd_14[vv,], q_redd_14[vv,])
            SDSS_name0 = SDSS_name_14[vv,]
            if bi_civ_14[vv] > 0:
                BALQSO = 'yes'
            fields = [SDSS_name0[0], int(plate_14[vv,]), int(mjd_14[vv,]), int(fiber_14[vv,]), float(log_mbh_14[vv,]),float(log_mbh_err_14[vv,]), float(log_lbol_14[vv,]), int(q_lbol_14[vv,]), float(log_redd_14[vv,]), int(q_redd_14[vv,]), BALQSO]
            append_row_to_csv(OUTFILE, fields)
            BALQSO ='no'
            #print(SDSS_name0)

        
#%%
infoParent2 = os.getcwd() + "/DR16parent_DR14INFO.csv"

#read EHVO csv file:
df = pd.read_csv(infoEHVO, header=None)
EHVO_plate = df[df.columns[3]].to_numpy()
EHVO_mdj = df[df.columns[4]].to_numpy()
EHVO_fiber = df[df.columns[5]].to_numpy()



#Read created ParentDR14info file:
df1 = pd.read_csv(infoParent2)
ParentinDR14_SDSSnames = df1[df1.columns[0]].to_numpy()
ParentinDR14_Plate = df1[df1.columns[1]].to_numpy()
ParentinDR14_mjd = df1[df1.columns[2]].to_numpy()
ParentinDR14_fiber = df1[df1.columns[3]].to_numpy()
ParentinDR14_mbh = df1[df1.columns[4]].to_numpy()
ParentinDR14_mbherr = df1[df1.columns[5]].to_numpy()
ParentinDR14_Lbol = df1[df1.columns[6]].to_numpy()
ParentinDR14_QLbol = df1[df1.columns[7]].to_numpy()
ParentinDR14_redd = df1[df1.columns[8]].to_numpy()
ParentinDR14_Qredd = df1[df1.columns[9]].to_numpy()
ParentinDR14_BALQSO = df1[df1.columns[10]].to_numpy()

#print(len(EHVO_plate))
MBH_parent = ParentinDR14_mbh
Lbol_parent = ParentinDR14_Lbol
Redd_parent = ParentinDR14_redd


OUTFILE2 = os.getcwd() + "/DR16EHVOs_inDR14INFO.csv"
clear_file(OUTFILE2)

MBH_EHVO =[] ; Lbol_EHVO = [] ; Redd_EHVO = []

#SDSS_Names_14, Plate_14, Mjd_14, Fiber_14, Mbh_14,Mbh_err_14, Lbol_14, Quality_Lbol_14, Redd_14, Quality_Redd_14 = ['SDSS_Names_14', 'Plate_14', 'Mjd_14', 'Fiber_14', 'Mbh_14','Mbh_err_14', 'Lbol_14', 'Quality_Lbol_14', 'Redd_14', 'Quality_Redd_14']
Parentheaders = ['SDSS_Names_14', 'Plate_14', 'Mjd_14', 'Fiber_14', 'Mbh_14','Mbh_err_14', 'Lbol_14', 'Quality_Lbol_14', 'Redd_14', 'Quality_Redd_14', 'BALQSO_14']


for jj in range(specnumEHVO):
    ww, = np.where((ParentinDR14_Plate == EHVO_plate[jj]) & (ParentinDR14_mjd == EHVO_mdj[jj]) & (ParentinDR14_fiber == EHVO_fiber[jj]))
    # print(EHVO_plate[jj])
    # print(ww)
    if (len(ww) != 0):
        Lbol_EHVO.append(ParentinDR14_Lbol[ww,])
        MBH_EHVO.append(ParentinDR14_mbh[ww,])
        Redd_EHVO.append(ParentinDR14_redd[ww,])
        #print(ParentinDR14_SDSSnames[ww]) 
        SDSS_name1 = ParentinDR14_SDSSnames[0]
        #using balQSO flags to get positions of Mass, Lbol, Redd
        fields2 = [SDSS_name1[0], int(ParentinDR14_Plate[ww,]), int(ParentinDR14_mjd[ww,]), int(ParentinDR14_fiber[ww,]), float(ParentinDR14_mbh[ww,]),float(ParentinDR14_mbherr[ww,]), float(ParentinDR14_Lbol[ww,]), int(ParentinDR14_QLbol[ww,]), float(ParentinDR14_redd[ww,]), int(ParentinDR14_Qredd[ww,]), ParentinDR14_BALQSO[ww,]]
        append_row_to_csv(OUTFILE2, fields2)
        
#%% BALQSO info
#infoParent2's last column has BALQSO flags, df1 = reading infoParent2
balinfo = df1[df1.columns[9]].to_numpy()
pos = np.where(balinfo < 0)

MBH_bal = ParentinDR14_mbh[pos]
Lbol_bal = ParentinDR14_Lbol[pos]
Redd_bal = ParentinDR14_redd[pos]
           
#%% bins width calc
# comb_mbhdata = np.concatenate([MBH_EHVOR , MBH_parentR])
# comb_edddata = np.concatenate([Redd_EHVOR , Redd_parentR])

op_binmbh = (stats.knuth_bin_width(MBH_EHVO)+stats.knuth_bin_width(MBH_bal))/2
op_binedd = (stats.knuth_bin_width(Redd_EHVO)+stats.knuth_bin_width(Redd_bal))/2
op_binlbol = (stats.knuth_bin_width(Lbol_EHVO)+stats.knuth_bin_width(Lbol_bal))/2

#%%
#Histogram:


def scatter_hist2(x, y, ax, ax_histx, ax_histy, color, area, mult, factor, ax_set, limx, limy):
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
    binwidth_x = op_binmbh
    binwidth_y = op_binlbol
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



# PLOTS

left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.005
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]


#%%
# Plot of  redd vs mass

fig = plt.figure(1)
fig = plt.figure(figsize=(8, 8))

x1m = MBH_parent
x2m = MBH_EHVO

y1 = Redd_parent
y2 = Redd_EHVO

x1 = Lbol_parent
x2 = Lbol_EHVO


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
scatter_hist2(x3m, y3, ax, ax_histx, ax_histy,'dodgerblue', 5, 'yes', 2, ax_set_2, limx=11.25, limy=48.3)
scatter_hist2(x4m, y4, ax, ax_histx, ax_histy ,'darkturquoise', 5, 'yes', 2, ax_set_2,limx=11.25, limy=48.3)
scatter_hist2(x2m, y2, ax, ax_histx, ax_histy ,'magenta', 60, 'yes', 10, ax_set_2,limx=11.25, limy=48.3)

# ax.add_artist(ax.legend(title='test'))
ax.legend(['Bal10k','Bal25k','EHVO'],loc='upper right')
ax.set_title("Title for second plot")
# plt.savefig(PLOT_DIREC + 'Redd_Lbol.png') 
plt.show()
plt.close()

#%%
# Plot of  redd vs Lbol

fig = plt.figure(1)
fig = plt.figure(figsize=(8, 8))

x1m = MBH_parent
x2m = MBH_EHVO

y1 = Redd_parent
y2 = Redd_EHVO

x1 = Lbol_parent
x2 = Lbol_EHVO

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
scatter_hist2(x3, y3, ax, ax_histx, ax_histy,'dodgerblue', 5, 'yes', 2, ax_set_2,limx=48.3, limy=48.3)
scatter_hist2(x4, y4, ax, ax_histx, ax_histy ,'darkturquoise', 5, 'yes', 2, ax_set_2,limx=48.3, limy=48.3)
scatter_hist2(x2, y2, ax, ax_histx, ax_histy ,'magenta', 60, 'yes', 10, ax_set_2,limx=48.3, limy=48.3)

# ax.add_artist(ax.legend(title='test'))
ax.legend(['Bal10k','Bal25k','EHVO'],loc='upper right')
ax.set_title("Title for second plot")
# plt.savefig(PLOT_DIREC + 'Redd_Lbol.png') 
plt.show()
plt.close()



