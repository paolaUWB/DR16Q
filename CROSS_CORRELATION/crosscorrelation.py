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
from utility_functions import read_spectra, append_row_to_csv, clear_file, print_to_file, read_list_spectra
import os
from os.path import exists
import csv
import pandas as pd

#specnumDR16=XXX
#specnumDR14=XXX
specnumEHVO=98 
specnumPARENT=int(18165)

wnewfile = 'no'
#wnewfile = 'yes'

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

#%%

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

#%%

hdu_14 = fits.open(infoDR14)
data_14 = hdu_14[1].data
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

#%%
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
    
    OUTFILE = os.getcwd() + "/DR16PA_DR14INFO.csv"
    clear_file(OUTFILE)
    
    MBH_parent=[] ; errMBH_parent=[] ; Redd_parent=[] ; Lbol_parent= [] ; 
    #SDSS_Names_14, Plate_14, Mjd_14, Fiber_14, Mbh_14,Mbh_err_14, Lbol_14, Quality_Lbol_14, Redd_14, Quality_Redd_14 = ['SDSS_Names_14', 'Plate_14', 'Mjd_14', 'Fiber_14', 'Mbh_14','Mbh_err_14', 'Lbol_14', 'Quality_Lbol_14', 'Redd_14', 'Quality_Redd_14']
    Parentheaders = ['SDSS_Names_14', 'Plate_14', 'Mjd_14', 'Fiber_14', 'Mbh_14','Mbh_err_14', 'Lbol_14', 'Quality_Lbol_14', 'Redd_14', 'Quality_Redd_14']
    
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
            fields = [SDSS_name0[0], int(plate_14[vv,]), int(mjd_14[vv,]), int(fiber_14[vv,]), float(log_mbh_14[vv,]),float(log_mbh_err_14[vv,]), float(log_lbol_14[vv,]), int(q_lbol_14[vv,]), float(log_redd_14[vv,]), int(q_redd_14[vv,])]
            append_row_to_csv(OUTFILE, fields)
            #print(SDSS_name0)

        
#%%
infoParent2 = os.getcwd() + "/DR16PA_DR14INFO.csv"

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

#print(len(EHVO_plate))
MBH_parent = ParentinDR14_mbh
Lbol_parent = ParentinDR14_Lbol
Redd_parent = ParentinDR14_redd

OUTFILE2 = os.getcwd() + "/DR16EHVOs_inDR14INFO.csv"
clear_file(OUTFILE2)

MBH_EHVO =[] ; Lbol_EHVO = [] ; Redd_EHVO = []

#SDSS_Names_14, Plate_14, Mjd_14, Fiber_14, Mbh_14,Mbh_err_14, Lbol_14, Quality_Lbol_14, Redd_14, Quality_Redd_14 = ['SDSS_Names_14', 'Plate_14', 'Mjd_14', 'Fiber_14', 'Mbh_14','Mbh_err_14', 'Lbol_14', 'Quality_Lbol_14', 'Redd_14', 'Quality_Redd_14']
Parentheaders = ['SDSS_Names_14', 'Plate_14', 'Mjd_14', 'Fiber_14', 'Mbh_14','Mbh_err_14', 'Lbol_14', 'Quality_Lbol_14', 'Redd_14', 'Quality_Redd_14']

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
        fields2 = [SDSS_name1[0], int(ParentinDR14_Plate[ww,]), int(ParentinDR14_mjd[ww,]), int(ParentinDR14_fiber[ww,]), float(ParentinDR14_mbh[ww,]),float(ParentinDR14_mbherr[ww,]), float(ParentinDR14_Lbol[ww,]), int(ParentinDR14_QLbol[ww,]), float(ParentinDR14_redd[ww,]), int(ParentinDR14_Qredd[ww,])]
        append_row_to_csv(OUTFILE2, fields2)
        
#%%
#Histogram codes:
# MBH_parent=np.hstack(MBH_parent)
# MBH_parent_toplot=MBH_parent[np.hstack(np.where(MBH_parent != 0.))]

# MBH_EHVO=np.hstack(MBH_EHVO)
# MBH_EHVO_toplot=MBH_EHVO[np.hstack(np.where(MBH_EHVO != 0.))]

# fig=plt.figure(1)
 
# iqr = np.subtract(*np.percentile(MBH_parent_toplot, [75, 25]))
# nhist=(max(MBH_parent_toplot)-min(MBH_parent_toplot))/(2*iqr*(len(MBH_parent_toplot)**(-1/3)))
 
# bins=np.linspace(min(MBH_parent_toplot),max(MBH_parent_toplot),int(nhist))
 
# plt.hist([MBH_parent_toplot],bins,color=['black'],label=['parent'],histtype='step')
# # This for later when you have the 3 samples: plt.hist([MBH,MBH_BAL_toplot,MBH_EHVO_toplot],bins,color=['black','blue','red'],label=['parent','BAL','EHVO'],histtype='step')
 
# plt.xlabel(r'log($M_{\mathrm{BH}}/M_{\odot})$')
# plt.ylabel('$n/N_{tot}$')
# plt.ylabel('Number')
# plt.legend(loc='upper left')
# plt.show()

fig = plt.figure(1)
x = MBH_parent
x2 = MBH_EHVO

# y = Lbol_parent
# y2 = Lbol_EHVO

y = Redd_parent
y2 = Redd_EHVO

y_= Lbol_parent
y2_ = Lbol_EHVO

#MBH RANGE:(8.972, 10.2896) Redd RANGE:(-0.7445, -0.175) Lbol RANGE: (46.937, 47.6849)

def scatter_hist2(x, y, ax, ax_histx, ax_histy, color, area, mult, factor, ax_set):
    # no labels
    ax_histx.tick_params(axis='x', labelbottom=False)
    ax_histy.tick_params(axis='y', labelleft=False)
    # the scatter plot:
    #ax.set_xlim([7.5, 11.25]) #MBH x limits
    #ax.set_ylim([-2.0, 1.0]) #Redd y limits
    ax.set_xlim([-2.0, 1.0]) #Redd x limits
    ax.set_ylim([46, 48.3]) #Lbol y limits
    ax.scatter(x, y, s = area, color = color)
    ax.text(10.5,0.5,'')
    # now determine nice limits by hand:
    binwidth_x = 0.2
    binwidth_y = 0.1
    limx = (int(11.25/binwidth_x) + 1) * binwidth_x
    limy = (int(48.3/binwidth_y) + 1) * binwidth_y
    binsx = np.arange(-limx, limx + binwidth_x, binwidth_x)
    binsy = np.arange(-limy, limy + binwidth_y, binwidth_y)
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
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.005
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]
# start with a square Figure
fig = plt.figure(figsize=(8, 8))
ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

#ax_set_1 = ax.set(xlabel=r'log($M_{\mathrm{BH}}/M_{\odot})$', ylabel='Eddington Ratio (log($L_{bol}$/$L_{Edd}$)')
#ax_set_2 = ax.set(xlabel=r'log($M_{\mathrm{BH}}/M_{\odot})$', ylabel='log($L_{bol}$/erg s$^{-1}$)')
ax_set_3 = ax.set(xlabel='Eddington Ratio (log($L_{bol}$/$L_{Edd}$)', ylabel='log($L_{bol}$/erg s$^{-1}$)')

# ax_set_xlim_MBH = ax.set_xlim([7.5, 11.25])
# ax_set_ylim_Redd = ax.set_ylim([-2.0, 1.0])
# ax_set_ylim_Lbol = ax.set_ylim([46, 48.3])
# ax_set_xlim_Redd = ax.set_xlim([-2.0, 1.0])

# use the previously defined function
scatter_hist2(y, y_, ax, ax_histx, ax_histy,'mediumblue', 2, 'yes', 3, ax_set_3)
scatter_hist2(y2, y2_, ax, ax_histx, ax_histy ,'purple', 50, 'yes', 10, ax_set_3)
plt.savefig(PLOT_DIREC + 'Redd_Lbol.png') 
plt.show()
plt.close()

# scatter_hist2(x, y_, ax, ax_histx, ax_histy,'mediumblue', 2, 'yes', 3, ax_set_2, ax_set_xlim_MBH, ax_set_ylim_Lbol)
# scatter_hist2(x2, y2_, ax, ax_histx, ax_histy ,'purple', 50, 'yes', 10, ax_set_2, ax_set_xlim_MBH, ax_set_ylim_Lbol)
# plt.savefig(PLOT_DIREC + 'MBH_Lbol_final.pdf')
# plt.show()
# plt.close()

# scatter_hist2(y, y_, ax, ax_histx, ax_histy,'mediumblue', 2, 'yes', 3, ax_set_3, ax_set_xlim_Redd, ax_set_ylim_Lbol)
# scatter_hist2(y2, y2_, ax, ax_histx, ax_histy ,'purple', 50, 'yes', 10, ax_set_3, ax_set_xlim_Redd, ax_set_ylim_Lbol)
# plt.savefig(PLOT_DIREC + 'Redd_Lbol_final.pdf')
# plt.show()
# plt.close()

