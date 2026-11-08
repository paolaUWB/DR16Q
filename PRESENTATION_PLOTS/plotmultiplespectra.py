# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 11:33:16 2017

This code plots pretty spectra for posters/presentations/papers/etc.

@author: Paola
"""

import numpy as np
import matplotlib.pyplot as plt 
import os
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd


##------ Inputs/Outputs to change
specdirec = os.getcwd() + '/./' + 'Temp_NV_Data/Normalized/'

save_format = 'pdf' # 'pdf' to save as pdf file, 'png' to save as png file

#-- TO CHANGE EVERY TIME:
csv_file = 'spectra_list_norm.csv'  # CSV file containing norm_spectra, quasar_name, zem columns

# Read the CSV
df = pd.read_csv(os.path.join(specdirec, csv_file))

# Initialize PdfPages for saving all plots to one PDF
output_pdf = os.getcwd() + '/PRESENTATION_PLOTS/OUTPUT_FILES/multi_spectra_plots.pdf'
pp = PdfPages(output_pdf)

#topylim = 2.7 xxx
topylim = 13

topemlabel = topylim - 0.03 # where you want to place the ion labels

zem_label_x = 1450 # x-coordinate for zem label (in restframe)
zem_label_y = 1.9 # y-coordinate for zem label

wavelength_emit1_initial = 940.  # left xlim in restframe
wavelength_emit2_initial = 1600.  # right xlim in restframe

#-- absorption shading: 'yes' to include
NVabs = 'yes'
OVIabs = 'no'
SiIVabs = 'yes'
Lyaabs = 'no'

OVIem = 'yes' # OVI emission label: 'yes' to include

n = 3 # smooth box car

#------ location (wavelength) for doublets
c = 300000. # speed of light

### WHY ARE THESE DEFINED TWICE?? (as CIVll and avr_CIV_doublet)??
CIVll = 1549.0524 # avr_CIV_doublet[weighted avg]; individuals: 1550.7700, 1548.1950
SiIVll = 1396.747 # avr_SiIV_doublet[weighted avg]; individuals: 1402.770, 1393.755
NVll = 1240.15  # avr_NV_doublet[weighted avg]; individuals: 1242.80, 1238.82
OVIll = 1033.8160 # avr_OVI_doublet[weighted avg]; individuals: 1037.6167, 1031.9261

CIVllred = 1550.7700 # red line of doublet for right (vmin)
CIVllblue = 1548.1950 # blue line of doublet for left (vmax)
# avr_CIV_doublet = 1549.0524 #weighted average ### DELETE?

SiIVllred = 1402.770
SiIVllblue = 1393.755
# avr_SiIV_doublet = 1396.747 # weighted average; individuals: 1402.770, 1393.755 ### DELETE?

CII_emitted = 1335.313 # (weighted average); individuals:
OI_emitted = 1303.4951 # weighted average; individuals page 20 in Verner Table

NVllred = 1242.80
NVllblue = 1238.82
# avr_NV_doublet = 1240.15 # weighted average; individuals: 1242.80, 1238.82 ### DELETE?

Lya = 1215.6700

OVIllred = 1037.6167
OVIllblue = 1031.9261
# avr_OVI_doublet=1033.8160 # weighted average; individuals: 1037.6167, 1031.9261 ### DELETE?

Lyb = 1025.7222

##----- functions
def smooth(norm_flux, box_pts):   
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(norm_flux, box, mode='same')
    return y_smooth

def zabs(v, zem):
    beta = -v/c
    Rc = np.sqrt((beta+1.)/(1.-beta))
    za = ((1.+zem)/Rc)-1.
    return za

# Pre-pass: compute max topylim per quasar_name
max_topylim_per_quasar = {}
for index, row in df.iterrows():
    quasar_name = row['quasar_name'].split('+')[0]
    zem = row['zem']
    
    data = np.loadtxt(specdirec + row['norm_spectra'])
    wavelength_observe1 = (zem+1.)*wavelength_emit1_initial
    wavelength_observe2 = (zem+1.)*wavelength_emit2_initial
    
    wavelength_lower_limit = np.where(data[:,0] > wavelength_observe1)
    wavelength_upper_limit = np.where(data[:,0] < wavelength_observe2)
    
    flux = data[np.min(wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0]),1]
    error = data[np.min(wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0]),2]
    
    aa = np.where(error > 2)
    if len(aa) > 0:
        error[aa] = 0
    
    max_flux = np.nanmax(np.concatenate([flux, error]))
    topylim_temp = max_flux * 0.8

    csv_topylim = pd.to_numeric(row.get('topylim', np.nan), errors='coerce')
    if not np.isnan(csv_topylim):
        topylim_temp = csv_topylim
    
    if quasar_name not in max_topylim_per_quasar:
        max_topylim_per_quasar[quasar_name] = topylim_temp
    else:
        max_topylim_per_quasar[quasar_name] = max(max_topylim_per_quasar[quasar_name], topylim_temp)
#________________________________________________________________________________________________________________

# Loop over each spectrum in the csv
for index, row in df.iterrows():
    norm_spectra = row['norm_spectra']
    quasar_name = row['quasar_name'].split('+')[0]
    zem = row['zem']
    vmin = np.fromstring(row['vmin'].strip('[]'), sep=' ')
    vmax = np.fromstring(row['vmax'].strip('[]'), sep=' ')

    vmin_temp = vmin
    vmax_temp = vmax

    vmin = np.zeros(1)
    vmax = np.zeros(1)
    
    vmin[0] = np.max(np.abs(vmin_temp)*-1)
    vmax[0] = np.max(np.abs(vmax_temp)*-1)
    
    norm_spectra_mjd = norm_spectra.split('-')[-2]
    
    spectra_count = index + 1
    print('Processing spectrum ' + str(spectra_count) + ': ' + norm_spectra + 'zem=', zem)
    
    data='' 
    data = np.loadtxt(specdirec+norm_spectra) 

    wavelength_observe1 = (zem+1.)*wavelength_emit1_initial #Shift start wavelength into frame |<--This makes our wavelength range for
    wavelength_observe2 = (zem+1.)*wavelength_emit2_initial #Shift end wavelength into frame   |     the region we want to look at
   
    wavelength_lower_limit = np.where(data[:,0] > wavelength_observe1)
    wavelength_upper_limit = np.where(data[:,0] < wavelength_observe2)

    minwave= np.min(wavelength_lower_limit[0])
    maxwave= np.max(wavelength_upper_limit[0])

    wavelength = data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit [0]),0] #Get wavelengths in our data set that fall into our region of study
    actual_wavelength= wavelength

    flux = data[np.min(wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),1] #Get flux values in our region
    actual_flux = flux

    error= data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2] #Get error values in our region
    messed_up_error = np.where ( data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2]  >3) #Get inexes of points with error > 3

    plerror=error
   
    wavelength_emit = wavelength/(zem+1) #Unshift(?) the wavelength, back to a rest frame

    normflux=flux
    error_normflux=error

    aa= np.where(error_normflux > 2)
    if len(aa) > 0:
        error_normflux[aa]=0

    zabs_min=vmin
    zabs_max=vmax


    for i in range(0,len(vmin)):
        zabs_min[i]=zabs(vmin[i], zem)
        zabs_max[i]=zabs(vmax[i], zem)

    # Use pre-computed max topylim for this quasar to match all observations of same quasar
    topylim = max_topylim_per_quasar[quasar_name]
    topemlabel = topylim - 0.03 #ion labels

# ------ 
    #SOMETIMES, THERE ARE PIXEL PROBLEMS, AND WE MIGHT GET AN ERROR OF 30 IN FLUX. TO AVOID THAT, WE HAVE DONE THIS. MESSED UP ERROR IS 
    ####       DEFINED ABOVE.
    #if len (messed_up_error[0]) > 0:######################################################original
        #plerror[messed_up_error[0]]=0####################################################
	#flux[messed_up_error[0]]=0

    fig, ay1 = plt.subplots()

# ay1 = fig.add_subplot(1, 1, 1)

    plt.title(quasar_name + ' (MJD: ' + norm_spectra_mjd + ')' + ' | ' + 'Normalized', fontsize=14, pad=13)
    ay1.set_xlabel(r"Observed Wavelength [$\rm \AA$]")
    ay1.set_ylabel(r"Flux")
     
    ay1.plot (wavelength, smooth(normflux,n),'k-', linewidth=0.5)
    ay1.plot (wavelength, error_normflux,'k--') 

    plt.plot([wavelength_observe1,wavelength_observe2],[1,1],'r--')
    color = ['xkcd:shocking pink', 'black', 'xkcd:purpleish blue']
    color = ['xkcd:shocking pink', 'xkcd:azure', 'blue', 'xkcd:purpleish blue', 'xkcd:slate']
# color = ['red', 'green', 'blue', 'orange', 'purple']

    for k in range(0,len(vmin)):
        plt.axvspan(CIVll*(1.+zabs_max[k]),CIVll*(1.+zabs_min[k]), alpha=0.2, color=color[0])

        if NVabs == 'yes':
            plt.axvspan(NVllblue*(1.+zabs_max[k]),NVllred*(1.+zabs_min[k]), alpha=0.2, color=color[1])
    
        if OVIabs == 'yes':
            plt.axvspan(OVIllblue*(1.+zabs_max[k]),OVIllred*(1.+zabs_min[k]), alpha=0.2, color=color[2])

        if SiIVabs == 'yes':
            plt.axvspan(SiIVllblue*(1.+zabs_max[k]),SiIVllred*(1.+zabs_min[k]), alpha=0.2, color=color[3])

        if Lyaabs == 'yes':
            plt.axvspan(Lya*(1.+zabs_max[k]),Lya*(1.+zabs_min[k]), alpha=0.2, color=color[4])
# Plot only the last colored text label for each ion
    k = len(vmin) - 1
    plt.text(CIVll*(1.+zabs_min[k])-30., topylim * 0.2 - 0.1*k, 'CIV', color=color[0], fontname='serif', weight='bold')
    
    if NVabs == 'yes':
        plt.text(NVll*(1.+zabs_min[k]), topylim * 0.8 - 0.1*k, 'NV', color='xkcd:azure', fontname='serif', weight='bold')
    
    if OVIabs == 'yes':
        plt.text(OVIll*(1.+zabs_min[k]), topylim * 0.6 - 0.1*k, 'OVI', color=color[2], fontname='serif', weight='bold')
    
    if SiIVabs == 'yes':
        plt.text(SiIVll*(1.+zabs_min[k]), topylim * 0.3 - 0.1*k, 'SiIV', color=color[3], fontname='serif', weight='bold')
    
    if Lyaabs == 'yes':
        plt.text(Lya*(1.+zabs_min[k])-20., topylim * 0.6 - 0.1*k, 'Lya', color=color[4], fontname='serif', weight='bold')

#matplotlib.rcParams['font.sans-serif'] = ['Source Han Sans TW', 'sans-serif']

    plt.text(1549.0*(1+zem)-30,topemlabel ,'CIV',color='black',rotation = 90,fontname='serif', verticalalignment = 'top')
    plt.text(1402.770*(1+zem)-40.,topemlabel,'SiIV+OIV]',color='black',rotation = 90,fontname='serif', verticalalignment = 'top')

    alpha = 'Ly' + chr(945)
    plt.text(1242.804*(1+zem)+30.,topemlabel, alpha + '+NV' ,color='black',rotation = 90,fontname='serif', verticalalignment = 'top')
    plt.text(1304.8576*(1+zem)-35.,topemlabel,'OI',color='black',rotation=90,fontname='serif', verticalalignment = 'top')
    plt.text(1334.5323*(1+zem)-30.,topemlabel,'CII',color='black',rotation=90,fontname='serif', verticalalignment = 'top')

    if OVIem == 'yes':
        plt.text(OVIll*(1+zem)-30.,topemlabel,'OVI',color='black',rotation=90,fontname='serif', verticalalignment = 'top')

    coem=zem+1.
    plt.xlim(wavelength_observe1,wavelength_observe2) 

    ay2 = ay1.twiny()

    ay2.plot(wavelength/coem,1000*smooth(flux,n)+1000.)
    ay2.set_xlabel(r"Restframe Wavelength [$\rm \AA$]")

    ay1.xaxis.set_label_coords(0.48, -0.08)
    ay2.xaxis.set_label_coords(0.48, 1.11)
    ay2.xaxis.set_major_locator(plt.MaxNLocator(5))

    zem_plot = "z = " + str(zem)
   # plt.text(zem_label_x, zem_label_y, zem_plot, bbox=dict(facecolor='none', edgecolor='black', pad=7.0))
    
    ay1.text(
    0.94, 0.12, zem_plot,
    transform=ay1.transAxes,
    horizontalalignment='right',
    verticalalignment='bottom',
    bbox=dict(facecolor='none', edgecolor='black', pad=7.0)
    )
    
    plt.xlim(wavelength_observe1/coem,wavelength_observe2/coem) 
    plt.ylim(0,topylim)
    ay1.yaxis.set_major_locator(plt.MaxNLocator(8))  # Consistent y-axis tick count across all plots

    fig.subplots_adjust(top=0.88, bottom=0.14)
    pp.savefig(fig, bbox_inches='tight')

# Close the PDF after processing all spectra
pp.close()
