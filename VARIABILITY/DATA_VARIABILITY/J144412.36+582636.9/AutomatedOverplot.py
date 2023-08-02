# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 19:04:47 2023

@author: Paola Rodriguez Hidalgo, Easton Pierce
"""

import re
import numpy as np
import matplotlib.pyplot as plt 
import os

##------ Inputs/Outputs to change

quasardirec = 'J144412.36+582636.9/' # Name of the object you are observing (same name as the directory folder it's contained in). If you don't have these directory names,
# Then set your specdirec to where your spec files are
quasardirecname = quasardirec[:-1]

specdirec = os.getcwd() + '/Documents/GitHub/DR16Q/VARIABILITY/DATA_VARIABILITY/' + quasardirec

save_format = 'png' # 'pdf' to save as pdf file, 'png' to save as png file

#-- TO CHANGE EVERY TIME:

save_file_name = quasardirecname

manual_ylim = 'no' # This allows you to manually set the ylim on the graphs. Default behavior on 'no' is 1.1x max flux.

zem =2.37# redshift of the quasar
coem=zem+1.

wavelength_emit1_initial = 1225.  # left xlim in restframe
wavelength_emit2_initial = 1600.  # right xlim in restframe

vmin = [-36900,-31500] # make smaller to move right line right (bigger to move right line left) For absorption shading
vmax = [-41000,-35000] # make bigger to move left line left (smaller to move left line right)

colors_norm = ['#2495DF','#C7301E', '#df9424', '#7b03fc', '#1417d9'] #Colors the program will run through
colors_dered = ['#2495DF','#C7301E', '#df9424', '#7b03fc', '#1417d9']

#-- absorption shading: 'yes' to include. CURRENTLY DEPRECIATED, DOESN'T WORK
CIV_abs = 'no'
NVabs = 'no'
OVIabs = 'no'
SiIVabs = 'no'
Lyaabs = 'no'

#Emission labels: 'yes' to include
CIVem = 'yes'
SiIVem = 'yes'
CIIem = 'yes'
OIem = 'yes'
LyNVem = 'yes'
OVIem = 'no' 

n_list_norm = [3,3,3,3,3] # smooth box car
n_list_dered = [3,3,3,3,3]
norm_error_threshold = 0.2
dered_error_threshold = 3
if manual_ylim == 'yes': # Manually change ylim here
    top_ylim_norm = 2.5
    top_ylim_dered = 37
    topemlabel = top_ylim_norm - 0.03 # where you want to place the ion labels
    topemlabel_dered = top_ylim_dered - 0.03
    zem_label_y_norm = top_ylim_norm/1.15 # y-coordinate for zem label
    zem_label_y_dered = top_ylim_dered/1.5
    
# Parameters to change the legend's location and size. Default is x = 0.58, y = 0.05 (norm and dered)
legend_fontsize = 8
norm_legend_xloc = 0.3
dered_legend_xloc = 0.3
norm_legend_yloc = 0.6
dered_legend_yloc = 0.6

###########################################################################################################################################################################################################################################

file_list = os.listdir(specdirec) # List of all the file names in our pointed directory

dered_list = []
norm_list = []

for i in range(len(file_list)): # Grab our dered and norm file names
    if file_list[i].endswith('dered.dr16') or file_list[i].endswith('dered.txt'):
       dered_list.append(file_list[i])
        
    if file_list[i].endswith('norm.dr16') or file_list[i].endswith('norm.txt'):
       norm_list.append(file_list[i])
    
mjd_norm_list = np.zeros_like(norm_list) 
mjd_dered_list = np.zeros_like(dered_list)    

for i in range(len(norm_list)): # Splitting files in order to grab our MJD
    mjd_norm_list[i] = re.split('-', norm_list[i])[-2]

for i in range(len(dered_list)):
    mjd_dered_list[i] = re.split('-', dered_list[i])[-3]


max_list_norm = np.zeros_like(norm_list)
max_list_dered = np.zeros_like(dered_list)


wavelength_observe1 = (zem+1.)*wavelength_emit1_initial #Shift start wavelength into frame |<--This makes our wavelength range for
wavelength_observe2 = (zem+1.)*wavelength_emit2_initial #Shift end wavelength into frame   |     the region we want to look at

zem_label_x_norm = (wavelength_observe1/coem + wavelength_observe2/coem)/2# x-coordinate for zem label (in restframe)
zem_label_x_dered = (wavelength_observe1 + wavelength_observe2)/2



#------ saving files
#-- save png
if save_format == 'png': 
    pp2 = '.png'

#-- save pdf
if save_format == 'pdf':
    pp2 = '.pdf'


#------ location (wavelength) for doublets
c = 300000. # speed of light


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

def zabs(v):
    beta = -v/c
    Rc = np.sqrt((beta+1.)/(1.-beta))
    za = ((1.+zem)/Rc)-1.
    return za

zabs_min=vmin
zabs_max=vmax

for i in range(0,len(vmin)):
    zabs_min[i]=zabs(vmin[i])
    zabs_max[i]=zabs(vmax[i])


fig, ay1 = plt.subplots()
ay2 = plt.twiny(ay1)

color = ['xkcd:shocking pink', 'black', 'xkcd:purpleish blue']
color = ['xkcd:shocking pink', 'xkcd:azure', 'blue', 'xkcd:purpleish blue', 'xkcd:slate']

for k in range(0,len(vmin)):
    if CIV_abs == 'yes':
        plt.axvspan(CIVll*(1.+zabs_max[k]),CIVll*(1.+zabs_min[k]), alpha=0.2, color=color[0])
        plt.text(CIVll*(1.+zabs_min[k])-30.,0.5-0.1*k,'CIV',color=color[0],fontname='serif',weight='bold')

    if NVabs == 'yes':
        plt.axvspan(NVllblue*(1.+zabs_max[k]),NVllred*(1.+zabs_min[k]), alpha=0.2, color=color[1])
        plt.text(NVll*(1.+zabs_min[k]),1.3-0.1*k,'NV',color='xkcd:azure',fontname='serif',weight='bold')
    
    if OVIabs == 'yes':
        plt.axvspan(OVIllblue*(1.+zabs_max[k]),OVIllred*(1.+zabs_min[k]), alpha=0.2, color=color[2])
        plt.text(OVIll*(1.+zabs_min[k]),1.4-0.1*k,'OVI',color=color[2],fontname='serif',weight='bold')

    if SiIVabs == 'yes':
        plt.axvspan(SiIVllblue*(1.+zabs_max[k]),SiIVllred*(1.+zabs_min[k]), alpha=0.2, color=color[3])
        plt.text(SiIVll*(1.+zabs_min[k]),0.25*k+2,'SiIV',color=color[3],fontname='serif',weight='bold')

    if Lyaabs == 'yes':
        plt.axvspan(Lya*(1.+zabs_max[k]),Lya*(1.+zabs_min[k]), alpha=0.2, color=color[4])
        plt.text(Lya*(1.+zabs_min[k])-20.,1.45-0.1*k,'Lya',color=color[4],fontname='serif',weight='bold')


ay1.set_title(quasardirecname)
ay1.set_xlabel(r"Observed Wavelength [$\rm \AA$]")
ay1.set_ylabel(r"Normalized Flux")
ay1.set_xlim(wavelength_observe1,wavelength_observe2)
ay1.xaxis.set_label_coords(0.48, -0.08)
ay2.xaxis.set_label_coords(0.48, 1.11)
ay2.xaxis.set_major_locator(plt.MaxNLocator(5))
ay2.set_xlim(wavelength_observe1/coem,wavelength_observe2/coem)
ay2.set_xlabel(r"Restframe Wavelength [$\rm \AA$]")

zem_plot = "z = " + str(zem)
if manual_ylim == 'yes':
    ay2.text(zem_label_x_norm, zem_label_y_norm, zem_plot, bbox=dict(facecolor='none', edgecolor='black', pad=7.0))
    plt.ylim(0,top_ylim_norm)



for i in range(len(norm_list)):
    n = n_list_norm[i]
       
    data = np.loadtxt(specdirec + norm_list[i])
    
       
    wavelength_lower_limit = np.where(data[:,0] > wavelength_observe1)
    wavelength_upper_limit = np.where(data[:,0] < wavelength_observe2)
    
    minwave= np.min(wavelength_lower_limit[0])
    maxwave= np.max(wavelength_upper_limit[0])
    
    wavelength = data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit [0]),0] #Get wavelengths in our data set that fall into our region of study
    actual_wavelength= wavelength
    
    flux = data[np.min(wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),1] #Get flux values in our region
    actual_flux = flux
    if manual_ylim == 'no':
        max_list_norm[i] = np.max(actual_flux)
        
    
    error= data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2] #Get error values in our region
    messed_up_error = np.where ( data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2]  > norm_error_threshold) #Get inexes of points with error > 3
    
    plerror=error
       
    wavelength_emit = wavelength/(zem+1) #Unshift(?) the wavelength, back to a rest frame
    
    normflux=flux
    error_normflux=error
    
    aa= np.where(error_normflux > 2)
    if len(aa) > 0:
        error_normflux[aa]=0
    
    # ------ 
    #SOMETIMES, THERE ARE PIXEL PROBLEMS, AND WE MIGHT GET AN ERROR OF 30 IN FLUX. TO AVOID THAT, WE HAVE DONE THIS. MESSED UP ERROR IS 
    ####       DEFINED ABOVE.
    if len (messed_up_error[0]) > 0:######################################################original
        plerror[messed_up_error[0]]=0###################################################
        flux[messed_up_error[0]]=0
    
    ay1.plot (wavelength, smooth(normflux,n),'-', color = colors_norm[i], label = 'MJD ' + mjd_norm_list[i], linewidth = 0.75)
    ay1.plot (wavelength, error_normflux,'--') 
    ay1.plot([wavelength_observe1,wavelength_observe2],[1,1],'r--')
    ay1.legend(bbox_to_anchor = (norm_legend_xloc,norm_legend_yloc),loc = 'lower center',fontsize = legend_fontsize)

if manual_ylim == 'no':
    top_ylim_norm = np.max(max_list_norm.astype(float))*1.15
    topemlabel = top_ylim_norm *0.99
    zem_label_y_norm = top_ylim_norm / 1.2
    ay2.text(zem_label_x_norm, zem_label_y_norm  , zem_plot, bbox=dict(facecolor='none', edgecolor='black', pad=5.0))
    plt.ylim(0,top_ylim_norm)

if CIVem == 'yes':
    ay1.text(1549.0*(1+zem)-30,topemlabel ,'CIV',color='black',rotation = 90,fontname='serif', verticalalignment = 'top')

if SiIVem == 'yes':
    ay1.text(1402.770*(1+zem)-40.,topemlabel,'SiIV+OIV]',color='black',rotation = 90,fontname='serif', verticalalignment = 'top')

if CIIem == 'yes' :
    ay1.text(1334.5323*(1+zem)-30.,topemlabel,'CII',color='black',rotation=90,fontname='serif', verticalalignment = 'top')

if OIem == 'yes':
    ay1.text(1304.8576*(1+zem)-35.,topemlabel,'OI',color='black',rotation=90,fontname='serif', verticalalignment = 'top')

if LyNVem == 'yes':
    alpha = 'Ly' + chr(945)
    ay1.text(1242.804*(1+zem)+30.,topemlabel, alpha + '+NV' ,color='black',rotation = 90,fontname='serif', verticalalignment = 'top')

if OVIem == 'yes':
    ay1.text(OVIll*(1+zem)-30.,topemlabel,'OVI',color='black',rotation=90,fontname='serif', verticalalignment = 'top')
    
fig.tight_layout() 
plt.savefig(os.getcwd() +  '/Documents/GitHub/DR16Q/VARIABILITY/DATA_VARIABILITY/' + quasardirec + quasardirecname + ' (Normalized)' + pp2, dpi=200)
plt.show()

########################################################################################################################################################################################################################################

fig2, ay3 = plt.subplots()
ay4 = plt.twiny(ay3)


ay3.set_title(quasardirecname)
ay3.set_xlabel(r"Observed Wavelength [$\rm \AA$]")
ay3.set_ylabel(r"Flux [$10^{-17}$ erg/s/cm$^{2}/{\rm \AA}$]")
ay3.set_xlim(wavelength_observe1,wavelength_observe2)
ay3.xaxis.set_label_coords(0.48, -0.08)
ay4.xaxis.set_label_coords(0.48, 1.11)
ay4.xaxis.set_major_locator(plt.MaxNLocator(5))
ay4.set_xlim(wavelength_observe1/coem,wavelength_observe2/coem)
ay4.set_xlabel(r"Restframe Wavelength [$\rm \AA$]")

zem_plot = "z = " + str(zem)

if manual_ylim == 'yes':
    ay3.text(zem_label_x_dered, zem_label_y_dered, zem_plot, bbox=dict(facecolor='none', edgecolor='black', pad=5.0))
    plt.ylim(0,top_ylim_dered)

for i in range(len(dered_list)):
    n = n_list_dered[i]
      
    data = np.loadtxt(specdirec + dered_list[i])
    
    
    wavelength_lower_limit = np.where(data[:,0] > wavelength_observe1)
    wavelength_upper_limit = np.where(data[:,0] < wavelength_observe2)
    
    minwave= np.min(wavelength_lower_limit[0])
    maxwave= np.max(wavelength_upper_limit[0])
    
    wavelength = data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit [0]),0] #Get wavelengths in our data set that fall into our region of study
    actual_wavelength= wavelength
    
    flux = data[np.min(wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),1] #Get flux values in our region
    actual_flux = flux
    
    if manual_ylim == 'no':
        max_list_dered[i] = np.max(actual_flux)
    
    error= data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2] #Get error values in our region
    messed_up_error = np.where ( data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2]  > dered_error_threshold) #Get inexes of points with error > 3
    
    plerror=error
      
    wavelength_emit = wavelength/(zem+1) #Unshift(?) the wavelength, back to a rest frame
    
    normflux=flux
    error_normflux=error
    
    aa= np.where(error_normflux > 2)
    if len(aa) > 0:
        error_normflux[aa]=0
    
    # ------ 
    #SOMETIMES, THERE ARE PIXEL PROBLEMS, AND WE MIGHT GET AN ERROR OF 30 IN FLUX. TO AVOID THAT, WE HAVE DONE THIS. MESSED UP ERROR IS 
    ####       DEFINED ABOVE.
    if len (messed_up_error[0]) > 0:######################################################original
        plerror[messed_up_error[0]]=0###################################################
        flux[messed_up_error[0]]=0
    
    ay3.plot(wavelength, smooth(normflux,n),color = colors_norm[i], label = 'MJD ' + mjd_dered_list[i], linewidth = 0.75)
    ay3.plot(wavelength, error_normflux, '--')
    ay3.legend(bbox_to_anchor = (dered_legend_xloc,dered_legend_yloc),loc = 'lower center', fontsize = legend_fontsize)

if manual_ylim == 'no':
    top_ylim_dered = np.max(max_list_dered.astype(float))*1.15
    topemlabel_dered = top_ylim_dered*0.99
    zem_label_y_dered = top_ylim_dered / 1.2
    ay3.text(zem_label_x_dered, zem_label_y_dered  , zem_plot, bbox=dict(facecolor='none', edgecolor='black', pad=5.0))
    plt.ylim(0,top_ylim_dered)

if CIVem == 'yes':
    ay3.text(1549.0*(1+zem)-30,topemlabel_dered ,'CIV',color='black',rotation = 90,fontname='serif', verticalalignment = 'top')

if SiIVem == 'yes':
    ay3.text(1402.770*(1+zem)-40.,topemlabel_dered,'SiIV+OIV]',color='black',rotation = 90,fontname='serif', verticalalignment = 'top')

if CIIem == 'yes' :
    ay3.text(1334.5323*(1+zem)-30.,topemlabel_dered,'CII',color='black',rotation=90,fontname='serif', verticalalignment = 'top')

if OIem == 'yes':
    ay3.text(1304.8576*(1+zem)-35.,topemlabel_dered,'OI',color='black',rotation=90,fontname='serif', verticalalignment = 'top')

if LyNVem == 'yes':
    alpha = 'Ly' + chr(945)
    ay3.text(1242.804*(1+zem)+30.,topemlabel_dered, alpha + '+NV' ,color='black',rotation = 90,fontname='serif', verticalalignment = 'top')

if OVIem == 'yes':
    ay3.text(OVIll*(1+zem)-30.,topemlabel_dered,'OVI',color='black',rotation=90,fontname='serif', verticalalignment = 'top')

fig2.tight_layout()

plt.savefig(os.getcwd() +  '/Documents/GitHub/DR16Q/VARIABILITY/DATA_VARIABILITY/' + quasardirec + quasardirecname +' (Dered)' + pp2, dpi=200)
