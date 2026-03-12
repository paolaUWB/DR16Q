# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 19:04:47 2023

@author: Paola Rodriguez Hidalgo, Easton Pierce, Morrigan Kalet, Maddi Gregory, Liliana Flores

This code includes the AutomatedOverplotFinal written by Easton from GitHub, url: https://github.com/paolaUWB/DR16Q/blob/Easton/VARIABILITY/AutomatedOverplotFinal.py

"""
import re
import numpy as np
import matplotlib.pyplot as plt 
import os

##------ Inputs/Outputs to change

quasardirec = '/J2318/' # Name of the object you are observing (same name as the directory folder it's contained in). If you don't have these directory names,
# Then set your specdirec to where your spec files are
quasardirecname = quasardirec[1:-1]

specdirec = os.getcwd() + quasardirec
# run this code one directory above where J2318 folder is located.

save_format = 'png' # 'pdf' to save as pdf file, 'png' to save as png file

#-- TO CHANGE EVERY TIME:

save_file_name = quasardirecname


use_norm = True # True to plot normalized spectra
use_div = True # True to plot divided spectra
reverse_division = False # True to divide stronger absorption by weaker absorption, False to divide weaker absorption by stronger absorption
use_smooth = True 


manual_ylim = 'yes' # This allows you to manually set the ylim on the graphs. Default behavior on 'no' is 1.1x max flux.

zem = 0 # redshift of the quasar, 0 in this case since it is in the rest frame.
coem=zem+1.

wavelength_emit1_initial = 1000.  # left xlim in restframe (set to 1000 to see carbon IV)
wavelength_emit2_initial = 1200.  # right xlim in restframe (set to 1200 to see carbon IV)

vmin = [-80000] # make smaller to move right line right (bigger to move right line left) For absorption shading
vmax = [-100000] # make bigger to move left line left (smaller to move left line right)

colors_norm = ['#2495DF','#C7301E', '#df9424', '#7b03fc', '#1417d9'] #Colors the program will run through
colors_dered = ['#2495DF','#C7301E', '#df9424', '#7b03fc', '#1417d9']
colors_div = ['#df9424','#2495DF']


#-- absorption shading: 'yes' to include. CURRENTLY DEPRECIATED, DOESN'T WORK
CIV_abs = 'yes'
NVabs = 'no'
OVIabs = 'no'
SiIVabs = 'yes'
Lyaabs = 'no'

# Emission labels: 'yes' to include
CIVem = 'no'
SiIVem = 'no'
CIIem = 'no'
OIem = 'no'
LyNVem = 'no'
OVIem = 'yes' 

n_list_norm = [3,3,3,3,3] # smooth box car
n_list_dered = [3,3,3,3,3]

error_diagnostics = False # Prints out the max error for each epoch, True or False

norm_error_threshold = 2025
dered_error_threshold = 4

if manual_ylim == 'yes': # Manually change ylim here
    if use_norm:
        top_ylim_norm = 1.3
        top_ylim_dered = 37
        topemlabel = top_ylim_norm - 0.03 # where you want to place the ion labels
        topemlabel_dered = top_ylim_dered - 0.03
        zem_label_y_norm = top_ylim_norm/1.15 # y-coordinate for zem label
        zem_label_y_dered = top_ylim_dered/1.5
    else:
        top_ylim_norm = 29.5
        top_ylim_dered = 37
        topemlabel = top_ylim_norm - 0.03 # where you want to place the ion labels
        topemlabel_dered = top_ylim_dered - 0.03
        zem_label_y_norm = top_ylim_norm/1.15 # y-coordinate for zem label
        zem_label_y_dered = top_ylim_dered/1.5
    
# Parameters to change the legend's location and size. Default is x = 0.58, y = 0.05 (norm and dered)
legend_fontsize = 12
norm_legend_xloc = 0.2
dered_legend_xloc = 0.2
norm_legend_yloc = 0.62
dered_legend_yloc = 0.65

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


if use_norm:
    list_to_use = norm_list
    mjd_list_to_use = mjd_norm_list
    n_list = n_list_norm
    max_list = max_list_norm
    mjd_list = mjd_norm_list
else:
    list_to_use = dered_list
    mjd_list_to_use = mjd_dered_list
    n_list = n_list_dered
    max_list = max_list_dered
    mjd_list = mjd_dered_list


# Added for dividing spectra: We want to divide the spectrum with weaker absorption (higher MJD) by the spectrum with stronger absorption (lower MJD)
for i in range(len(mjd_list_to_use)):
    div_index_weaker = np.argmax(mjd_list_to_use.astype(int))
    div_index_stronger = np.argmin(mjd_list_to_use.astype(int))
    print("index of spectrum with weaker absorption (higher MJD): " + str(div_index_weaker))
    print(mjd_list_to_use.astype(int))
# ------


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


##----- Functions
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

if use_smooth:
    ay1.set_title(quasardirecname + " | " + "Smoothed | RP2")
if use_smooth & use_norm:
    ay1.set_title(quasardirecname + " | " + "Smoothed | Normalized | RP2")
if not use_smooth & use_norm:
    ay1.set_title(quasardirecname + " | " + "Not Smoothed | Normalized | RP2")
if not use_smooth and not use_norm:
    ay1.set_title(quasardirecname + " | " + "Not Smoothed | Unnormalized | RP2")

if not use_div:
    ay1.set_xlabel(r"Observed Wavelength [$\rm \AA$]") #Only set x label for top panel if not dividing, since X-axis will be shared with bottom panel if dividing

ay1.set_ylabel("Flux", fontsize=12, labelpad=8)

ay1.set_xlim(wavelength_observe1,wavelength_observe2)
ay1.xaxis.set_label_coords(0.48, -0.08)
ay2.xaxis.set_label_coords(0.48, 1.11)
ay2.xaxis.set_major_locator(plt.MaxNLocator(5))
ay2.set_xlim(wavelength_observe1/coem,wavelength_observe2/coem)
ay2.set_xlabel(r"Restframe Wavelength [$\rm \AA$]", fontsize=12)

zem_plot = "z = " + str(zem)
if manual_ylim == 'yes':
    ay2.text(zem_label_x_norm, zem_label_y_norm, zem_plot, bbox=dict(facecolor='none', edgecolor='black', pad=7.0))
    plt.ylim(0,top_ylim_norm)


for i, file in enumerate(list_to_use):
    n = n_list[i]

    data = np.loadtxt(specdirec + file)
       
    wavelength_lower_limit = np.where(data[:,0] > wavelength_observe1)
    wavelength_upper_limit = np.where(data[:,0] < wavelength_observe2)
    
    minwave= np.min(wavelength_lower_limit[0])
    maxwave= np.max(wavelength_upper_limit[0])
    
    wavelength = data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit [0]),0] #Get wavelengths in our data set that fall into our region of study
    
    flux = data[np.min(wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),1] #Get flux values in our region
    
    # ------ Added for dividing spectra
    
    if i == div_index_stronger :
        flux1 = flux
        wavelength1 = wavelength 
        error1 = data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2] 
      
    elif i == div_index_weaker :
        flux2 = flux
        wavelength2 = wavelength
        error2 =  data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2]


    # Print files to check their order in directory
    print('File ' + str(i) + ': ' + list_to_use[i])
        
    # ------ 
    
    if manual_ylim == 'no':
        max_list[i] = np.max(flux)
        
    
    error = data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2] #Get error values in our region
    
    if error_diagnostics == True:
        
        print(mjd_dered_list[i] + ' Norm error max: ' + str(np.max(error)))
        
    messed_up_error = np.where ( data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2]  > norm_error_threshold) #Get inexes of points with error > 3
    
    plerror = error
       
    wavelength_emit = wavelength/(zem+1) #Unshift(?) the wavelength, back to a rest frame
    
    aa = np.where(error > 2)
    if len(aa) > 0:
        error[aa] = 0
    
    # ------ 
    # SOMETIMES, THERE ARE PIXEL PROBLEMS, AND WE MIGHT GET AN ERROR OF 30 IN FLUX. TO AVOID THAT, WE HAVE DONE THIS. MESSED UP ERROR IS 
    ####       DEFINED ABOVE.
    if len (messed_up_error[0]) > 0:######################################################original
        plerror[messed_up_error[0]]=0###################################################
        flux[messed_up_error[0]]=0
    print(flux.shape)
    print(error.shape)
    print(wavelength.shape)
    print(np.min(flux), np.max(flux))
    print(np.min(error), np.max(error))
    print(np.min(wavelength), np.max(wavelength))

    if use_smooth:
        if use_div:
            ay1.plot (wavelength, smooth(flux, n),'-', color = colors_div[i], label = 'MJD ' + mjd_list[i], linewidth = 0.75)

        if not use_div:
            ay1.plot (wavelength, smooth(flux,n),'-', color = colors_dered[i], label = 'MJD ' + mjd_list[i], linewidth = 0.75)
            ay1.plot([wavelength_observe1,wavelength_observe2],[1,1],'r--')

        ay1.plot (wavelength, (smooth(error, n))/np.sqrt(n),'--') 
        ay1.legend(fontsize = legend_fontsize)

    else:
        if use_div:
            ay1.plot (wavelength, flux,'-', color = colors_div[i], label = 'MJD ' + mjd_list[i], linewidth = 0.75)

        if not use_div:
            ay1.plot (wavelength, flux,'-', color = colors_dered[i], label = 'MJD ' + mjd_list[i], linewidth = 0.75)
            ay1.plot([wavelength_observe1,wavelength_observe2],[1,1],'r--')

        ay1.plot (wavelength, (smooth(error, n))/np.sqrt(n),'--') 
        ay1.legend(fontsize = legend_fontsize)


###########################################################################################################################################################################################################################################

# Dividing spectra by interpolating spectra with weaker absorption (flux1) onto wavelength2 grid.

if use_div:
    common_wavelength = wavelength2
    interp_flux1 = np.interp(wavelength2, wavelength1, flux1)
    divided_spectra =  flux2 / interp_flux1
    if reverse_division:
        divided_spectra =  interp_flux1 / flux2

    # For error on division, we can use error propagation for division
    interp_error1 = np.interp(wavelength2, wavelength1, error1)

    divided_error = np.abs(divided_spectra) * np.sqrt((interp_error1/interp_flux1)**2 + (error2/flux2)**2)

    # Add a residual panel below the main plot. Relative residual = ratio - 1
    rel_res = divided_spectra - 1.0 
    res_smooth = smooth(rel_res, 13)

    # Add small residual axis under ay1 and create it
    pos = ay1.get_position()  # Bbox in figure coordinates
    res_height = pos.height*0.30  # 25% of original axes height for residual
    ay1.set_position([pos.x0, pos.y0 + res_height, pos.width, pos.height - res_height])
    ax_res = fig.add_axes([pos.x0, pos.y0, pos.width, res_height], sharex=ay1)

    if reverse_division:
        ax_res.plot(common_wavelength, -res_smooth, color='green', linewidth=1.0, label=(mjd_list[div_index_stronger] + ' / ' + mjd_list[div_index_weaker]))
    else:
        ax_res.plot(common_wavelength, res_smooth, color='green', linewidth=1.0, label=(mjd_list[div_index_weaker] + ' / ' + mjd_list[div_index_stronger]))

    ax_res.legend(fontsize=legend_fontsize)
    ax_res.axhline(0.0, color='darkgray', linestyle='--', linewidth=0.8)
    # Symmetric y-limits for clearer visualization
    ylim_val = np.max(np.abs(res_smooth)) if np.max(np.abs(res_smooth)) > 0 else 1e-3
    ax_res.set_ylim(-1.1*ylim_val, 1.1*ylim_val)
    ax_res.set_ylabel('Division Ratio', fontsize=12)
    ax_res.tick_params(axis='both', which='major', labelsize=10)
    ax_res.set_xlabel(r"Observed Wavelength [$\rm \AA$]", fontsize = 12,labelpad = 8)  # Shared x-label for the bottom panel
    ay1.xaxis.set_visible(False)     # Hide x-axis labels and ticks on ay1 since they're shared with ax_res

#Output divided spectra to run in absorption.py
division_output = np.column_stack((common_wavelength, divided_spectra, divided_error))
headers = '\tWavelength\t Divided_Spectra\t Divided_Error'

if reverse_division:
    if use_norm:
        np.savetxt(specdirec + str(mjd_list[div_index_stronger]) + '_over_' + str(mjd_list[div_index_weaker]) + '_norm.txt', division_output, header=headers)
        np.save
    else:
        np.savetxt(specdirec + str(mjd_list[div_index_stronger]) + '_over_' + str(mjd_list[div_index_weaker]) + '_unnorm.txt', division_output, header=headers)
else:
    if use_norm:
        np.savetxt(specdirec + str(mjd_list[div_index_weaker]) + '_over_' + str(mjd_list[div_index_stronger]) + '_norm.txt', division_output, header=headers)
    else:
        np.savetxt(specdirec + str(mjd_list[div_index_weaker]) + '_over_' + str(mjd_list[div_index_stronger]) + '_unnorm.txt', division_output, header=headers)


        


###########################################################################################################################################################################################################################################


if manual_ylim == 'no':
    top_ylim_norm = np.max(max_list.astype(float))*1.15
    topemlabel = top_ylim_norm *0.99
    zem_label_y_norm = top_ylim_norm / 1.2
    ay2.text(zem_label_x_norm, zem_label_y_norm  , zem_plot, bbox=dict(facecolor='none', edgecolor='black', pad=5.0))
    plt.ylim(0,top_ylim_norm)
    
for k in range(0,len(vmin)):
    if CIV_abs == 'yes':
        ay1.axvspan(CIVll*(1.+zabs_max[k]),CIVll*(1.+zabs_min[k]), alpha=0.2, color=color[0])
        ay1.text(CIVll*(1.+zabs_min[k])-30.,0.5-0.1*k,'CIV',color=color[0],fontname='serif',weight='bold')

    if NVabs == 'yes':
        ay1.axvspan(NVllblue*(1.+zabs_max[k]),NVllred*(1.+zabs_min[k]), alpha=0.2, color=color[1])
        ay1.text(NVll*(1.+zabs_min[k]),1.3-0.1*k,'NV',color='xkcd:azure',fontname='serif',weight='bold')
    
    if OVIabs == 'yes':
        ay1.axvspan(OVIllblue*(1.+zabs_max[k]),OVIllred*(1.+zabs_min[k]), alpha=0.2, color=color[2])
        ay1.text(OVIll*(1.+zabs_min[k]),1.4-0.1*k,'OVI',color=color[2],fontname='serif',weight='bold')

    if SiIVabs == 'yes':
        ay1.axvspan(SiIVllblue*(1.+zabs_max[k]),SiIVllred*(1.+zabs_min[k]), alpha=0.2, color=color[3])
        # ay1.text(SiIVll*(1.+zabs_min[k]),0.25*k+2,'SiIV',color=color[3],fontname='serif',weight='bold')
# Delete 279 later
    if Lyaabs == 'yes':
        ay1.axvspan(Lya*(1.+zabs_max[k]),Lya*(1.+zabs_min[k]), alpha=0.2, color=color[4])
        ay1.text(Lya*(1.+zabs_min[k])-20.,1.45-0.1*k,'Lya',color=color[4],fontname='serif',weight='bold')

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

if use_div:   
    fig.subplots_adjust(top=0.86, bottom=0.36)
else:
    fig.tight_layout()  

plt.savefig(specdirec + quasardirecname + 'norm' + pp2, dpi=200)
plt.show()

########################################################################################################################################################################################################################################

# ------ Plot dereddened spectra

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
ay4.set_xlabel(r"Restframe Wavelength [$\rm \AA$]", )

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
    
    flux = data[np.min(wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),1] #Get flux values in our region
    actual_flux = flux
    
    if manual_ylim == 'no':
        max_list_dered[i] = np.max(actual_flux)
    
    error= data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2] #Get error values in our region
    if error_diagnostics == True:
        
        print( mjd_dered_list[i] + ' Dered error max: ' + str(np.max(error)))
    messed_up_error = np.where ( data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2]  > dered_error_threshold) #Get indexes of points with error > 3
    
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
    ay3.plot(wavelength, (smooth(error_normflux, n))/np.sqrt(n), '--')
    ay3.legend(bbox_to_anchor = (dered_legend_xloc,dered_legend_yloc),loc = 'lower center', fontsize = legend_fontsize)

if manual_ylim == 'no':
    top_ylim_dered = np.max(max_list_dered.astype(float))*1.15
    topemlabel_dered = top_ylim_dered*0.99
    zem_label_y_dered = top_ylim_dered / 1.2
    ay3.text(zem_label_x_dered, zem_label_y_dered  , zem_plot, bbox=dict(facecolor='none', edgecolor='black', pad=5.0))
    plt.ylim(0,top_ylim_dered)
    
for k in range(0,len(vmin)):
    if CIV_abs == 'yes':
        ay3.axvspan(CIVll*(1.+zabs_max[k]),CIVll*(1.+zabs_min[k]), alpha=0.2, color=color[0])
        ay3.text(CIVll*(1.+zabs_min[k])-63.,top_ylim_dered/8-0.1*k,'CIV',color=color[0],fontname='serif',weight='bold')

    if NVabs == 'yes':
        ay3.axvspan(NVllblue*(1.+zabs_max[k]),NVllred*(1.+zabs_min[k]), alpha=0.2, color=color[1])
        ay3.text(NVll*(1.+zabs_min[k]),1.3-0.1*k,'NV',color='xkcd:azure',fontname='serif',weight='bold')
    
    if OVIabs == 'yes':
        ay3.axvspan(OVIllblue*(1.+zabs_max[k]),OVIllred*(1.+zabs_min[k]), alpha=0.2, color=color[2])
        ay3.text(OVIll*(1.+zabs_min[k]),1.4-0.1*k,'OVI',color=color[2],fontname='serif',weight='bold')

    if SiIVabs == 'yes':
        ay3.axvspan(SiIVllblue*(1.+zabs_max[k]),SiIVllred*(1.+zabs_min[k]), alpha=0.2, color=color[3])
        ay3.text(SiIVll*(1.+zabs_min[k]),0.25*k+2,'SiIV',color=color[3],fontname='serif',weight='bold')

    if Lyaabs == 'yes':
        ay3.axvspan(Lya*(1.+zabs_max[k]),Lya*(1.+zabs_min[k]), alpha=0.2, color=color[4])
        ay3.text(Lya*(1.+zabs_min[k])-20.,1.45-0.1*k,'Lya',color=color[4],fontname='serif',weight='bold')

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
plt.savefig(os.getcwd() + quasardirec + quasardirecname +' (Dered)Expanded') # '/Documents/GitHub/DR16Q/VARIABILITY/DATA_VARIABILITY/' + quasardirec + quasardirecname +' (Dered)Expanded' + pp2, dpi=200)

