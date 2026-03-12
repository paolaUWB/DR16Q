# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 19:04:47 2023

@author: Paola Rodriguez Hidalgo, Easton Pierce, Morrigan Kalet, Maddi Gregory, Liliana Flores

This code includes the AutomatedOverplotFinal written by Easton from GitHub, url: https://github.com/paolaUWB/DR16Q/blob/Easton/VARIABILITY/AutomatedOverplotFinal.py

"""
import re
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
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
div_spectra_name = "57328" # MJD of the spectrum to be used as the numerator for division (the spectrum with a weaker absorption feature)


manual_ylim = 'yes' # This allows you to manually set the ylim on the graphs. Default behavior on 'no' is 1.1x max flux.

zem = 0 # redshift of the quasar, 0 in this case since it is in the rest frame.
coem=zem+1.

wavelength_emit1_initial = 990.  # left xlim in restframe (set to 1000 to see carbon IV)
wavelength_emit2_initial = 6000.  # right xlim in restframe (set to 1200 to see carbon IV) #CHANGE THIS LATER!!! XXX

vmin = [-80000] # make smaller to move right line right (bigger to move right line left) For absorption shading
vmax = [-100000] # make bigger to move left line left (smaller to move left line right)

colors_norm = ['#2495DF','#C7301E', '#df9424', '#7b03fc', '#1417d9'] #Colors the program will run through
colors_dered = ['#2495DF','#C7301E', '#df9424', '#7b03fc', '#1417d9']
colors_div = ['#2495DF','#2495DF', '#C7301E', '#7b03fc', '#1417d9']


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
legend_fontsize = 13
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

if use_norm:
    for i in range(len(norm_list)): # Splitting files in order to grab our MJD
        mjd_norm_list[i] = re.split('-', norm_list[i])[-2]
else:
    for i in range(len(dered_list)):
        mjd_dered_list[i] = re.split('-', dered_list[i])[-3]


max_list_norm = np.zeros_like(norm_list)
max_list_dered = np.zeros_like(dered_list)


if use_norm:
    list_to_use = norm_list
    n_list = n_list_norm
    max_list = max_list_norm
    mjd_list = mjd_norm_list
else:
    list_to_use = dered_list
    n_list = n_list_dered
    max_list = max_list_dered
    mjd_list = mjd_dered_list


# Set the index of the reference spectrum for division
div_index_weaker = np.where(mjd_list.astype(str) == div_spectra_name)[0][0]
print("div_index_weaker:", div_index_weaker)
print("Initial mjd_list: ", mjd_list)
print("Initial list_to_use: ", list_to_use)
print("Initial n_list: ", n_list)
print("Initial max_list: ", max_list)


# Make the reference spectrum the first element of the working lists, so that it gets plotted first and used as the numerator for division
if use_div and div_index_weaker != 0:
    # Swap using temporary variables for numpy arrays
    temp_file = list_to_use[0]
    list_to_use[0] = list_to_use[div_index_weaker]
    list_to_use[div_index_weaker] = temp_file
    
    temp_file = mjd_list[0]
    mjd_list[0] = mjd_list[div_index_weaker]
    mjd_list[div_index_weaker] = temp_file
    
    # Also keep the smoothing & maximum flux lists synced
    temp_file = n_list[0]
    n_list[0] = n_list[div_index_weaker]
    n_list[div_index_weaker] = temp_file
    
    temp_file = max_list[0]
    max_list[0] = max_list[div_index_weaker]
    max_list[div_index_weaker] = temp_file
    
print("After swapping indices: ")
print("Updated mjd_list: ", mjd_list)
print("Updated list_to_use: ", list_to_use)



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

plot_title = quasardirecname + ' | '

if use_smooth:
    plot_title += 'Smoothed | '
else:
    plot_title += 'Not Smoothed | '
if use_norm:
    plot_title += 'Normalized | '
else:
    plot_title += 'Unnormalized | '

plot_title += 'RP2'
ay1.set_title(plot_title, fontsize=20, pad=28)

if not use_div:
    ay1.set_xlabel(r"Observed Wavelength [$\rm \AA$]", fontsize=16) #Only set x label for top panel if not dividing, since X-axis will be shared with bottom panel if dividing
ay1.set_ylabel("Flux", fontsize=16)

ay1.set_xlim(wavelength_observe1,wavelength_observe2)
ay1.xaxis.set_label_coords(0.48, -0.08)
ay2.xaxis.set_label_coords(0.48, 1.11)
ay2.xaxis.set_major_locator(plt.MaxNLocator(5))
ay2.set_xlim(wavelength_observe1/coem,wavelength_observe2/coem)
ay2.set_xlabel(r"Restframe Wavelength [$\rm \AA$]", fontsize=16)

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
    
    if manual_ylim == 'no':
        max_list[i] = np.max(flux)

    error = data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2] #Get error values in our region
    
    if error_diagnostics == True:
        
        print(mjd_dered_list[i] + ' Norm error max: ' + str(np.max(error)))
        
    messed_up_error = np.where ( data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2]  > norm_error_threshold) #Get indices of points with error > 3
    
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
 
    if use_smooth:
        ay1.plot (wavelength, smooth(flux,n),'-', color = colors_dered[i], label = 'MJD ' + mjd_list[i], linewidth = 0.75)
        ay1.plot([wavelength_observe1,wavelength_observe2],[1,1],'r--')
        ay1.plot (wavelength, (smooth(error, n))/np.sqrt(n),'--') 
        ay1.legend(fontsize = legend_fontsize)

    else:
        ay1.plot (wavelength, flux,'-', color = colors_dered[i], label = 'MJD ' + mjd_list[i], linewidth = 0.75)
        ay1.plot([wavelength_observe1,wavelength_observe2],[1,1],'r--')
        ay1.plot (wavelength, (smooth(error, n))/np.sqrt(n),'--') 
        ay1.legend(fontsize = legend_fontsize)

###########################################################################################################################################################################################################################################

# Dividing spectra by interpolating spectra with weaker absorption (flux1) onto wavelength2 grid.

if use_div:
    for i, file in enumerate(list_to_use):
        if i == 0:
            # Store the reference (weaker absorption) spectrum data
            n_ref = n_list[i]
            data_ref = np.loadtxt(specdirec + file)
            
            wavelength_lower_limit_ref = np.where(data_ref[:,0] > wavelength_observe1)
            wavelength_upper_limit_ref = np.where(data_ref[:,0] < wavelength_observe2)
            
            wavelength_ref = data_ref[np.min(wavelength_lower_limit_ref[0]) : np.max(wavelength_upper_limit_ref[0]), 0]
            flux_ref = data_ref[np.min(wavelength_lower_limit_ref[0]) : np.max(wavelength_upper_limit_ref[0]), 1]
            error_ref = data_ref[np.min(wavelength_lower_limit_ref[0]) : np.max(wavelength_upper_limit_ref[0]), 2]
            continue  # Skip to next iteration; don't plot yet
        
        # For all other comparison spectra (i > 0), create individual plots
        fig, ay1 = plt.subplots()
        ay2 = plt.twiny(ay1)
        
        n = n_list[i]
        data = np.loadtxt(specdirec + file)
        
        wavelength_lower_limit = np.where(data[:,0] > wavelength_observe1)
        wavelength_upper_limit = np.where(data[:,0] < wavelength_observe2)
        
        wavelength = data[np.min(wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0]), 0]
        flux = data[np.min(wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0]), 1]
        error = data[np.min(wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0]), 2]
        
        if manual_ylim == 'no':
            max_list[i] = np.max(flux)

        # ------ Error diiagnostics
        if error_diagnostics == True:
            print(mjd_dered_list[i] + ' Norm error max: ' + str(np.max(error)))
            
        messed_up_error = np.where ( data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2]  > norm_error_threshold) #Get indices of points with error > 3
        plerror = error
        wavelength_emit = wavelength/(zem+1) #Unshift(?) the wavelength, back to a rest frame
        
        aa = np.where(error > 2)
        if len(aa) > 0:
            error[aa] = 0
        
        # SOMETIMES, THERE ARE PIXEL PROBLEMS, AND WE MIGHT GET AN ERROR OF 30 IN FLUX. TO AVOID THAT, WE HAVE DONE THIS. MESSED UP ERROR IS DEFINED ABOVE.
        if len (messed_up_error[0]) > 0: # original
            plerror[messed_up_error[0]]=0
            flux[messed_up_error[0]]=0
        # ------
        
        # Plot title
        plot_title = quasardirecname + ' | ' + mjd_list[0] + ' / ' + mjd_list[i]

        if use_norm:
            plot_title += ' | Normalized'
        else:
            plot_title += ' | Unnormalized'

        if use_smooth:
            plot_title += ' | Smoothed'
        else:
            plot_title += ' | Not Smoothed'
        ay1.set_title(plot_title, fontsize=22, pad=28)
        
        ay1.set_ylabel("Flux", fontsize=18, labelpad=17)
        ay1.set_xlim(wavelength_observe1, wavelength_observe2)
        ay1.xaxis.set_label_coords(0.48, -0.08)
        ay2.xaxis.set_label_coords(0.48, 1.11)
        ay2.xaxis.set_major_locator(plt.MaxNLocator(5))
        ay2.set_xlim(wavelength_observe1/coem, wavelength_observe2/coem)
        ay2.set_xlabel(r"Restframe Wavelength [$\rm \AA$]", fontsize=18)
        
        zem_plot = "z = " + str(zem)
        if manual_ylim == 'yes':
            ay2.text(zem_label_x_norm, zem_label_y_norm, zem_plot, bbox=dict(facecolor='none', edgecolor='black', pad=7.0))
            plt.ylim(0, top_ylim_norm)
        
        # Plot both spectra
        if i != 0:
            if use_smooth:
                ay1.plot(wavelength_ref, smooth(flux_ref, n_ref), '-', color=colors_div[0], label='MJD ' + mjd_list[0], linewidth=0.75)
                ay1.plot(wavelength, smooth(flux, n), '-', color=colors_div[i], label='MJD ' + mjd_list[i], linewidth=0.75)
                ay1.plot(wavelength_ref, smooth(error_ref, n_ref)/np.sqrt(n_ref), '--', color=colors_div[0], alpha=0.5)
                ay1.plot(wavelength, smooth(error, n)/np.sqrt(n), '--', color=colors_div[i], alpha=0.5)
            else:
                ay1.plot(wavelength_ref, flux_ref, '-', color=colors_div[0], label='MJD ' + mjd_list[0], linewidth=0.75)
                ay1.plot(wavelength, flux, '-', color=colors_div[i], label='MJD ' + mjd_list[i], linewidth=0.75)
                ay1.plot(wavelength_ref, error, '--', color=colors_div[0], alpha=0.5)
                ay1.plot(wavelength, error, '--', color=colors_div[i], alpha=0.5)
            
            ay1.legend(fontsize=legend_fontsize, loc='upper right')
            
            # Interpolate flux onto ref wavelength and divide
            interp_flux = np.interp(wavelength_ref, wavelength, flux)
            divided_spectra = interp_flux / flux_ref # ref (0) / current (i)
            if reverse_division:
                divided_spectra = flux_ref / interp_flux
            
            # Error propagation for division
            interp_error = np.interp(wavelength_ref, wavelength, error)
            divided_error = np.abs(divided_spectra) * np.sqrt((interp_error/interp_flux)**2 + (error_ref/flux_ref)**2)

            print(wavelength_ref)
            
            # ------ Plot ratio panel
            rel_ratio= divided_spectra - 1.0 # Relative (subtract 1 to center around 0)
            ratio_smooth = smooth(rel_ratio, n_ref)
            
            pos = ay1.get_position()
            res_height = pos.height * 0.30
            ay1.set_position([pos.x0, pos.y0 + res_height, pos.width, pos.height - res_height])
            ax_res = fig.add_axes([pos.x0, pos.y0, pos.width, res_height], sharex=ay1)
            
            if reverse_division:
                ax_res.plot(wavelength, -ratio_smooth, color='green', linewidth=1.0, label=(mjd_list[i] + ' / ' + mjd_list[0]))
            else:
                ax_res.plot(wavelength, ratio_smooth, color='green', linewidth=1.0, label=(mjd_list[0] + ' / ' + mjd_list[i]))
            
            ax_res.legend(fontsize=legend_fontsize, loc='upper right')
            ax_res.axhline(0.0, color='darkgray', linestyle='--', linewidth=0.8)
            # Symmetric y-limits for clearer visualization
            ylim_val = np.max(np.abs(ratio_smooth)) if np.max(np.abs(ratio_smooth)) > 0 else 1.0
            ax_res.set_ylim(-1.1*ylim_val, 1.1*ylim_val)
            #ax_res.set_ylim(-2.5, 2.5) # Fixed y-limits for better comparison across panels
            ax_res.set_ylabel ('Division Ratio', fontsize=18, labelpad=9)
            ax_res.tick_params(axis='both', which='major', labelsize=10)
            ax_res.set_xlabel(r"Observed Wavelength [$\rm \AA$]", fontsize=18, labelpad=10) # Shared x-label for the bottom panel
            ay1.xaxis.set_visible(False) # Hide x-axis labels and ticks on ay1 since they're shared with ax_res
            # ------

            # ------ Add absorption shading and emission line labels to the division plot
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
            # ------
            
            # Save division output
            div_output = np.column_stack((wavelength, divided_spectra, divided_error))
            headers = '\tWavelength\t Divided_Spectra\t Divided_Error'
            
            if reverse_division:
                if use_norm:
                    np.savetxt(specdirec + str(mjd_list[i]) + '_over_' + str(mjd_list[0]) + '_normalized.txt', div_output, header=headers)
                else:
                    np.savetxt(specdirec + str(mjd_list[i]) + '_over_' + str(mjd_list[0]) + '_unnormalized.txt', div_output, header=headers)
            else:
                if use_norm:
                    np.savetxt(specdirec + str(mjd_list[0]) + '_over_' + str(mjd_list[i]) + '_normalized.txt', div_output, header=headers)
                else:
                    np.savetxt(specdirec + str(mjd_list[0]) + '_over_' + str(mjd_list[i]) + '_unnormalized.txt', div_output, header=headers)
            
            # Save figure
            fig.subplots_adjust(top=0.81, bottom=0.36)
            outputdirec = '/Output_Plots/'
            div_fig_filename = specdirec + quasardirecname + mjd_list[0] + '_over_' + mjd_list[i] + pp2
            plt.savefig(div_fig_filename, dpi=300, bbox_inches='tight')
    plt.show()

    # Create CSV file for divided spectra
    divided_spectra_files = []
    for file in os.listdir(specdirec):
        # Match files with pattern like "57328_over_57329_normalized.txt" or similar
        if '_over_' in file and file.endswith('.txt'):
            divided_spectra_files.append(file)

    data = {
        'NORM SPECTRA FILE NAME': divided_spectra_files,
        'REDSHIFT': [0] * len(divided_spectra_files),
        'CALCULATED SNR': [0] * len(divided_spectra_files),
        'NEEDS RECALCULATION': ['N'] * len(divided_spectra_files),
        'Masked Regions': ['[]'] * len(divided_spectra_files)
    }
    df = pd.DataFrame(data)

    # Save to current directory
    csv_filename = specdirec + quasardirecname + '_divided.csv'
    print(csv_filename)
    print(specdirec)
    df.to_csv(csv_filename, index=False)
###########################################################################################################################################################################################################################################

# Add absorption shading and emission line labels to the original plot (not the division plots)
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

fig.subplots_adjust(left=0.1, right=2.5, top=0.9, bottom=0.6)
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
