# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 11:33:16 2017

This code plots pretty spectra for posters/presentations/papers/etc.

@author: Paola
"""

import numpy as np
import matplotlib.pyplot as plt 
import os
from matplotlib.backends.backend_pdf import PdfPages ###Not used in code


##------ Inputs/Outputs to change
specdirec = "/home/jmaharaj/"

save_format = 'pdf' # 'pdf' to save as pdf file, 'png' to save as png file

#-- TO CHANGE EVERY TIME:
save_file_name = 'spec001' 
norm_spectra = 'spec-10431-58137-0135norm.dr16' ###Spectra file with all the data (3 columns: Observed wavelength, normalized flux, error)

zem = 1.91713968767695 # redshift of norm_spectra

topylim = 2.7 ###Max Norm flux on the graph
topemlabel = topylim - 0.03 # where you want to place the ion labels ###Like OI, CII, and CIV

zem_label_x = 1450 # x-coordinate for zem label (in restframe)
zem_label_y = 1.9 # y-coordinate for zem label

wavelength_emit1_initial = 1000.  # left xlim in restframe
wavelength_emit2_initial = 1600.  # right xlim in restframe

vmin = [-37190.96072] # make smaller to move right line right (bigger to move right line left)
vmax = [-40381.02541] # make bigger to move left line left (smaller to move left line right)

#-- absorption shading: 'yes' to include
NVabs = 'yes'
OVIabs = 'no'
SiIVabs = 'yes'
Lyaabs = 'no'

OVIem = 'no' # OVI emission label: 'yes' to include

n = 3 # smooth box car

#------ saving files
#-- save png
if save_format == 'png': 
    pp2 = save_file_name + '.png'

#-- save pdf
if save_format == 'pdf':
    pp2 = save_file_name + '.pdf'


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
def smooth(norm_flux, box_pts):   #<--- reduces noise ###norm flux = og flux val, box_pts = smooth box car -> n=3
    box = np.ones(box_pts)/box_pts ###np.ones makes array with ones, box_pts = 3 so np.ones(box_pts) is [1, 1, 1] then [1, 1, 1]/3 = [0.33, 0.33, 0.33]
    y_smooth = np.convolve(norm_flux, box, mode='same') ###np.convolve does moving average smoothing, np.convolve(data, smoothing filter, mode), mode='same' includes data in box while mode='full' does all data included outside of box bounds  
    return y_smooth

def zabs(v): ###calculates absorption redshift of absorp line corresponding to the outflow velocity v
    beta = -v/c ###v=-37190.96072 so -v=37190.96072 and c=300000 so beta=0.123969 which mean gas is moving 0.123969 times speed of light
    Rc = np.sqrt((beta+1.)/(1.-beta)) ###Rc is relativistic Doppler factor which is 1.132706
    za = ((1.+zem)/Rc)-1. ###za is absorption redshift of 1.57537, Quasar emission lines are at z=zem, Absorbing gas produces lines at z=a
    return za

spectra_count = 1 ###Number of spectra being processed, right now just 1 graph so = 1
print('spec_name' + ": " + norm_spectra)
print(spectra_count)
print('zem=',zem) ###All 3 will print spec_name: spec-10431-58137-0135norm.dr16, 1, and zem= 1.91713968767695 <- prints in terminal when code is ran
    
data='' ###sets data as an empty string but gets overwritten in next line
data = np.loadtxt(specdirec+norm_spectra) ###specdirec = /home/jmaharaj/, norm_spectra = "spec-10431-58137-0135norm.dr16", so full path becomes /home/jmaharaj/spec-10431-58137-0135norm.dr16
                                          ###np.loadtxt() reads txt files numbers and stores it as numpy array -> 1, 2, 3 becomes [1, 2, 3]
                                          ###data becomes array so data = [1, 2, 3]

wavelength_observe1 = (zem+1.)*wavelength_emit1_initial #Shift start wavelength into frame |<--This makes our wavelength range for
wavelength_observe2 = (zem+1.)*wavelength_emit2_initial #Shift end wavelength into frame   |     the region we want to look at
                                                        ###We only want to look at the spectra from wave_obs1 to wave_obs2 like 2900 to 4600
                                                        ###2900 to 4600 is 1000 to 1600 in restframe like seen on the graph

wavelength_lower_limit = np.where(data[:,0] > wavelength_observe1) ###in data theres column 0, 1, and 2 -> column 0 = wavelength vals
                                                                   ###data[:,0] pulls all wavelength vals
                                                                   ###wave_obs1=2917 so wave_low_lim is positions where wavelength is > 2917
wavelength_upper_limit = np.where(data[:,0] < wavelength_observe2) ###same thing but with less than
                                                                   ###so wave_up_lim is positions where wavelength is < 4667

minwave= np.min(wavelength_lower_limit[0]) ###finds the smallest index from wave_low_lim which is in row 0 -> wavelength is 3594
maxwave= np.max(wavelength_upper_limit[0]) ###finds the biggest index from wave_up_lim which is in row 1134 -> wavelength is 4666
                                          
wavelength = data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit [0]),0] #Get wavelengths in our data set that fall into our region of study
                              ###sets wavelength equal to the range of data we need which is from row 0 to 1134 or wavelength vals 3594 to 4666
actual_wavelength= wavelength ###renames variable

flux = data[np.min(wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),1] #Get flux values in our region
###does same as last two lines but for flux over the range we are looking at, uses 1 for the 2nd column since 2 is flux (0 is wave in 1st column)
actual_flux = flux ###renames

error= data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2] #Get error values in our region
###gives error in same format over range of what we are observing and uses 2 for column 3 which is error
messed_up_error = np.where ( data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2]  >3) #Get inexes of points with error > 3
###messed up error makes any val in error array ignored if it is greater than 3, this seems to not be used as mentioned later
plerror=error ###renames
   
wavelength_emit = wavelength/(zem+1) #Unshift(?) the wavelength, back to a rest frame ###observed to restframe wavelength

normflux=flux ###renamed
error_normflux=error ###renamed

aa= np.where(error_normflux > 2) ###Finds where in the array error is greater than 2
if len(aa) > 0: ###if there are any error values
    error_normflux[aa]=0 ###replace the error vals found with val of 0

zabs_min=vmin ###renamed, vmin=outlfow vel
zabs_max=vmax ###renamed

for i in range(0,len(vmin)): ###saying to repeat code for all vals of vmin (only 1 val)
    zabs_min[i]=zabs(vmin[i]) ###velocity limits converted to absorption redshifts, vmin to zabs_min, so absorption features are drawn on graph
    zabs_max[i]=zabs(vmax[i])
    

# ------ 
    #SOMETIMES, THERE ARE PIXEL PROBLEMS, AND WE MIGHT GET AN ERROR OF 30 IN FLUX. TO AVOID THAT, WE HAVE DONE THIS. MESSED UP ERROR IS 
    ####       DEFINED ABOVE.
    #if len (messed_up_error[0]) > 0:######################################################original
        #plerror[messed_up_error[0]]=0####################################################
	#flux[messed_up_error[0]]=0

fig, ay1 = plt.subplots() ###makes figure and axes 

# ay1 = fig.add_subplot(1, 1, 1)

# plt.title(spectrum)
ay1.set_xlabel(r"Observed Wavelength [$\rm \AA$]") ###x axis label
ay1.set_ylabel(r"Normalized Flux") ###y axis label
     
ay1.plot (wavelength, smooth(normflux,n),'k-') ###plots (x,y) which is (wave,flux) and in smooth(normflux,n), n=3 like before, style is k- (normal)
ay1.plot (wavelength, error_normflux,'k--') ###plots the error line, style is k-- which is dashed

plt.plot([wavelength_observe1,wavelength_observe2],[1,1],'r--') ###draws normalized flux line at y=1
color = ['xkcd:shocking pink', 'black', 'xkcd:purpleish blue'] ###list of colors used later on
color = ['xkcd:shocking pink', 'xkcd:azure', 'blue', 'xkcd:purpleish blue', 'xkcd:slate']
# color = ['red', 'green', 'blue', 'orange', 'purple']



for k in range(0,len(vmin)): ###loops for each absorption system which is each pair of vmin/vmax
                             ###converts velocity limits into absorption redshift vals
                             ###uses redshifts to get observed wavelength ranges
    plt.axvspan(CIVll*(1.+zabs_max[k]),CIVll*(1.+zabs_min[k]), alpha=0.2, color=color[0]) ###shades CIV region in pink
    plt.text(CIVll*(1.+zabs_min[k])-30.,0.5-0.1*k,'CIV',color=color[0],fontname='serif',weight='bold') ###text for CIV

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



#matplotlib.rcParams['font.sans-serif'] = ['Source Han Sans TW', 'sans-serif']

plt.text(1549.0*(1+zem)-30,topemlabel ,'CIV',color='black',rotation = 90,fontname='serif', verticalalignment = 'top') ###Names CIV at rest wavelength
plt.text(1402.770*(1+zem)-40.,topemlabel,'SiIV+OIV]',color='black',rotation = 90,fontname='serif', verticalalignment = 'top') ###Same for SilV, OIV

alpha = 'Ly' + chr(945) ###For Greek text
plt.text(1242.804*(1+zem)+30.,topemlabel, alpha + '+NV' ,color='black',rotation = 90,fontname='serif', verticalalignment = 'top') ###rest label for NV
plt.text(1304.8576*(1+zem)-35.,topemlabel,'OI',color='black',rotation=90,fontname='serif', verticalalignment = 'top') ###rest label for OI
plt.text(1334.5323*(1+zem)-30.,topemlabel,'CII',color='black',rotation=90,fontname='serif', verticalalignment = 'top') ###rest label for CII

if OVIem == 'yes': ###does same only if OVI is enabled
    plt.text(OVIll*(1+zem)-30.,topemlabel,'OVI',color='black',rotation=90,fontname='serif', verticalalignment = 'top') ###rest label for OVI

coem=zem+1. ###stores redshift + 1 for later
plt.xlim(wavelength_observe1,wavelength_observe2) ###sets the x axis limits from 2917 to 4667

ay2 = ay1.twiny() ###creates 2nd x axis for the restframe wavelength at the top

ay2.plot(wavelength/coem,1000*smooth(flux,n)+1000.) ###converts observed to restframe for the top x axis
ay2.set_xlabel(r"Restframe Wavelength [$\rm \AA$]") ###adds label for restframe wavelegth in top x axis

ay1.xaxis.set_label_coords(0.48, -0.08) ###location of bottom x axis label
ay2.xaxis.set_label_coords(0.48, 1.11) ###location of top x axis label
ay2.xaxis.set_major_locator(plt.MaxNLocator(5)) ###limits top x axis to 5 tick marks

zem_plot = "z = " + str(zem) ###converts z val into text for redshift label
plt.text(zem_label_x, zem_label_y, zem_plot, bbox=dict(facecolor='none', edgecolor='black', pad=7.0)) ###places box around redshoft label
plt.xlim(wavelength_observe1/coem,wavelength_observe2/coem) ###sets x axis limits in range of the restframe from 1000 to 1600
plt.ylim(0,topylim) ###sets y axis range of normalized flux from 0 to 2.7

fig.tight_layout() ###adjusts spacing and fit to ensure everything is properly visible

outdir = os.path.join(os.getcwd(), "OUTPUT_FILES") ###creates path for where output files will be saved
os.makedirs(outdir, exist_ok=True) ###creates output folder if not already exisitng

plt.savefig(os.path.join(outdir, pp2), dpi=100) ###saves spectrum plot