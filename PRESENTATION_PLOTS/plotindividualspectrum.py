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


##------ Inputs/Outputs to change
specdirec = os.getcwd() + '/Documents/GitHub/DR16Q/DR16Q_EHVO/' + 'NORM_DR16Q_EHVO/'

save_format = 'pdf' # 'pdf' to save as pdf file, 'png' to save as png file

#-- TO CHANGE EVERY TIME:
save_file_name = 'spec-5372-55978-0742' 
norm_spectra = 'spec-5372-55978-0742norm.dr16'

zem =2.382# redshift of norm_spectra

topylim = 2.5
topemlabel = topylim - 0.03 # where you want to place the ion labels

zem_label_x = 1400# x-coordinate for zem label (in restframe)
zem_label_y = 1.5# y-coordinate for zem label

wavelength_emit1_initial = 1060.  # left xlim in restframe
wavelength_emit2_initial = 1600.  # right xlim in restframe

vmin = [-36900,-31500] # make smaller to move right line right (bigger to move right line left)
vmax = [-41000,-35000] # make bigger to move left line left (smaller to move left line right)

#-- absorption shading: 'yes' to include
NVabs = 'yes'
OVIabs = 'no'
SiIVabs = 'yes'
Lyaabs = 'yes'

#Emission labels: 'yes' to include
CIVem = 'yes'
SiIVem = 'yes'
CIIem = 'yes'
OIem = 'yes'
LyNVem = 'yes'
OVIem = 'no' 

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
def smooth(norm_flux, box_pts):   
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(norm_flux, box, mode='same')
    return y_smooth

def zabs(v):
    beta = -v/c
    Rc = np.sqrt((beta+1.)/(1.-beta))
    za = ((1.+zem)/Rc)-1.
    return za


spectra_count = 1
print('spec_name' + ": " + norm_spectra)
print(spectra_count)
print('zem=',zem)
    
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
    zabs_min[i]=zabs(vmin[i])
    zabs_max[i]=zabs(vmax[i])
    

# ------ 
    #SOMETIMES, THERE ARE PIXEL PROBLEMS, AND WE MIGHT GET AN ERROR OF 30 IN FLUX. TO AVOID THAT, WE HAVE DONE THIS. MESSED UP ERROR IS 
    ####       DEFINED ABOVE.
    #if len (messed_up_error[0]) > 0:######################################################original
        #plerror[messed_up_error[0]]=0####################################################
	#flux[messed_up_error[0]]=0

fig, ay1 = plt.subplots()

# ay1 = fig.add_subplot(1, 1, 1)

# plt.title(spectrum)
ay1.set_xlabel(r"Observed Wavelength [$\rm \AA$]")
ay1.set_ylabel(r"Normalized Flux")
     
ay1.plot (wavelength, smooth(normflux,n),'k-')
ay1.plot (wavelength, error_normflux,'k--') 

plt.plot([wavelength_observe1,wavelength_observe2],[1,1],'r--')
color = ['xkcd:shocking pink', 'black', 'xkcd:purpleish blue']
color = ['xkcd:shocking pink', 'xkcd:azure', 'blue', 'xkcd:purpleish blue', 'xkcd:slate']
# color = ['red', 'green', 'blue', 'orange', 'purple']

for k in range(0,len(vmin)):
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



#matplotlib.rcParams['font.sans-serif'] = ['Source Han Sans TW', 'sans-serif']

if CIVem == 'yes':
    plt.text(1549.0*(1+zem)-30,topemlabel ,'CIV',color='black',rotation = 90,fontname='serif', verticalalignment = 'top')

if SiIVem == 'yes':
    plt.text(1402.770*(1+zem)-40.,topemlabel,'SiIV+OIV]',color='black',rotation = 90,fontname='serif', verticalalignment = 'top')

if CIIem == 'yes' :
    plt.text(1334.5323*(1+zem)-30.,topemlabel,'CII',color='black',rotation=90,fontname='serif', verticalalignment = 'top')

if OIem == 'yes':
    plt.text(1304.8576*(1+zem)-35.,topemlabel,'OI',color='black',rotation=90,fontname='serif', verticalalignment = 'top')

if LyNVem == 'yes':
    alpha = 'Ly' + chr(945)
    plt.text(1242.804*(1+zem)+30.,topemlabel, alpha + '+NV' ,color='black',rotation = 90,fontname='serif', verticalalignment = 'top')

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
plt.text(zem_label_x, zem_label_y, zem_plot, bbox=dict(facecolor='none', edgecolor='black', pad=7.0))
plt.xlim(wavelength_observe1/coem,wavelength_observe2/coem) 
plt.ylim(0,topylim)

fig.tight_layout() 

plt.savefig(os.getcwd() + '/Documents/GitHub/DR16Q/PRESENTATION_PLOTS/OUTPUT_FILES/' + pp2, dpi=200)
