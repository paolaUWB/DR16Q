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
specdirec = os.getcwd() + '/../' #+ 'EHVO_NORM_DR16Q/'

save_format = 'pdf' # 'pdf' to save as pdf file, 'png' to save as png file

#-- TO CHANGE EVERY TIME:
save_file_name = 'spec002' 
norm_spectra = 'COPY_bh9_b364c3_l-2_al075_mw70_wa200_sorted72.dat'

zem = 0 #4.479 # redshift of norm_spectra

topylim = 3.0
topemlabel = topylim - 0.03 # where you want to place the ion labels

plot_z = 'yes'
zem_label_x = 1450 # x-coordinate for zem label (in restframe)
zem_label_y = 1.5 # y-coordinate for zem label


wavelength_emit1_initial = 400.  # left xlim in restframe
wavelength_emit2_initial = 1600.  # right xlim in restframe

vmin = [-56000] # make smaller to move right line right (bigger to move right line left)
vmax = [-52200] # make bigger to move left line left (smaller to move left line right)

#-- absorption shading: 'yes' to include

#Averages

CIIIabs = 'no'
CIVabs = 'yes'
SiIVabs = 'yes'
NVabs = 'yes'
OVIabs = 'no'
Lyaabs = 'yes'
Lybabs = 'no'
AlIII = 'no'
MgXabs = 'no'

#Doublets 

CIVllredabs = 'no'
CIVllblueabs = 'no'

SiIVllredabs = 'no'
SiIVllblueabs = 'no'

NVllredabs = 'no'
NVllblueabs = 'no'

OIVllredabs = 'no'
OIVllblueabs = 'no'

OIV2ndredabs = 'no'
OIV2ndblueabs = 'no'

OVIllredabs = 'no'
OVIllblueabs = 'no'

PVllredabs = 'no'
PVllblueabs = 'no'    

NeVIIIredabs = 'no'
NeVIIIblueabs = 'no'

NaIXllredabs = 'no'
NaIXllblueabs = 'no'

MgXllblueabs = 'no'
MgXllredabs = 'no'

Al3llblueabs = 'no'
Al3llredabs = 'no'


plot_em = 'no' # plot emission text and lines labels: 'yes' to include
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


#------ Weighted Averages
    
#Verner Tables:
CIVll = 1549.0524  #Carbon (C)
SiIVll = 1396.7470 #Silicon (Si)
NVll = 1240.1500  #Nitrogen (N)   
OVIll = 1033.8160 #Oxygen (O)
OIV2nd = 609.3506

#1600+ Angstroms, Verner Tables:
AlIII = 1857.4000   #Aluminum (Al)

#Below 1000 Angstroms, Verner Tables:
MgX = 614.7600        

#------ Weighted Doublets, Verner tables
CIVllred = 1550.7700 # red line of doublet for right (vmin) #Carbon
CIVllblue = 1548.1950 # blue line of doublet for left (vmax)

SiIVllred = 1402.770 #Silicon 
SiIVllblue = 1393.755

NVllred = 1242.80 #Nitrogen
NVllblue = 1238.82

OVIllred = 1037.6167 #Oxygen
OVIllblue = 1031.9261

PVllred = 1128.28   #Phosphorous - No weighted average
PVllblue = 1117.57

#1600+ Angstroms, Verner Tables:
Al3llred = 1862.7900  #Aluminum
Al3llblue = 1854.7160

#Below 1000 Angstroms, Verner Tables:
MgXllred = 624.9410  #Magnisium
MgXllblue = 609.7930

OIV2ndred = 609.8286   #Oxygen  
OIV2ndblue = 608.3968 

#------ From url link https://arxiv.org/pdf/1303.0043
NeVIIIllred = 780.3  #Neon
NeVIIIllblue = 770.4

NaIXllred = 694.3  #Sodium
NaIXllblue = 681.7

OIVllred =  553.3  #Oxygen
OIVllblue = 554.1 

    

#Emission wavelengths

CII_emitted = 1335.313 # (weighted average); individuals:
OI_emitted = 1303.4951 # weighted average; individuals page 20 in Verner Table

#Ly lines

Lya = 1215.6700
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
color = ['xkcd:shocking pink', 'xkcd:azure', 'blue', 'xkcd:purpleish blue', 'xkcd:slate', 'xkcd:shocking pink', \
         'black']
# color = ['red', 'green', 'blue', 'orange', 'purple']

for k in range(0,len(vmin)):
    
#Plot Averages 
    
    if CIVabs == 'yes':
        plt.axvspan(CIVllblue*(1.+zabs_max[k]),CIVllred*(1.+zabs_min[k]), alpha=0.2, color=color[5])
        plt.text(CIVll*(1.+zabs_min[k])-30.,1.135-0.1*k,'CIV',color=color[0],fontname='serif',weight='bold')

    if SiIVabs == 'yes':
        plt.axvspan(SiIVllblue*(1.+zabs_max[k]),SiIVllred*(1.+zabs_min[k]), alpha=0.2, color=color[3])
        plt.text(SiIVll*(1.+zabs_min[k])-30.,1.135-0.1*k,'SiIV',color=color[3],fontname='serif',weight='bold')
    
    if NVabs == 'yes':
        plt.axvspan(NVllblue*(1.+zabs_max[k]),NVllred*(1.+zabs_min[k]), alpha=0.2, color=color[1])
        plt.text(NVll*(1.+zabs_min[k]),1.8-0.1*k,'NV',color='xkcd:azure',fontname='serif',weight='bold')
    
    if OVIabs == 'yes':
        plt.axvspan(OVIllblue*(1.+zabs_max[k]),OVIllred*(1.+zabs_min[k]), alpha=0.2, color=color[2])
        plt.text(OVIll*(1.+zabs_min[k]),1.14-0.1*k,'OVI',color=color[2],fontname='serif',weight='bold')

    if AlIII == 'yes':
        plt.axvspan(Al3llblue*(1.+zabs_max[k]),Al3llred*(1.+zabs_min[k]), alpha=0.2, color=color[1])
        plt.text(AlIII*(1.+zabs_min[k])-20.,1.45-0.1*k,'AlIII',color=color[4],fontname='serif',weight='bold')

    if MgXabs == 'yes':
      plt.axvspan(MgXllblue*(1.+zabs_max[k]),MgXllred*(1.+zabs_min[k]), alpha=0.2, color=color[1])
      plt.text(MgX*(1.+zabs_min[k])-20.,1.45-0.1*k,'MgX',color=color[0],fontname='serif',weight='bold')

    if Lyaabs == 'yes':
        plt.axvspan(Lya*(1.+zabs_max[k]),Lya*(1.+zabs_min[k]), alpha=0.2, color=color[4])
        plt.text(Lya*(1.+zabs_min[k])-20.,1.45-0.1*k,'Lya',color=color[4],fontname='serif',weight='bold')

    if Lybabs == 'yes':
        plt.axvspan(Lyb*(1.+zabs_max[k]),Lyb*(1.+zabs_min[k]), alpha=0.2, color=color[4])
        plt.text(Lyb*(1.+zabs_min[k])-20.,1.45-0.1*k,'Lyb',color=color[4],fontname='serif',weight='bold')


#Plotting Doublets

#Carbon
    if CIVllredabs == 'yes' and CIVllblueabs == 'yes' :
        plt.axvspan(CIVllblue*(1.+zabs_max[k]),CIVllblue*(1.+zabs_min[k]), alpha=0.2, color=color[0])
        plt.axvspan(CIVllred*(1.+zabs_max[k]),CIVllred*(1.+zabs_min[k]), alpha=0.2, color=color[0])
        plt.text(CIVllblue*(1.+zabs_min[k])-30.,1.135-0.1*k,'CIV',color=color[0],fontname='serif',weight='bold')
  #     plt.text(CIVllred*(1.+zabs_min[k])-30.,1.135-0.1*k,'CIV',color=color[0],fontname='serif',weight='bold') 

#Silicon
    if SiIVllredabs == 'yes' and SiIVllblueabs == 'yes' :
        plt.axvspan(SiIVllblue*(1.+zabs_max[k]),SiIVllblue*(1.+zabs_min[k]), alpha=0.2, color=color[0])
        plt.axvspan(SiIVllred*(1.+zabs_max[k]),SiIVllred*(1.+zabs_min[k]), alpha=0.2, color=color[0])
        plt.text(SiIVllblue*(1.+zabs_min[k])-30.,1.135-0.1*k,'CIV',color=color[0],fontname='serif',weight='bold')
  #     plt.text(SiIVllred*(1.+zabs_min[k])-30.,1.135-0.1*k,'CIV',color=color[0],fontname='serif',weight='bold')

#Oxygen 4
    if OIVllredabs == 'yes' and OIVllblueabs == 'yes' :
        plt.axvspan(OIVllblue*(1.+zabs_max[k]),OIVllblue*(1.+zabs_min[k]), alpha=0.2, color=color[2])
        plt.axvspan(OIVllred*(1.+zabs_max[k]),OIVllred*(1.+zabs_min[k]), alpha=0.2, color=color[2])
        plt.text(OIVllblue*(1.+zabs_min[k]),2.5-0.1*k,'OIV',color=color[2],fontname='serif',weight='bold')
        # plt.text(OIVllred*(1.+zabs_min[k])-20.,1.45-0.1*k,'OIV',color=color[2],fontname='serif',weight='bold')


#Oxygen 4 2nd value
    if OIV2ndredabs == 'yes' and OIV2ndblueabs == 'yes' :
        plt.axvspan(OIV2ndblue*(1.+zabs_max[k]),OIV2ndblue*(1.+zabs_min[k]), alpha=0.2, color=color[2])
        plt.axvspan(OIV2ndred*(1.+zabs_max[k]),OIV2ndred*(1.+zabs_min[k]), alpha=0.2, color=color[2])
        plt.text(OIV2ndblue*(1.+zabs_min[k]),0.80-0.1*k,'OIV 2nd',color=color[2],fontname='serif',weight='bold')
        # plt.text(OIV2ndred*(1.+zabs_min[k]),1.14-0.1*k,'OIV 2nd',color=color[2],fontname='serif',weight='bold')

#Oxygen 6
    if OVIllredabs == 'yes' and OIVllblueabs == 'yes' :
        plt.axvspan(OVIllblue*(1.+zabs_max[k]),OVIllblue*(1.+zabs_min[k]), alpha=0.2, color=color[2])
        plt.axvspan(OVIllred*(1.+zabs_max[k]),OVIllred*(1.+zabs_min[k]), alpha=0.2, color=color[2])
        plt.text(OVIllblue*(1.+zabs_min[k]),1.14-0.1*k,'OVI',color=color[2],fontname='serif',weight='bold')
        # plt.text(OVIllred*(1.+zabs_min[k]),1.14-0.1*k,'OVI',color=color[2],fontname='serif',weight='bold')


#Neon 8
    if NeVIIIredabs == 'yes' and NeVIIIblueabs == 'yes' :
        plt.axvspan(NeVIIIllblue*(1.+zabs_max[k]),NeVIIIllblue*(1.+zabs_min[k]), alpha=0.2, color=color[0])
        plt.axvspan(NeVIIIllred*(1.+zabs_max[k]),NeVIIIllred*(1.+zabs_min[k]), alpha=0.2, color=color[0])
        plt.text(NeVIIIllblue*(1.+zabs_min[k]),2.7-0.1*k,'NeVIII',color=color[0],fontname='serif',weight='bold')
        # plt.text(NeVIIIllred*(1.+zabs_min[k])-20.,1.45-0.1*k,'NeVIII',color=color[0],fontname='serif',weight='bold')

#Sodium 9
    if NaIXllredabs == 'yes' and NaIXllblueabs == 'yes' :
        plt.axvspan(NaIXllblue*(1.+zabs_max[k]),NaIXllblue*(1.+zabs_min[k]), alpha=0.2, color=color[2])
        plt.axvspan(NaIXllred*(1.+zabs_max[k]),NaIXllred*(1.+zabs_min[k]), alpha=0.2, color=color[2])
        plt.text(NaIXllblue*(1.+zabs_min[k]),1.19-0.1*k,'NaIX',color=color[2],fontname='serif',weight='bold')
        # plt.text(NaIXllred*(1.+zabs_min[k])-20.,1.45-0.1*k,'NaIX',color=color[2],fontname='serif',weight='bold')
        
#Magnisium 10
    if MgXllredabs == 'yes' and MgXllblueabs == 'yes' :
        plt.axvspan(MgXllblue*(1.+zabs_max[k]),MgXllblue*(1.+zabs_min[k]), alpha=0.2, color=color[0])
        plt.axvspan(MgXllred*(1.+zabs_max[k]),MgXllred*(1.+zabs_min[k]), alpha=0.2, color=color[0])
        plt.text(MgXllblue*(1.+zabs_min[k]),2.7-0.1*k,'MgX',color=color[0],fontname='serif',weight='bold')
        # plt.text(MgXllred*(1.+zabs_min[k])-20.,1.45-0.1*k,'MgX',color=color[0],fontname='serif',weight='bold')

# NV 
    if NVllredabs == 'yes' and NVllblueabs == 'yes' :
        plt.axvspan(NVllblue*(1.+zabs_max[k]),NVllblue*(1.+zabs_min[k]), alpha=0.2, color=color[1])
        plt.axvspan(NVllred*(1.+zabs_max[k]),NVllred*(1.+zabs_min[k]), alpha=0.2, color=color[1])
        plt.text(NVllblue*(1.+zabs_min[k]),2.8-0.1*k,'NV',color='xkcd:azure',fontname='serif',weight='bold')
        plt.text(NVllred*(1.+zabs_min[k]),2.8-0.1*k,'NV',color='xkcd:azure',fontname='serif',weight='bold')
  
# PV 
    if PVllredabs == 'yes' and PVllblueabs == 'yes' :
        plt.axvspan(PVllblue*(1.+zabs_max[k]),PVllblue*(1.+zabs_min[k]), alpha=0.2, color=color[1])
        plt.axvspan(PVllred*(1.+zabs_max[k]),PVllred*(1.+zabs_min[k]), alpha=0.2, color=color[1])
        plt.text(PVllblue*(1.+zabs_min[k]),2.16-0.1*k,'PV',color=color[1],fontname='serif',weight='bold')
        plt.text(PVllred*(1.+zabs_min[k]),2.16-0.1*k,'PV',color=color[1],fontname='serif',weight='bold')


# AlIII
    if Al3llblueabs == 'yes' and Al3llredabs == 'yes' :
        plt.axvspan(Al3llblue*(1.+zabs_max[k]),Al3llblue*(1.+zabs_min[k]), alpha=0.2, color=color[1])
        plt.axvspan(Al3llred*(1.+zabs_max[k]),Al3llred*(1.+zabs_min[k]), alpha=0.2, color=color[1])
        plt.text(Al3llblue*(1.+zabs_min[k]),1.19-0.1*k,'AlIII',color=color[4],fontname='serif',weight='bold')
#        plt.text(Al3llred*(1.+zabs_min[k])-20.,1.45-0.1*k,'AlIII',color=color[4],fontname='serif',weight='bold')

    
#matplotlib.rcParams['font.sans-serif'] = ['Source Han Sans TW', 'sans-serif']

#Text for emission lines

if plot_em == 'yes' :

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

if plot_z == 'yes': 
    zem_plot = "z = " + str(zem)
    plt.text(zem_label_x, zem_label_y, zem_plot, bbox=dict(facecolor='none', edgecolor='black', pad=7.0))
 
plt.xlim(wavelength_observe1/coem,wavelength_observe2/coem) 
plt.ylim(0,topylim)

fig.tight_layout() 

plt.savefig (pp2)





# plt.savefig(os.getcwd() + '/PRESENTATION_PLOTS/OUTPUT_FILES/' + pp2, dpi=100)
