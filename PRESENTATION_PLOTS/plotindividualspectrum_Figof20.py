# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 11:33:16 2017

This code plots spectra

@author: Paola
"""

import numpy as np
import matplotlib.pyplot as plt 
import os

from matplotlib.backends.backend_pdf import PdfPages



## Inputs/Outputs to change -------------------
specnum=1 
# specdirec='/Users/Paola/SDSSIII/DR9Q/FILES_FOR_ABDUL/' #Set location of spectrum files
# specdirec='/Users/Paola/QUASAR/Work_EHVO/ROUTINES/NORM/'
specdirec = os.getcwd() + '/../' + "EHVO_NORM_DR16Q/"
pp2= PdfPages('EHVO1.pdf') # normalized pdf

spectra='spec-6806-56419-0302norm.dr16'


zem=2.671
topy=2
topemlabel=1.8 # This is the location where you want to place the ion labels

n=3 # smooth box car
wavelength_emit1_initial=1290.  # in restframe
wavelength_emit2_initial=1500.  # in restframe

vmin = [-30701.4298, -39761.2702]
vmax = [-31931.3069, -41456.5876]

OVIem='no'

NVabs='no'
OVIabs='no'
SiIVabs='no'
Lyaabs='no'

#-------------------

c = 300000.
CIVll = 1549.0524 # avr_CIV_doublet 
NVll = 1240.15  # avr_NV_doublet
OVIll = 1033.8160 # avr_OVI_doublet
SiIVll = 1396.747 # avr_SiIV_doublet

# I am going to use blue line of the doublet for left (vmax)
# and red line for right (vmin) 
CIVllred = 1550.7700
CIVllblue = 1548.1950
avr_CIV_doublet = 1549.0524 #weighted average

avr_SiIV_doublet = 1396.747 # weighted average; individuals: 1402.770, 1393.755
SiIVllred = 1402.770
SiIVllblue = 1393.755

CII_emitted = 1335.313 # (weighted average); individuals:
OI_emitted = 1303.4951 # weighted average; individuals pag 20 in Verner Table

avr_NV_doublet = 1240.15 # weighted average; individuals: 1242.80, 1238.82
NVllred = 1242.80
NVllblue = 1238.82

avr_OVI_doublet=1033.8160 # weighted average; individuals: 1037.6167, 1031.9261
OVIllred = 1037.6167
OVIllblue = 1031.9261

Lya=1215.6700
Lyb=1025.7222


def smooth(norm_flux, box_pts):   
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(norm_flux, box, mode='same')
    return y_smooth

def zabs(v):
    beta=-v/c
    Rc=np.sqrt((beta+1.)/(1.-beta))
    za=((1.+zem)/Rc)-1.
    return za
    
spectra_count = 1
print('spec_name' + ": " + spectra)
print(spectra_count)
print('zem=',zem)
    
data='' 
data = np.loadtxt(specdirec+spectra) 
  
wavelength_observe1 =(zem+1.)*wavelength_emit1_initial #Shift start wavelength into frame |<--This makes our wavelength range for
wavelength_observe2 =(zem+1.)*wavelength_emit2_initial #Shift end wavelength into frame   |     the region we want to look at
   
wavelength_lower_limit = np.where(data[:,0] >wavelength_observe1)
wavelength_upper_limit = np.where(data[:,0] <wavelength_observe2)

minwave= np.min(wavelength_lower_limit[0])
maxwave= np.max(wavelength_upper_limit[0])
wavelength = data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit [0]),0] #Get wavelengths in out data set that fall into our region of study
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
    


# ----------------------------------------------------------------------
  


    #SOMETIMES, THERE ARE PIXEL PROBLEMS, AND WE MIGHT GET AN ERROR OF 30 IN FLUX. TO AVOID THAT, WE HAVE DONE THIS. MESSED UP ERROR IS 
    ####       DEFINED ABOVE.
    #if len (messed_up_error[0]) > 0:######################################################original
        #plerror[messed_up_error[0]]=0####################################################
	#flux[messed_up_error[0]]=0

fig, ay1 = plt.subplots()

#ay1 = fig.add_subplot(1, 1, 1)

#plt.title(spectrum)
ay1.set_xlabel(r"Observed Wavelength [$\rm \AA$]")
ay1.set_ylabel(r"Flux [10$^{-17}$ erg/cm$^2$/$\rm \AA$]")
     
ay1.plot (wavelength, smooth(normflux,n),'k-')
ay1.plot (wavelength, error_normflux,'k--') 

plt.plot([wavelength_observe1,wavelength_observe2],[1,1],'r--')
#plt.plot([wavelength_observe1,wavelength_observe2],[0.9,0.9],'r--')


for k in range(0,len(vmin)):

    plt.axvspan(CIVll*(1.+zabs_max[k]),CIVll*(1.+zabs_min[k]), alpha=0.2, color='red')
    plt.text(CIVll*(1.+zabs_min[k])-30.,0.5-0.1*k,'CIV',color='red',fontname='serif',weight='bold')

    if NVabs == 'yes':
        plt.axvspan(NVllblue*(1.+zabs_max[k]),NVllred*(1.+zabs_min[k]), alpha=0.2, color='green')
        plt.text(NVll*(1.+zabs_min[k]),1.3-0.1*k,'NV',color='green',fontname='serif',weight='bold')
    
    if (OVIabs == 'yes'): # and k == 0):
        plt.axvspan(OVIllblue*(1.+zabs_max[k]),OVIllred*(1.+zabs_min[k]), alpha=0.2, color='blue')
        plt.text(OVIll*(1.+zabs_min[k]),1.4-0.1*k,'OVI',color='blue',fontname='serif',weight='bold')

    if SiIVabs == 'yes':
        plt.axvspan(SiIVllblue*(1.+zabs_max[k]),SiIVllred*(1.+zabs_min[k]), alpha=0.2, color='orange')
        plt.text(SiIVll*(1.+zabs_min[k]),0.25,'SiIV',color='orange',fontname='serif',weight='bold')
        
    if (Lyaabs == 'yes'): # and k == 0):
        plt.axvspan(Lya*(1.+zabs_max[k]),Lya*(1.+zabs_min[k]), alpha=0.2, color='purple')
        plt.text(Lya*(1.+zabs_min[k])-20.,1.45-0.1*k,'Lya',color='purple',fontname='serif',weight='bold')



#matplotlib.rcParams['font.sans-serif'] = ['Source Han Sans TW', 'sans-serif']

plt.text(1549.0*(1+zem)-30,topemlabel ,'CIV',color='black',rotation = 90,fontname='serif')  #weight='bold'
plt.text(1402.770*(1+zem)-40.,topemlabel,'SiIV+OIV]',color='black',rotation = 90,fontname='serif')
#text(1242.804*(1+zem)+10.,topemlabel,'NV',color='black',rotation=90)
alpha = 'Ly' + chr(945)
plt.text(1242.804*(1+zem)+30.,topemlabel, alpha + '+NV' ,color='black',rotation = 90,fontname='serif')
plt.text(1304.8576*(1+zem)-35.,topemlabel,'OI',color='black',rotation=90,fontname='serif')
plt.text(1334.5323*(1+zem)-30.,topemlabel,'CII',color='black',rotation=90,fontname='serif')

if OVIem == 'yes':
    plt.text(OVIll*(1+zem)-30.,topemlabel,'OVI',color='black',rotation=90,fontname='serif')

## How to get to have the right values of restframe??? Relate to em redshift, how?? 

coem=zem+1.
plt.xlim(wavelength_observe1,wavelength_observe2) 

ay2 = ay1.twiny()

ay2.plot(wavelength/coem,1000*smooth(flux,n)+1000.)
ay2.set_xlabel(r"Restframe Wavelength [$\rm \AA$]")

ay1.xaxis.set_label_coords(0.48, -0.08)
ay2.xaxis.set_label_coords(0.48, 1.11)
ay2.xaxis.set_major_locator(plt.MaxNLocator(5))

plt.xlim(wavelength_observe1/coem,wavelength_observe2/coem) 
plt.ylim(0,topy)

fig.tight_layout()

#plot((1393.755*(1+zabs),1393.755*(1+zabs)),(0,max(actual_flux)),'g--')
#plot((1402.770*(1+zabs),1402.770*(1+zabs)),(0,max(actual_flux)),'g--')
#text(1402.770*(1+zabs)+10.,topemlabel+2,'SiIV+OIV]',color='green')

#plot((1302.1685*(1+zabs),1302.1685*(1+zabs)),(0,max(actual_flux)),'b--')
#plot((1304.8576*(1+zabs),1304.8576*(1+zabs)),(0,max(actual_flux)),'b--')
#text(1304.8576*(1+zabs)+10.,topemlabel,'OI',color='blue')
#plot((1334.5323*(1+zabs),1334.5323*(1+zabs)),(0,max(actual_flux)),'b--')
#text(1334.5323*(1+zabs)+10.,topemlabel,'CII',color='blue')

#plot((1548.195*(1+zabs),1548.195*(1+zabs)),(0,max(actual_flux)),'r--')
#plot((1550.770*(1+zabs),1550.770*(1+zabs)),(0,max(actual_flux)),'r--')
#text(1548.195*(1+zabs)+10.,topemlabel,'CIV',color='red')

#plot((1215.6701*(1+zabs),1215.6701*(1+zabs)),(0,max(actual_flux)),'c--')
#text(1215.6701*(1+zabs)+10.,topemlabel+2,'Lya',color='cyan')

#plot((1238.821*(1+zabs),1238.821*(1+zabs)),(0,max(actual_flux)),'b--')
#plot((1242.804*(1+zabs),1242.804*(1+zabs)),(0,max(actual_flux)),'b--')
#text(1242.804*(1+zabs)+10.,topemlabel,'NV',color='blue')

#plot((1031.927*(1+zabs),1031.927*(1+zabs)),(0,max(actual_flux)),'m--')
#plot((1037.616*(1+zabs),1037.616*(1+zabs)),(0,max(actual_flux)),'m--')
#text(1037.616*(1+zabs)+10.,topemlabel,'OVI',color='magenta')

pp2.savefig()
 
pp2.close()

    
