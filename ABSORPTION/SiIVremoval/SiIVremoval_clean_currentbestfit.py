"""
Created on Tue Feb 25 17:47:32 2025

@author: lilianaflores
"""

#############################################################

from astropy.io import fits
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from functools import partial
#from functions import wavelength_to_velocity_adj ----- I don't have functions file with these functions. I made an adjustable wavelength to velocity function for myself. - LEF
from SiV_functions import wavelength_to_velocity_adj #Paola could comment this out and uncomment above ^^^^


#__________________________________________________________
# INPUTS & VARIABLES TO CHANGE:

# Location of DATA: (User will need to adjust for file path?)
path = os.getcwd()
spec_direc = path + '/DATA/'

filename = 'J2157_Zappacosta_Onken_median_VIS_NIR_normalized.fits' #currently (1/21) used normalized data

# redshift value
redshift = 4.692 

#--------------------------------
#  INPUTS & VARIABLES NOT TO CHANGE:
# Necessary data from Verner table
wavelength_CIV_emit1 = 1548.1950
wavelength_CIV_emit2 = 1550.7700
wavelength_SiIV_emit1 = 1393.755
wavelength_SiIV_emit2 = 1402.770
avr_CIV_doublet = 1549.0524 # weighted average
avr_SiIV_doublet = 1396.747 # weighted average

#_______________________________________________________
# FUNCTIONS: 

def tau_v(v,v0,b,tau0):
    return tau0*np.exp(-(v-v0)**2./(b**2.))

#function for curve fitting down where used

#-------------------------------------
# MAIN CODE:
    
    
#...................................................
# Read data information from file:

file = fits.open(spec_direc + filename, memmap=True)
data = file[1].data
#print(file[1].columns)  # Do this first time to check if you have the names right
# letting you know what the column headers are 'Restframe wavelength' and 'Flux normalized'


wavelength = data['Restframe wavelength']
flux = data['Restframe Flux normalized']


file.close()


##############################################################################

########### Drawing SiIV doublets in SiIV velocity in CIV EHVO ##################

#...................................................
# Convert wavelength to velocity:
    
redshift = 0 # Because data is already restframed 

velo_SiIV = wavelength_to_velocity_adj(redshift, wavelength, wavelength_SiIV_emit1)

# .........................
# Select region for SiIV 

x1 = np.min(np.where(velo_SiIV > 1000))
x2 = np.max(np.where(velo_SiIV < -8800))

xfit = velo_SiIV[x2:x1]
yfit = flux[x2:x1]

w1 = np.min(np.where(wavelength > 1370))
w2 = np.max(np.where(wavelength < 1400))

#...................................................
# calculate vdiff

wavelengthSi = np.array([np.mean(wavelength[w1:w2])])

vdiff_SiIV = wavelength_to_velocity_adj(redshift, wavelengthSi,
                    wavelength_SiIV_emit1) - wavelength_to_velocity_adj(
                    redshift, wavelengthSi, wavelength_SiIV_emit2)
             
                     
vdiff =  vdiff_SiIV

# Plot the original spectrum
plt.xlim(-15000,1000)
plt.ylim(0,1.25)
plt.plot(velo_SiIV, flux, 'k')

v=xfit
Cf=0.96
I0=1
vdiff=vdiff_SiIV
tau_ratio=2


tau0 = 0.68
v0 = -7573.098
b = 85.784
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

tau02 = 0.3
v02 = -4900.
b2 = 200.
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v02,b2,tau02)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

tau03 = 0.85
v03 = -4620.
b3 = 200.
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v03,b3,tau03)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v03+vdiff),b3,tau03/tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

tau04 = 0.05
v04 = -4150.
b4 = 380.
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v04,b4,tau04)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v04+vdiff),b4,tau04/tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

tau05 = 0.5
v05 = -3750.
b5 = 190.
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v05,b5,tau05)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v05+vdiff),b5,tau05/tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

tau06 = 0.8
v06 = -3400.
b6 = 220.
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v06,b6,tau06)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v06+vdiff),b6,tau06/tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

tau07 = 0.15
v07 = -3200.
b7 = 200.
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v07,b7,tau07)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v07+vdiff),b7,tau07/tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   



#combined curve
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)) * np.exp(-tau_v(v,v02,b2,tau02)) * np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio))* 
                                           np.exp(-tau_v(v,v03,b3,tau03)) * np.exp(-tau_v(v,(v03+vdiff),b3,tau03/tau_ratio))* np.exp(-tau_v(v,v04,b4,tau04)) * np.exp(-tau_v(v,(v04+vdiff),b4,tau04/tau_ratio))* 
                                           np.exp(-tau_v(v,v05,b5,tau05)) * np.exp(-tau_v(v,(v05+vdiff),b5,tau05/tau_ratio)) * np.exp(-tau_v(v,v06,b6,tau06)) * np.exp(-tau_v(v,(v06+vdiff),b6,tau06/tau_ratio))* 
                                           np.exp(-tau_v(v,v07,b7,tau07)) * np.exp(-tau_v(v,(v07+vdiff),b7,tau07/tau_ratio)))), color='purple' )

plt.title('Drawing SiIV Doublets')
plt.xlabel('Velocity (km/s)')
plt.ylabel('Normalized Flux')

plt.figure()
#stop

##############################################################################

########### Drawing SiIV doublets in SiIV velocity in CIV EHVO zoomed ##################

plt.xlim(-8500,1000)
plt.ylim(0,1.25)
plt.plot(velo_SiIV, flux, 'k')


Cf=1
I0=1
vdiff=vdiff_SiIV
tau_ratio=2

plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v0, b, tau0)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v0 + vdiff, b, tau0 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v02, b2, tau02)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v02 + vdiff, b2, tau02 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v03, b3, tau03)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v03 + vdiff, b3, tau03 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v04, b4, tau04)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v04 + vdiff, b4, tau04 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v05, b5, tau05)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v05 + vdiff, b5, tau05 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v06, b6, tau06)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v06 + vdiff, b6, tau06 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v07, b7, tau07)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v07 + vdiff, b7, tau07 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

#combined SiIV curve
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)) * np.exp(-tau_v(v,v02,b2,tau02)) * np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio))* 
                                           np.exp(-tau_v(v,v03,b3,tau03)) * np.exp(-tau_v(v,(v03+vdiff),b3,tau03/tau_ratio))* np.exp(-tau_v(v,v04,b4,tau04)) * np.exp(-tau_v(v,(v04+vdiff),b4,tau04/tau_ratio))* 
                                           np.exp(-tau_v(v,v05,b5,tau05)) * np.exp(-tau_v(v,(v05+vdiff),b5,tau05/tau_ratio)) * np.exp(-tau_v(v,v06,b6,tau06)) * np.exp(-tau_v(v,(v06+vdiff),b6,tau06/tau_ratio))* 
                                           np.exp(-tau_v(v,v07,b7,tau07)) * np.exp(-tau_v(v,(v07+vdiff),b7,tau07/tau_ratio)))), color='purple' )

plt.title('Drawing SiIV Doublets')
plt.xlabel('Velocity (km/s)')
plt.ylabel('Normalized Flux')

plt.figure()
#stop

########################################################################

################### plotting overlay SiIV doublets in CIV in CIV velocity ####################################

#...................................................
redshift = 0 # Because data is already restframed 

velo_CIV = wavelength_to_velocity_adj(redshift, wavelength, wavelength_CIV_emit1)

# .........................
# Select region for CIV 

x1 = np.min(np.where(velo_CIV > -1000))
x2 = np.max(np.where(velo_CIV < -8000))

xfit = velo_CIV[x2:x1]
yfit = flux[x2:x1]

w1 = np.min(np.where(wavelength > 1535))
w2 = np.max(np.where(wavelength < 1520))

#...................................................
# calculate vdiff

wavelengthC = np.array([np.mean(wavelength[w2:w1])])

vdiff_CIV = wavelength_to_velocity_adj(redshift, wavelengthC,
                    wavelength_CIV_emit1) - wavelength_to_velocity_adj(
                    redshift, wavelengthC, wavelength_CIV_emit2)
             
                     
plt.xlim(-8000, -1000)
plt.ylim(0,1.25)
plt.plot(velo_CIV, flux, 'k')

minflux=np.min(yfit)
print(f'Min flux: {minflux}')
v=xfit

Cf=  1-minflux   # was 0.98
I0=1
vdiff=vdiff_CIV
tau_ratio=2

tau0 = 0.15
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v0, b, tau0)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v0 + vdiff, b, tau0 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

tau02 = 4.5
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v02, b2, tau02)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v02 + vdiff, b2, tau02 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

tau03 = 2.8
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v03, b3, tau03)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v03 + vdiff, b3, tau03 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

tau04 = 2.3
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v04, b4, tau04)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v04 + vdiff, b4, tau04 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

tau05 = 2.1
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v05, b5, tau05)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v05 + vdiff, b5, tau05 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

tau06 = 2.5
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v06, b6, tau06)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v06 + vdiff, b6, tau06 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

tau07 = 2.5
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v07, b7, tau07)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v07 + vdiff, b7, tau07 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   


#combined curve
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)) * np.exp(-tau_v(v,v02,b2,tau02)) * np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio))* 
                                           np.exp(-tau_v(v,v03,b3,tau03)) * np.exp(-tau_v(v,(v03+vdiff),b3,tau03/tau_ratio))* np.exp(-tau_v(v,v04,b4,tau04)) * np.exp(-tau_v(v,(v04+vdiff),b4,tau04/tau_ratio))* 
                                           np.exp(-tau_v(v,v05,b5,tau05)) * np.exp(-tau_v(v,(v05+vdiff),b5,tau05/tau_ratio)) * np.exp(-tau_v(v,v06,b6,tau06)) * np.exp(-tau_v(v,(v06+vdiff),b6,tau06/tau_ratio))* 
                                           np.exp(-tau_v(v,v07,b7,tau07)) * np.exp(-tau_v(v,(v07+vdiff),b7,tau07/tau_ratio)))), color='purple' )


plt.title('Plotting Overlay SiIV Doublets in CIV')
plt.xlabel('Velocity (km/s)')
plt.ylabel('Normalized Flux')

plt.figure()

#stop

######################################################################

############ fitting entire CIV EHVO with SiIV in SiIV velocity

redshift = 0 # Because data is already restframed 

### whole CIV trough

x1 = np.min(np.where(velo_SiIV > 2000))
x2 = np.max(np.where(velo_SiIV < -11800))

xfit = velo_SiIV[x2:x1]
yfit = flux[x2:x1]

w1 = np.min(np.where(wavelength > 1370))
w2 = np.max(np.where(wavelength < 1400))

plt.xlim(-15000,2000)
plt.ylim(0,1.25)
plt.plot(velo_SiIV, flux, 'k')


#fixed parameters
Cf=1
I0=1
tau_ratio=2
v=xfit

#SiIV doublet tau guesses (other parameters, v0 & b, are fixed and called from above sections ^^)
tau0 = 0.4
tau02= 0.15 
tau03= 0.55 
tau04= 0.2  
tau05= 0.5
tau06= 0.40 
tau07= 0.01 


#CIV doublet parameter guesses
v0_CIV = -9400
tau0_CIV = 0.17
b_CIV = 750

v02_CIV=-7000
tau02_CIV=0.35
b2_CIV=2700

v03_CIV = -1900
tau03_CIV = 0.25
b3_CIV = 2100

#function is defined down here since I need to be able to reference the ordr of parameters to be able to fix values
def curve_func_all(v, tau0, tau02, tau03, tau04, tau05, tau06, tau07, tau0_CIV, tau02_CIV, tau03_CIV, 
                   v0_CIV, v02_CIV, v03_CIV, b_CIV, b2_CIV, b3_CIV, Cf, v0, v02, v03, v04, v05, v06, v07, 
                   b, b2, b3, b4, b5, b6, b7, vdiff_SiIV, vdiff_CIV, I0, tau_ratio):
    return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff_SiIV),b,tau0/tau_ratio)) * np.exp(-tau_v(v,v02,b2,tau02)) * np.exp(-tau_v(v,(v02+vdiff_SiIV),b2,tau02/tau_ratio))* 
                                       np.exp(-tau_v(v,v03,b3,tau03)) * np.exp(-tau_v(v,(v03+vdiff_SiIV),b3,tau03/tau_ratio)) * np.exp(-tau_v(v,v04,b4,tau04)) * np.exp(-tau_v(v,(v04+vdiff_SiIV),b4,tau04/tau_ratio)) * 
                                       np.exp(-tau_v(v,v05,b5,tau05)) * np.exp(-tau_v(v,(v05+vdiff_SiIV),b5,tau05/tau_ratio)) * np.exp(-tau_v(v,v06,b6,tau06)) * np.exp(-tau_v(v,(v06+vdiff_SiIV),b6,tau06/tau_ratio)) * 
                                       np.exp(-tau_v(v,v07,b7,tau07)) * np.exp(-tau_v(v,(v07+vdiff_SiIV),b7,tau07/tau_ratio)) * np.exp(-tau_v(v,v0_CIV,b_CIV,tau0_CIV)) * np.exp(-tau_v(v,(v0_CIV+vdiff_CIV),b_CIV,tau0_CIV/tau_ratio)) * 
                                       np.exp(-tau_v(v,v02_CIV,b2_CIV,tau02_CIV)) * np.exp(-tau_v(v,(v02_CIV+vdiff_CIV),b2_CIV,tau02_CIV/tau_ratio))* np.exp(-tau_v(v,v03_CIV,b3_CIV,tau03_CIV)) * np.exp(-tau_v(v,(v03_CIV+vdiff_CIV),b3_CIV,tau03_CIV/tau_ratio))))



#Fitting with free tau, CIV: v0 and b  (Cf = 1)
(tau0, tau02, tau03, tau04, tau05, tau06, tau07, tau0_CIV, tau02_CIV, tau03_CIV, v0_CIV, v02_CIV, v03_CIV, b_CIV, b2_CIV, b3_CIV), covar = curve_fit(partial(curve_func_all, Cf=Cf, v0=v0, v02=v02, v03=v03, v04=v04, v05=v05, v06=v06, v07=v07, b=b, b2=b2, b3=b3, b4=b4, b5=b5, b6=b6, b7=b7,  
vdiff_SiIV=vdiff_SiIV, vdiff_CIV=vdiff_CIV, I0=I0, tau_ratio=tau_ratio), xfit, yfit, p0=[tau0, tau02, tau03, tau04, tau05, tau06, tau07, tau0_CIV, tau02_CIV, tau03_CIV, v0_CIV, v02_CIV, v03_CIV, b_CIV, b2_CIV, b3_CIV], 
maxfev=10000, bounds=((0.,0.,0.,0.,0.,0.,0.,0.,0.,0., -np.inf, -np.inf, -np.inf,0.,0.,0.), (10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,0., 0., 0.,np.inf,np.inf,np.inf)))
plt.plot(xfit, (curve_func_all(v, tau0, tau02, tau03, tau04, tau05, tau06, tau07, tau0_CIV, tau02_CIV, tau03_CIV, 
                    v0_CIV, v02_CIV, v03_CIV, b_CIV, b2_CIV, b3_CIV, Cf, v0, v02, v03, v04, v05, v06, v07, 
                   b, b2, b3, b4, b5, b6, b7,  
                   vdiff_SiIV, vdiff_CIV, I0, tau_ratio)), label='Curve Fit', color='purple')

pcov= np.sqrt(np.diag(covar))


print()
print('Fitting all doublets:')
print(f"tau0_SiIV: {tau0} +/- {pcov[0]}")
print(f"tau02_SiIV: {tau02} +/- {pcov[1]}")
print(f"tau03_SiIV: {tau03} +/- {pcov[2]}")
print(f"tau04_SiIV: {tau04} +/- {pcov[3]}")
print(f"tau05_SiIV: {tau05} +/- {pcov[4]}")
print(f"tau06_SiIV: {tau06} +/- {pcov[5]}")
print(f"tau07_SiIV: {tau07} +/- {pcov[6]}")
print()
print(f"tau0_CIV: {tau0_CIV} +/- {pcov[7]}")
print(f"tau02_CIV: {tau02_CIV} +/- {pcov[8]}")
print(f"tau03_CIV: {tau03_CIV} +/- {pcov[9]}")
print()


print(f"Cf: {Cf}")
print()

print(f"v0_CIV: {v0_CIV} +/- {pcov[10]}")
print(f"v02_CIV: {v02_CIV} +/- {pcov[11]}")
print(f"v03_CIV: {v03_CIV} +/- {pcov[12]}")
print()
print(f"v0_SiIV: {v0}")
print(f"v02_SiIV: {v02}")
print(f"v03_SiIV: {v03}")
print(f"v04_SiIV: {v04}")
print(f"v05_SiIV: {v05}")
print(f"v06_SiIV: {v06}")
print(f"v07_SiIV: {v07}")
print()

print(f"b_CIV: {b_CIV} +/- {pcov[13]}")
print(f"b2_CIV: {b2_CIV} +/- {pcov[14]}")
print(f"b3_CIV: {b3_CIV} +/- {pcov[15]}")

print()
print(f"b_SiIV: {b}")
print(f"b_SiIV: {b2}")
print(f"b_SiIV: {b3}")
print(f"b_SiIV: {b4}")
print(f"b_SiIV: {b5}")
print(f"b_SiIV: {b6}")
print(f"b_SiIV: {b7}")

#plotting CIV doublets
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0_CIV,b_CIV,tau0_CIV)))), ':', label = 'Gaussian One', color = 'blue')
plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v0_CIV+vdiff_CIV),b_CIV,tau0_CIV/tau_ratio)))), ':', label = 'Gaussian Two', color = 'red')       

plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v02_CIV,b2_CIV,tau02_CIV)))), ':', label = 'Gaussian Three', color = 'blue')       
plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v02_CIV+vdiff_CIV),b2_CIV,tau02_CIV/tau_ratio)))), ':', label = 'Gaussian Four', color = 'red')       

plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v03_CIV,b3_CIV,tau03_CIV)))), ':', label = 'Gaussian Three', color = 'blue')       
plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v03_CIV+vdiff_CIV),b3_CIV,tau03_CIV/tau_ratio)))), ':', label = 'Gaussian Four', color = 'red')

#plotting SiIV doublets
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v0, b, tau0)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v0 + vdiff_SiIV, b, tau0 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v02, b2, tau02)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v02 + vdiff_SiIV, b2, tau02 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v03, b3, tau03)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v03 + vdiff_SiIV, b3, tau03 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v04, b4, tau04)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v04 + vdiff_SiIV, b4, tau04 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v05, b5, tau05)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v05 + vdiff_SiIV, b5, tau05 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v06, b6, tau06)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v06 + vdiff_SiIV, b6, tau06 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   

plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v07, b7, tau07)))), '--', label = 'Gaussian One', color = 'b')
plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v, v07 + vdiff_SiIV, b7, tau07 / tau_ratio)))), '--', label = 'Gaussian Two', color = 'r')   


plt.title('Fitting: free tau, CIV v0, & CIV b (Cf = 1)')
plt.xlabel('Velocity (km/s)')
plt.ylabel('Normalized Flux')

plt.figure()

stop
#Column density calculations
from column_density_calcs import column_density

print()
print('Column Density Calculations:')

tau_values_SiIV = [tau0, tau02, tau03, tau04, tau05, tau06, tau07]  
tau_values_CIV = [tau0_CIV, tau02_CIV, tau03_CIV]
Cf_value = Cf  
v0_values_CIV = [v0_CIV, v02_CIV, v03_CIV] 
v0_values_SiIV = [v0, v02, v03, v04, v05, v06, v07]  
b_values_CIV = [b_CIV, b2_CIV, b3_CIV] 
b_values_SiIV = [b, b2, b3, b4, b5, b6, b7]  

values_SiIV = []
for i in range(len(tau_values_SiIV)):
    values_SiIV.append([tau_values_SiIV[i], Cf_value, v0_values_SiIV[i], b_values_SiIV[i]])

values_CIV = []
for i in range(len(tau_values_CIV)):
    values_CIV.append([tau_values_CIV[i], Cf_value, v0_values_CIV[i], b_values_CIV[i]])


headers1 = ['tau_SiIV', 'Cf_SiIV', 'v0_SiIV', 'b_SiIV']
headers2 = ['tau_CIV', 'Cf_CIV', 'v0_CIV', 'b_CIV']

np.savetxt('SiIV_parameters.txt', values_SiIV, fmt='%f', delimiter='\t', header='\t'.join(headers1))
np.savetxt('CIV_parameters.txt', values_CIV, fmt='%f', delimiter='\t', header='\t'.join(headers2))


col_den_SiIV =[]

for i in range(len(tau_values_SiIV)):
    v0, Cf, b, tau0, logNion_m1, logNion_m2, logNH_m1, logNH_m2 = column_density(v0_values_SiIV[i], b_values_SiIV[i], tau_values_SiIV[i], Cf_value, xfit)
    col_den_SiIV.append(logNH_m1)

col_den_CIV =[]
for i in range(len(tau_values_CIV)):
    v0, Cf, b, tau0, logNion_m1, logNion_m2, logNH_m1, logNH_m2 = column_density(v0_values_CIV[i], b_values_CIV[i], tau_values_CIV[i], Cf_value, xfit)
    col_den_CIV.append(logNH_m1)

logNH_SiIV_data = np.array(col_den_SiIV)
logNH_SiIV = np.log10(sum(10**logNH_SiIV_data))
     
logNH_CIV_data = np.array(col_den_CIV)
logNH_CIV = np.log10(sum(10**logNH_CIV_data))

#subtracted = np.log10(sum(10**logNH_CIV_data)-sum(10**logNH_SiIV_data))


#print(f'Individual doublet log(NH) SiIV: {logNH_SiIV_data}')
#print(f'Total log(NH) for SiIV: {logNH_SiIV}')
print()
print(f'Individual doublet log(NH) CIV: {logNH_CIV_data}')
print(f'Total log(NH) for CIV: {logNH_CIV}')
print()
#print(f'log(NH) of CIV EHVO without SiIV: {subtracted}')





