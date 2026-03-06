#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 17:16:10 2024

@author: lilianaflores


Key Features:
    - Fits CIV EHVO absorption troughs for all three observations of J2318.
    - Applies sigma clipping function and additional spectral masking and re-addition of data points to
      improve fit accuracy.
    - Generates and displays plots for each observation.
    - Allows toggling between multiple spectral data types, including
      smoothed, unsmoothed, and alternative power-law.

Inputs:
    - Twelve .dat files enabling selection among four spectral variants
      for each of the three observations.

Outputs:
    - Fitting plots for each observation, saved to a specified
      output directory.
    - .txt files containing fitted spectra (wavelength, flux, error),
      formatted for direct input into absorption.py to calculate absorption measurements.
    - .csv files containing CIV fitting parameters, to be used as initial
      guesses and fixed values when fitting the SiIV EHVOs.

Note: To obtain results from the f-test you must rerun the program changing 'doublets' in the 
changeable variables as 1 then again as 2 then 3 then 4 so that fittings done using 
1, 2, 3, and 4 doublets can be compared.
"""

# Imports
import os
import sys
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
sys.path.insert(0,os.getcwd()+'/../')
from abs_function_module import wavelength_to_velocity
sys.path.insert(0,os.getcwd()+'/../../')
from utility_functions import read_spectra
from sigma_clipping_functions import ChoiClip
from Fastest_CIV_EHVO_Functions import plot_masked_regions, f_test 
from curve_fit_function import plotting_curvefit_test, curve_fit_area, tau_v


# Defining fitting equations:
def curve_func( v, tau0, v0, b, Cf, I0, vdiff, tau_ratio): 
    return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))) 
    
    
def curve_func2( v, tau0, tau02, v0, v02, b, b2, Cf, I0, vdiff, tau_ratio): 
    return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)) * np.exp(-tau_v(v,v02,b2,tau02)) * np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)))) 


def curve_func3( v, tau0, tau02, tau03, v0, v02, v03, b, b2, b3, Cf, I0, vdiff, tau_ratio): 
    return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)) * np.exp(-tau_v(v,v02,b2,tau02)) * np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)) *np.exp(-tau_v(v,v03,b3,tau03)) * np.exp(-tau_v(v,(v03+vdiff),b3,tau03/tau_ratio))))
     
def curve_func4( v, tau0, tau02, tau03, tau04, v0, v02, v03, v04, b, b2, b3, b4, Cf, I0, vdiff, tau_ratio): 
    return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)) * np.exp(-tau_v(v,v02,b2,tau02)) * np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)) *np.exp(-tau_v(v,v03,b3,tau03)) * np.exp(-tau_v(v,(v03+vdiff),b3,tau03/tau_ratio)) *np.exp(-tau_v(v,v04,b4,tau04)) * np.exp(-tau_v(v,(v04+vdiff),b4,tau04/tau_ratio))))


#.....................................................................................................................
#.....................................................................................................................
# CHANGEABLE VARIABLES

# Establishing directories
OUTDIREC = os.getcwd() + '/output_files/'

DATADIREC = os.getcwd() + '/data/'

f_test_OUTDIREC = os.getcwd() + '/output_files/f_test'

# save figures and .txt files of spectra
want_save_fig = 'y'
want_txt_save = 'y'

#want_save_fig = 'n'
#want_txt_save = 'n'

# adding back in points for obs 59188
#readd = 'no'
readd = 'yes'

# manual masking some points in obs 57328
manual_mask = True
#manual_mask = False

# run and plot in wavelength or velocity
wavelength_or_velocity = 'v'
#wavelength_or_velocity = 'W'

# Redshift set to zero as spectra are provided in restframe
redshift = 0

# establishing set plot limits/settings
x_min = -130000
x_max = -50000
y_min = 0#0.2
y_max = 1.3#1.1
xlim =(-105000,-60000)
xtick = 10000


# Different plot tile choices for figures
title = 'yr' #CIV EHVO: year
title = 'mjd' #CIV EHVO (mjd:#, year:# ) Sigma Clip
title = 'org' #J2318_obs#_DZ/UZ_RP/RP2norm (Sigma Clip std=X, GK =X)
title = 'paper' #J2318 MJD:# CIV EHVO

# singma clipping function settings
sigma_smooth = 3
sigma_clip = 0.85

# fitting attempts with differnt numbers of doublets for f-test
# when fitting attempts with other numbers of doublets set the default_fit to False
# when wanting to plot default fittings, fitted with ideal number of doublets, set default_fit to True
doublets = 4 #1-4 # to get f-test results 
default_fit = True
#default_fit = False


# Necessary data from Verner table to convert wavelength to velocity
wavelength_CIV_emit1 = 1548.1950
wavelength_CIV_emit2 = 1550.7700


#grabbing spectra data files
currentRP_57 = DATADIREC + 'J2318_57328_UZ1_RPnorm.dat'
currentRP2_57 = DATADIREC + 'J2318_57328_UZ1_RP2norm.dat'
currentRP_smooth_57 = DATADIREC + 'J2318_57328_DZ1_RPnorm.dat'
currentRP2_smooth_57 = DATADIREC + 'J2318_57328_DZ1_RP2norm.dat'


currentRP_59 = DATADIREC + 'J2318_59188_UZ3_RPnorm.dat'
currentRP2_59 = DATADIREC + 'J2318_59188_UZ3_RP2norm.dat'
currentRP_smooth_59 = DATADIREC + 'J2318_59188_DZ3_RPnorm.dat'
currentRP2_smooth_59 = DATADIREC + 'J2318_59188_DZ3_RP2norm.dat'


currentRP_60 = DATADIREC + 'J2318_60251_UZ2_RPnorm.dat'
currentRP2_60 = DATADIREC + 'J2318_60251_UZ2_RP2norm.dat'
currentRP_smooth_60 = DATADIREC + 'J2318_60251_DZ2_RPnorm.dat'
currentRP2_smooth_60 = DATADIREC + 'J2318_60251_DZ2_RP2norm.dat'

#choose data to use
data = [currentRP_57, currentRP_59, currentRP_60]
data = [currentRP2_57, currentRP2_59, currentRP2_60]
data = [currentRP_smooth_57, currentRP_smooth_59, currentRP_smooth_60]
data = [currentRP2_smooth_57, currentRP2_smooth_59, currentRP2_smooth_60]

#.....................................................................................................................
#.....................................................................................................................
if default_fit == False: #we don't want to save .txt files or figures for unfinalized fittings (default)
    want_save_fig = 'n'
    want_txt_save = 'n'
    
############################################################################################################################################################################

################60251 observation#########################################################################################
print('Observation 60251')
print()


if data[2] == currentRP_60:
    name1 = "J2318_60251_UZ2_RPnorm"
    smoothed_or_unsmoothed = 'not smoothed'
    power = 'RP'
elif data[2] == currentRP2_60:
    name1 = "J2318_60251_UZ2_RP2norm"
    power = 'RP2'
    smoothed_or_unsmoothed = 'not smoothed'
elif data[2] == currentRP_smooth_60:
    name1 = "J2318_60251_DZ2_RPnorm"
    smoothed_or_unsmoothed = 'smoothed'
    power = 'RP'
elif data[2] == currentRP2_smooth_60:
    name1 = "J2318_60251_DZ2_RP2norm"
    smoothed_or_unsmoothed = 'smoothed'
    power = 'RP2'

norm_spectra_data = np.loadtxt(data[2])

#reading in spectra
wavelength, norm_flux, norm_error = read_spectra(norm_spectra_data)

#converting wavelength to velocity
beta = wavelength_to_velocity(redshift, wavelength, wavelength_CIV_emit1)

if wavelength_or_velocity == 'v':
    x = beta
    x_label = 'Velocity km/s'
elif wavelength_or_velocity == 'w':
    x = wavelength
    x_label = 'Wavelength'

#...................................................
# plotting spectra

plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.axhline(y = 1.0)
plt.plot(x, norm_flux, color = 'k')
plt.xlabel(str(x_label))
plt.ylabel('Normalized Flux')
plt.xticks(np.arange(xlim[1], xlim[0] - 1, -xtick))


if title == 'yr':
    plt.title('CIV EHVO: 2023')
elif title == 'mjd':
    plt.title('CIV EHVO (mjd:60251, year: 2023) Sigma Clip')
elif title == 'paper':
    plt.title('J2318 MJD:60251 CIV EHVO')
else:
    plt.title(str(name1)+ f"  (Sigma Clip std={sigma_clip}, GK={sigma_smooth})")



#...................................................
# calculate vdiff

w1 = np.min(np.where(wavelength > 1500))
w2 = np.max(np.where(wavelength < 1490))

wavelengthC = np.array([np.mean(wavelength[w2:w1])]) 

vdiff_CIV = wavelength_to_velocity(redshift, wavelengthC,
                    wavelength_CIV_emit1) - wavelength_to_velocity(
                    redshift, wavelengthC, wavelength_CIV_emit2)
             
                     
vdiff =  vdiff_CIV #497.56973


#...................................................
# sigma clipping

FLAGGUE=(ChoiClip(wavelength, norm_flux, norm_error, sigma_smooth, sigma_clip))[0]

#......................................................
#Defining area of the spectra we will be curve fitting  
(xfit, yfit, errfit) = curve_fit_area(xstart = -72000, xend = -92500, xvalues = x[FLAGGUE==1], flux = norm_flux[FLAGGUE==1], error = norm_error[FLAGGUE==1])


plot_masked_regions(beta, norm_flux, xfit) #plotting masked regions as grey

#plt.scatter(x[FLAGGUE==0],norm_flux[FLAGGUE==0],c='orange',marker='o', s=12) #plotting masked data points
#plt.scatter(x[FLAGGUE==1],norm_flux[FLAGGUE==1],c='green',marker='o', s=5) #plotting points used in fitting

plt.plot(x, norm_error, color = 'lightgrey', label = 'Error')
plt.xlim(xlim)

# if statement for redefining number of doulbets used to fit if the ideal/default number of doublets should be used
if default_fit == True:
    doublets = 2

else:
    None

# testing one doublet 
if doublets == 1:
    
    # Define the variables
    v0 = -84435.3649881088
    
    tau0 = 0.550913156287789
    
    b = 5674.607897367083
    
    Cf = 1
    I0 = 1
    vdiff = vdiff_CIV
    tau_ratio = 2
    
    #plt.plot(x[FLAGGUE==0],norm_flux[FLAGGUE==0],c='b')#FLAGGUE==0 means flagged/clipped
    
    # implementing curvefitting function 
    v0_60, tau0_60, b_60, Cf_60, v02_60, tau02_60, b2_60, v03_60, b3_60, tau03_60, v04_60, b4_60, tau04_60  = plotting_curvefit_test(1, 'CIV', x, xfit, yfit, 'Cf', Cf, I0, vdiff, tau_ratio, v0, tau0, b)
    #v0, tau0, b, Cf, v02, tau02, b2, v03, b3, tau03, v04, b4, tau04
                                                                                                                                    #d, ion, x, xfit, yfit, fixed, Cf, I0, vdiff, tau_ratio, v0, tau0, b
                                
                                
    #saving the fitting array for the purpose of completing the f-test to evaluate the ideal number of doublets
    f_test_ydata = np.array(curve_func( xfit, tau0_60, v0_60, b_60, Cf_60, I0, vdiff, tau_ratio))
    np.savetxt(os.path.join(f_test_OUTDIREC, f'60fit_{doublets}doublets' ), f_test_ydata, comments='')

    
# fitting with 2 doublets (this is the default number of doublets for this observation)
elif doublets == 2:

    # Define the variables
    v0 = -84435.3649881088
    v02 = -75736.57431789643
    tau0 = 0.550913156287789
    tau02 = 0.08630589893577298
    b = 5674.607897367083
    b2 = 1669.790877039717
    Cf = 1
    I0 = 1
    vdiff = vdiff_CIV
    tau_ratio = 2
    
    #plt.plot(x[FLAGGUE==0],norm_flux[FLAGGUE==0],c='b')#FLAGGUE==0 means flagged/clipped
    
    # implementing curvefitting function 
    v0_60, tau0_60, b_60, Cf_60, v02_60, tau02_60, b2_60, v03_60, b3_60, tau03_60, v04_60, b4_60, tau04_60  = plotting_curvefit_test(2, 'CIV', x, xfit, yfit, 'Cf', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2)
    
    # used later to save doublet fitting parameters in a csv file
    parameters_60_1 = np.array([v0_60, tau0_60, b_60])
    parameters_60_2 = np.array([v02_60, tau02_60, b2_60])

    #saving the fitting array for the purpose of completing the f-test to evaluate the ideal number of doublets
    f_test_ydata = np.array(curve_func2( xfit, tau0_60, tau02_60, v0_60, v02_60, b_60, b2_60, Cf_60, I0, vdiff, tau_ratio))
    np.savetxt(os.path.join(f_test_OUTDIREC, f'60fit_{doublets}doublets' ), f_test_ydata, comments='')


    #saving fit as spectra to be run in absorption.py as this is the default/ideal fit 
    ydata = np.array(curve_func2( x, tau0_60, tau02_60, v0_60, v02_60, b_60, b2_60, Cf, I0, vdiff, tau_ratio))
    
    newy = np.ones_like(norm_flux)
    
    for i, x in enumerate(newy):
        if x in xfit:  
            newy[i] = ydata[i]
    
    error = norm_error
    
    output_data = np.column_stack((wavelength, ydata, error))
    
    if want_txt_save == 'y':
        np.savetxt(os.path.join(OUTDIREC, str(name1) + "_fit_spectra.txt"), output_data, header="#wavelengths flux error", comments='')
    
    

# fitting with 3 doublets
elif doublets == 3:
    
     # Define the variables
     v0 = -84435.3649881088
     v02 = -75736.57431789643
     v03 = -79000
     
     tau0 = 0.550913156287789
     tau02 = 0.08630589893577298
     tau03 = 0.2
     
     b = 5674.607897367083
     b2 = 1669.790877039717
     b3 = 2500
     
     Cf = 1
     I0 = 1
     vdiff = vdiff_CIV
     tau_ratio = 2
     
     #plt.plot(x[FLAGGUE==0],norm_flux[FLAGGUE==0],c='b')#FLAGGUE==0 means flagged/clipped
     
     # implementing curvefitting function
     v0_60, tau0_60, b_60, Cf_60, v02_60, tau02_60, b2_60, v03_60, b3_60, tau03_60, v04_60, b4_60, tau04_60  = plotting_curvefit_test(3, 'CIV', x, xfit, yfit, 'Cf', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2, v03, tau03, b3)
    
    #saving the fitting array for the purpose of completing the f-test to evaluate the ideal number of doublets
     f_test_ydata = np.array(curve_func3( xfit, tau0_60, tau02_60, tau03_60, v0_60, v02_60, v03_60, b_60, b2_60, b3_60, Cf_60, I0, vdiff, tau_ratio))
     np.savetxt(os.path.join(f_test_OUTDIREC, f'60fit_{doublets}doublets' ), f_test_ydata, comments='')

# we are not testing a fit of 4 doublets as 2 is ideal for this observation
elif doublets == 4:
    print('WARNING: Not fitting 4 doublets. Will not perform F-test with 4 doublet fitting.')



# final settings on plot and saving image
plt.legend(loc='lower right')

if want_save_fig == 'y':
    plt.savefig(os.path.join(OUTDIREC, f'{name1}_plottedfit.png'))
plt.show()
plt.close()



# F-test ------------------
#compare 1 doublets to 2 doublets
try:
    fitting_old = np.loadtxt(f_test_OUTDIREC + '/60fit_1doublets')
    fitting_new = np.loadtxt(f_test_OUTDIREC + '/60fit_2doublets')
    
    param_old = 3
    param_new = 6


    f_test_1_2, v, f_crit = f_test(yfit, param_new, fitting_new, param_old, fitting_old, 0.001)

    print()
    print(f'The resulting value of the f-test between 1 and 2 doublet fits is: F = {f_test_1_2} & DOF = {v} : {f_crit}')

    
except:
    print()
    print('Fittings for obs 60251 with 1 or 2 doublets have not been run!')


#compare 2 doublets to 3 doublets
try:
    fitting_old = np.loadtxt(f_test_OUTDIREC + '/60fit_2doublets')
    fitting_new = np.loadtxt(f_test_OUTDIREC + '/60fit_3doublets')

    param_old = 6
    param_new = 9

    
    f_test_2_3, v, f_crit = f_test(yfit, param_new, fitting_new, param_old, fitting_old, 0.001)

    print()
    print(f'The resulting value of the f-test between 2 and 3 doublet fits is: F = {f_test_2_3} & DOF = {v} : {f_crit}')


except:
    print()
    print('Fittings for obs 59188 with 2 or 3 doublets have not been run!')
    

############################################################################################################################################################################

################59188 observation#########################################################################################
print()
print()
print('Observation 59188')
print()

# assign string variables for file naming
if data[1] == currentRP_59:
    name3 = "J2318_59188_UZ3_RPnorm"
    smoothed_or_unsmoothed = 'not smoothed'
    power = 'RP'
elif data[1] == currentRP2_59:
    name3 = "J2318_59188_UZ3_RP2norm"
    power = 'RP2'
    smoothed_or_unsmoothed = 'not smoothed'
elif data[1] == currentRP_smooth_59:
    name3 = "J2318_59188_DZ3_RPnorm"
    smoothed_or_unsmoothed = 'smoothed'
    power = 'RP'
elif data[1] == currentRP2_smooth_59:
    name3 = "J2318_59188_DZ3_RP2norm"
    smoothed_or_unsmoothed = 'smoothed'
    power = 'RP2'

# loading data chosen in changeable variables
norm_spectra_data = np.loadtxt(data[1])

#reading in spectra
wavelength, norm_flux, norm_error = read_spectra(norm_spectra_data)

# convert wavelength to velocity
beta = wavelength_to_velocity(redshift, wavelength, wavelength_CIV_emit1)

# use wavelength or velocity depending on choice in changeable variables
if wavelength_or_velocity == 'v':
    x = beta
    x_label = 'Velocity km/s'
elif wavelength_or_velocity == 'w':
    x = wavelength
    x_label = 'Wavelength'

# establishing plot parameters
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.axhline(y = 1.0)
plt.plot(x, norm_flux, color = 'k')
plt.xlabel(str(x_label))
plt.ylabel('Normalized Flux')
plt.xticks(np.arange(xlim[1], xlim[0] - 1, -xtick))

if title == 'yr':
    plt.title('CIV EHVO: 2020')
elif title == 'mjd':
    plt.title('CIV EHVO (mjd:59188, year:2020) Sigma Clip')
elif title == 'paper':
    plt.title('J2318 MJD:59188 CIV EHVO')
else:
    plt.title(str(name3)+ f"  (Sigma Clip std={sigma_clip}, GK={sigma_smooth})")



#...................................................
# calculate vdiff

w1 = np.min(np.where(wavelength > 1500))
w2 = np.max(np.where(wavelength < 1490))

wavelengthC = np.array([np.mean(wavelength[w2:w1])]) #made np.array bc my version of w to v itterates over a list of values - LEF

vdiff_CIV = wavelength_to_velocity(redshift, wavelengthC,
                    wavelength_CIV_emit1) - wavelength_to_velocity(
                    redshift, wavelengthC, wavelength_CIV_emit2)
             
                     
vdiff =  vdiff_CIV #497.56973

#...................................................
# sigma clipping

FLAGGUE=(ChoiClip(wavelength, norm_flux, norm_error, sigma_smooth, sigma_clip))[0]

#......................................................
#Defining area of the spectra we will be curve fitting  
(xfit, yfit, errfit) = curve_fit_area(xstart = -78000, xend = -97500, xvalues = x[FLAGGUE==1], flux = norm_flux[FLAGGUE==1], error = norm_error)

# manually adding back in points
if readd == 'yes':
    left = -87200  # left limit of range to add
    right = -86910  # right limit of range to add
    oldxfit = xfit
    oldyfit = yfit
    olderrfit = errfit
    add_indices = np.where((x >= left) & (x <= right))[0]
    newx = x[add_indices]
    newy = norm_flux[add_indices]
    newerr = norm_error[add_indices]
    xfit = np.hstack((oldxfit, newx))
    yfit = np.hstack((oldyfit, newy))
    errfit = np.hstack((olderrfit, newerr))
    sort_indices = np.argsort(xfit)
    xfit = xfit[sort_indices]
    yfit = yfit[sort_indices]
    errfit = errfit[sort_indices]
      
if readd == 'yes':
    left = -85920  # left limit of range to add
    right = -85600  # right limit of range to add
    oldxfit = xfit
    oldyfit = yfit
    olderrfit = errfit
    add_indices = np.where((x >= left) & (x <= right))[0]
    newx = x[add_indices]
    newy = norm_flux[add_indices]
    newerr = norm_error[add_indices]
    xfit = np.hstack((oldxfit, newx))
    yfit = np.hstack((oldyfit, newy))
    errfit = np.hstack((olderrfit, newerr))
    sort_indices = np.argsort(xfit)
    xfit = xfit[sort_indices]
    yfit = yfit[sort_indices]
    errfit = errfit[sort_indices]
    
# plotting the masked regions in grey
plot_masked_regions(beta, norm_flux, xfit)

# plotting error and establishing xlimits
plt.plot(x, norm_error, color = 'lightgrey', label = 'Error')
plt.xlim(xlim)

# redefining the number of doublets in the default fit is selected
if default_fit == True:
    doublets = 3
else:
    None

# fitting with 1 doublet
if doublets == 1:
    
    # Define the variables
    v0 = -90000
    tau0 = 0.2
    b = 10000

    Cf = 1
    I0 = 1
    vdiff = vdiff_CIV
    tau_ratio = 2

    # implementing curvefitting function
    v0_59, tau0_59, b_59, Cf_59, v02_59, tau02_59, b2_59, v03_59, b3_59, tau03_59, v04_59, b4_59, tau04_59 = plotting_curvefit_test(1, 'CIV', x, xfit, yfit, 'Cf', Cf, I0, vdiff, tau_ratio, v0, tau0, b)

    #saving the fitting array for the purpose of completing the f-test to evaluate the ideal number of doublets
    f_test_ydata = np.array(curve_func( xfit, tau0_59, v0_59, b_59, Cf_59, I0, vdiff, tau_ratio))
    np.savetxt(os.path.join(f_test_OUTDIREC, f'59fit_{doublets}doublets' ), f_test_ydata, comments='')

# fitting with 2 doublets
elif doublets == 2:
    
    # Define the variables
    v0 = -89000
    v02 = -83000

    tau0 = 0.14
    tau02 = 0.55
   
    b = 5300
    b2 = 1600
   

    Cf = 1
    I0 = 1
    vdiff = vdiff_CIV
    tau_ratio = 2

    # implementing curve fitting function
    v0_59, tau0_59, b_59, Cf_59, v02_59, tau02_59, b2_59, v03_59, b3_59, tau03_59, v04_59, b4_59, tau04_59 = plotting_curvefit_test(2, 'CIV', x, xfit, yfit, 'Cf', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2)

    #saving the fitting array for the purpose of completing the f-test to evaluate the ideal number of doublets
    f_test_ydata = np.array(curve_func2( xfit, tau0_59, tau02_59, v0_59, v02_59, b_59, b2_59, Cf_59, I0, vdiff, tau_ratio))
    np.savetxt(os.path.join(f_test_OUTDIREC, f'59fit_{doublets}doublets' ), f_test_ydata, comments='')

# fitting with 3 doublets (this is the default/ideal fit)
elif doublets == 3:
    # Define the variables
    v0 = -94000
    v02 = -87000
    v03 = -83000
    
    tau0 = 0.2
    tau02 = 0.5
    tau03 = 1
    
    b = 900
    b2 = 2000
    b3 = 1500
    
    
    Cf = 1
    I0 = 1
    vdiff = vdiff_CIV
    tau_ratio = 2
    
    
    # implementing curve fitting function    
    v0_59, tau0_59, b_59, Cf_59, v02_59, tau02_59, b2_59, v03_59, b3_59, tau03_59, v04_59, b4_59, tau04_59 = plotting_curvefit_test(3, 'CIV', x, xfit, yfit, 'Cf', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2, v03, tau03, b3)
    
    #saving the fitting array for the purpose of completing the f-test to evaluate the ideal number of doublets    
    f_test_ydata = np.array(curve_func3( xfit, tau0_59, tau02_59, tau03_59, v0_59, v02_59, v03_59, b_59, b2_59, b3_59, Cf_59, I0, vdiff, tau_ratio))
    np.savetxt(os.path.join(f_test_OUTDIREC, f'59fit_{doublets}doublets' ), f_test_ydata, comments='')
    
    # saving doublet parameters for csv file-----------------
    parameters_59_3 = np.array([v0_59, tau0_59, b_59])
    parameters_59_2 = np.array([v02_59, tau02_59, b2_59])
    parameters_59_1 = np.array([v03_59, tau03_59, b3_59])
                                  
    # saving default/ideal fitting to .txt file to run in absorption.py
    ydata = np.array(curve_func3( x, tau0_59, tau02_59, tau03_59, v0_59, v02_59, v03_59, b_59, b2_59, b3_59, Cf_59, I0, vdiff, tau_ratio))

    newy = np.ones_like(norm_flux)

    for i, x in enumerate(newy):
        if x in xfit:  
            newy[i] = ydata[i]
            
    error = norm_error

    output_data = np.column_stack((wavelength, ydata, error))

    if want_txt_save == 'y':
        np.savetxt(os.path.join(OUTDIREC, str(name3) + "_fit_spectra.txt"), output_data, header="#wavelengths flux error", comments='')


# fitting with 4 doublets
elif doublets == 4:
    # Define the variables
    v0 = -94000
    v02 = -87000
    v03 = -83000
    v04 = -86000

    
    tau0 = 0.2
    tau02 = 0.5
    tau03 = 1
    tau04 = 0.7

    
    b = 900
    b2 = 2000
    b3 = 1500
    b4 = 500
    
    
    Cf = 1
    I0 = 1
    vdiff = vdiff_CIV
    tau_ratio = 2
    
    # implementing curvefitting function    
    v0_59, tau0_59, b_59, Cf_59, v02_59, tau02_59, b2_59, v03_59, b3_59, tau03_59, v04_59, b4_59, tau04_59 = plotting_curvefit_test(4, 'CIV', x, xfit, yfit, 'Cf', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2, v03, tau03, b3, v04, tau04, b4)
    
    #saving the fitting array for the purpose of completing the f-test to evaluate the ideal number of doublets
    f_test_ydata = np.array(curve_func4( xfit, tau0_59, tau02_59, tau03_59, tau04_59, v0_59, v02_59, v03_59, v04_59, b_59, b2_59, b3_59, b4_59, Cf_59, I0, vdiff, tau_ratio))
    np.savetxt(os.path.join(f_test_OUTDIREC, f'59fit_{doublets}doublets' ), f_test_ydata, comments='')
    
    
    
# last settings for plot and saving figure
plt.legend(loc='lower right')

if want_save_fig == 'y':
    plt.savefig(os.path.join(OUTDIREC, f'{name3}_plottedfit.png'))
    
plt.show()
plt.close()


#F-test ------------------
# compare 1 doublet to 2 doublets
try:
    fitting_old = np.loadtxt(f_test_OUTDIREC + '/59fit_1doublets')
    fitting_new = np.loadtxt(f_test_OUTDIREC + '/59fit_2doublets')
    
    param_old = 3
    param_new = 6
        
    f_test_1_2, v, f_crit = f_test(yfit, param_new, fitting_new, param_old, fitting_old, 0.001)

    print()
    print(f'The resulting value of the f-test between 1 and 2 doublet fits is: F = {f_test_1_2} & DOF = {v} : {f_crit}')
 
except:
    print()
    print('Fittings for obs 59188 with 1 or 2 doublets have not been run!')
    

# compare 2 doublets to 3 doublets
try:
    fitting_old = np.loadtxt(f_test_OUTDIREC + '/59fit_2doublets')
    fitting_new = np.loadtxt(f_test_OUTDIREC + '/59fit_3doublets')
    
    param_old = 6
    param_new = 9


    f_test_2_3, v, f_crit = f_test(yfit, param_new, fitting_new, param_old, fitting_old, 0.001)

    print()
    print(f'The resulting value of the f-test between 2 and 3 doublet fits is: F = {f_test_2_3} & DOF = {v} : {f_crit}')

except:
    print()
    print('Fittings for obs 59188 with 2 or 3 doublets have not been run!')
    

    
# compare 3 doublets to 4 doublets 
try:
    fitting_old = np.loadtxt(f_test_OUTDIREC + '/59fit_3doublets')
    fitting_new = np.loadtxt(f_test_OUTDIREC + '/59fit_4doublets')
        
    param_old = 9
    param_new = 12


    f_test_3_4, v, f_crit = f_test(yfit, param_new, fitting_new, param_old, fitting_old, 0.001)

    print()
    print(f'The resulting value of the f-test between 3 and 4 doublet fits is: F = {f_test_3_4} & DOF = {v} : {f_crit}')
    
except:
    print()
    print('Fittings for obs 59188 with 3 or 4 doublets have not been run!')
    



############################################################################################################################################################################

######## 57328 Observation ######################################################################################
print()
print()
print('Observation 57328')
print()

# assign string variables for file naming
if data[0] == currentRP_57:
    name2 = "J2318_57328_UZ1_RPnorm"
    smoothed_or_unsmoothed = 'not smoothed'
    power = 'RP'
elif data[0] == currentRP2_57:
    name2 = "J2318_57328_UZ1_RP2norm"
    power = 'RP2'
    smoothed_or_unsmoothed = 'not smoothed'
elif data[0] == currentRP_smooth_57:
    name2 = "J2318_57328_DZ1_RPnorm"
    smoothed_or_unsmoothed = 'smoothed'
    power = 'RP'
elif data[0] == currentRP2_smooth_57:
    name2 = "J2318_57328_DZ1_RP2norm"
    smoothed_or_unsmoothed = 'smoothed'
    power = 'RP2'

# loading data chosen in changeable variables
norm_spectra_data = np.loadtxt(data[0])


#reading in spectra
wavelength, norm_flux, norm_error = read_spectra(norm_spectra_data)

# converting wavelength to velocity
beta = wavelength_to_velocity(redshift, wavelength, wavelength_CIV_emit1)

# plotting with wavelength or velocity using selection in chageable variables
if wavelength_or_velocity == 'v':
    x = beta
    x_label = 'Velocity km/s'
elif wavelength_or_velocity == 'w':
    x = wavelength
    x_label = 'Wavelength'

# establishing plot parameters
plt.ylim(y_min, y_max)
plt.xlim(x_min, x_max)
plt.axhline(y = 1.0)
plt.plot(x, norm_flux, color = 'k')
plt.xlabel(str(x_label))
plt.ylabel('Normalized Flux')
plt.xticks(np.arange(xlim[1], xlim[0] - 1, -xtick))
plt.xlim(xlim)


if title == 'yr':
    plt.title('CIV EHVO: 2015')
elif title == 'mjd':
    plt.title('CIV EHVO (mjd:57328, year:2015) Sigma Clip')
elif title == 'paper':
    plt.title('J2318 MJD:57328 CIV EHVO')
else:
    plt.title(str(name2)+ f"  (Sigma Clip std={sigma_clip}, GK={sigma_smooth})")


#...................................................
# sigma clipping

FLAGGUE=(ChoiClip(wavelength, norm_flux, norm_error, sigma_smooth, sigma_clip))[0]

#......................................................
#Defining area of the spectra we will be curve fitting  
(xfit, yfit, errfit) = curve_fit_area(xstart = -78600, xend = -92500, xvalues = x[FLAGGUE==1], flux = norm_flux[FLAGGUE==1], error = norm_error)

# manual masking regions
if manual_mask == True:
    left = -84300
    right = -83000#-87200
    oldxfit = xfit
    oldyfit = yfit
    xfit1 = xfit2 = 0
    xfit1=xfit[0:np.max(np.where(xfit < left))]
    xfit2=xfit[np.min(np.where(xfit > right)):]
    xfit=np.hstack((xfit1,xfit2))
    yfit1=yfit[0:np.max(np.where(oldxfit < left))]
    yfit2=yfit[np.min(np.where(oldxfit > right)):]
    yfit=np.hstack((yfit1,yfit2))
    errfit1=errfit[0:np.max(np.where(oldxfit < left))]
    errfit2=errfit[np.min(np.where(oldxfit > right)):]
    errfit=np.hstack((errfit1,errfit2))

# plotting masked regions in grey
plot_masked_regions(beta, norm_flux, xfit)

# plotting error
plt.plot(x, norm_error, color = 'lightgrey', label = 'Error')


# redefining the number of doublets if the default fit is selected in the changeable variables section
if default_fit == True:
    doublets = 2

else:
    None  
    
    
# fitting with 1 doublet
if doublets == 1:
    
    
    tau0 = 0.187
    v0 = -86191
    b = 5888
    #b = 2000
    
    Cf = 1
    I0 = 1
    vdiff = vdiff_CIV
    tau_ratio = 2
    
    # implementing curve fitting function
    v0_57, tau0_57, b_57, Cf_57, v02_57, tau02_57, b2_57, v03_57, b3_57, tau03_57, v04_57, b4_57, tau04_57 = plotting_curvefit_test(1, 'CIV', x, xfit, yfit, 'Cf', Cf, I0, vdiff, tau_ratio, v0, tau0, b)
    
    #saving the fitting array for the purpose of completing the f-test to evaluate the ideal number of doublets
    f_test_ydata = np.array(curve_func( xfit, tau0_57, v0_57, b_57, Cf_57, I0, vdiff, tau_ratio))
    np.savetxt(os.path.join(f_test_OUTDIREC, f'57fit_{doublets}doublets' ), f_test_ydata, comments='')

    
# fitting with 2 doublets (this is the default/ideal fit)
elif doublets == 2:

    tau0 = 0.12
    v0 = -90000
    b = 5000


    tau0 = 0.5
    v02 = -83000
    b = 3000


    # implementing curve fitting function
    v0_57, tau0_57, b_57, Cf_57, v02_57, tau02_57, b2_57,  v03_57, b3_57, tau03_57, v04_57, b4_57, tau04_57 = plotting_curvefit_test(2, 'CIV', x, xfit, yfit, 'Cf', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2)

    
    #saving the fitting array for the purpose of completing the f-test to evaluate the ideal number of doublets
    f_test_ydata = np.array(curve_func2( xfit, tau0_57, tau02_57, v0_57, v02_57, b_57, b2_57, Cf_57, I0, vdiff, tau_ratio))
    np.savetxt(os.path.join(f_test_OUTDIREC, f'57fit_{doublets}doublets' ), f_test_ydata, comments='')
    
    # saving doublet parameters to be saved into a csv file
    parameters_57_1 = np.array([v0_57, tau0_57, b_57])
    parameters_57_2 = np.array([v02_57, tau02_57, b2_57])

  
    # saving fitting to .txt file to be run through absorption.py
    ydata = np.array(curve_func2( x, tau0_57, tau02_57, v0_57, v02_57, b_57, b2_57, Cf_57, I0, vdiff, tau_ratio))

    newy = np.ones_like(norm_flux)

    for i, x in enumerate(newy):
        if x in xfit:  
            newy[i] = ydata[i]

    error = norm_error

    output_data = np.column_stack((wavelength, ydata, error))

    if want_txt_save == 'y':
        np.savetxt(os.path.join(OUTDIREC, str(name2) + "_fit_spectra.txt"), output_data, header="#wavelengths flux error", comments='')

# fitting with 3 doublets   
elif doublets == 3:
    
    print(xfit)

    v0 = -90000
    v02 = -81000
    v03 = -80000
    
    tau0 = 0.550913156287789
    tau02 = 0.08630589893577298
    tau03 = 0.05
    
    b = 5674.607897367083
    b2 = 1669.790877039717
    b3 = 1000
    
    Cf = 1
    I0 = 1
    vdiff = vdiff_CIV
    tau_ratio = 2
    
    #plt.plot(x[FLAGGUE==0],norm_flux[FLAGGUE==0],c='b')#FLAGGUE==0 means flagged/clipped
    
    # implementing curve fitting function
    v0_57, tau0_57, b_57, Cf_57, v02_57, tau02_57, b2_57, v03_57, b3_57, tau03_57, v04_57, b4_57, tau04_57  = plotting_curvefit_test(3, 'CIV', x, xfit, yfit, 'Cf', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2, v03, tau03, b3)
    
    #saving the fitting array for the purpose of completing the f-test to evaluate the ideal number of doublets
    f_test_ydata = np.array(curve_func3( xfit, tau0_57, tau02_57, tau03_57, v0_57, v02_57, v03_57, b_57, b2_57, b3_57, Cf_57, I0, vdiff, tau_ratio))
    np.savetxt(os.path.join(f_test_OUTDIREC, f'57fit_{doublets}doublets' ), f_test_ydata, comments='')

# Not fitting with 4 doublets as 2 doublets was determined ideal with the f-test
elif doublets == 4:
    print('WARNING: Not fitting 4 doublets. Will not perform F-test with 4 doublet fitting.')

# establishing last plot settings and saving the figures
plt.legend(loc='lower right')
if want_save_fig == 'y':
    plt.savefig(os.path.join(OUTDIREC, f'{name2}_plottedfit.png'))
plt.show()
plt.close()

##################

# F-test ----------------------------------------------------------------------
#compare 1 doublet to 2 doublets
try:
    fitting_old = np.loadtxt(f_test_OUTDIREC + '/57fit_1doublets')
    fitting_new = np.loadtxt(f_test_OUTDIREC + '/57fit_2doublets')
    
    param_old = 3
    param_new = 6
    
    f_test_1_2, v, f_crit = f_test(yfit, param_new, fitting_new, param_old, fitting_old, 0.001)

    print()
    print(f'The resulting value of the f-test between 1 and 2 doublet fits is: F = {f_test_1_2} & DOF = {v} : {f_crit}')
except:
    print()
    print('Fittings for obs 59188 with 1 or 2 doublets have not been run!')

# compare 2 doublets to 3 doublets
try:
    #compare 2 doublets to 3 doublets
    fitting_old = np.loadtxt(f_test_OUTDIREC + '/57fit_2doublets')
    fitting_new = np.loadtxt(f_test_OUTDIREC + '/57fit_3doublets')
    
    param_old = 6
    param_new = 9
    
    
    f_test_2_3, v, f_crit = f_test(yfit, param_new, fitting_new, param_old, fitting_old, 0.001)

    print()
    print(f'The resulting value of the f-test between 2 and 3 doublet fits is: F = {f_test_2_3} & DOF = {v} : {f_crit}')

except:
    print()
    print('Fittings for obs 59188 with 2 or 3 doublets have not been run!')



#----------------------------------------------------------------------------------------------
#################################################################################################################################################################################################
 ###### export csv of doublet parameters to be used in SiIV fitting
try:
    parameters_dic = {'60_1': parameters_60_1,
                       '60_2': parameters_60_2,
                       '59_1': parameters_59_1,
                       '59_2': parameters_59_2,
                       '59_3': parameters_59_3,
                       '57_1': parameters_57_1,
                       '57_2': parameters_57_2}
     
    df_params = pd.DataFrame(parameters_dic, index=['v0', 'tau0', 'b'])
    print(OUTDIREC)
    df_params.to_csv(os.path.join(OUTDIREC, f'doublet_params_{power}_{smoothed_or_unsmoothed}.csv'), index=True)
    
except:
    print('Run "default_fit = True" in changable variables to save paramters of true fittings.')




