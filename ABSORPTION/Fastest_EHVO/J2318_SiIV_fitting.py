#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 17:16:10 2024

@author: lilianaflores


Key Features:
    - Fits SiIV EHVO absorption troughs for all three observations of J2318.
    - Applies a sigma clipping function and additional spectral masking to 
      improve fit accuracy.
    - Generates and displays plots for each observation.
    - Allows toggling between multiple spectral data types, including
      smoothed, unsmoothed, and alternative power-law.
    - Performs SiIV fits on error-shifted spectra that will be used to constrain 
      absorption measurement uncertainties.

Inputs:
    - Twelve .dat files enabling selection among four spectral variants
      for each of the three observations.
    - CSV files containing CIV fitting parameters used as initial guesses
      and fixed values when fitting the SiIV EHVOs.

Outputs:
    - Fitting plots for each observation, saved to a specified output directory.
    - Optional plots generated for constraining absorption measurements
      by ±0.5σ spectral shifts. Corresponding .txt files of the shifted
      spectral fits are saved to the output directory.
    - .txt files containing fitted spectra (wavelength, flux, error),
      formatted for direct input into absorption.py for to get measuremnets.


"""

# imports ................................................................................................
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
from Fastest_CIV_EHVO_Functions import plot_masked_regions 
from curve_fit_function import plotting_curvefit_test, curve_fit_area, tau_v

# curve function for exporting fitting data
def curve_func( v, tau0, tau02, tau03, v0, v02, v03, b, b2, b3, Cf, I0, vdiff, tau_ratio): 
    return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)) * np.exp(-tau_v(v,v02,b2,tau02)) * np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)) * np.exp(-tau_v(v,v03,b3,tau03)) * np.exp(-tau_v(v,(v03+vdiff),b3,tau03/tau_ratio)))) 


def curve_func2( v, tau0, tau02, v0, v02, b, b2, Cf, I0, vdiff, tau_ratio): 
    return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)) * np.exp(-tau_v(v,v02,b2,tau02)) * np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)))) 
#.....................................................................................................................
###########################################################################################################################################################################################################
#CHANGEABLE VARIABLES
###########################################################################################################################################################################################################

# establish directory pathways
OUTDIREC = os.getcwd() + '/output_files/'

DATADIREC = os.getcwd() + '/data/'


# save figures and .txt files of spectra
want_save_fig = 'y'
want_txt_save = 'y'

#want_save_fig = 'n'
#want_txt_save = 'n'


# manual masking some points in obs
manual_mask = True
#manual_mask = False

# run and plot in wavelength or velocity
wavelength_or_velocity = 'v'
#wavelength_or_velocity = 'W'

# redshift set to zero as spectra are provided in restframe
redshift = 0

# different plot tile choices for figures
title = 'yr' #SiIV EHVO: year
title = 'mjd' #SiIV EHVO (mjd:#, year:# ) Sigma Clip
title = 'org' #J2318_obs#_DZ/UZ_RP/RP2norm (Sigma Clip std=X, GK =X)
title = 'paper' #J2318 MJD:# SiIV EHVO

# sigma clipping function settings
sigma_smooth = 3
sigma_clip = 0.85

# fit with free b, not a fixed value from CIV EHVO fitting.
free_b = False
#free_b = True


# number of doublets to fit obs 59188
D59=2 # fitting with 2 doublets as when fitting with 3 there is no absorption fitted for the first doublet
#D59=3

# fit settings for obs 57328
# find an upper limit fit by shifting the spectra down by subtracting by 0.5*error
attempt57_fit = 'og'
#attempt57_fit = 'shifted'

# additional spectra fittings for shifted spectra to constrain errors for absorption measurements
#shifted_spec_fit = True
shifted_spec_fit = False

# establishing set plot limits
x_min = -130000
x_max = -50000
y_min = 0#0.2
y_max = 1.3#1.1
xlim =(-105000,-60000)
xtick = 10000

# V3 spectra data used
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

# necessary data from Verner table
wavelength_SiIV_emit1 = 1393.755 # we use this wavelength
wavelength_SiIV_emit2 = 1402.770

#selecting data files to use
data = [currentRP_57, currentRP_59, currentRP_60]
#data = [currentRP2_57, currentRP2_59, currentRP2_60]
#data = [currentRP_smooth_57, currentRP_smooth_59, currentRP_smooth_60]
#data = [currentRP2_smooth_57, currentRP2_smooth_59, currentRP2_smooth_60]

###########################################################################################################################################################################################################
###########################################################################################################################################################################################################


################## 60251 ###################################################################################################################################################################################################################################################################
###########################################################################################################################################################################################################
print()
print()
print('**** FITTING MJD 60251 ****')
print()

# conditional statements based on data choice. Reading in doublet parameters for CIV fittings to be used for SiIV fittings
if data[2] == currentRP_60:
    name2 = "J2318_60251_UZ2_RPnorm"
    smoothed_or_unsmoothed = 'not smoothed'
    power = 'RP'
    
    parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP_not smoothed.csv'))
    df_parameters = parameters.set_index(parameters.columns[0])
    
    tau0 = df_parameters['60_1']['tau0']
    v0 = df_parameters['60_1']['v0']
    b = df_parameters['60_1']['b']

    tau02 = df_parameters['60_2']['tau0']
    v02 = df_parameters['60_2']['v0']
    b2 = df_parameters['60_2']['b']
    
elif data[2] == currentRP2_60:
    name2 = "J2318_60251_UZ2_RP2norm"
    power = 'RP2'
    smoothed_or_unsmoothed = 'not smoothed' 
    
    parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP2_not smoothed.csv'))
    df_parameters = parameters.set_index(parameters.columns[0])
    
    tau0 = df_parameters['60_1']['tau0']
    v0 = df_parameters['60_1']['v0']
    b = df_parameters['60_1']['b']

    tau02 = df_parameters['60_2']['tau0']
    v02 = df_parameters['60_2']['v0']
    b2 = df_parameters['60_2']['b']
    
elif data[2] == currentRP_smooth_60:
    name2 = "J2318_60251_DZ2_RPnorm"
    smoothed_or_unsmoothed = 'smoothed'
    power = 'RP'
    
    parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP_smoothed.csv'))
    df_parameters = parameters.set_index(parameters.columns[0])
    
    tau0 = df_parameters['60_1']['tau0']
    v0 = df_parameters['60_1']['v0']
    b = df_parameters['60_1']['b']

    tau02 = df_parameters['60_2']['tau0']
    v02 = df_parameters['60_2']['v0']
    b2 = df_parameters['60_2']['b']
    
elif data[2] == currentRP2_smooth_60:
    name2 = "J2318_60251_DZ2_RP2norm"
    smoothed_or_unsmoothed = 'smoothed'
    power = 'RP2'
    
    parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP2_smoothed.csv'))
    df_parameters = parameters.set_index(parameters.columns[0])
    
    tau0 = df_parameters['60_1']['tau0']
    v0 = df_parameters['60_1']['v0']
    b = df_parameters['60_1']['b']

    tau02 = df_parameters['60_2']['tau0']
    v02 = df_parameters['60_2']['v0']
    b2 = df_parameters['60_2']['b']
    
#reading in spectra
norm_spectra_data = np.loadtxt(data[2])

wavelength, norm_flux, norm_error = read_spectra(norm_spectra_data)

# converting wavelength to velocity
beta = wavelength_to_velocity(redshift, wavelength, wavelength_SiIV_emit1) 


# selecting to plot in wavelength or velocity based on changeable variable section selection
if wavelength_or_velocity == 'v':
    x = beta
    x_label = 'Velocity km/s'
elif wavelength_or_velocity == 'w':
    x = wavelength
    x_label = 'Wavelength'


# calculate vdiff............................................................................
w1 = np.min(np.where(wavelength > 1500))
w2 = np.max(np.where(wavelength < 1490))

wave = np.array([np.mean(wavelength[w2:w1])]) 


vdiff_SiIV = wavelength_to_velocity(redshift, wave,
                    wavelength_SiIV_emit1) - wavelength_to_velocity(
                    redshift, wave, wavelength_SiIV_emit2)
             
                     
vdiff =  vdiff_SiIV 
 #................................................................................................

# plot settings
plt.xlim(xlim)
plt.xticks(np.arange(xlim[1], xlim[0] - 1, -xtick))
plt.ylim(y_min, y_max)
plt.axhline(y = 1.0)
plt.plot(x, norm_flux, color = 'k')   
plt.xlabel(str(x_label))
plt.ylabel('Normalized Flux')


# establishing plot title with coice made in changeable variable section ......................................................
if title == 'yr':
    plt.title('SiIV: 2023')
elif title == 'mjd':
    plt.title('SiIV  (mjd:60251, year: 2023) Sigma Clip')
elif title == 'paper':
    plt.title('J2318 MJD:60251 SiIV EHVO')
else:
    plt.title(str(name2)+ f"  (Sigma Clip std={sigma_clip}, GK={sigma_smooth})")
#................................................................................................



#...................................................
#...................................................
# sigma clipping

start=np.ones(len(wavelength))
totale=len(np.where(wavelength<1215.7)[0])

FLAGGUE=(ChoiClip(wavelength, norm_flux, norm_error, sigma_smooth, sigma_clip))[0]
    
#......................................................
#Defining area of the spectra we will be curve fitting  
(xfit, yfit, errfit) = curve_fit_area(xstart = -71000, xend = -94000, xvalues = x[FLAGGUE==1], flux = norm_flux[FLAGGUE==1], error = norm_error)


# additional fitting data adjustments
if manual_mask == True:
    
    left = -87700
    right = -86500#-87200
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
    

    left = -84500
    right = -84200
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


    left = -83750#-83700
    right = -82950#-83100
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
    

    left = -82300
    right = -81700
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
    

    left = -81400
    right = -81000
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


# establishing curvefitting parameters
Cf = 1
I0 = 1
vdiff = vdiff_SiIV
tau_ratio = 2

# implementing curvefitting function
if free_b == True:
    b2 =800
    b=500
    v0_60, tau0_60, b_60, Cf, v02_60, tau02_60, b2_60, v03, b3, tau03, v04, b4, tau04 = plotting_curvefit_test(2, 'SiIV', x, xfit, yfit, 'v0, Cf', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2)
else:
    v0_60, tau0_60, b_60, Cf, v02_60, tau02_60, b2_60, v03, b3, tau03, v04, b4, tau04 = plotting_curvefit_test(2, 'SiIV', x, xfit, yfit, 'Cf, v0, b', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2)

# plot settings
plt.plot(x, norm_error, color = 'lightgrey', label = 'Error')
plt.legend(loc='lower right')
plt.savefig(os.path.join(OUTDIREC, f'SiIV_60251_{smoothed_or_unsmoothed}_{power}.png'))
plt.figure()
plt.show()
plt.close()



ydata = np.array(curve_func2( x, tau0_60, tau02_60, v0_60, v02_60, b_60, b2_60, Cf, I0, vdiff, tau_ratio))

newy = np.ones_like(norm_flux)

for i, x in enumerate(newy):
    if x in xfit:  
        newy[i] = ydata[i]

#error = np.zeros_like(wavelength)  
error = norm_error

output_data = np.column_stack((wavelength, ydata, error))

if want_txt_save == 'y':
    np.savetxt(os.path.join(OUTDIREC, str(name2) + "_SiIVfit.txt"), output_data, header="#wavelengths flux error", comments='')
###########################################################################################################################################################################################################

################## 59188 ###################################################################################################################################################################################################################################################################
print()
print()
print('**** FITTING MJD 59188 ****')
print()
# conditional statements based on data choice. Reading in doublet parameters for CIV fittings to be used for SiIV fittings
if data[1] == currentRP_59:
    name1 = "J2318_59188_UZ3_RPnorm"
    smoothed_or_unsmoothed = 'not smoothed'
    power = 'RP'
    
    parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP_not smoothed.csv'))
    df_parameters = parameters.set_index(parameters.columns[0])
    
    tau0 = df_parameters['59_1']['tau0']
    v0 = df_parameters['59_1']['v0']
    b = df_parameters['59_1']['b']

    tau02 = df_parameters['59_2']['tau0']
    v02 = df_parameters['59_2']['v0']
    b2 = df_parameters['59_2']['b']
    
    tau03 = df_parameters['59_3']['tau0']
    v03 = df_parameters['59_3']['v0']
    b3 = df_parameters['59_3']['b']
    
elif data[1] == currentRP2_59:
    name1 = "J2318_59188_UZ3_RP2norm"
    power = 'RP2'
    smoothed_or_unsmoothed = 'not smoothed'
    
    parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP2_not smoothed.csv'))
    df_parameters = parameters.set_index(parameters.columns[0])
    
    tau0 = df_parameters['59_1']['tau0']
    v0 = df_parameters['59_1']['v0']
    b = df_parameters['59_1']['b']

    tau02 = df_parameters['59_2']['tau0']
    v02 = df_parameters['59_2']['v0']
    b2 = df_parameters['59_2']['b']
    
    tau03 = df_parameters['59_3']['tau0']
    v03 = df_parameters['59_3']['v0']
    b3 = df_parameters['59_3']['b']
    
elif data[1] == currentRP_smooth_59:
    name1 = "J2318_59188_DZ3_RPnorm"
    smoothed_or_unsmoothed = 'smoothed'
    power = 'RP'
    
    parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP_smoothed.csv'))
    df_parameters = parameters.set_index(parameters.columns[0])
    
    tau0 = df_parameters['59_1']['tau0']
    v0 = df_parameters['59_1']['v0']
    b = df_parameters['59_1']['b']

    tau02 = df_parameters['59_2']['tau0']
    v02 = df_parameters['59_2']['v0']
    b2 = df_parameters['59_2']['b']
    
    tau03 = df_parameters['59_3']['tau0']
    v03 = df_parameters['59_3']['v0']
    b3 = df_parameters['59_3']['b']
    
elif data[1] == currentRP2_smooth_59:
    name1 = "J2318_59188_DZ3_RP2norm"
    smoothed_or_unsmoothed = 'smoothed'
    power = 'RP2'
    
    parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP2_smoothed.csv'))
    df_parameters = parameters.set_index(parameters.columns[0])
    
    tau0 = df_parameters['59_1']['tau0']
    v0 = df_parameters['59_1']['v0']
    b = df_parameters['59_1']['b']

    tau02 = df_parameters['59_2']['tau0']
    v02 = df_parameters['59_2']['v0']
    b2 = df_parameters['59_2']['b']
    
    tau03 = df_parameters['59_3']['tau0']
    v03 = df_parameters['59_3']['v0']
    b3 = df_parameters['59_3']['b']

    
#reading in spectra
norm_spectra_data = np.loadtxt(data[1])

wavelength, norm_flux, norm_error = read_spectra(norm_spectra_data)

# converting wavelength to velocity
beta = wavelength_to_velocity(redshift, wavelength, wavelength_SiIV_emit1)   


# selecting to plot in wavelength or velocity based on changeable variable section selection
if wavelength_or_velocity == 'v':
    x = beta
    x_label = 'Velocity km/s'
elif wavelength_or_velocity == 'w':
    x = wavelength
    x_label = 'Wavelength'


# plot settings
plt.xlim(xlim)
plt.xticks(np.arange(xlim[1], xlim[0] - 1, -xtick))
plt.ylim(y_min, y_max)
plt.axhline(y = 1.0)
plt.plot(x, norm_flux, color = 'k')
plt.xlabel(str(x_label))
plt.ylabel('Normalized Flux')

# establishing plot title with coice made in changeable variable section ......................................................
if title == 'yr':
    plt.title('SiIV:2020')
elif title == 'mjd':
    plt.title('SiIV EHVO (mjd:59188, year:2020) Sigma Clip')
elif title == 'paper':
    plt.title('J2318 MJD:59188 SiIV EHVO')
else:
    plt.title(str(name1)+ f"  (Sigma Clip std={sigma_clip}, GK={sigma_smooth})")

#...................................................
#...................................................
# sigma clipping

start=np.ones(len(wavelength))
totale=len(np.where(wavelength<1215.7)[0])

FLAGGUE=(ChoiClip(wavelength, norm_flux, norm_error, sigma_smooth, sigma_clip))[0]
#......................................................

#Defining area of the spectra we will be curve fitting  
(xfit, yfit, errfit) = curve_fit_area(xstart  = -78000, xend = -92000, xvalues = x[FLAGGUE==1], flux = norm_flux[FLAGGUE==1], error = norm_error)

# additional fitting data adjustments
if manual_mask == True:
   
    left = -89000
    right = -88300#-88500
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

    left = -87100
    right = -86800 #tested  -86400 to raise peak# was -86800
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

# establishing curvefitting parameters
Cf = 1
I0 = 1
vdiff = vdiff_SiIV
tau_ratio = 2

# implementing curve fitting function depending on changeable variables
if D59 == 3:
    
    if free_b == True:
        b=1200
        b2=900
        v0_59, tau0_59, b_59, Cf, v02_59, tau02_59, b2_59, v03_59, b3_59, tau03_59, v04, b4, tau04 = plotting_curvefit_test(3, 'SiIV', x, xfit, yfit, 'Cf, v0', Cf, I0, vdiff, tau_ratio,v0, tau0, b, v02, tau02, b2, v03, tau03, b3)
    else:
        v0_59, tau0_59, b_59, Cf, v02_59, tau02_59, b2_59, v03_59, b3_59, tau03_59, v04, b4, tau04 = plotting_curvefit_test(3, 'SiIV', x, xfit, yfit, 'Cf, v0, b', Cf, I0, vdiff, tau_ratio, 
                                                                 v0, tau0, b, v02, tau02, b2, v03, tau03, b3)
                                                                

    # saving fitting flux data
    ydata = np.array(curve_func( x, tau0_59, tau02_59, tau03_59, v0_59, v02_59, v03_59, b_59, b2_59, b3_59, Cf, I0, vdiff, tau_ratio))


else:
    if free_b == True:
        b=1200
        b2=900
        v0_59, tau0_59, b_59, Cf, v02_59, tau02_59, b2_59, v03, b3, tau03, v04, b4, tau04 = plotting_curvefit_test(2, 'SiIV', x, xfit, yfit, 'v0, Cf', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2)
    else:            
        v0_59, tau0_59, b_59, Cf, v02_59, tau02_59, b2_59, v03, b3, tau03, v04, b4, tau04 = plotting_curvefit_test(2, 'SiIV', x, xfit, yfit, 'Cf, v0, b', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2)

    # saving fitting flux data
    ydata = np.array(curve_func2( x, tau0_59, tau02_59, v0_59, v02_59, b_59, b2_59, Cf, I0, vdiff, tau_ratio))


# plot settings
plt.plot(x, norm_error, color = 'lightgrey', label = 'Error')
plt.legend(loc='lower right')
plt.savefig(os.path.join(OUTDIREC, f'SiIV_59188_{smoothed_or_unsmoothed}_{power}.png'))
plt.figure()
plt.show()
plt.close()


# saving fitting data
newy = np.ones_like(norm_flux)

for i, x in enumerate(newy):
    if x in xfit:  
        newy[i] = ydata[i]

error = norm_error

output_data = np.column_stack((wavelength, ydata, error))

if want_txt_save == 'y':
    np.savetxt(os.path.join(OUTDIREC, str(name1) + "_SiIVfit.txt"), output_data, header="#wavelengths flux error", comments='')
###########################################################################################################################################################################################################


################## 57328 ##############################################################################################################################################
print()
print()
print('**** FITTING MJD 57328 ****')
print()

# conditional statements based on data choice. Reading in doublet parameters for CIV fittings to be used for SiIV fittings
if data[0] == currentRP_57:
    name2 = "J2318_57328_UZ1_RPnorm"
    smoothed_or_unsmoothed = 'not smoothed'
    power = 'RP'
    
    parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP_not smoothed.csv'))
    df_parameters = parameters.set_index(parameters.columns[0])
    
    tau0 = df_parameters['57_1']['tau0']
    v0 = df_parameters['57_1']['v0']
    b = df_parameters['57_1']['b']

    tau02 = df_parameters['57_2']['tau0']
    v02 = df_parameters['57_2']['v0']
    b2 = df_parameters['57_2']['b']
    
elif data[0] == currentRP2_57:
    name2 = "J2318_57328_UZ1_RP2norm"
    power = 'RP2'
    smoothed_or_unsmoothed = 'not smoothed'
    
    parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP2_not smoothed.csv'))
    df_parameters = parameters.set_index(parameters.columns[0])
    
    tau0 = df_parameters['57_1']['tau0']
    v0 = df_parameters['57_1']['v0']
    b = df_parameters['57_1']['b']

    tau02 = df_parameters['57_2']['tau0']
    v02 = df_parameters['57_2']['v0']
    b2 = df_parameters['57_2']['b']
    
elif data[0] == currentRP_smooth_57:
    name2 = "J2318_57328_DZ1_RPnorm"
    smoothed_or_unsmoothed = 'smoothed'
    power = 'RP'
    
    parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP_smoothed.csv'))
    df_parameters = parameters.set_index(parameters.columns[0])
    
    tau0 = df_parameters['57_1']['tau0']
    v0 = df_parameters['57_1']['v0']
    b = df_parameters['57_1']['b']

    tau02 = df_parameters['57_2']['tau0']
    v02 = df_parameters['57_2']['v0']
    b2 = df_parameters['57_2']['b']
    
elif data[0] == currentRP2_smooth_57:
    name2 = "J2318_57328_DZ1_RP2norm"
    smoothed_or_unsmoothed = 'smoothed'
    power = 'RP2'
    
    parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP2_smoothed.csv'))
    df_parameters = parameters.set_index(parameters.columns[0])
    
    tau0 = df_parameters['57_1']['tau0']
    v0 = df_parameters['57_1']['v0']
    b = df_parameters['57_1']['b']

    tau02 = df_parameters['57_2']['tau0']
    v02 = df_parameters['57_2']['v0']
    b2 = df_parameters['57_2']['b']

#reading in spectra
norm_spectra_data = np.loadtxt(data[0])

wavelength, norm_flux, norm_error = read_spectra(norm_spectra_data)

# converting wavelength to velocity
beta = wavelength_to_velocity(redshift, wavelength, wavelength_SiIV_emit1)   


# selecting to plot in wavelength or velocity based on changeable variable section selection
if wavelength_or_velocity == 'v':
    x = beta
    x_label = 'Velocity km/s'
elif wavelength_or_velocity == 'w':
    x = wavelength
    x_label = 'Wavelength'



# calculate vdiff............................................................................
w1 = np.min(np.where(wavelength > 1500))
w2 = np.max(np.where(wavelength < 1490))

wave = np.array([np.mean(wavelength[w2:w1])]) 


vdiff_SiIV = wavelength_to_velocity(redshift, wave,
                    wavelength_SiIV_emit1) - wavelength_to_velocity(
                    redshift, wave, wavelength_SiIV_emit2)
                           
vdiff =  vdiff_SiIV 
 #..................................................................................................

# plot settings
plt.xlim(xlim)
plt.xticks(np.arange(xlim[1], xlim[0] - 1, -xtick))
plt.ylim(y_min, y_max)
plt.axhline(y = 1.0)
plt.xlabel(str(x_label))


if attempt57_fit == 'shifted':
    
    '''
    Here I am attempting to find an upper limit fit by shifting the spectra down by subtracting 
    by 0.5*error.
    
    '''
    # plot settings
    plt.ylabel('Normalized Flux - 0.5(error)')
    norm_flux = norm_flux - 0.5 *norm_error
    plt.plot(x, norm_flux, color = 'k')
    #...................................................
    #...................................................
    # sigma clipping

    start=np.ones(len(wavelength))
    totale=len(np.where(wavelength<1215.7)[0])

    FLAGGUE=(ChoiClip(wavelength, norm_flux, norm_error, sigma_smooth, sigma_clip))[0]
    #......................................................

    #Defining area of the spectra we will be curve fitting  
    (xfit, yfit, errfit) = curve_fit_area(xstart  = -75000, xend = -98500, xvalues = x[FLAGGUE==1], flux = norm_flux[FLAGGUE==1], error = norm_error)

    
    # plotting masked regions in grey
    plot_masked_regions(beta, norm_flux, xfit)

    # setting curve fitting parameters
    Cf = 1
    I0 = 1
    vdiff = vdiff_SiIV
    tau_ratio = 2

    # implementing curve fitting function based on changeable variable settings
    if free_b == True:
        b2 =800 #not been changed for 57
        b=500
        v0_57, tau0_57, b_57, Cf, v02_57, tau02_57, b2_57, v03, b3, tau03, v04, b4, tau04 = plotting_curvefit_test(2, 'SiIV', x, xfit, yfit, 'v0, Cf', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2)
    else:
        v0_57, tau0_57, b_57, Cf, v02_57, tau02_57, b2_57, v03, b3, tau03, v04, b4, tau04 = plotting_curvefit_test(2, 'SiIV', x, xfit, yfit, 'Cf, v0, b', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2)


    # saving fit data
    ydata = np.array(curve_func2( x, tau0_57, tau02_57, v0_57, v02_57, b_57, b2_57, Cf, I0, vdiff, tau_ratio))
    
    newy = np.ones_like(norm_flux)
    
    for i, x in enumerate(newy):
        if x in xfit:  
            newy[i] = ydata[i]
    
    error = norm_error
    
    output_data = np.column_stack((wavelength, ydata, error))
    
    if want_txt_save == 'y':
        np.savetxt(os.path.join(OUTDIREC, str(name2) + "_SiIVfit_spectra_down.txt"), output_data, header="#wavelengths flux error", comments='')
    

elif attempt57_fit == 'og':
    
    '''
    Here I am attempting to find if there is any SiIV absorption with no shifting.
    
    '''
    plt.ylabel('Normalized Flux')
    norm_flux = norm_flux
    plt.plot(x, norm_flux, color = 'k')
    #...................................................
    #...................................................
    # sigma clipping

    start=np.ones(len(wavelength))
    totale=len(np.where(wavelength<1215.7)[0])

    FLAGGUE=(ChoiClip(wavelength, norm_flux, norm_error, sigma_smooth, sigma_clip))[0]
    #......................................................

    #Defining area of the spectra we will be curve fitting  
    (xfit, yfit, errfit) = curve_fit_area(xstart  = -75000, xend = -98500, xvalues = x[FLAGGUE==1], flux = norm_flux[FLAGGUE==1], error = norm_error)

    # additional fitting data adjustments
    if manual_mask == True:
       
        left = -77400
        right = -77050
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
       
        left = -82780
        right = -82000
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

        left = -89500
        right = -88500 #tested  -86400 to raise peak# was -86800
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

    
    # setting curve fitting parameters
    Cf = 1
    I0 = 1
    vdiff = vdiff_SiIV
    tau_ratio = 2

    # implementing curve fitting function based on changeable variables
    if free_b == True:
        print('Warning the b guesses have not been adjusted to fit with free b for obs 57328')
        b2 =800 #not been changed for 57
        b=500
        v0_57, tau0_57, b_57, Cf, v02_57, tau02_57, b2_57, v03, b3, tau03, v04, b4, tau04 = plotting_curvefit_test(2, 'SiIV', x, xfit, yfit, 'v0, Cf', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2)
    else:
        v0_57, tau0_57, b_57, Cf, v02_57, tau02_57, b2_57, v03, b3, tau03, v04, b4, tau04 = plotting_curvefit_test(2, 'SiIV', x, xfit, yfit, 'Cf, v0, b', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2)

    # plotting masked regions as grey
    plt.plot(x, norm_error, color = 'lightgrey', label = 'Error')

    # saving fitting data
    ydata = np.array(curve_func2( x, tau0_57, tau02_57, v0_57, v02_57, b_57, b2_57, Cf, I0, vdiff, tau_ratio))
    
    newy = np.ones_like(norm_flux)
    
    for i, x in enumerate(newy):
        if x in xfit:  
            newy[i] = ydata[i]
    
    error = norm_error
    
    output_data = np.column_stack((wavelength, ydata, error))
    
    if want_txt_save == 'y':
        np.savetxt(os.path.join(OUTDIREC, str(name2) + "_SiIVfit.txt"), output_data, header="#wavelengths flux error", comments='')
    

else:
    plt.plot(x, norm_flux, color = 'k')   
    plt.ylabel('Normalized Flux')

# establishing plot title with coice made in changeable variable section ......................................................
if title == 'yr':
    plt.title('SiIV: 2015')
elif title == 'mjd':
    plt.title('SiIV  (mjd:57328, year: 2015) Sigma Clip')
elif title == 'paper':
    plt.title('J2318 MJD:57328 SiIV EHVO')
else:
    plt.title(str(name2)+ f"  (Sigma Clip std={sigma_clip}, GK={sigma_smooth})")

# plot settings
plt.legend(loc='lower right')
plt.savefig(os.path.join(OUTDIREC, f'SiIV_57328_{smoothed_or_unsmoothed}_{power}.png'))
plt.figure()
plt.show()
plt.close()

##############################################################################################################################################
##############################################################################################################################################
    
#Running shifted spectra fittings for constraining absorption.py measurement errors ------------------------------------------------------------------------------------------------------------
###########################################################################################################################################################################################################
if shifted_spec_fit == True:
    for i in range(0,2):
        # conditional statements based on data choice. Reading in doublet parameters for CIV fittings to be used for SiIV fittings
        if data[2] == currentRP_60:
            name2 = "J2318_60251_UZ2_RPnorm"
            smoothed_or_unsmoothed = 'not smoothed'
            power = 'RP'
            
            parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP_not smoothed.csv'))
            df_parameters = parameters.set_index(parameters.columns[0])
            
            tau0 = df_parameters['60_1']['tau0']
            v0 = df_parameters['60_1']['v0']
            b = df_parameters['60_1']['b']
        
            tau02 = df_parameters['60_2']['tau0']
            v02 = df_parameters['60_2']['v0']
            b2 = df_parameters['60_2']['b']
            
        elif data[2] == currentRP2_60:
            name2 = "J2318_60251_UZ2_RP2norm"
            power = 'RP2'
            smoothed_or_unsmoothed = 'not smoothed'
            
            parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP2_not smoothed.csv'))
            df_parameters = parameters.set_index(parameters.columns[0])
            
            tau0 = df_parameters['60_1']['tau0']
            v0 = df_parameters['60_1']['v0']
            b = df_parameters['60_1']['b']
        
            tau02 = df_parameters['60_2']['tau0']
            v02 = df_parameters['60_2']['v0']
            b2 = df_parameters['60_2']['b']
            
        elif data[2] == currentRP_smooth_60:
            name2 = "J2318_60251_DZ2_RPnorm"
            smoothed_or_unsmoothed = 'smoothed'
            power = 'RP'
            
            parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP_smoothed.csv'))
            df_parameters = parameters.set_index(parameters.columns[0])
            
            tau0 = df_parameters['60_1']['tau0']
            v0 = df_parameters['60_1']['v0']
            b = df_parameters['60_1']['b']
        
            tau02 = df_parameters['60_2']['tau0']
            v02 = df_parameters['60_2']['v0']
            b2 = df_parameters['60_2']['b']
            
        elif data[2] == currentRP2_smooth_60:
            name2 = "J2318_60251_DZ2_RP2norm"
            smoothed_or_unsmoothed = 'smoothed'
            power = 'RP2'
            
            parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP2_smoothed.csv'))
            df_parameters = parameters.set_index(parameters.columns[0])
            
            tau0 = df_parameters['60_1']['tau0']
            v0 = df_parameters['60_1']['v0']
            b = df_parameters['60_1']['b']
        
            tau02 = df_parameters['60_2']['tau0']
            v02 = df_parameters['60_2']['v0']
            b2 = df_parameters['60_2']['b']
            
        #reading in spectra
        norm_spectra_data = np.loadtxt(data[2])
        
        wavelength, norm_flux, norm_error = read_spectra(norm_spectra_data)
        
        # shifting flux by percentage of error
        if i == 0:
            shift = 0.5
            shiftname = 'up0.5err'
            pltname = ' + 0.5err'
            
        if i == 1:
            shift = - 0.5
            shiftname = 'down0.5err'
            pltname = ' - 0.5err'
    
        
        norm_flux = norm_flux + shift * norm_error
        
        # converting wavelength to velocity
        beta = wavelength_to_velocity(redshift, wavelength, wavelength_SiIV_emit1)   
        
        
    
        # selecting to plot in wavelength or velocity based on changeable variable section selection
        if wavelength_or_velocity == 'v':
            x = beta
            x_label = 'Velocity km/s'
        elif wavelength_or_velocity == 'w':
            x = wavelength
            x_label = 'Wavelength'
        
        
        
        ################## 60251 ########################################################
        
        # calculate vdiff............................................................................
        w1 = np.min(np.where(wavelength > 1500))
        w2 = np.max(np.where(wavelength < 1490))
        
        wave = np.array([np.mean(wavelength[w2:w1])]) #made np.array bc my version of w to v itterates over a list of values - LEF
        
        
        vdiff_SiIV = wavelength_to_velocity(redshift, wave,
                            wavelength_SiIV_emit1) - wavelength_to_velocity(
                            redshift, wave, wavelength_SiIV_emit2)
                     
                             
        vdiff =  vdiff_SiIV 
        #.............................................................................................
        
        # plot settings
        plt.xlim(xlim)
        plt.xticks(np.arange(xlim[1], xlim[0] - 1, -xtick))
        plt.ylim(y_min, y_max)
        plt.axhline(y = 1.0)
        plt.plot(x, norm_flux, color = 'k')
        plt.xlabel(str(x_label))
        plt.ylabel(f'Normalized Flux{pltname}')
        
        # setting plot title based on changeable variable selection
        if title == 'yr':
            plt.title('SiIV: 2023')
        elif title == 'mjd':
            plt.title('SiIV  (mjd:60251, year: 2023) Sigma Clip')
        elif title == 'paper':
            plt.title('J2318 MJD:60251 SiIV EHVO')
        else:
            plt.title(str(name2)+ f"  (Sigma Clip std={sigma_clip}, GK={sigma_smooth})")
        
        #...................................................
        #...................................................
        # sigma clipping
        
        start=np.ones(len(wavelength))
        totale=len(np.where(wavelength<1215.7)[0])
        
        FLAGGUE=(ChoiClip(wavelength, norm_flux, norm_error, sigma_smooth, sigma_clip))[0]
        
        #......................................................
        #Defining area of the spectra we will be curve fitting  
        (xfit, yfit, errfit) = curve_fit_area(xstart = -71000, xend = -94000, xvalues = x[FLAGGUE==1], flux = norm_flux[FLAGGUE==1], error = norm_error)
        
        # additional fitting data adjustments
        if manual_mask == True:
            
            left = -87700
            right = -86500#-87200
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
            
         
            left = -84500
            right = -84200
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
    
         
            left = -83750#-83700
            right = -82950#-83100
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
            
         
            left = -82300
            right = -81700
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
            
         
            left = -81400
            right = -81000
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
        
        # setting parameters for curvefitting function
        Cf = 1
        I0 = 1
        vdiff = vdiff_SiIV
        tau_ratio = 2
        
        if free_b == True:
            b2 =800
            b=500
            v0_60, tau0_60, b_60, Cf, v02_60, tau02_60, b2_60, v03, b3, tau03, v04, b4, tau04 = plotting_curvefit_test(2, 'SiIV', x, xfit, yfit, 'v0, Cf', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2)
        else:
            v0_60, tau0_60, b_60, Cf, v02_60, tau02_60, b2_60, v03, b3, tau03, v04, b4, tau04 = plotting_curvefit_test(2, 'SiIV', x, xfit, yfit, 'Cf, v0, b', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2)
        
        # plot settings
        plt.plot(x, norm_error, color = 'lightgrey', label = 'Error')
        plt.legend(loc='lower right')
        plt.figure()
        plt.show()
        plt.close()
        
        # saving fitting data
        ydata = np.array(curve_func2( x, tau0_60, tau02_60, v0_60, v02_60, b_60, b2_60, Cf, I0, vdiff, tau_ratio))
        
        newy = np.ones_like(norm_flux)
        
        for i, x in enumerate(newy):
            if x in xfit:  
                newy[i] = ydata[i]
        
        error = norm_error
        
        output_data = np.column_stack((wavelength, ydata, error))
        
        if want_txt_save == 'y':
            np.savetxt(os.path.join(OUTDIREC, str(name2) + f"_SiIVfit_spectra_{shiftname}.txt"), output_data, header="#wavelengths flux error", comments='')
        
        
        
        ################## 59188 ########################################################
        # conditional statements based on data choice. Reading in doublet parameters for CIV fittings to be used for SiIV fittings
        if data[1] == currentRP_59:
            name1 = "J2318_59188_UZ3_RPnorm"
            smoothed_or_unsmoothed = 'not smoothed'
            power = 'RP'
            
            parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP_not smoothed.csv'))
            df_parameters = parameters.set_index(parameters.columns[0])
            
            tau0 = df_parameters['59_1']['tau0']
            v0 = df_parameters['59_1']['v0']
            b = df_parameters['59_1']['b']
        
            tau02 = df_parameters['59_2']['tau0']
            v02 = df_parameters['59_2']['v0']
            b2 = df_parameters['59_2']['b']
            
            tau03 = df_parameters['59_3']['tau0']
            v03 = df_parameters['59_3']['v0']
            b3 = df_parameters['59_3']['b']
            
        elif data[1] == currentRP2_59:
            name1 = "J2318_59188_UZ3_RP2norm"
            power = 'RP2'
            smoothed_or_unsmoothed = 'not smoothed'
            
            parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP2_not smoothed.csv'))
            df_parameters = parameters.set_index(parameters.columns[0])
            
            tau0 = df_parameters['59_1']['tau0']
            v0 = df_parameters['59_1']['v0']
            b = df_parameters['59_1']['b']
        
            tau02 = df_parameters['59_2']['tau0']
            v02 = df_parameters['59_2']['v0']
            b2 = df_parameters['59_2']['b']
            
            tau03 = df_parameters['59_3']['tau0']
            v03 = df_parameters['59_3']['v0']
            b3 = df_parameters['59_3']['b']
            
        elif data[1] == currentRP_smooth_59:
            name1 = "J2318_59188_DZ3_RPnorm"
            smoothed_or_unsmoothed = 'smoothed'
            power = 'RP'
            
            parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP_smoothed.csv'))
            df_parameters = parameters.set_index(parameters.columns[0])
            
            tau0 = df_parameters['59_1']['tau0']
            v0 = df_parameters['59_1']['v0']
            b = df_parameters['59_1']['b']
        
            tau02 = df_parameters['59_2']['tau0']
            v02 = df_parameters['59_2']['v0']
            b2 = df_parameters['59_2']['b']
            
            tau03 = df_parameters['59_3']['tau0']
            v03 = df_parameters['59_3']['v0']
            b3 = df_parameters['59_3']['b']
            
        elif data[1] == currentRP2_smooth_59:
            name1 = "J2318_59188_DZ3_RP2norm"
            smoothed_or_unsmoothed = 'smoothed'
            power = 'RP2'
            
            parameters = pd.read_csv(os.path.join(OUTDIREC,'doublet_params_RP2_smoothed.csv'))
            df_parameters = parameters.set_index(parameters.columns[0])
            
            tau0 = df_parameters['59_1']['tau0']
            v0 = df_parameters['59_1']['v0']
            b = df_parameters['59_1']['b']
        
            tau02 = df_parameters['59_2']['tau0']
            v02 = df_parameters['59_2']['v0']
            b2 = df_parameters['59_2']['b']
            
            tau03 = df_parameters['59_3']['tau0']
            v03 = df_parameters['59_3']['v0']
            b3 = df_parameters['59_3']['b']
        
            
        #reading in spectra
        norm_spectra_data = np.loadtxt(data[1])
        
        wavelength, norm_flux, norm_error = read_spectra(norm_spectra_data)
        
        # shifting flux with percentage of error
        if i == 0:
            shift = 0.5
            shiftname = 'up0.5err'
            pltname = ' + 0.5err'
            
        if i == 1:
            shift = - 0.5
            shiftname = 'down0.5err'
            pltname = ' - 0.5err'
    
        
        norm_flux = norm_flux + shift * norm_error
        
        # converting wavelength to velocity
        beta = wavelength_to_velocity(redshift, wavelength, wavelength_SiIV_emit1)   
        
        
        # selecting to plot in wavelength or velocity based on changeable variable section selection
        if wavelength_or_velocity == 'v':
            x = beta
            x_label = 'Velocity km/s'
        elif wavelength_or_velocity == 'w':
            x = wavelength
            x_label = 'Wavelength'
        
        
        
        
        # establishing plot title with coice made in changeable variable section ......................................................
        plt.xlim(xlim)
        plt.xticks(np.arange(xlim[1], xlim[0] - 1, -xtick))
        plt.ylim(y_min, y_max)
        plt.axhline(y = 1.0)
        plt.plot(x, norm_flux, color = 'k')
        plt.xlabel(str(x_label))
        plt.ylabel(f'Normalized Flux{pltname}')
        
        # establishing plot title with coice made in changeable variable section ......................................................
        if title == 'yr':
            plt.title('SiIV:2020')
        elif title == 'mjd':
            plt.title('SiIV EHVO (mjd:59188, year:2020) Sigma Clip')
        elif title == 'paper':
            plt.title('J2318 MJD:59188 SiIV EHVO')
        else:
            plt.title(str(name1)+ f"  (Sigma Clip std={sigma_clip}, GK={sigma_smooth})")
        
        #...................................................
        #...................................................
        # sigma clipping
        
        start=np.ones(len(wavelength))
        totale=len(np.where(wavelength<1215.7)[0])
        
        FLAGGUE=(ChoiClip(wavelength, norm_flux, norm_error, sigma_smooth, sigma_clip))[0]
        
        #......................................................
        #Defining area of the spectra we will be curve fitting  
        (xfit, yfit, errfit) = curve_fit_area(xstart  = -78000, xend = -92000, xvalues = x[FLAGGUE==1], flux = norm_flux[FLAGGUE==1], error = norm_error)
        
        # additional fitting data adjustments
        if manual_mask == True:
             
            left = -89000
            right = -88300#-88500
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
        
         
            left = -87100
            right = -86800 #tested  -86400 to raise peak# was -86800
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
        
        # setting curve fitting parameters
        Cf = 1
        I0 = 1
        vdiff = vdiff_SiIV
        tau_ratio = 2
        
        # implementing curve fitting function based on changeable variables
        if D59 == 3:
            
            if free_b == True:
                b=1200
                b2=900
                v0_59, tau0_59, b_59, Cf, v02_59, tau02_59, b2_59, v03_59, b3_59, tau03_59, v04, b4, tau04 = plotting_curvefit_test(3, 'SiIV', x, xfit, yfit, 'Cf, v0', Cf, I0, vdiff, tau_ratio,v0, tau0, b, v02, tau02, b2, v03, tau03, b3)
            else:
                v0_59, tau0_59, b_59, Cf, v02_59, tau02_59, b2_59, v03_59, b3_59, tau03_59, v04, b4, tau04 = plotting_curvefit_test(3, 'SiIV', x, xfit, yfit, 'Cf, v0, b', Cf, I0, vdiff, tau_ratio, 
                                                                         v0, tau0, b, v02, tau02, b2, v03, tau03, b3)
                                                                        
        
            # saving fit flux data
            ydata = np.array(curve_func( x, tau0_59, tau02_59, tau03_59, v0_59, v02_59, v03_59, b_59, b2_59, b3_59, Cf, I0, vdiff, tau_ratio))
        
        
        else:
            if free_b == True:
                b=1200
                b2=900
                v0_59, tau0_59, b_59, Cf, v02_59, tau02_59, b2_59, v03, b3, tau03, v04, b4, tau04 = plotting_curvefit_test(2, 'SiIV', x, xfit, yfit, 'v0, Cf', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2)
            else:            
                v0_59, tau0_59, b_59, Cf, v02_59, tau02_59, b2_59, v03, b3, tau03, v04, b4, tau04 = plotting_curvefit_test(2, 'SiIV', x, xfit, yfit, 'Cf, v0, b', Cf, I0, vdiff, tau_ratio, v0, tau0, b, v02, tau02, b2)
        
            def curve_func( v, tau0, tau02, v0, v02, b, b2, Cf, I0, vdiff, tau_ratio): 
                return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)) * np.exp(-tau_v(v,v02,b2,tau02)) * np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)))) 
        
            # saving fit flux data
            ydata = np.array(curve_func( x, tau0_59, tau02_59, v0_59, v02_59, b_59, b2_59, Cf, I0, vdiff, tau_ratio))
        
        # plot settings        
        plt.plot(x, norm_error, color = 'lightgrey', label = 'Error')
        plt.legend(loc='lower right')
        plt.figure()
        plt.show()
        plt.close()
        
        
        # saving fit data
        newy = np.ones_like(norm_flux)
        
        for i, x in enumerate(newy):
            if x in xfit:  
                newy[i] = ydata[i]
        
        error = norm_error
        
        output_data = np.column_stack((wavelength, ydata, error))
        
        if want_txt_save == 'y':
            np.savetxt(os.path.join(OUTDIREC, str(name1) + f"_SiIVfit_spectra_{shiftname}.txt"), output_data, header="#wavelengths flux error", comments='')
        
else:
    None
