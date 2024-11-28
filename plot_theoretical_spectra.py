#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 15:19:10 2024

@author: lilianaflores
"""

import os
import sys
import numpy as np 
import pandas as pd
import math
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.getcwd()))
from utility_functions import read_spectra
sys.path.insert(0, os.path.dirname(os.getcwd()) + '/ABSORPTION')
from ABSORPTION.abs_function_module import wavelength_to_velocity, smooth
'''
spectra = os.getcwd() + "/bh9_b364c3_l-2_al075_mw90_wa423.spec"


df = pd.read_csv(spectra, delim_whitespace=True, on_bad_lines='warn')
wavelength = df[df.columns[1]].to_numpy()
flux_20 = df[df.columns[10]].to_numpy()
flux_40 = df[df.columns[11]].to_numpy()
flux_60 = df[df.columns[12]].to_numpy()
flux_67 = df[df.columns[13]].to_numpy()
flux_69 = df[df.columns[14]].to_numpy()
flux_70 = df[df.columns[15]].to_numpy()
flux_71 = df[df.columns[16]].to_numpy()
flux_72 = df[df.columns[17]].to_numpy()
flux_73 = df[df.columns[18]].to_numpy()
flux_74 = df[df.columns[19]].to_numpy()
flux_75 = df[df.columns[20]].to_numpy()
flux_77 = df[df.columns[21]].to_numpy()
flux_80 = df[df.columns[22]].to_numpy()
flux_82 = df[df.columns[23]].to_numpy()
flux_85 = df[df.columns[24]].to_numpy()
flux_89 = df[df.columns[25]].to_numpy()
'''


def plot_spectra( q, wavelength_or_velocity, want_to_smooth, boxcar_size, save_plot=False, plot_filename='plot.png', x_min=None, x_max=None, y_min=None, y_max=None):
    
    
        
    spectra = os.getcwd() + "/bh9_b364c3_l-2_al075_mw90_wa423.spec"
    #Need to delete commented header in .spec file for this code to work


    df = pd.read_csv(spectra, delim_whitespace=True, on_bad_lines='warn')
    wavelength = df[df.columns[1]].to_numpy()
    redshift = 0
    
    
    if q == 20:
        flux = df[df.columns[10]].to_numpy()
    elif q == 40:
        flux = df[df.columns[11]].to_numpy()
    elif q == 60:
        flux = df[df.columns[12]].to_numpy()
    elif q == 67:
        flux = df[df.columns[13]].to_numpy()
    elif q == 69:
        flux = df[df.columns[14]].to_numpy()
    elif q == 70:
        flux = df[df.columns[15]].to_numpy()
    elif q == 71:
        flux = df[df.columns[16]].to_numpy()
    elif q == 72:
        flux = df[df.columns[17]].to_numpy()
    elif q == 73:
        flux = df[df.columns[18]].to_numpy()
    elif q == 74:
        flux = df[df.columns[19]].to_numpy()
    elif q == 75:
        flux = df[df.columns[20]].to_numpy()
    elif q == 77:
        flux = df[df.columns[21]].to_numpy()
    elif q == 80:
        flux = df[df.columns[22]].to_numpy()
    elif q == 82:
        flux = df[df.columns[23]].to_numpy()
    elif q == 85:
        flux = df[df.columns[24]].to_numpy
    elif q == 89:
        flux = df[df.columns[25]].to_numpy
   
    redshift = 0
    velocity = wavelength_to_velocity(redshift, wavelength)
        
    if want_to_smooth == 'yes':
        flux = smooth(flux, boxcar_size)
        #error = smooth(error, boxcar_size) / math.sqrt(boxcar_size)
        
    if wavelength_or_velocity == 'v':
        n = velocity
        x_label = 'Velocity km/s'
    elif wavelength_or_velocity == 'w':
        n = wavelength
        x_label = 'Wavelength'
        
    if x_min is not None and x_max is not None:
        plt.xlim(x_min, x_max)
        
    else:
        plt.xlim(min(n) - 3000, max(n) + 3000)
        
    if y_min is not None and y_max is not None:
        plt.ylim(y_min, y_max)
    else:
        plt.ylim(min(flux) - 1, max(flux) + 1)
        
    plt.plot(n, flux, color = 'k')
    #plt.plot(n, error, color = 'grey')
    plt.xlabel(str(x_label))
    plt.ylabel('Flux')
    plt.title('LOS Spectra '+ str(q)+'\u00b0')
    
   
    if save_plot:
        plt.savefig(plot_filename)

    #plt.show()
    

plot_spectra(70, 'v', 'no', 0, x_min=-200000, x_max=0, y_min=0, y_max=1.0)
       
