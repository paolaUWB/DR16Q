#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 15:20:10 2026

@author: lilianaflores
"""

import matplotlib.pyplot as plt
import os
import sys
import numpy as np
from utility_functions import clear_file, read_list_spectra, read_spectra, append_row_to_csv
from abs_function_module import wavelength_to_velocity
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)

file = os.getcwd() + '/../DR16Q_EHVO/NORM_DR16Q_EHVO/spec-5212-56016-0348norm.dr16'
z=2.138

def scat_plot_spec(file,z, name):
    norm_spectra_data = np.loadtxt(file)
    
    # setting a variable for each of those values from the spectra
    wavelength, normalized_flux, normalized_error = read_spectra(norm_spectra_data)
    
    beta = wavelength_to_velocity(z, wavelength)
    
    plt.scatter(beta, normalized_flux, s=10, color='blue')
    #plt.plot(beta, normalized_flux, color='blue',linestyle='--')
    plt.plot(beta, np.ones_like(beta), linestyle='--', color='k')
    plt.xlim(-50000,-30000)
    #plt.xlim(-70000,0)
    plt.ylim(0,1.75)
    plt.yticks(np.arange(0, 1.76, 0.25))
    plt.ylabel('Normalized Flux')
    plt.xlabel('Velocity km/s')
    plt.title(name)
    plt.show()
    plt.close()
    
file = os.getcwd() + '/../DR16Q_EHVO/NORM_DR16Q_EHVO/spec-6199-56220-0264norm.dr16'
z=2.762
name = 'spec-6199-56220-0264norm.dr16'

#scat_plot_spec(file, z, name)


file = os.getcwd() + '/../DR16Q_EHVO/NORM_DR16Q_EHVO/spec-5896-56047-0658norm.dr16'
z=2.892
name = 'spec-5896-56047-0658norm.dr16'
scat_plot_spec(file, z, name)




