#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 12:21:18 2026

@author: lilianaflores
"""

import os
import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.getcwd() + '/../')
from abs_function_module import wavelength_to_velocity


#spectra files for the 3 preliminary EHVOs out of the x-ray selected sample of quasars from Hiremath+2025
file = os.getcwd() + '/Recons_Hiremath2025/spec-allepoch-59381-4350968445.fits'
file2 = os.getcwd() + '/Recons_Hiremath2025/spec-allepoch-60071-27021598160032302.fits'
file3 = os.getcwd() + '/Recons_Hiremath2025/spec-allepoch-59318-4570380920.fits'



def Plot_spec_compare_morphed(file, xlims = None, ylims = None):
    '''
    Parameters
    ----------
    file : str
        Provide a file pathway.
    xlims: tuple
        Optionally specify x limits when plotting.
    ylims: tuple
        Optionally specify y limits when plotting.

    Returns
    -------
    Plot of x-ray selected quasar spectra comparing the Morphed, Normalized, and Reconstruction spectra.

    '''
    
    data = fits.open(file)
    
    #print(data[1].columns) #run to print column names
    
    data_tab = data[1].data
    
    wavelength = data_tab['wave']
    flux = data_tab['flux']
    noise = data_tab['noise']
    mask = data_tab['mask']
    morph = data_tab['morph']
    recon = data_tab['recon']
    
    
    beta2 = wavelength_to_velocity(0, wavelength)
    
    plt.plot(beta2, flux, color = 'red', label='Morphed')
    plt.plot(beta2, flux/recon, color = 'purple', label = 'Normalized')
    plt.plot(beta2, recon, color = 'orange', label = 'Recon')
    plt.plot(beta2,np.ones_like(beta2), color='k', linestyle='--')
    plt.plot(beta2,noise, color='grey', label='Noise')

    plt.ylabel('Normalized Flux')
    plt.xlabel('Velocity km/s')
    plt.ylim(ylims)
    plt.xlim(xlims)
    plt.legend(loc='upper right')
    plt.title(file[70:])
#    plt.show()
#    plt.close()
    
    
def Plot_spec_compare_og(file, xlims = None, ylims = None):
    '''
    Parameters
    ----------
    file : str
        Provide a file pathway.
    xlims: tuple
        Optionally specify x limits when plotting.
    ylims: tuple
        Optionally specify y limits when plotting.

    Returns
    -------
    Plot of x-ray selected quasar spectra comparing the original spectrum and reconstruction.

    '''
    
    data = fits.open(file)
    
    #print(data[1].columns) #run to print column names
    
    data_tab = data[1].data
    
    wavelength = data_tab['wave']
    flux = data_tab['flux']
    noise = data_tab['noise']
    mask = data_tab['mask']
    morph = data_tab['morph']
    recon = data_tab['recon']
    
    
    beta2 = wavelength_to_velocity(0, wavelength)
    
    plt.plot(beta2, flux*morph, color = 'red', label='Original')
    plt.plot(beta2, recon*morph, color = 'orange', label = 'Recon')
    plt.plot(beta2,noise*morph, color='grey', label='Noise')

    plt.ylabel('Flux')
    plt.xlabel('Velocity km/s')
    plt.ylim(ylims)
    plt.xlim(xlims)
    plt.legend(loc='upper right')
    plt.title(file[70:])
#   plt.show()
#   plt.close()
    
    
    
#Plotting spectra comparing plots for all three spectra files using a for loop
files = [file, file2, file3]

xlims = -70000,20000
ylims = 0,3


for i in np.arange(len(files)):
    Plot_spec_compare_morphed(files[i])
    Plot_spec_compare_og(files[i])





