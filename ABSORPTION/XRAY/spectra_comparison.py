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


# take in two files, one morphed/og and the other is the SDSS DR-19
#     Add parameter for the file
#     Convert wavelength to velocity
# Morphed = normalized
def Plot_spec_compare_morphed(recon_file, xlims = None, ylims = None):
    '''
    Parameters
    ----------
    file : str
        Provide a file pathway.
    xlims: tuple
        Optionally specify x limits when plotting.
    ylims: tuple
        Optionally specify y limits when plotting.
        If not, it will automatically take the range of the mask of the given parameters in the flux and find the max value and at 5% to it

    Returns
    -------
    Plot of x-ray selected quasar spectra comparing the Morphed, Normalized, and Reconstruction spectra.

    '''
    

    data = fits.open(recon_file)
    
    #print(data[1].columns) #run to print column names
    
    data_tab = data[1].data

    # Use these and whatever column names I need
    # Run columns line

    wavelength = data_tab['wave']
    flux = data_tab['flux']
    noise = data_tab['noise']
#   mask = data_tab['mask']
#   morph = data_tab['morph']
    recon = data_tab['recon']
    
    
    beta2 = wavelength_to_velocity(0, wavelength)
    
    plt.plot(beta2, flux, color = 'red', label='Morphed')
    plt.plot(beta2, flux/recon, color = 'purple', label = 'Normalized')
    plt.plot(beta2, recon, color = 'orange', label = 'Recon')
    plt.plot(beta2,np.ones_like(beta2), color='k', linestyle='--')
    plt.plot(beta2,noise, color='grey', label='Noise')

    plt.ylabel('Normalized Flux')
    plt.xlabel('Velocity km/s')
    
    #Check if ylims were defined
    if ylims is not None:
        plt.ylim(ylims)

    #If ylims were not defines, check if xlims were defined if it was we will find the ylims in the range of xlims
    elif xlims is not None:
        xmin, xmax = xlims
        
        #Make sure beta2 is in the range of the velocites we care about
        mask_range = (beta2 >= xmin) & (beta2 <= xmax)
    
        if np.any(mask_range):
            y_candidates = np.concatenate([
                flux[mask_range],
                (flux/recon)[mask_range],
                recon[mask_range],
                noise[mask_range],
            ])
            ymax = np.max(y_candidates)
            plt.ylim(top=ymax + (0.05*ymax))
    
    plt.xlim(xlims)
    plt.legend(loc='upper right')
    plt.title(file[70:])
#    plt.show()
#    plt.close()


# Flux * Morph = og
# Take that out and plot SDSS file as it is not normalized, keep it the same
# May need to convert wavelength to velocity
"""
def Plot_spec_compare_full_sdss(og_file, sdss_file, xlims = None, ylims = None):
    '''
    Parameters
    ----------
    file : str
        Provide a file pathway.
    xlims: tuple
        Optionally specify x limits when plotting.
    ylims: tuple
        Optionally specify y limits when plotting.
        If not, it will automatically take the range of the mask of the given parameters in the flux and find the max value and at 5% to it

    Returns
    -------
    Plot of x-ray selected quasar spectra comparing the original spectrum and reconstruction.

    '''
    
    # ---------------- OG DATA ----------------
    data = fits.open(og_file)
    
    #print(data[1].columns) #run to print column names
    
    data_tab = data[1].data
    
    wavelength = data_tab['wave']
    flux = data_tab['flux']
    noise = data_tab['noise']
#   mask = data_tab['mask']
    morph = data_tab['morph']
    recon = data_tab['recon']
    
    data.close()
    
    
    # ------------------SDSS DATA -------------------------
    sdss = fits.open(sdss_file)
    sdss_tab = sdss[1].data
    '''
    flux = data_tab['flux']
    noise = data_tab['noise']
#   mask = data_tab['mask']
    morph = data_tab['morph']
    recon = data_tab['recon']
    '''
    sdss_wavelength = 10**sdss_tab['LOGLAM']
    sdss_flux = sdss_tab['FLUX']
    sdss_ivar = sdss_tab['IVAR']
    sdss_noise = np.zeros_like(sdss_flux)
    sdss_redshift = np.zeros_like(sdss_flux)
    
    sdss.close()
    # ----------------- PLOTTING ----------------
    beta = wavelength_to_velocity(0, wavelength)
    beta2 = wavelength_to_velocity(0, sdss_wavelength)
    
    
    plt.plot(beta2, sdss_flux, color = 'red', label='sdss full')
    plt.plot(beta, recon*morph, color = 'blue', label = 'Recon')
    plt.plot(beta, noise*morph, color='grey', label='Noise')
    

    plt.ylabel('og Flux')
    plt.xlabel('Velocity km/s')
    
    #Check if ylims were defined
    if ylims is not None:
        plt.ylim(ylims)

    #If ylims were not defines, check if xlims were defined if it was we will find the ylims in the range of xlims
    elif xlims is not None:
        xmin, xmax = xlims
        
        #Make sure beta2 is in the range of the velocites we care about
        mask_range = (beta2 >= xmin) & (beta2 <= xmax)
    
        if np.any(mask_range):
            y_candidates = np.concatenate([
                flux[mask_range],
                (flux/recon)[mask_range],
                recon[mask_range],
                noise[mask_range],
            ])
            ymax = np.max(y_candidates)
            plt.ylim(top=ymax + (0.05*ymax))
            

    
    plt.xlim(xlims)
    plt.legend(loc='upper right')
    plt.title(file[70:])
"""
def Plot_spec_compare_full_sdss(og_file, sdss_file, xlims=None, ylims=None):

    # -------- OG DATA --------
    data = fits.open(og_file)
    data_tab = data[1].data

    wavelength = data_tab['wave']
    flux = data_tab['flux']
    morph = data_tab['morph']

    data.close()

    # -------- SDSS DATA --------
    sdss = fits.open(sdss_file)
    sdss_tab = sdss[1].data

    sdss_wavelength = 10**sdss_tab['LOGLAM']
    sdss_flux = sdss_tab['FLUX']

    sdss.close()

    # --- SCALE SDSS TO MATCH OG ---
    scale = np.nanmedian(flux * morph) / np.nanmedian(sdss_flux)
    sdss_flux_scaled = sdss_flux * scale
    
    # -------- VELOCITY --------
    beta  = wavelength_to_velocity(0, wavelength)
    beta2 = wavelength_to_velocity(0, sdss_wavelength)

    # -------- PLOT --------
    fig, ax = plt.subplots()

    ax.plot(beta2, sdss_flux_scaled, color='red', label='SDSS Full')
    ax.plot(beta, flux*morph, color='blue', label='Original')

    ax.set_ylabel('Flux')
    ax.set_xlabel('Velocity km/s')
    ax.legend(loc='upper right')
    ax.set_title(os.path.basename(og_file))

    # -------- Y-LIMITS --------
    if xlims is not None:
        xmin, xmax = xlims

        mask_sdss = (beta2 >= xmin) & (beta2 <= xmax)
        mask_og   = (beta  >= xmin) & (beta  <= xmax)

        y_values = []

        if np.any(mask_sdss):
            y_values.append(sdss_flux_scaled[mask_sdss])

        if np.any(mask_og):
            y_values.append((flux*morph)[mask_og])

        if len(y_values) > 0:
            y_all = np.concatenate(y_values)
            ymin = np.min(y_all)
            ymax = np.max(y_all)

            ax.set_ylim(ymin * 0.95, ymax * 1.05)

    if xlims is not None:
        ax.set_xlim(xlims)

    return fig


"""    
#Plotting spectra comparing plots for all three spectra files using a for loop
files = [file, file2, file3]

xlims = -70000,20000
ylims = 0,3


for i in np.arange(len(files)):
    Plot_spec_compare_morphed(files[i])
    Plot_spec_compare_og(files[i])
"""

def Plot_spec_compare_og(og_file, xlims = None, ylims = None):
    '''
    Parameters
    ----------
    file : str
        Provide a file pathway.
    xlims: tuple
        Optionally specify x limits when plotting.
    ylims: tuple
        Optionally specify y limits when plotting.
        If not, it will automatically take the range of the mask of the given parameters in the flux and find the max value and at 5% to it

    Returns
    -------
    Plot of x-ray selected quasar spectra comparing the original spectrum and reconstruction.

    '''
    
 
    data = fits.open(og_file)
    
    #print(data[1].columns) #run to print column names
    
    data_tab = data[1].data
    
    wavelength = data_tab['wave']
    flux = data_tab['flux']
    noise = data_tab['noise']
#   mask = data_tab['mask']
    morph = data_tab['morph']
    recon = data_tab['recon']
    
    
    beta2 = wavelength_to_velocity(0, wavelength)
    
    plt.plot(beta2, flux*morph, color = 'red', label='Original')
    plt.plot(beta2, recon*morph, color = 'blue', label = 'Recon')
    plt.plot(beta2,noise*morph, color='grey', label='Noise')

    plt.ylabel('Flux')
    plt.xlabel('Velocity km/s')
    
    #Check if ylims were defined
    if ylims is not None:
        plt.ylim(ylims)

    #If ylims were not defines, check if xlims were defined if it was we will find the ylims in the range of xlims
    elif xlims is not None:
        xmin, xmax = xlims
        
        #Make sure beta2 is in the range of the velocites we care about
        mask_range = (beta2 >= xmin) & (beta2 <= xmax)
        if np.any(mask_range):
            y_candidates = np.concatenate([
            (flux*morph)[mask_range],
            (recon*morph)[mask_range],
            (noise*morph)[mask_range],
            ])
    
            ymax = np.max(y_candidates)
            plt.ylim(top=ymax + (0.05*ymax))
            
    plt.xlim(xlims)
    plt.legend(loc='upper right')
    plt.title(file[70:])
#   plt.show()
#   plt.close()





