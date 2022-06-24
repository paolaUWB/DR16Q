# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 12:00:53 2022

@author: bunge
"""

import os
import numpy as np
import normalization_module as norm
import matplotlib.pyplot as plt
# from spec_pca_module import plot_sdss_lines
from scipy.optimize import curve_fit
from utility_functions import read_file

config_path = "C:/Users/Dakota/Documents/GitHub/DR16Q/DR9_sorted_norm.csv"
fontsize = 20
figsize = (12,6)

# Configure parameters
plt.rcParams.update({'font.size': fontsize, 'figure.figsize': figsize})

# Default tick label size
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2

plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['axes.linewidth'] = 2


def load_config(config_path, data_name):
    """
    Load redshifts and snr from DR config file.

    Parameters
    ----------
    config_path : str
        Location of config file.
    data_names : arr
        Array of specnames in the form 'plate-mjd-fiber'.

    Returns
    -------
    parent_zs : arr
        Redshifts in the order provided by data_names.
    data_names : arr
        SNRs in the order provided by data_names.

    """
    #load data from config
    
    parent_zs_config, parent_snrs_config, parent_names_config_raw = read_file(config_path)
    parent_names_config = []
    for file in parent_names_config_raw:
        specname = file[:19]
        parent_names_config.append(specname)
    parent_names_config = np.asarray(parent_names_config_raw)
    parent_zs_config = np.asarray(parent_zs_config)
    parent_snrs_config = np.asarray(parent_snrs_config)
    
    #match redshifts and snrs to config file
    
    parent_z = []
    parent_snr = []
    for i in range(len(parent_names_config)):
        if parent_names_config[i] == data_name:
            parent_z = parent_zs_config[i]
            print(parent_z)
            parent_snr = parent_snrs_config[i]
    
    
    return(parent_z, parent_snr)

def normalize_spec(wave, flux, err, z=0):
    """
    Normalize a spectrum. Based on the normalization module. Just for purposes
    of this script - not for general use. Anchor points, delta, and initial
    powerlaw parameters are hardcoded.
    
    Parameters
    ----------
    wave : arr
        Wavelenght bins.
    flux : arr
        Flux values.
    err : arr
        Error values.
    z : float, optional
        Redshift (0 for composite).

    Returns
    -------
    bf : float
        PL norm.
    cf : float
        PL slope.
    normed_flux : arr
        Normalized spectrum.

    """
    # user_anchors = [1285,1425,1680,1820] #anchor points (restframe) for fitting PL
    user_anchors = [1285,1425,1680]
    user_delta = 5 #distance (A) around points for fitting
    b = 1200 # initial PL norm
    c = -0.5 # initial PL slope
    
    current_spectra_data = np.array((wave, flux, err)).T
    anchor_pts = norm.dynamic_find_anchor_points(current_spectra_data, z, \
                                                 user_anchors, user_delta, verbose=False)

    power_law_wave = []
    power_law_flux = []
    
    for point in anchor_pts:
                
        power_law_wave.append(point[0])
        power_law_flux.append(point[1])
    
    
    pars, covar = curve_fit(norm.powerlaw, power_law_wave, power_law_flux, p0=[b, c], maxfev=10000)
        
    bf, cf = pars[0], pars[1]
    
    normed_flux = flux / norm.powerlaw(wave, bf, cf)
    
    return(bf, cf, normed_flux)

def plot_spectra(directory, lines=False, error=False, wavemin=1200, wavemax=1500):
    """
    
    Parameters
    ----------
    spectra : str
        Directory in which the files are contained, or in which the
        directories are contained.
    lines : bool, optional
        Whether or not you want to plot common emission
        lines. The default is False.

    Returns
    -------
    None.

    """
   
    wave = []
    flux = []
    err = []
    name = []
    label = []
    colors = ['b', 'darkgreen', 'purple', 'cyan', 'magenta']
    for f in os.listdir(directory):
        if f.startswith('J'):
            z = np.loadtxt(directory + f + '/info.txt', skiprows=1,delimiter = ',',dtype=str)[1].astype(np.float64)
            for i, file in enumerate(os.listdir(directory + f)):
                if file.startswith('spec'):
                    data = np.loadtxt(directory + f + '/' + file, dtype=np.float64).T
                    wave = data[0]
                    flux = data[1]
                    err = data[2]
                    name = f
                    try:
                        wave = wave / (1 + z)
                    except UnboundLocalError:
                        pass
                    if file[-3:] == 'txt':
                        try:
                            bf, cf, norm_flux = normalize_spec(wave, flux, err)
                            flux = norm_flux
                        except ValueError:
                            pass
                        label = 'SDSS-I'
                    if file[-3:] == 'dr9':
                        label = 'SDSS-II'
                    if file[-3:] == 'r16':
                        label = 'SDSS-III/IV'
                    # else:
                    #     label = '?'
                    mask = np.where((wave > wavemin) & (wave < wavemax))
                    plt.plot(wave[mask], flux[mask], label=label, alpha=0.75, color = colors[i])
                    plt.axhline(1, color = 'r')
            if lines:
                plot_sdss_lines(wavemin, wavemax)
            if error:
                plt.plot(wave, err, 'grey')
            plt.title(name)
            plt.xlabel("Restframe Wavelength (A)")
            plt.ylabel("Normalized Flux Density")
            plt.legend()
            # plt.ylim((0,2))
            plt.show()
            plt.clf()
            
plot_spectra('DATA_VARIABILITY/', lines=False)