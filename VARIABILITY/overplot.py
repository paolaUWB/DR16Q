# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 12:00:53 2022

@author: bunge
"""

import os
import sys
import numpy as np
from scipy.optimize import curve_fit
# from spec_pca_module import plot_sdss_lines
sys.path.insert(0, os.getcwd() + '/../' + 'DR16Q') # changes the directory to the DR16Q --> all paths after this will need to be written as if this was in the top level of the DR16Q
from useful_wavelength_flux_error_modules import wavelength_flux_error_for_points, wavelength_flux_error_in_range, calculate_snr
# from draw_figures import powerlaw, draw_dynamic, draw_dynamic_points, draw_original_figure, draw_normalized_figure
# import normalization_module as norm
import draw_figures as norm
from utility_functions import read_file
import matplotlib.pyplot as plt
# from spec_pca_module import plot_sdss_lines


config_path = os.getcwd() + '/DR16_sorted_norm.csv' #"C:/Users/Dakota/Documents/GitHub/DR16Q/DR9_sorted_norm.csv"
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



def plot_orig_spectra(directory, lines=False, error=False, wavemin=1200, wavemax=1500):
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
            # if lines:
                # plot_sdss_lines(wavemin, wavemax)
            if error:
                plt.plot(wave, err, 'grey')
            plt.title(name)
            plt.xlabel("Restframe Wavelength (A)")
            plt.ylabel("Normalized Flux Density")
            plt.legend()
            # plt.ylim((0,2))
            plt.show()
            plt.clf()
            


def plot_norm_spectra(directory, lines=False, error=False, wavemin=1200, wavemax=1500):
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
            # if lines:
                # plot_sdss_lines(wavemin, wavemax)
            if error:
                plt.plot(wave, err, 'grey')
            plt.title(name)
            plt.xlabel("Restframe Wavelength (A)")
            plt.ylabel("Normalized Flux Density")
            plt.legend()
            # plt.ylim((0,2))
            plt.show()
            plt.clf()
            
plot_spectra(os.getcwd() + '/VARIABILITY/DATA_VARIABILITY/', lines=False)