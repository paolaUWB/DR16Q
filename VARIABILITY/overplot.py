# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 12:00:53 2022

@author: bunge
"""

import os
import sys
import numpy as np
#import utility_functions
#from scipy.optimize import curve_fit
#from spec_pca_module import plot_sdss_lines
sys.path.insert(0, os.getcwd() + '/../' + 'DR16Q') # changes the directory to the DR16Q --> all paths after this will need to be written as if this was in the top level of the DR16Q
#from useful_wavelength_flux_error_modules import wavelength_flux_error_for_points, wavelength_flux_error_in_range, calculate_snr
#from draw_figures import powerlaw, draw_dynamic, draw_dynamic_points, draw_original_figure, draw_normalized_figure
#import normalization_module as norm
#import draw_figures as norm
#from utilities import read_file
from utility_functions import read_spectra, read_list_spectra
from data_types import ColumnIndexes
import matplotlib.pyplot as plt
#from spec_pca_module import plot_sdss_lines


config_path = os.getcwd() + '/DR16_EHVO_sorted_norm.csv' #"C:/Users/Dakota/Documents/GitHub/DR16Q/DR9_sorted_norm.csv"
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

def plot_norm_spectra(directory, lines=False, error=False): # Normalized Spectra
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

    # wave = []
    # flux = []
    # err = []
    # name = []
    #label = []
    #colors = ['b', 'darkgreen', 'purple', 'cyan', 'magenta']
    #file = ['directory']

   # data_columns = np.loadtxt('/Users/rachelfulda/Desktop/Galaxies/DR16Q/DR16Q_EHVO/NORM_DR16Q_EHVO/spec-11634-58484-0334norm.dr16', dtype=str).T
   # print(data_columns[1])
     
    # z = np.loadtxt(directory, skiprows=1, delimiter = ',',dtype=str)[1].astype(np.float64)
    
    #if directory.startswith('spec') and directory.endswith('norm.dr16' or 'norm.txt'):
    #data_columns = np.loadtxt('/Users/rachelfulda/Desktop/Galaxies/DR16Q/DR16Q_EHVO/NORM_DR16Q_EHVO/spec-11634-58484-0334norm.dr16', dtype=float)
    #path = '/Users/rachelfulda/Desktop/Galaxies/DR16Q/DR16Q_EHVO/NORM_DR16Q_EHVO/spec-4485-55836-0092norm.dr16'
    
    norm_direc = '/Users/rachelfulda/Desktop/Galaxies/DR16Q/DR16Q_EHVO/NORM_DR16Q_EHVO'
    norm_spectrum_file_name = '/spec-4485-55836-0092norm.dr16'
    norm_spectra_data = np.loadtxt(norm_direc + norm_spectrum_file_name)
    wave, flux, err = read_spectra(norm_spectra_data)
    #data = np.loadtxt(directory + '/' + file, dtype=np.float64).T
    #wave = data_columns[0]
    #flux = data_columns[1]
    #err =  data_columns[2]
     
    # name = file 
    # try:
    #     bf, cf, norm_flux = file(wave, flux, err)
    #     flux = norm_flux
    # except ValueError:
    #     pass
    #label = 'SDSS-I'
    # if file[-3:] == 'dr9':
    #     label = 'SDSS-II'
    # if file[-3:] == 'dr16':
    #     label = 'SDSS-III/IV'
         
    #ax1 = plt.subplots(1,1)

    # ay1 = fig.add_subplot(1, 1, 1)
    # plt.title(spectrum)
    plt.plot(wave, flux,'k-')
    plt.plot(wave, err, 'r--') 
    plt.xlabel(r"Observed Wavelength [$\rm \AA$]")
    plt.ylabel(r"Normalized Flux")
    
    plt.xlim(7000,10000)
    plt.ylim(0,1.5)

    #plt.plot([wavelength_observe1,wavelength_observe2],[1,1],'r--')
    #color = ['xkcd:shocking pink', 'black', 'xkcd:purpleish blue']
    #color = ['xkcd:shocking pink', 'xkcd:azure', 'blue', 'xkcd:purpleish blue', 'xkcd:slate']
    #color = ['red', 'green', 'blue', 'orange', 'purple']
    
    #     else:
    #         label = '?'
    #         mask = np.where((wave > wavemin) & (wave < wavemax))
    #         plt.plot(wave[mask], flux[mask], label=label, alpha=0.75, color = colors[i])
    #         plt.axhline(1, color = 'r')
    # if lines:
    #     plot_sdss_lines(wavemin, wavemax)
    # if error:
    #     plt.plot(wave, err, 'grey')
    #ax1.title(name)
    # ax1.set_xlabel("Restframe Wavelength (A)")
    # ax1.set_ylabel("Normalized Flux Density")
    # ax1.legend()
    # plt.ylim((0,2))
    #plt.subplots()
    # plt.show()
    #plt.clf()
        
plot_norm_spectra('/Users/rachelfulda/Desktop/Galaxies/DR16Q/DR16Q_EHVO/NORM_DR16Q_EHVO/spec-11634-58484-0334norm.dr16', lines=False) #, wavemin=1200, wavemax=1600

###
# def plot_orig_spectra(directory, lines=False, error=False, wavemin=1200, wavemax=1500):
#     """
    
#     Parameters
#     ----------
#     spectra : str
#         Directory in which the files are contained, or in which the
#         directories are contained.
#     lines : bool, optional
#         Whether or not you want to plot common emission
#         lines. The default is False.

#     Returns
#     -------
#     None.

#     """
   
#     wave = []
#     flux = []
#     err = []
#     name = []
#     label = []
#     colors = ['b', 'darkgreen', 'purple', 'cyan', 'magenta']
#     for f in os.listdir(directory):
#         if f.startswith('J'):
#             z = np.loadtxt(directory + f + '/info.txt', skiprows=1,delimiter = ',',dtype=str)[1].astype(np.float64)
#             for i, file in enumerate(os.listdir(directory + f)):
#                 if file.startswith('spec') and file.endswith('norm.dr16' or 'norm.txt'):
#                     data = np.loadtxt(file, dtype=np.float64).T
#                     wave = data[0]
#                     flux = data[1]
#                     err = data[2]
#                     name = file
#                     try:
#                         wave = wave / (1 + z)
#                     except UnboundLocalError:
#                         pass
#                     if file[-3:] == 'txt':
#                         try:
#                             bf, cf, norm_flux = normalize_spec(wave, flux, err)
#                             flux = norm_flux
#                         except ValueError:
#                             pass
#                         label = 'SDSS-I'
#                     if file[-3:] == 'dr9':
#                         label = 'SDSS-II'
#                     if file[-3:] == 'r16':
#                         label = 'SDSS-III/IV'
#                     # else:
#                     #     label = '?'
#                     mask = np.where((wave > wavemin) & (wave < wavemax))
#                     plt.plot(wave[mask], flux[mask], label=label, alpha=0.75, color = colors[i])
#                     plt.axhline(1, color = 'r')
#             # if lines:
#                 # plot_sdss_lines(wavemin, wavemax)
#             if error:
#                 plt.plot(wave, err, 'grey')
#             plt.title(file)
#             plt.xlabel("Restframe Wavelength (A)")
#             plt.ylabel("Normalized Flux Density")
#             plt.legend()
#             # plt.ylim((0,2))
#             plt.show()
#             plt.clf()
            
# plot_spectra(os.getcwd() + '/VARIABILITY/DATA_VARIABILITY/', lines=False)

###