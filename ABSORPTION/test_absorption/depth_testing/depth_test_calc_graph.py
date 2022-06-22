"""
=============
depth_test_calc_graph.py
=============

@author Wendy Garcia Naranjo  
based on code prepared by Abdul Khatri and Paola Rodriguez Hidalgo, Mikel Charles, Nathnael Kahassai, Michael Parker

Short description:
    Writing code to draw 10 different spectra with 10 different boxcar values (in the smooth function) on one graph.
    Will use depth_test.csv.
"""

########################################### IMPORTS ############################################################################
from matplotlib import pyplot as plt
import numpy as np
from data_types import Range
from abs_function_module import smooth
################################################################################################################################

def depth_figure(spectra_index, velocity, flux_normalized, error, savefile_name, spectra_name, redshift, snr, max_peak, new_flux, number):
    """ Makes a flux vs velocity graph, that also has the error vs velocity on the same graph. Has text that identifies 
    what the graph number is, what the spectra name is, the signal to noise ratio, and the redshift value used.
    
    Parameters
    ----------
    spectra_index: int or list
        The index number of the graph you are plotting.

    velocity: array or list
        The velocity values to plot on the x-axis.

    flux_normalized: array or list
        The flux values that will be plotted on the y-axis.

    error: array or list
        The error values to be plotted.

    savefile_name: string
        The name that you want the plot to be saved as.

    spectra_name: string
        The full name of the spectra. 

    redshift: int or list or array
        The redshift value used for calculations.

    snr: int or list or array
        The signal to noise ratio value.

    max_peak: int or array
        The max peak value which is used to scale the y-axis. 

    smooth_values: list
        The list of the vlaues for the boxcar size in the smooth function.

    Returns
    -------
    None.
    """

    plt.plot(velocity, flux_normalized, color = 'k', linewidth = 0.7) # original spectra

    plt.plot(velocity, new_flux, color = 'r', label = 'smooth value = ' + str(number))

    #plt.plot(velocity, flux_normalized, color = 'k', linestyle = '-', linewidth = 0.7) # original spectra over top

    plt.plot(velocity, error, color = 'grey')
    plt.xlabel("Velocity (km/s)")
    plt.ylabel("Normalized Flux")
    plt.legend(loc = 'lower right', prop = {"size":5})
    plt.xlim(-70000, 0)
    snr = round(snr, 2)
    min_peak = -0.1
    plt.title(str(spectra_index) + ' tot: ' + str(spectra_name) + ', z=' + str(redshift) + ' snr=' + str(snr))
    plt.axhline(y = 0.9, color='r', linestyle = '--')    
    plt.axhline(y = 1.0)
    plt.ylim(min_peak, max_peak + (max_peak / 4))
    savefile_name.savefig()
    plt.close()

def new_smooth_normalzied_flux_value(flux_normalized, snr, smooth_values):
    if snr < 20:
        new_smooth_normalzied_flux_value = smooth(flux_normalized, smooth_values[0])
        number = smooth_values[0]
    elif snr < 30:
        new_smooth_normalzied_flux_value = smooth(flux_normalized, smooth_values[1])
        number = smooth_values[1]
    else: 
        new_smooth_normalzied_flux_value = smooth(flux_normalized, smooth_values[2])
        number = smooth_values[2]

    return new_smooth_normalzied_flux_value, number

def depth_testing_calculation(normalized_flux,vmins_index, vmaxs_index):

    # using 7 pixel average method to calculate depth value ###########################################
    final_depth = np.min(normalized_flux[vmaxs_index:vmins_index])
    final_depth_index = list(np.where(normalized_flux == final_depth))
    
    for i in range(len(final_depth_index)):

        number_minus_range = []
        appx_minus_values = []
        for j in range(2): 
            number_minus = final_depth_index[i] - (j + 1)
            number_minus_range.append(number_minus)

            appx_minus = normalized_flux[number_minus]
            appx_minus_values.append(appx_minus)

        number_plus_range = []
        appx_plus_values = []
        for k in reversed(range(0, 2)):
            number_plus = final_depth_index[i] + (k + 1)
            number_plus_range.append(number_plus)

            appx_plus = normalized_flux[number_plus]
            appx_plus_values.append(appx_plus)

        average = (sum(appx_minus_values) + sum(appx_plus_values) + final_depth) / (1 + len(appx_minus_values) + len(appx_plus_values))

        final_depth = np.round((1. - average), 2)
        final_depth_individual = []     
        final_depth_individual.append(final_depth)

    return final_depth_individual