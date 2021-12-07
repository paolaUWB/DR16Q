"""
=============
depth_graph_test.py
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

def depth_figure(spectra_index, velocity, flux_normalized, error, savefile_name, spectra_name, redshift, snr, max_peak, smooth_values):
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

    if snr < 20:
        new_flux = smooth(flux_normalized, smooth_values[0])
        number = smooth_values[0]
        plt.plot(velocity, new_flux, color = 'b', label = 'smooth value = ' + str(number))
    elif snr < 30:
        new_flux = smooth(flux_normalized, smooth_values[1])
        number = smooth_values[1]
        plt.plot(velocity, new_flux, color = 'g', label = 'smooth value = ' + str(number))
    else: 
        new_flux = smooth(flux_normalized, smooth_values[1])
        number = smooth_values[2]
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

def depth_testing_calculation(normalized_flux,vmins_index, vmaxs_index):

    final_depth = np.min(normalized_flux[vmaxs_index:vmins_index])
    final_depth_index = np.where(normalized_flux == final_depth)
    final_depth_index = final_depth_index.astype(int)

    number1 = final_depth_index - 3
    number2 = final_depth_index - 2
    number3 = final_depth_index - 1

    number5 = final_depth_index + 3
    number6 = final_depth_index + 2
    number7 = final_depth_index + 1

    seven_value_apprx_1 = normalized_flux[number1]
    seven_value_apprx_2 = normalized_flux[number2]
    seven_value_apprx_3 = normalized_flux[number3]

    seven_value_apprx_5 = normalized_flux[number5]
    seven_value_apprx_6 = normalized_flux[number6]
    seven_value_apprx_7 = normalized_flux[number7]

    average = (seven_value_apprx_1 + seven_value_apprx_2 + seven_value_apprx_3 + final_depth + seven_value_apprx_5 + seven_value_apprx_6 + seven_value_apprx_7) / 7

    final_depth = round((1. - average), 2)
    final_depth_individual = []     
    final_depth_individual.append(average)

    return final_depth_individual