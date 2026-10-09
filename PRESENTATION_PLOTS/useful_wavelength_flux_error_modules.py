#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 17:10:41 2021

@author: wendygarcia mikelcharles
"""
import numpy as np 
from data_types import ColumnIndexes, RangesData, PointData
from scipy import signal

######################################### shpinx ######################################### 
"""
useful_wavelength_flux_error_modules
====================================
Useful functions. 
"""
#############################################################################################

## COLUMN INDEX FOR WAVELENGTH, FLUX, ERROR IN SPECTRA (DR) FILES 
column_index = ColumnIndexes(0, 1, 2)
#wavelength: list,
def wavelength_flux_error_for_points( starting_point: float, ending_point: float, z: float, spectra_data) -> PointData: 
    """Returns one point of wavelength, flux and error for a range of wavelengths.
    ** OBSERVED WAVELENGTH NEEDED - WILL NOT WORK PROPERLY FOR RESTFRAME
    Records the wavelengths, flux, and error within the range of wavelengths provided. Using the observed wavelengths, it finds the average 
    wavelength, median flux and median error for each of the anchor points.

    Parameters
    ----------
    wavelength: list
        
    starting_point : float
        Uses the range defined by the following variables: WAVELENGTH_RESTFRAME_FOR_LEFT_POINT, 
        WAVELENGTH_RESTFRAME_FOR_RIGHT_POINT, WAVELENGTH_RESTFRAME_FOR_MIDDLE_POINT.
    ending_point: float
        Also, uses the range defined by the following variables: WAVELENGTH_RESTFRAME_FOR_LEFT_POINT
        WAVELENGTH_RESTFRAME_FOR_RIGHT_POINT, WAVELENGTH_RESTFRAME_FOR_MIDDLE_POINT.
    z: float
        Values from the data base of the redshift, DR16Q (for now..)
    spectra_data: list
        Current spectra data from files, DR16Q (for now...)

    Returns
    -------
    PointData.

    Examples
    --------
    wavelength, flux, and error would be replaced with data points.
    [(wavelength, flux, error),
    (wavelength, flux, error),
    (wavelength, flux, error)]
    """
    
    ###### CHECK OVER DOCUMENTATION FOR THIS ######
    
    wavelength_column = spectra_data[:, column_index.wavelength]
    point_from = np.max(np.where(wavelength_column <= starting_point))
    point_to = np.min(np.where(wavelength_column >= ending_point))

    wavelength = spectra_data[point_from:point_to, column_index.wavelength]
    flux = spectra_data[point_from:point_to, column_index.flux] 
    error = spectra_data[point_from:point_to, column_index.error] 
    
    point = PointData(
        np.average(wavelength),
        np.median(flux),
        np.median(error))

    return point

def wavelength_flux_error_in_range(starting_point: float, ending_point: float, z: float, spectra_data) -> RangesData:
    """Returns a range of a wavelength, flux and error defined by starting and ending points.

    Parameters
    ----------
    starting_point : float
        Uses the range defined by the following variable: WAVELENGTH_RESTFRAME.start,     
    ending_point: float
        Also, uses the range defined by the following variable: WAVELENGTH_RESTFRAME.end
    z: float
        Values from the data base of the redshift, DR16Q (for now..)
    spectra_data: list
        Current spectra data from files, DR16Q (for now...)

    Returns
    -------
    RangesData.

    Examples
    --------
    RangesData(wavelength, flux, error). List of tuples of arrays.
    RangesData creates a list of tuples, within each tuple it stores arrays for the ranges of 
    wavelength, flux, and error values.
    """
    wavelength_column = spectra_data[:, column_index.wavelength]

    wavelength_observed_from = (z + 1) * starting_point
    wavelength_observed_to = (z + 1) * ending_point

    wavelength_lower_limit = np.where(wavelength_column > wavelength_observed_from)
    wavelength_upper_limit = np.where(wavelength_column < wavelength_observed_to)
    
    wavelength = spectra_data[np.min(wavelength_lower_limit[column_index.wavelength]):np.max(wavelength_upper_limit[column_index.wavelength]), column_index.wavelength]
    flux = spectra_data[np.min(wavelength_lower_limit[column_index.wavelength]): np.max(wavelength_upper_limit[column_index.wavelength]), column_index.flux]
    error = spectra_data[np.min(wavelength_lower_limit[column_index.wavelength]): np.max(wavelength_upper_limit[column_index.wavelength]), column_index.error]
    
    return RangesData(wavelength, flux, error)
 
#def calculate_snr(wavelength, z: float, WAVELENGTH_FOR_SNR: range, error_normalized):
def calculate_snr(wavelength, z: float, WAVELENGTH_FOR_SNR: range, flux, error):
    """ Calculates the snr (signal to noise ratio). [Want a high SNR value].

    Parameters
    ----------
    wavelength: array
        Comes from RangesData().    
    z: float
        Values from the data base of the redshift, DR9Q (for now..).
    WAVELENGTH_FOR_SNR: range
        A range defined in the beginning of the normalization code. Can 
        be changed by user.
    error_normalized: array
        The error from RangesData() divided by the power law. 

    Returns
    -------
    Float or Int.
        snr_mean_in_ehvo().
    """
    wavelengths_for_snr_lower = np.where (wavelength/(z + 1.) < WAVELENGTH_FOR_SNR.start)
    wavelengths_for_snr_upper = np.where (wavelength/(z + 1.) > WAVELENGTH_FOR_SNR.end)
    snr_mean_in_ehvo = round(np.mean(flux[np.max(wavelengths_for_snr_lower):np.min(wavelengths_for_snr_upper)]/error[np.max(wavelengths_for_snr_lower):np.min(wavelengths_for_snr_upper)]), 5)
    #snr_mean_in_ehvo = round(np.mean(1./error_normalized[np.max(wavelengths_for_snr_lower[0]):np.min(wavelengths_for_snr_upper)]), 5)
    return snr_mean_in_ehvo 