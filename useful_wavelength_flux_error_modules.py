#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 17:10:41 2021

@author: wendygarcia mikelcharles
"""
import numpy as np 
from data_types import ColumnIndexes, RangesData, PointData
from scipy import signal

column_index = ColumnIndexes(0, 1, 2)


def wavelength_flux_error_for_points(starting_point: float, ending_point: float, z: float, spectra_data) -> PointData: 
    """Returns one point of wavelength, flux and error based on a range of values in a good defined range.

    Uses the red shift to find the observed wavelenghts, and between those two wavelengths records all of the 
    wavelengths, all of the flux, and all the error. Using the observed wavelengths finds the average wavelength, 
    median flux and median error for the right, left and middle point.

    Positional Input Parameter:
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

    Keyword Input Parameters:
        None.

    Returns:
        Point. Tuple.
        Ex: wavelength, flux, and error would be replaced with data points.
        [(wavelength, flux, error),
        (wavelength, flux, error),
        (wavelength, flux, error)]
    """
    wavelength_column = spectra_data[:, column_index.wavelength]

    wavelength_observed_start = (z + 1) * starting_point
    wavelength_observed_end = (z + 1) * ending_point

    point_from = np.max(np.where(wavelength_column < wavelength_observed_start))
    point_to = np.min(np.where(wavelength_column > wavelength_observed_end))

    wavelength = spectra_data[point_from:point_to, column_index.wavelength]
    flux = spectra_data[point_from:point_to, column_index.flux] 
    error = spectra_data[point_from:point_to, column_index.error] 
    
    point = PointData(
        np.average(wavelength),
        np.median(flux),
        np.median(error))

   
    return point

def wavelength_flux_error_in_range(starting_point: float, ending_point: float, z: float, spectra_data) -> RangesData:
    """Function: returns a range of a wavelength, flux and error defined by starting and ending points
        TBD. #############################################

    Positional Input Parameter:
        starting_point : float
            Uses the range defined by the following variable: WAVELENGTH_RESTFRAME.start,     
        ending_point: float
            Also, uses the range defined by the following variable: WAVELENGTH_RESTFRAME.end
        z: float
            Values from the data base of the redshift, DR16Q (for now..)
        spectra_data: list
            Current spectra data from files, DR16Q (for now...)

    Keyword Input Parameters:
        None.

    Returns:
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
 
def calculate_snr(wavelength, z: float, WAVELENGTH_FOR_SNR: range, error_normalized):
    """ Calculates the snr (signal to noise ratio). [Want a high SNR value]

    Positional Input Parameter:
        wavelength: array
            Comes from RangesData().    
        z: float
            Values from the data base of the redshift, DR16Q (for now..)
        WAVELENGTH_FOR_SNR: range
            A range defined in the beginning of the normalization code. Can 
            be changed by user.
        error_normalized: array
            The error from RangesData() divided by the power law. 

    Keyword Input Parameters:
        None.

    Returns:
        snr_mean_in_ehvo(). Float or Int.
    """
    wavelengths_for_snr_lower = np.where (wavelength/(z + 1.) < WAVELENGTH_FOR_SNR.start)
    wavelengths_for_snr_upper = np.where (wavelength/(z + 1.) > WAVELENGTH_FOR_SNR.end)
    snr_mean_in_ehvo = round(np.mean(1./error_normalized[np.max(wavelengths_for_snr_lower[0]):np.min(wavelengths_for_snr_upper)]), 5)
    return snr_mean_in_ehvo 


def smooth(norm_flux, box_size):
    y_smooth = signal.savgol_filter(norm_flux,box_size,2)  #linear
    return y_smooth