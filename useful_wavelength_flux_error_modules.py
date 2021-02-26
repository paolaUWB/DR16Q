#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 17:10:41 2021

@author: wendygarcia mikelcharles
"""
import numpy as np 
from data_types import ColumnIndexes, RangesData, PointData

column_index = ColumnIndexes(0, 1, 2)

###  Function: returns one point of wavelength, flux and error based on a range of values in a good defined range
def wavelength_flux_error_for_points(starting_point: float, ending_point: float, z: float, spectra_data) -> PointData: 
    wavelength_column = spectra_data[:, column_index.wavelength]

    wavelength_observed_start = (z + 1) * starting_point
    wavelength_observed_end = (z + 1) * ending_point

    point_from = np.max(np.where(wavelength_column < wavelength_observed_start))
    point_to = np.min(np.where(wavelength_column > wavelength_observed_end))

    wavelength = spectra_data[point_from:point_to, column_index.wavelength]
    flux = spectra_data[point_from:point_to, column_index.flux] 
    error = spectra_data[point_from:point_to, column_index.error] 
    
   point = PointData(
        np.average(point_ranges.wavelength),
        np.median(point_ranges.flux),
        np.median(point_ranges.error))

   
    return PointData(wavelength, flux, error)

###  Function: returns a range of of wavelength, flux and error defined by starting and ending points
def wavelength_flux_error_in_range(starting_point: float, ending_point: float, z: float, spectra_data) -> RangesData:
    wavelength_column = spectra_data[:, column_index.wavelength]

    wavelength_observed_from = (z + 1) * starting_point
    wavelength_observed_to = (z + 1) * ending_point

    wavelength_lower_limit = np.where(wavelength_column > wavelength_observed_from)
    wavelength_upper_limit = np.where(wavelength_column < wavelength_observed_to)
    
    wavelength = spectra_data[np.min(wavelength_lower_limit[column_index.wavelength]):np.max(wavelength_upper_limit[column_index.wavelength]), column_index.wavelength]
    flux = spectra_data[np.min(wavelength_lower_limit[column_index.wavelength]): np.max(wavelength_upper_limit[column_index.wavelength]), column_index.flux]
    error = spectra_data[np.min(wavelength_lower_limit[column_index.wavelength]): np.max(wavelength_upper_limit[column_index.wavelength]), column_index.error]
    
    return RangesData(wavelength, flux, error)

### then calculates the snr with this and "flag" cases with SNR lower a value will be excluded later. 
def calculate_snr(wavelength, z: float, WAVELENGTH_FOR_SNR: range, error_normalized):
    wavelengths_for_snr_lower = np.where (wavelength/(z + 1.) < WAVELENGTH_FOR_SNR.start)
    wavelengths_for_snr_upper = np.where (wavelength/(z + 1.) > WAVELENGTH_FOR_SNR.end)
    snr_mean_in_ehvo = round(np.mean(1./error_normalized[np.max(wavelengths_for_snr_lower[0]):np.min(wavelengths_for_snr_upper)]), 5)
    return snr_mean_in_ehvo 