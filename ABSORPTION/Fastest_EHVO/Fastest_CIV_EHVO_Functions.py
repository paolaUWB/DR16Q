#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 13:11:31 2025

@author: lilianaflores
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import f


def round_custom(value):
    '''
    Round a numeric value to a single significant leading digit
    using custom rules based on the first digit.

    The function converts the input to an integer. 
    It then rounds the number depending on its first digit:

    - If the first digit is greater than 1, the value is rounded 
      to one significant digit (e.g., 345 -> 300).
    - If the first digit is equal to 1, the value is rounded to 
      two significant digits (e.g., 145 -> 150).
    - If the first digit is 0 (or cannot be processed), the value 
      is returned unchanged and a warning is printed.
    - If the value is 0 or NaN, NaN is returned.

    Parameters
    ----------
    value : int or float
        Numeric value that needs to be rounded.

    Returns
    -------
    value_new : int or float
        Rounded value according to the custom significant-digit rule.
        Returns np.nan if the input is 0 or NaN.

    '''
    try:
        value = int(value)
    except:
        value = value
    n = len(str(value)) #getting the length of the error
    d1 = str(value)[0] #getting the first digit of the error
    
    if value == 0 or np.isnan(value):
        print(f'Warning: cannot round NaN or None value ({value})')
        return np.nan
    
    if int(d1) > 1:
        value_new = round(value, -(n-1))
        
    elif int(d1) == 1:
        value_new = round(value, -(n-2))
        
    else:
        value_new = value
        print(f'Error value {value} not rounded')
        
    return value_new

def round_by_error(value, error):
    '''
    Round a value based on the magnitude of its associated error.

    The function first rounds the provided error using ``round_custom``.
    The number of significant digits in the rounded error is then used
    to determine the decimal position to which the value should be rounded.
    The value is rounded to match the precision implied by the error.

    Parameters
    ----------
    value : int or float
        The numeric value to be rounded.

    error : int or float
        The associated uncertainty/error used to determine the
        rounding precision.

    Returns
    -------
    value_new : float
        The input value rounded to the same significant digit
        position as the rounded error.
    '''
    round_err = round_custom(error)
    n = len(str(round_err))
    sigfigs = len(str(round_err).strip('0'))
    
    value_new = np.round(value, -(n-sigfigs))
    
    return value_new

     
def get_values_str(string):
    '''
    Extract numeric values from a string representation of a list.

    The function removes square brackets and commas, then splits
    the string into individual elements based on whitespace.
    
    This may be used when reading in the result csv from absorption.py 
    which may produce multiple values in one table cell such as [70000,40000].

    Parameters
    ----------
    string : str
        String containing numeric values, typically formatted like
        a list (e.g., "[1, 2, 3]" or "1 2 3").

    Returns
    -------
    values : list of str
        List of string elements representing the extracted values.
    '''
    string = string.strip('[]').replace(',', ' ')
    values = string.split()
    return values

def get_max(str_values):
    '''
    Compute the maximum absolute value from a string of numbers.

    The input string is made into numeric values using
    ``get_values_str``, converted to floats in a list, and the maximum
    of their absolute values is returned.

    Parameters
    ----------
    str_values : str
        String containing numeric values (e.g., "[1, -2, 3]").

    Returns
    -------
    max_value : float
        Maximum absolute value among the numbers.
    '''
    values = get_values_str(str_values)
        
    numbers = []
    for x in values:
        numbers.append(float(x))
    
    max_value = max(np.abs(numbers))
    return max_value

def get_min(str_values):
    '''
    Compute the minimum absolute value from a string of numbers.

    The input string is made into numeric values using
    ``get_values_str``, converted to floats in a list, and the minimum
    of their absolute values is returned.

    Parameters
    ----------
    str_values : str
        String containing numeric values (e.g., "[1, -2, 3]").

    Returns
    -------
    min_value : float
        Minimum absolute value among the numbers.
    '''
    values = get_values_str(str_values)
        
    numbers = []
    for x in values:
        numbers.append(float(x))
    
    min_value = min(np.abs(numbers))
    return min_value

def get_sum(str_values):
    '''
    Compute the sum of the value from a string of numbers.

    The input string is made into numeric values using
    ``get_values_str``, converted to floats in a list, and the sum
    of their absolute values is returned.

    Parameters
    ----------
    str_values : str
        String containing numeric values (e.g., "[1, -2, 3]").

    Returns
    -------
    sum_values : float
        Sum of the absolute values of the numbers.
    '''
    values = get_values_str(str_values)
        
    numbers = []
    for x in values:
        numbers.append(float(x))
    
    sum_values = sum(np.abs(numbers))
    return sum_values       


def plot_masked_regions(x,y,xfit, want_mask_beyond_xfit = True):
    '''
    Plot masked regions of a spectrum relative to a fitted x-range.

    This function identifies regions of `x` that are not included in
    the fitted x-values (`xfit`) and plots those regions in grey.
    Masked regions are constructed to excluded segments (grey overplot)
    that would be visually continuous when plotted between unmasked and 
    masked points. This avoids a false display of data points used.

    By default, values outside the range of `xfit` are also treated
    as masked regions. This behavior can be disabled.

    Parameters
    ----------
    x : array-like
        Full x-axis values (e.g., wavelength or velocity).

    y : array-like
        Corresponding y-axis values (e.g., flux values).

    xfit : array-like
        Subset of x-values used for fitting. Points in `x` that are
        not approximately equal (via `np.isclose`) to any value in
        `xfit` are considered masked.

    want_mask_beyond_xfit : bool, optional
        If True (default), regions outside the minimum and maximum
        of `xfit` are plotted as masked regions.
        If False, those regions are ignored (set to NaN).

    Returns
    -------
    None
        The function does not return any values. It adds a line to the
        current matplotlib axes representing masked regions.
    '''
    norm_flux = y
    start = np.min(xfit)
    end = np.max(xfit)
        
    xm = [] #x values masked
    new_flux = [] #new flux aligning x masked values
    
    in_mask = False #tracking whether in a masked region
    prev_flux = np.nan #tracking previous flux value
    prev_x = np.nan #tracking previous x value
    
    for i in range(len(x)):
        
        if x[i] < start or x[i] > end: #assigning nan to x and flux values outside of fit range
            if want_mask_beyond_xfit == True: #edited to make masking outside the xfit range optional.Default is to mask outside the xfit range.
                repx = x[i]
                repy = norm_flux[i]
                
            else:
                repx = np.nan
                repy = np.nan
        
            xm.append(repx)
            new_flux.append(repy)
            continue
    
        in_xfit = np.any(np.isclose(xfit, x[i])) #checks if x[i] is close to any value in xfit to determine whether in xfit (floating point numbers->close to)
    
        if not in_xfit: # Entering a masked region
            if not in_mask: #then add last value to plot full masked region
                xm.append(prev_x)
                new_flux.append(prev_flux)
            xm.append(x[i])
            new_flux.append(norm_flux[i])
            in_mask = True #are in masked region add current points to plot

        else: #in xfit, leaving masked region
            if in_mask: #if was previously just in masked region then append current points to plot complete masked region
                xm.append(x[i])
                new_flux.append(norm_flux[i])
                in_mask = False #esablish that no longer in masked region
            else: # not in masked region add nan
                xm.append(np.nan)
                new_flux.append(np.nan)
    
            prev_flux = norm_flux[i] 
            prev_x = x[i]
    
    xm = np.array(xm)
    new_flux = np.array(new_flux)
    
    plt.plot(xm, new_flux, color='grey', label='Masked Regions')
    
    

def f_test(yfit, param_new, fitting_new, param_old, fitting_old, alpha):
    """Preforms f-test to quantify whether a fitting is improved. The idea is that you compare the Chi^2 of fittings as you add numbers of doublets.
    For example you would compare a fitting done with one doublet (or gaussian, etc.) to a fitting done with 2 doublets.

    Parameters
    ----------
    yfit: array
        Array of y data values used in fitting. (Likely flux values)
    param_new: int
        The number of parameters that are being fitted for in the new fit.
    fitting_new: Array
        The resulting y values of your new fit. You would provide this by using the function you used for the new fit with the parameter values produced through fitting.
    
    param_old: int
        The number of parameters that are being fitted for in the old fit.
    fitting_old: Array
        The resulting y values of your old fit. You would provide this by using the function you used for the old fit with the parameter values produced through fitting.
    
    
    Returns
    -------
    F_x: float
        The F value of the f-test.
        
    v: int
        The degrees of freedom (DOF) of the new fit
    """
    
    num_points = np.size(yfit) #Finding the number of data points used in the fitting
    
    deltap = param_new - param_old #change in the number of parameters being fitted for
    
    v = num_points - param_new #degrees of freedom: number of data points in the fitted area - number of parameters fitted for in the new fit
    
    chi2_new = np.sum((yfit - fitting_new)**2)/np.sum((yfit - np.mean(yfit))**2) #chi^2 of new fitting
    
    chi2_old = np.sum((yfit - fitting_old)**2)/np.sum((yfit - np.mean(yfit))**2) #chi^2 of old fitting

    Rchi2_new = chi2_new/v #new reduced chi^2

    delta_chi2 = chi2_old - chi2_new #change of chi^2 between fittings
    
    F_x = (delta_chi2/deltap)/Rchi2_new #f-test: (change of chi^2/change in # of parameters fitted for)/new reduced chi^2
    
    F_crit = f.ppf(0.999, deltap, v)
    
    
    return F_x, v, F_crit



