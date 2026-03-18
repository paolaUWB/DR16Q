#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 13:58:56 2026

@author: lilianaflores
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def linear_model(x, m, b):
    return m * x + b


def bisector(x, y, ax = None, plot=True):
    """
    Fits a line to (x, y),
    draws a perpendicular line through the data center,
    and returns masks for points above and below the line.
    """
    # fitting the data points of sample provided
    (m, b), _ = curve_fit(linear_model, x, y)

    # getting the mid point of the sample to place bisector
    x_mid = np.mean(x)
    y_mid = np.mean(y)

    # getting the slope of the perpendicular bisector
    m_perp = -1 / m

    # function for perpendicular line
    def perp_line(x_val):
        return m_perp * (x_val - x_mid) + y_mid

    # split the sample above and below the perpendicular bisector
    above_mask = y > perp_line(x)
    below_mask = y < perp_line(x)

    # optionally plot the fit line and perpendicular bisecting line
    if plot:
        
        if plot and ax is not None:

            # fit line
            x_fit = np.linspace(min(x), max(x), 200)
            ax.plot(x_fit, linear_model(x_fit, m, b), color='k', linewidth=1, zorder=1000)
    
            # perpendicular bisector
            x_range = np.linspace(x_mid - 1, x_mid + 1, 100)
            y_perp = perp_line(x_range)
            ax.plot(x_range, y_perp, 'r--', zorder=1000)
        
        else:
            # plot fit line
            x_fit = np.linspace(min(x), max(x), 200)
            plt.plot(x_fit, linear_model(x_fit, m, b), color='k', linewidth=1, label='Fit', zorder=1000)
    
            # plot perpendicular bisector
            x_range = np.linspace(x_mid - 1, x_mid + 1, 100)
            y_perp = perp_line(x_range)
            plt.plot(x_range, y_perp, 'r--', label='Perpendicular Bisector', zorder=1000)
    

    return above_mask, below_mask