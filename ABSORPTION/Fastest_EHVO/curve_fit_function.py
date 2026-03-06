#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 17:05:36 2024

@author: lilianaflores
"""
import numpy as np
import matplotlib.pyplot as plt
from functools import partial
from scipy.optimize import curve_fit


###############################

def curve_fit_area(xstart, xend, xvalues, flux, error):
    '''
    Select a subsection of data around a specified x-range for curve fitting.

    The function determines index boundaries based on the maximum
    value of `xstart` and `xend`, then extracts a slice of the
    input arrays (`xvalues`, `flux`, `error`). An additional 10 
    indices on both sides is included to extend the fitting region 
    slightly beyond the boundaries.

    Parameters
    ----------
    xstart : int or flt
        Starting x-value that define the region of interest.

    xend : int or flt
        Ending x-values that define the region of interest.

    xvalues : array
        Full x-axis data (e.g., wavelength or velocity).

    flux : array
        Flux (or y-values) corresponding to `xvalues`.

    error : array
        Error values corresponding to `flux`.

    Returns
    -------
    xfit : array
        Subset of `xvalues` within the selected region (including padding).

    yfit : array
        Corresponding subset of `flux`.

    errfit : array
        Corresponding subset of `error`.
    '''
    xmax = np.max(xstart)
    xmin = np.max(xend)
    index_xmax = np.min(np.where(xvalues > xmax))
    index_xmin = np.max(np.where(xvalues < xmin))
    xfit=xvalues[(index_xmin-10):(index_xmax+10)]
    yfit=flux[(index_xmin-10):(index_xmax+10)]
    errfit=error[(index_xmin-10):(index_xmax+10)]
    return (xfit, yfit, errfit)


def tau_v(v,v0,b,tau0):
    return tau0*np.exp(-(v-v0)**2./(b**2.))


##############################################################################################################
def plotting_curvefit_test(d, ion, x, xfit, yfit, fixed, Cf, I0, vdiff, tau_ratio, 
v0, tau0, b, 
v02=None, tau02=None, b2=None, 
v03=None, tau03=None, b3=None,
v04=None, tau04=None, b4=None,
v05=None, tau05=None, b5=None,
v06=None, tau06=None, b6=None,
v07=None, tau07=None, b7=None,
v08=None, tau08=None, b8=None,
v09=None, tau09=None, b9=None,
v10=None, tau10=None, b10=None):
    
    """
    Fit one or more absorption doublets to spectral data and plot the result.

    This function performs non-linear least-squares fitting (via ``curve_fit``)
    to model one or more absorption doublets using an exponential optical depth
    profile defined by ``tau_v``. The number of doublet components is set by `d`
    (currently supported up to 4). 

    The model includes:
        - Covering fraction (Cf)
        - Central velocity (v0, v02, v03, ...)
        - Width (b, b2, b3, ...)
        - Optical depth (tau0, tau02, tau03, ...)
        - Continuum level (I0)
        - Doublet velocity separation (vdiff)
        - Optical depth ratio between doublet lines (tau_ratio)

    The parameter(s) specified in `fixed` are held constant during the fit,
    while the remaining parameters are optimized. After fitting, the function:

        - Plots the combined model profile
        - Plots individual doublet components
        - Prints best-fit values and 1σ uncertainties
        - Returns the fitted parameters

    Parameters
    ----------
    d : int
        Number of doublet components to fit (1–4 supported).

    ion : str
        Name of the ion (used for plot labeling Ex: 'CIV').

    x : array
        Full x-axis array used for plotting the final model.

    xfit : array
        Subset of x-values used for fitting.

    yfit : array
        Flux values corresponding to `xfit`.

    fixed : str
        Comma-separated string specifying which parameters to fix
        during fitting (e.g., 'Cf', 'v0', 'b', 'v0, b', etc.).
        Available options depend on the number of doublets.

    Cf : float
        Covering fraction.

    I0 : float
        Continuum intensity level.

    vdiff : float
        Velocity separation between the two lines of each doublet.

    tau_ratio : float
        Ratio of optical depths between doublet lines.

    v0, v02, v03, v04 : float, more than 1 is optional
        Central velocities of each doublet component.

    tau0, tau02, tau03, tau04 : float, more than 1 is optional
        Optical depths of each doublet component.

    b, b2, b3, b4 : float, more than 1 is optional
        Width parameters of each component.

    Returns
    -------
    tuple
        Fitted parameters in the order:
        (v0, tau0, b, Cf,
         v02, tau02, b2,
         v03, tau03, b3,
         v04, tau04, b4)

    Notes
    -----
    - Requires `tau_v` to compute the optical depth profile.
    - Uses `scipy.optimize.curve_fit` for fitting.
    - Designed for absorption-line spectroscopy applications.
    - Not implemented for more than four doublets.
    """
    
    v = xfit
    if d==1:
        if fixed == 'v0':
            print('***You chose one doublet and to fix v0***')
            print()
            def curve_func( v, tau0, Cf, b, v0, I0, vdiff, tau_ratio):
                return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))) 
            
            (tau0, Cf, b), covar = curve_fit(partial(curve_func, v0=v0, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio), xfit, yfit, p0=[tau0, Cf, b], maxfev=10000, bounds=((0.,0.,0.), (10.,1., np.inf)))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, tau0, Cf, b, v0, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio))), '--', color = 'r')                            
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))),'--', label = str(ion) + ' Doublet', color = 'b')
            plt.legend()
            pcov= np.sqrt(np.diag(covar))

            print()
            print('tau0 =',tau0,"+/-", pcov[0])
            print('Cf =',Cf,"+/-", pcov[1])
            print('b =',b,"+/-", pcov[2])
            
        elif fixed == 'tau0':
            
            print('***You chose one doublet and to fix tau0***')
            print()
            def curve_func(v, v0, Cf, b, tau0, I0, vdiff, tau_ratio):
                return I0 * (1. - Cf) + (Cf * I0 * 
                        (np.exp(-tau_v(v, v0, b, tau0)) * 
                         np.exp(-tau_v(v, (v0 + vdiff), b, tau0 / tau_ratio))))
            
            (v0, Cf, b), covar = curve_fit(partial(curve_func, tau0=tau0, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio),
                                           xfit, yfit, p0=[v0, Cf, b], maxfev=10000,
                                           bounds=((-np.inf, 0., 0.), (0., 1., np.inf)))
            xfit = x
            v = x
            plt.plot(xfit, curve_func(v, v0, Cf, b, tau0, I0, vdiff, tau_ratio),
                     label='Combined Fit', color='purple')
            plt.plot(xfit, I0 * (1. - Cf) + Cf * I0 * 
                     np.exp(-tau_v(v, (v0 + vdiff), b, tau0 / tau_ratio)),
                     '--', color = 'r')
            plt.plot(xfit, I0 * (1. - Cf) + Cf * I0 * 
                     np.exp(-tau_v(v, v0, b, tau0)),
                     '--', label=f'{ion} Doublet', color='b')
            plt.legend()
            pcov = np.sqrt(np.diag(covar))

            print()
            print('v0 =', v0, "+/-", pcov[0])
            print('Cf =', Cf, "+/-", pcov[1])
            print('b =', b, "+/-", pcov[2])
            
        elif fixed == 'b':
            
            print('***You chose one doublet and to fix b***')
            print()
            def curve_func( v, tau0, v0, Cf, b, I0, vdiff, tau_ratio):
                return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))) 
            
            (tau0, v0, Cf), covar = curve_fit(partial(curve_func, b=b, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio), xfit, yfit, p0=[tau0, v0, Cf], maxfev=10000, bounds=((0.,-np.inf,0.), (10.,0., 1.)))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, tau0, v0, Cf, b, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio))), '--', color = 'r')                            
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))),'--', label = str(ion)+ ' Doublet', color = 'b')
            plt.legend()
            pcov= np.sqrt(np.diag(covar))

            print()
            print('tau0 =',tau0,"+/-", pcov[0])
            print('velocity =',v0,"+/-", pcov[1])
            print('Cf =',Cf,"+/-", pcov[2])
            
        elif fixed == 'Cf':
            
            print('***You chose one doublet and to fix Cf***')
            print()
            def curve_func( v, tau0, v0, b, Cf, I0, vdiff, tau_ratio):
                return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))) 
           
            (tau0, v0, b), covar = curve_fit(partial(curve_func, Cf=Cf, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio ), xfit, yfit, p0=[tau0, v0, b], maxfev=10000, bounds=((0.,-np.inf,0), (10.,0.,np.inf)))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, tau0, v0, b, Cf, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio))), '--', color = 'r')                            
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))),'--', label = str(ion)+ ' Doublet', color = 'b')
            plt.legend()
            pcov= np.sqrt(np.diag(covar))

            v0, tau0, b, Cf = v0, tau0, b, Cf
            print()
            print('velocity =',v0,"+/-", pcov[1])
            print('b =',b,"+/-", pcov[2])
            print('tau0 =',tau0,"+/-", pcov[0])
            
        elif fixed == 'v0, tau0':
            
            print('***You chose one doublet and to fix v0 and tau0***')
            print()
            def curve_func(v, Cf, b, v0, tau0, I0, vdiff, tau_ratio):
                return I0 * (1. - Cf) + (Cf * I0 *
                        (np.exp(-tau_v(v, v0, b, tau0)) *
                         np.exp(-tau_v(v, (v0 + vdiff), b, tau0 / tau_ratio))))
            
            (Cf, b), covar = curve_fit(partial(curve_func, v0=v0, tau0=tau0, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio),
                                       xfit, yfit, p0=[Cf, b], maxfev=10000,
                                       bounds=((0., 0.), (1., np.inf)))
            xfit = x
            v = x
            plt.plot(xfit, curve_func(v, Cf, b, v0, tau0, I0, vdiff, tau_ratio),
                     label='Combined Fit', color='purple')
            plt.plot(xfit, I0 * (1. - Cf) + Cf * I0 *
                     np.exp(-tau_v(v, (v0 + vdiff), b, tau0 / tau_ratio)),
                     '--', color = 'r')
            plt.plot(xfit, I0 * (1. - Cf) + Cf * I0 *
                     np.exp(-tau_v(v, v0, b, tau0)),
                     '--', label=f'{ion} Doublet', color='b')
            plt.legend()
            pcov = np.sqrt(np.diag(covar))

            print()
            print('Cf =', Cf, "+/-", pcov[0])
            print('b =', b, "+/-", pcov[1])
            
            
        elif fixed == 'v0, b':
            
            print('***You chose one doublet and to fix v0 and b***')
            print()
            def curve_func( v, tau0, Cf, v0, b, I0, vdiff, tau_ratio):
                return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))) 
            
            (tau0, Cf), covar = curve_fit(partial(curve_func, v0=v0, b=b, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio), xfit, yfit, p0=[tau0, Cf], maxfev=10000, bounds=((0., 0.), (10., 1.)))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, tau0, Cf, v0, b, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple', linewidth=1)
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio))), '--',  color = 'r') 
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))),'--', label = str(ion)+ ' Doublet', color = 'b')
            plt.legend()
            pcov= np.sqrt(np.diag(covar))

            print()
            print('tau0 =',tau0,"+/-", pcov[0])
            print('Cf =',Cf,"+/-", pcov[1])
            
        elif fixed == 'v0, Cf':
            
            print('***You chose one doublet and to fix v0 and Cf***')
            print()
            def curve_func(v, tau0, b, v0, Cf, I0, vdiff, tau_ratio):
                return I0 * (1. - Cf) + (Cf * I0 *
                        (np.exp(-tau_v(v, v0, b, tau0)) *
                         np.exp(-tau_v(v, (v0 + vdiff), b, tau0 / tau_ratio))))
            
            (tau0, b), covar = curve_fit(partial(curve_func, v0=v0, Cf=Cf, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio),
                                         xfit, yfit, p0=[tau0, b], maxfev=10000,
                                         bounds=((0., 0.), (10., np.inf)))
            xfit = x
            v = x
            plt.plot(xfit, curve_func(v, tau0, b, v0, Cf, I0, vdiff, tau_ratio),
                     label='Combined Fit', color='purple')
            plt.plot(xfit, I0 * (1. - Cf) + Cf * I0 *
                     np.exp(-tau_v(v, (v0 + vdiff), b, tau0 / tau_ratio)),
                     '--', color = 'r')
            plt.plot(xfit, I0 * (1. - Cf) + Cf * I0 *
                     np.exp(-tau_v(v, v0, b, tau0)),
                     '--', label=f'{ion} Doublet', color='b')
            plt.legend()
            pcov = np.sqrt(np.diag(covar))

            print()
            print('tau0 =', tau0, "+/-", pcov[0])
            print('b =', b, "+/-", pcov[1])
            
            
        elif fixed == 'tau0, b':
            
            print('***You chose one doublet and to fix tau0 and b***')
            print()
            def curve_func( v, v0, Cf, tau0, b, I0, vdiff, tau_ratio):
                return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))) 
            
            (v0, Cf), covar = curve_fit(partial(curve_func, tau0=tau0, b=b, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio), xfit, yfit, p0=[v0, Cf], maxfev=10000, bounds=((0., 0.), (-np.inf, 1.)))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, v0, Cf, tau0, b, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple', linewidth=1)
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio))), '--',  color = 'r') 
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))),'--', label = str(ion)+ ' Doublet', color = 'b')
            plt.legend()
            pcov= np.sqrt(np.diag(covar))

            print()
            print('v0 =',v0,"+/-", pcov[0])
            print('Cf =',Cf,"+/-", pcov[1])
            
        elif fixed == 'tau0, Cf':
            
            print('***You chose one doublet and to fix tau0 and Cf***')
            print()
            def curve_func( v, v0, b, tau0, Cf, I0, vdiff, tau_ratio):
                return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))) 
            
            (v0, b), covar = curve_fit(partial(curve_func, tau0=tau0, Cf=Cf, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio), xfit, yfit, p0=[v0, b], maxfev=10000, bounds=((0., 0.), (-np.inf, np.inf)))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, v0, b, tau0, Cf, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple', linewidth=1)
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio))), '--',  color = 'r') 
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))),'--', label = str(ion)+ ' Doublet', color = 'b')
            plt.legend()
            pcov= np.sqrt(np.diag(covar))

            print()
            print('v0 =',v0,"+/-", pcov[0])
            print('b =',b,"+/-", pcov[1])
            
        elif fixed == 'b, Cf':
            
            print('***You chose one doublet and to fix b and Cf***')
            print()
            def curve_func( v, v0, tau0, b, Cf, I0, vdiff, tau_ratio):
                return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))) 
            
            (v0, tau0), covar = curve_fit(partial(curve_func, b=b, Cf=Cf, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio), xfit, yfit, p0=[v0, tau0], maxfev=10000, bounds=((0., 0.), (-np.inf, 10.)))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, v0, tau0, b, Cf, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple', linewidth=1)
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio))), '--',  color = 'r') 
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))),'--', label = str(ion)+ ' Doublet', color = 'b')
            plt.legend()
            pcov= np.sqrt(np.diag(covar))

            print()
            print('v0 =',v0,"+/-", pcov[0])
            print('tau0 =',tau0,"+/-", pcov[1])
            
          
        elif fixed == 'v0, tau0, b': 
            
            print('***You chose one doublet and to fix v0, tau0, and b***')
            print()
            def curve_func( v, Cf, v0, tau0, b, I0, vdiff, tau_ratio):
                return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))) 
       
            (tau0), covar = curve_fit(partial(curve_func, v0=v0, tau0=tau0, b=b, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio), xfit, yfit, p0=[Cf], maxfev=10000, bounds=(( 0.), (1.)))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, Cf, v0, tau0, b, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio))), '--',  color = 'r') 
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))),'--', label = str(ion)+ ' Doublet', color = 'b')
            plt.legend()
            pcov= np.sqrt(np.diag(covar))
            
            print()
            print('Cf =',Cf,"+/-", pcov[0])
            
        elif fixed == 'v0, tau0, Cf':
            
            print('***You chose one doublet and to fix v0, tau0, and Cf***')
            print()
            def curve_func( v, b, v0, tau0, Cf, I0, vdiff, tau_ratio):
                return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))) 
       
            (tau0), covar = curve_fit(partial(curve_func, tau0=tau0, Cf=Cf, b=b, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio), xfit, yfit, p0=[b], maxfev=10000, bounds=(( 0.), (np.inf)))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func( v, b, v0, tau0, Cf, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio))), '--',  color = 'r') 
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))),'--', label = str(ion)+ ' Doublet', color = 'b')
            plt.legend()
            pcov= np.sqrt(np.diag(covar))
            
            print()
            print('b =',b,"+/-", pcov[0])
            
        elif fixed == 'Cf, v0, b':
            
            print('***You chose one doublet and to fix Cf, v0, and b***')
            print()
            def curve_func( v, tau0, Cf, v0, b, I0, vdiff, tau_ratio):
                return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))) 
       
            (tau0), covar = curve_fit(partial(curve_func, Cf=Cf, v0=v0, b=b, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio), xfit, yfit, p0=[tau0], maxfev=10000, bounds=(( 0.), (10.)))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, tau0, Cf, v0, b, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio))), '--',  color = 'r') 
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))),'--', label = str(ion)+ ' Doublet', color = 'b')
            plt.legend()
            pcov= np.sqrt(np.diag(covar))

            print()
            print('tau0 =',tau0,"+/-", pcov[0])
            
        elif fixed == 'tau0, Cf, b':

            print('***You chose one doublet and to fix tau0, Cf, and b***')
            print()
            def curve_func( v, v0, tau0, Cf, b, I0, vdiff, tau_ratio):
                return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))) 
       
            (tau0), covar = curve_fit(partial(curve_func, tau0=tau0, Cf=Cf, b=b, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio), xfit, yfit, p0=[v0], maxfev=10000, bounds=(( -np.inf), (0.)))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, v0, tau0, Cf, b, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio))), '--',  color = 'r') 
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))),'--', label = str(ion)+ ' Doublet', color = 'b')
            plt.legend()
            pcov= np.sqrt(np.diag(covar))
            
            print()
            print('v0 =',v0,"+/-", pcov[0])

  
            
    elif d==2:
                
        if fixed == 'v0':
            print('***You chose one doublet and to fix v0***')
            print()
            def curve_func(v, b, b2, tau0, tau02, Cf, v0, v02, I0, vdiff, tau_ratio):
                return I0 * (1. - Cf) + (Cf * I0 *
                        (np.exp(-tau_v(v, v0, b, tau0)) *
                         np.exp(-tau_v(v, v0 + vdiff, b, tau0 / tau_ratio))) *
                        (np.exp(-tau_v(v, v02, b2, tau02)) *
                         np.exp(-tau_v(v, v02 + vdiff, b2, tau02 / tau_ratio))))
            
            (b, b2, tau0, tau02, Cf), covar = curve_fit(partial(curve_func, v0=v0, v02=v02, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio),xfit, yfit, 
            p0=[b, b2, tau0, tau02, Cf], maxfev=20000, bounds=([0., 0., 0., 0., 0.], [np.inf, np.inf, 10., 10., 1.]))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, b, b2, tau0, tau02, Cf, v0, v02, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))), '--', color = 'r')       
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))), '--', label = str(ion)+' Doublet 1', color = 'b')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)))), ':', color = 'r')       
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v02,b2,tau02)))), ':', label = str(ion)+' Doublet 2', color = 'b')    
            plt.legend()

            pcov = np.sqrt(np.diag(covar))
            print()
            print('b =', b, "+/-", pcov[0])
            print('tau0 =', tau0, "+/-", pcov[2])
            print('b2 =', b2, "+/-", pcov[1])
            print('tau02 =', tau02, "+/-", pcov[3])
            print('Cf=', Cf, "+/-", pcov[4])
            
        elif fixed == 'tau0':
            
            print('***You chose one doublet and to fix tau0***')
            print()
            def curve_func(v, b, b2, Cf, v0, v02, tau0, tau02, I0, vdiff, tau_ratio):
                return I0 * (1. - Cf) + (Cf * I0 *
                        (np.exp(-tau_v(v, v0, b, tau0)) *
                         np.exp(-tau_v(v, v0 + vdiff, b, tau0 / tau_ratio))) *
                        (np.exp(-tau_v(v, v02, b2, tau02)) *
                         np.exp(-tau_v(v, v02 + vdiff, b2, tau02 / tau_ratio))))
            
            (b, b2, Cf, v0, v02), covar = curve_fit(partial(curve_func, tau0=tau0, tau02=tau02, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio),xfit, yfit, 
            p0=[b, b2, Cf, v0, v02], maxfev=20000, bounds=([0., 0., 0., -np.inf, -np.inf], [np.inf, np.inf, 1., 0., 0.]))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, b, b2, Cf, v0, v02, tau0, tau02, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))), '--', color = 'r')       
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))), '--', label = str(ion)+' Doublet 1', color = 'b')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)))), ':', color = 'r')       
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v02,b2,tau02)))), ':', label = str(ion)+' Doublet 2', color = 'b')    
            plt.legend()

            pcov = np.sqrt(np.diag(covar))
            print()
            print('b =', b, "+/-", pcov[0])
            print('v0 =', v0, "+/-", pcov[3])
            print('b2 =', b2, "+/-", pcov[1])
            print('v02 =', v02, "+/-", pcov[4])
            print('Cf=', Cf, "+/-", pcov[2])
            
        elif fixed == 'b':
            
            def curve_func(v, Cf, v0, v02, tau0, tau02, b, b2, I0, vdiff, tau_ratio):
                return I0 * (1. - Cf) + (Cf * I0 *
                        (np.exp(-tau_v(v, v0, b, tau0)) *
                         np.exp(-tau_v(v, v0 + vdiff, b, tau0 / tau_ratio))) *
                        (np.exp(-tau_v(v, v02, b2, tau02)) *
                         np.exp(-tau_v(v, v02 + vdiff, b2, tau02 / tau_ratio))))
            
            (Cf, v0, v02, tau0, tau02), covar = curve_fit(partial(curve_func, b=b, b2=b2, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio),xfit, yfit, 
            p0=[Cf, v0, v02, tau0, tau02], maxfev=20000, bounds=([0., -np.inf, -np.inf, 0., 0.], [1., 0., 0., 10., 10.]))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, Cf, v0, v02, tau0, tau02, b, b2, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))), '--', color = 'r')       
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))), '--', label = str(ion)+' Doublet 1', color = 'b')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)))), ':', color = 'r')       
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v02,b2,tau02)))), ':', label = str(ion)+' Doublet 2', color = 'b')    
            plt.legend()

            pcov = np.sqrt(np.diag(covar))
            print()
            print('v0 =', v0, "+/-", pcov[1])
            print('tau0 =', tau0, "+/-", pcov[3])
            print('v02 =', v02, "+/-", pcov[2])
            print('tau02 =', tau02, "+/-", pcov[4])
            print('Cf=', Cf, "+/-", pcov[0])
            
        elif fixed == 'Cf':
            
            print('***You chose two doublets and to fix Cf***')
            print()
            def curve_func(v, v0, b, tau0, v02, b2, tau02, Cf, I0, vdiff, tau_ratio):
                return I0 * (1. - Cf) + (Cf * I0 *
                        (np.exp(-tau_v(v, v0, b, tau0)) *
                         np.exp(-tau_v(v, v0 + vdiff, b, tau0 / tau_ratio))) *
                        (np.exp(-tau_v(v, v02, b2, tau02)) *
                         np.exp(-tau_v(v, v02 + vdiff, b2, tau02 / tau_ratio))))
            
            (v0, b, tau0, v02, b2, tau02), covar = curve_fit(partial(curve_func, Cf=Cf, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio),xfit, yfit, 
            p0=[v0, b, tau0, v02, b2, tau02], maxfev=20000, bounds=([-np.inf, 0., 0., -np.inf, 0., 0.], [0., np.inf, 10., 0., np.inf, 10.]))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, v0, b, tau0, v02, b2, tau02, Cf, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))), '--', color = 'r')       
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))), '--', label = str(ion)+' Doublet 1', color = 'b')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)))), ':', color = 'r')       
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v02,b2,tau02)))), ':', label = str(ion)+' Doublet 2', color = 'b')    
            plt.legend()

            pcov = np.sqrt(np.diag(covar))
            print()
            print('v0 =', v0, "+/-", pcov[0])
            print('b =', b, "+/-", pcov[1])
            print('tau0 =', tau0, "+/-", pcov[2])
            print('v02 =', v02, "+/-", pcov[3])
            print('b2 =', b2, "+/-", pcov[4])
            print('tau02 =', tau02, "+/-", pcov[5])
            
        elif fixed == 'v0, tau0':
            
            print('***You chose one doublet and to fix v0 and tau0***')
            print()
            
            
            
        elif fixed == 'v0, b':
            
            print('***You chose one doublet and to fix v0 and b***')
            print()
            
            
        elif fixed == 'v0, Cf':
            
            print('***You chose one doublet and to fix b and Cf***') #has altered limits of possible b values to fit
            print()
            def curve_func( v, tau0, tau02, b, b2, v0, v02, Cf, I0, vdiff, tau_ratio):
                return I0 * (1. - Cf) + (Cf * I0 *
                        (np.exp(-tau_v(v, v0, b, tau0)) *
                         np.exp(-tau_v(v, v0 + vdiff, b, tau0 / tau_ratio))) *
                        (np.exp(-tau_v(v, v02, b2, tau02)) *
                         np.exp(-tau_v(v, v02 + vdiff, b2, tau02 / tau_ratio))))
            
            (tau0, tau02, b, b2), covar = curve_fit(partial(curve_func, v0=v0, v02=v02, Cf=Cf, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio), xfit, yfit, p0=[tau0, tau02, b, b2], maxfev=10000, bounds=((0., 0.,0.,0.), (10., 10., np.inf, np.inf)))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, tau0, tau02, b, b2, v0, v02, Cf, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple', linewidth=1)
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio))), '--', color = 'r') 
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))),'--', label = str(ion)+ ' Doublet 1', color = 'b')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)))), ':', color = 'r')       
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v02,b2,tau02)))), ':', label = str(ion)+' Doublet 2', color = 'b')    
            
            
            plt.legend()
            pcov= np.sqrt(np.diag(covar))

            print()
            print('b =',b,"+/-", pcov[2])
            print('b2 =',b2,"+/-", pcov[3])
            #print('**has altered limits of possible b values to fit**')

            print('tau0 =',tau0,"+/-", pcov[0])
            print('tau02 =',tau02,"+/-", pcov[1])
            
        elif fixed == 'tau0, b':
            
            print('***You chose one doublet and to fix tau0 and b***')
            print()
            
            
        elif fixed == 'tau0, Cf':
            
            print('***You chose one doublet and to fix tau0 and Cf***')
            print()
            
            
        elif fixed == 'b, Cf':
            
            print('***You chose one doublet and to fix b and Cf***')
            print()
            def curve_func( v, v0, v02, tau0, tau02, b, b2, Cf, I0, vdiff, tau_ratio):
                return I0 * (1. - Cf) + (Cf * I0 *
                        (np.exp(-tau_v(v, v0, b, tau0)) *
                         np.exp(-tau_v(v, v0 + vdiff, b, tau0 / tau_ratio))) *
                        (np.exp(-tau_v(v, v02, b2, tau02)) *
                         np.exp(-tau_v(v, v02 + vdiff, b2, tau02 / tau_ratio))))
            
            (v0, v02, tau0, tau02), covar = curve_fit(partial(curve_func, b=b, b2=b2, Cf=Cf, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio), xfit, yfit, p0=[v0, v02, tau0, tau02], maxfev=10000, bounds=((-np.inf, -np.inf,0., 0.), (0.,0.,10., 10.)))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, v0, v02, tau0, tau02, b, b2, Cf, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple', linewidth=1)
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio))), '--', color = 'r') 
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))),'--', label = str(ion)+ ' Doublet 1', color = 'b')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)))), ':', color = 'r')       
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v02,b2,tau02)))), ':', label = str(ion)+' Doublet 2', color = 'b')    
            
            
            plt.legend()
            pcov= np.sqrt(np.diag(covar))

            print()
            print('v0 =',v0,"+/-", pcov[0])
            print('v02 =',v02,"+/-", pcov[1])

            print('tau0 =',tau0,"+/-", pcov[2])
            print('tau02 =',tau02,"+/-", pcov[3])

            
            
          
        elif fixed == 'v0, tau0, b': 
            
            print('***You chose one doublet and to fix v0, tau0, and b***')
            print()
            
            
        elif fixed == 'v0, tau0, Cf':
            
            print('***You chose one doublet and to fix v0, tau0, and Cf***')
            print()
            
            
        elif fixed == 'Cf, v0, b':
            print()
            print('***You chose two doublet and to fix Cf, v0, and b***')
            print()
            def curve_func(v, tau0, tau02, v0, b, v02, b2, Cf, I0, vdiff, tau_ratio):
                return I0 * (1. - Cf) + (Cf * I0 *
                        (np.exp(-tau_v(v, v0, b, tau0)) *
                         np.exp(-tau_v(v, v0 + vdiff, b, tau0 / tau_ratio))) *
                        (np.exp(-tau_v(v, v02, b2, tau02)) *
                         np.exp(-tau_v(v, v02 + vdiff, b2, tau02 / tau_ratio))))
            
            (tau0, tau02), covar = curve_fit(partial(curve_func, v0=v0, b=b, v02=v02, b2=b2, Cf=Cf, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio),xfit, yfit, 
            p0=[tau0, tau02], maxfev=20000, bounds=([0., 0.], [10.,10.]))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, tau0, tau02, v0, b, v02, b2, Cf, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))), '--', color = 'r')       
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))), '--', label = str(ion)+' Doublet 1', color = 'b')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)))), ':', color = 'r')       
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v02,b2,tau02)))), ':', label = str(ion)+' Doublet 2', color = 'b')    
            plt.legend()
            pcov= np.sqrt(np.diag(covar))
            
            print()
            print('tau0 =',tau0,"+/-", pcov[0])
            print('tau02 =',tau02,"+/-", pcov[1])
            
            
        elif fixed == 'tau0, Cf, b':

            print('***You chose one doublet and to fix tau0, Cf, and b***')
            print()
            
            
        elif fixed == 'b, Cf':
            print('You chose two doublets and to fix b, b2, Cf. Parameters requiring guesses: v0, v02, tau0, tau02. Constant parameters: I0, vdiff, tau_ratio.')
            print()
            def curve_func( v, v0, v02, tau0, tau02, b, b2, Cf, I0, vdiff, tau_ratio): 
                return I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)) * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)) * np.exp(-tau_v(v,v02,b2,tau02)) * np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)))) 
            
            
            (v0, v02, tau0, tau02), covar = curve_fit(partial(curve_func, b=b, b2=b2, Cf=Cf, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio), xfit, yfit, p0=[v0, v02, tau0, tau02], maxfev=10000, bounds=((-np.inf,-np.inf,0.,0.), (0., 0.,10.,10.)))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, v0, v02, tau0, tau02, b, b2, Cf, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))), '--', color = 'r')       
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))), '--', label = str(ion)+' Doublet 1', color = 'b')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)))), ':', color = 'r')       
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v02,b2,tau02)))), ':', label = str(ion)+' Doublet 2', color = 'b')    
            plt.legend()
            pcov= np.sqrt(np.diag(covar))

            print()
            print('velocity =',v0,"+/-", pcov[0])
            print('velocity2 =',v02,"+/-", pcov[1])
            print('tau0 =',tau0,"+/-", pcov[2])
            print('tau02 =',tau02,"+/-", pcov[3])
        
    elif d == 3:
        if fixed == 'Cf':
                
            print('***You chose two doublets and to fix Cf***')
            print()
            def curve_func(v, v0, b, tau0, v02, b2, tau02, v03, b3, tau03, Cf, I0, vdiff, tau_ratio):
                return I0 * (1. - Cf) + (Cf * I0 *
                        (np.exp(-tau_v(v, v0, b, tau0)) *
                             np.exp(-tau_v(v, v0 + vdiff, b, tau0 / tau_ratio))) *
                            (np.exp(-tau_v(v, v02, b2, tau02)) *
                             np.exp(-tau_v(v, v02 + vdiff, b2, tau02 / tau_ratio)))*
                             (np.exp(-tau_v(v, v03, b3, tau03)) *
                              np.exp(-tau_v(v, v03 + vdiff, b3, tau03 / tau_ratio))))
                
            (v0, b, tau0, v02, b2, tau02, v03, b3, tau03), covar = curve_fit(partial(curve_func, Cf=Cf, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio),xfit, yfit, 
            p0=[v0, b, tau0, v02, b2, tau02, v03, b3, tau03], maxfev=20000, bounds=([-np.inf, 0., 0., -np.inf, 0., 0.,  -np.inf, 0., 0.], [0., np.inf, 10., 0., np.inf, 10.,  0., np.inf, 10.]))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, v0, b, tau0, v02, b2, tau02, v03, b3, tau03, Cf, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))), '--', color = 'r')       
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))), '--', label = str(ion)+' Doublet 1', color = 'b')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)))), ':', color = 'r')       
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v02,b2,tau02)))), ':', label = str(ion)+' Doublet 2', color = 'b')    
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v03+vdiff),b3,tau03/tau_ratio)))), '-.', color = 'r')       
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v03,b3,tau03)))), '-.', label = str(ion)+' Doublet 3', color = 'b')    
            pcov= np.sqrt(np.diag(covar))

            plt.legend()
            
            print()
            print('tau0 =',tau0,"+/-", pcov[2])
            print('tau02 =',tau02,"+/-", pcov[5])
            print('tau03 =',tau03,"+/-", pcov[8])
            print()
            print('v0 =',v0,"+/-", pcov[0])
            print('v02 =',v02,"+/-", pcov[3])
            print('v03 =',v03,"+/-", pcov[6])
            print()
            print('b =',b,"+/-", pcov[1])
            print('b2 =',b2,"+/-", pcov[4])
            print('b3 =',b3,"+/-", pcov[7])
                
        elif fixed == 'Cf, v0, b':
            print()
            print('***You chose one doublet and to fix Cf, v0, and b***')
            print()
            def curve_func(v, tau0, tau02, tau03, v0, b, v02, b2, v03, b3, Cf, I0, vdiff, tau_ratio):
                return I0 * (1. - Cf) + (Cf * I0 *
                        (np.exp(-tau_v(v, v0, b, tau0)) *
                         np.exp(-tau_v(v, v0 + vdiff, b, tau0 / tau_ratio))) *
                        (np.exp(-tau_v(v, v02, b2, tau02)) *
                         np.exp(-tau_v(v, v02 + vdiff, b2, tau02 / tau_ratio)))*
                         (np.exp(-tau_v(v, v03, b3, tau03)) *
                          np.exp(-tau_v(v, v03 + vdiff, b3, tau03 / tau_ratio))))
            
            (tau0, tau02, tau03), covar = curve_fit(partial(curve_func, v0=v0, b=b, v02=v02, b2=b2, v03=v03, b3=b3, Cf=Cf, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio),xfit, yfit, 
            p0=[tau0, tau02, tau03], maxfev=20000, bounds=([0., 0., 0.], [10.,10.,10.]))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, tau0, tau02, tau03, v0, b, v02, b2, v03, b3, Cf, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))), '--', color = 'r')       
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))), '--', label = str(ion)+' Doublet 1', color = 'b')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)))), ':', color = 'r')       
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v02,b2,tau02)))), ':', label = str(ion)+' Doublet 2', color = 'b')    
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v03+vdiff),b3,tau03/tau_ratio)))), '-.', color = 'r')       
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v03,b3,tau03)))), '-.', label = str(ion)+' Doublet 3', color = 'b')    
            plt.legend()
            pcov= np.sqrt(np.diag(covar))
            
            print()
            print('tau0 =',tau0,"+/-", pcov[0])
            print('tau02 =',tau02,"+/-", pcov[1])
            print('tau03 =',tau03,"+/-", pcov[2])
            
            
            
        elif fixed == 'Cf, v0': #In progress
            print('***You chose one doublet and to fix b and Cf***')
            print()
            def curve_func( v, tau0, tau02, tau03, b, b2, b3, v0, v02, v03, Cf, I0, vdiff, tau_ratio):
                return I0 * (1. - Cf) + (Cf * I0 *
                        (np.exp(-tau_v(v, v0, b, tau0)) *
                         np.exp(-tau_v(v, v0 + vdiff, b, tau0 / tau_ratio))) *
                        (np.exp(-tau_v(v, v02, b2, tau02)) *
                         np.exp(-tau_v(v, v02 + vdiff, b2, tau02 / tau_ratio))) *
                         (np.exp(-tau_v(v, v02, b2, tau02)) *
                          np.exp(-tau_v(v, v02 + vdiff, b2, tau02 / tau_ratio))))
            
            (tau0, tau02, tau03, b, b2, b3), covar = curve_fit(partial(curve_func, v0=v0, v02=v02, v03=v03, Cf=Cf, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio), xfit, yfit, p0=[tau0, tau02, tau03, b, b2, b3], maxfev=10000, bounds=((0., 0.,0.,0.,0., 0.), (10., 10., 10., np.inf, np.inf, np.inf)))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func( v, tau0, tau02, tau03, b, b2, b3, v0, v02, v03, Cf, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple', linewidth=1)
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio))), '--', color = 'r') 
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))),'--', label = str(ion)+ ' Doublet 1', color = 'b')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)))), ':', color = 'r')       
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v02,b2,tau02)))), ':', label = str(ion)+' Doublet 2', color = 'b')    
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v03+vdiff),b3,tau03/tau_ratio)))), '-.', color = 'r')       
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v03,b3,tau03)))), '-.', label = str(ion)+' Doublet 3', color = 'b')    
            
            
            plt.legend()
            pcov= np.sqrt(np.diag(covar))
    
            print()
            print('b =',b,"+/-", pcov[2])
            print('b2 =',b2,"+/-", pcov[3])
            print('b3 =',b3,"+/-", pcov[3])

    
            print('tau0 =',tau0,"+/-", pcov[0])
            print('tau02 =',tau02,"+/-", pcov[1])
            print('tau03 =',tau03,"+/-", pcov[2])
            
            
    elif d == 4:
        if fixed == 'Cf':
                    
            print('***You chose two doublets and to fix Cf***')
            print()
            def curve_func(v, v0, b, tau0, v02, b2, tau02, v03, b3, tau03, v04, b4, tau04, Cf, I0, vdiff, tau_ratio):
                return I0 * (1. - Cf) + (Cf * I0 *
                            (np.exp(-tau_v(v, v0, b, tau0)) *
                                 np.exp(-tau_v(v, v0 + vdiff, b, tau0 / tau_ratio))) *
                                (np.exp(-tau_v(v, v02, b2, tau02)) *
                                 np.exp(-tau_v(v, v02 + vdiff, b2, tau02 / tau_ratio)))*
                                 (np.exp(-tau_v(v, v03, b3, tau03)) *
                                  np.exp(-tau_v(v, v03 + vdiff, b3, tau03 / tau_ratio)))*
                                  (np.exp(-tau_v(v, v04, b4, tau04)) *
                                   np.exp(-tau_v(v, v04 + vdiff, b4, tau04 / tau_ratio))))
                    
            (v0, b, tau0, v02, b2, tau02, v03, b3, tau03, v04, b4, tau04), covar = curve_fit(partial(curve_func, Cf=Cf, I0=I0, vdiff=vdiff, tau_ratio=tau_ratio),xfit, yfit, 
            p0=[v0, b, tau0, v02, b2, tau02, v03, b3, tau03, v04, b4, tau04], maxfev=20000, bounds=([-np.inf, 0., 0., -np.inf, 0., 0.,  -np.inf, 0., 0.,  -np.inf, 0., 0.], [0., np.inf, 10., 0., np.inf, 10.,  0., np.inf, 10.,  0., np.inf, 10.]))
            xfit = x
            v = x
            plt.plot(xfit, (curve_func(v, v0, b, tau0, v02, b2, tau02, v03, b3, tau03, v04, b4, tau04, Cf, I0, vdiff, tau_ratio)), label='Combined Fit', color='purple')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v0+vdiff),b,tau0/tau_ratio)))), '--', color = 'r')       
            plt.plot(xfit, I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v0,b,tau0)))), '--', label = str(ion)+' Doublet 1', color = 'b')
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v02+vdiff),b2,tau02/tau_ratio)))), ':', color = 'r')       
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v02,b2,tau02)))), ':', label = str(ion)+' Doublet 2', color = 'b')    
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v03+vdiff),b3,tau03/tau_ratio)))), '-.', color = 'r')       
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v03,b3,tau03)))), '-.', label = str(ion)+' Doublet 3', color = 'b')    
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,(v04+vdiff),b4,tau04/tau_ratio)))), '-.', color = 'r')       
            plt.plot(xfit,I0 *(1. - Cf) + (Cf * I0 * (np.exp(-tau_v(v,v04,b4,tau04)))), '-.', label = str(ion)+' Doublet 4', color = 'b')    
                
            pcov= np.sqrt(np.diag(covar))

            plt.legend()
                
            print()
            print('tau0 =',tau0,"+/-", pcov[2])
            print('tau02 =',tau02,"+/-", pcov[5])
            print('tau03 =',tau03,"+/-", pcov[8])
            print('tau04 =',tau04,"+/-", pcov[11])
            print()
            print('v0 =',v0,"+/-", pcov[0])
            print('v02 =',v02,"+/-", pcov[3])
            print('v03 =',v03,"+/-", pcov[6])
            print('v04 =',v04,"+/-", pcov[9])

            print()
            print('b =',b,"+/-", pcov[1])
            print('b2 =',b2,"+/-", pcov[4])
            print('b3 =',b3,"+/-", pcov[7])
            print('b4 =',b4,"+/-", pcov[10])

    else:
        print('The curve fit function has not been finished to accomadate a number of doublets above 4.')

        

    return v0, tau0, b, Cf, v02, tau02, b2, v03, b3, tau03, v04, b4, tau04

############################################################################################################################






