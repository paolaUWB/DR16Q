#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 24 18:00:40 2025

@author: lilianaflores
"""

from astropy.convolution import Gaussian1DKernel, convolve
#from https://github.com/sbailey/empca
from empca import SavitzkyGolay as SG
import numpy as np
from matplotlib import pyplot as plt
import os
from utility_functions import read_spectra


g0 = Gaussian1DKernel(stddev=3)
g1 = Gaussian1DKernel(stddev=5)
g2 = Gaussian1DKernel(stddev=8)
SGsmooth1=SG(width=30)
SGsmooth2=SG(width=50)

'''
These functions were provided by Joseph Choi. To perform sigma clipping you must call ChoiClip function. 
The ChoiClip function then uses the GaussClip function. 
'''


def GaussClip(WL,INT,ERR,FLAG,std=4,err=1.5):
    '''
    

    Parameters
    ----------
    WL : array
        Wavelength values.
    INT : array
        Intensity (flux values).
    ERR : array
        Error values.
    FLAG : array
        Array of 0 or 1 values indicating whether they are flagged (from previous iteration). 
        FLAG == 1: Not flagged
        FLAG == 0: Flagged
    std : float, optional
        Standard deviation to decide how far bellow the smoothed spectra points will be flagged. 
        The default is 4.
    err : float, optional
        Used within the function as a multiplier to the errors. The ERR*err is used to decide what FLAG
        values copied as clip_logic wil be replaced with 0. This means those points are not flagged as 
        they are not impacted by the total establish error (ERR*err). The default is 1.5.

    Returns
    -------
    clip_logic : array
        Array of 1s and 0s to indicate points that are clipped.

    '''
    
    
    g_kernel = Gaussian1DKernel(stddev=std) #std is used for the Gaussian kernel to smooth the spectra while the err value is the sigma value that decides the cuttoff for values to flag
    temp_INT=np.interp(WL,WL[FLAG==1],INT[FLAG==1]) #placing an interpolated data points?
    temp=np.where(np.logical_and(WL<1215.7,convolve(temp_INT,g0)-temp_INT>ERR*err))[0]
    #print(np.where(np.logical_and(WL<1215.7,convolve(temp_INT,g0)-temp_INT>ERR*err)))
    #what is convolve? ^ temp values are where wavelength values are in Lya forest and convolution of the 
    #interpolated data points and established 1D gaussian kernel minus the interpolated data points is greater
    #than the array of ERR (Error) values multipled by a provided factor err. 

    clip_logic=np.copy(FLAG)
    clip_logic[temp]=0
    return clip_logic




def ChoiClip(WL,INT,ERR,std=4,err=1.5,maxiter=10):
    '''
    Sigma clipping function. This is the main function you will call to perform any Lya line clippings. 
    The funtion starts with a smoothed spectrum made with the Gaussian kernel of set width. any data points 
    that fall defined sigma (std) below the smoothed spectrum are flagged. The function then makes another smooth 
    spectrum version where those flagged points are are replaced with interpolated values based on the surrounding 
    unflagged points. The the cycle restarts until there are no more data points being flagged OR the maximum iteration
    limit (maxiter) is met.


    Parameters
    ----------
    WL : array
        Wavelength values.
    INT : array
        Intensity (flux values).
    ERR : array
        Error values.
    std : float, optional
        Standard deviation to decide how far bellow the smoothed spectra points will be flagged. The default is 4.
    err : TYPE, optional
        Error is used within the GaussClip function called within this one (see GaussClip function for more details). 
        The default is 1.5.
    maxiter : TYPE, optional
        Maximum iteration limit stops the function process if the limit is met and data points are still being flagged. 
        The default is 10.

    Returns
    -------
    clip_loop: array
        Array of float values representing the new spectrum without the clipped data points and with interpolated points
        to replace removed points.
    histoire_loop: array
        Array of the total data points clipped after each iteration

    '''
    
    start_loop=np.ones(len(WL))
    totale_loop=len(np.where(WL<1215.7)[0]) #Starting a wavelength where Lya forest starts?
    histoire_loop=[]
    for i in range(maxiter):
        if i==0:
            loopdeloop=np.ones(len(WL))
        elif i>0:
            loopdeloop=clip_loop
        clip_loop=GaussClip(WL,INT,ERR,loopdeloop,std,err)
        histoire_loop.append(sum(start_loop)-sum(clip_loop)) #Tracking the number of data points clipped in each iteration
        print('clipped: ',histoire_loop[-1],' ',round(histoire_loop[-1]/totale_loop*100,2),' %')
        if len(histoire_loop)>1:
            if histoire_loop[-1]-histoire_loop[-2]==0: #Why does it break if two iterations clipped the same amount of data points?
                break
        if i==maxiter-1:
            'Reached Maximum Iteration'
    return clip_loop,histoire_loop



