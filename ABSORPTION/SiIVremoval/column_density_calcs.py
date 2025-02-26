#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 23:54:32 2025

@author: lilianaflores
"""
import numpy as np
from scipy.integrate import quad
from astropy import constants as K
from astropy import units as u
from math import sqrt, pi
from curve_fit_functions import tau_v




def column_density(v0, b, tau0, Cf, xfit):

    ##### column density calcs################
    # Calculate Nion:
    me = 9.109e-28     # electron mass (g)
    qe = 4.8032e-10          #electron charge (esu)
    c = 2.99792458e10        #light speed (cm/s)
    osca = 0.190             # oscillator strength stronger line
    oscb = 0.952e-1          #       ”            weaker line
   
    l0a_cm = 1548.195 * 1e-8  # cm  # A and it needs to be cm --  that is why I added the 1.e8
    l0b_cm = 1550.770 * 1e-8  # cm
    k_N=(me*c)/(pi *(qe**2))
    integral_1a = quad(tau_v,np.min(xfit),np.max(xfit),args=(v0,b ,tau0))
    N1a = (k_N/(osca*l0a_cm)) * integral_1a[0] * 1e5 # 1e5 to convert km to cm
    integral_1b = quad(tau_v,np.min(xfit),np.max(xfit),args=(v0,b,tau0))
    N1b = (k_N/(oscb*l0b_cm)) * integral_1b[0] * 1e5 # 1e5 to convert km to cm
    #print(’log Nion = ’, np.log10(N1a+N1b))  #Can I add both lines? I don’t think so. Only one line (use right line)
    #print('method in paper -- log Nion = ', np.log10(N1b))  #Can I add both lines? I don’t think so. Only one line (use right line)
    e2_me_c = (K.e.gauss**2 / (K.m_e*K.c)).to(u.cm**2 / u.s).value
    b_cm_s = b * 1e5
    logN1a = np.log10(tau0 * b_cm_s / (sqrt(pi) * e2_me_c * osca * l0a_cm))  # this was copied from a Github -- why the sqrt(pi)?!?!
    logN1b = np.log10(tau0 * b_cm_s / (sqrt(pi) * e2_me_c * oscb * l0b_cm))
    #print(‘second method: log Nion = ’,np.log10(10**logN1a+10**logN1b))
    #print('second method from Github: log Nion = ',np.log10(10**logN1b))
    N_CIV1 = N1b
    N_CIV2 = 10**np.log10(10**logN1a+10**logN1b) # Not using this, both lines
    N_CIV2 = 10**np.log10(10**logN1b) # Only using red line
    AC=8.43 #Asplund+09
    NHC=10.**(12.-AC)
    # Ionization fractions f(Mi) at the top of their line in Figure A1 Hamann+11:
    NCCIV=1./10.**(-1.25) # 0.4 before
    NH_C1=NHC*NCCIV*N_CIV1
    NH_C2=NHC*NCCIV*N_CIV2
    #print('logNH derived from C (first method) =',np.log10(NH_C1))
    #print('logNH derived from C (second method from Github) =',np.log10(NH_C2))
    ###########################################
    logNion_m1 = np.log10(N1b)
    logNion_m2 = np.log10(10**logN1b)
    logNH_m1 = np.log10(NH_C1)
    logNH_m2 = np.log10(NH_C2)

    return v0, Cf, b, tau0, logNion_m1, logNion_m2, logNH_m1, logNH_m2
