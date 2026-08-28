#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 21:36:12 2026

@author: jmaharaj
"""

import matplotlib.pyplot as plt
from lines import (
    CIVll, SiIVll, NVll, OVIll,
    CIVllred, CIVllblue,
    SiIVllred, SiIVllblue,
    CII_emitted, OI_emitted,
    NVllred, NVllblue,
    Lya,
    OVIllred, OVIllblue,
    Lyb
)

def plot_emission_labels(zem, topemlabel,                            ###function for plotting emission line labels on the spectrum
                         CIVem, SiIVem, LyaNVem, OIem, CIIem, OVIem,
                         CIVll, SiIVll, NVll, OVIll,
                         CIVllred, CIVllblue,
                         SiIVllred, SiIVllblue,
                         CII_emitted, OI_emitted,
                         NVllred, NVllblue,
                         Lya,
                         OVIllred, OVIllblue,
                         Lyb):

    if CIVem:
        plt.text(CIVll*(1+zem)-30,topemlabel,'CIV',color='black',rotation=90,fontname='serif',verticalalignment='top') ###plots emission label

    if SiIVem:
        plt.text(SiIVllred*(1+zem)-40.,topemlabel,'SiIV+OIV]',color='black',rotation=90,fontname='serif',verticalalignment='top')

    if LyaNVem:
        alpha = 'Ly' + chr(945)
        plt.text(NVllred*(1+zem)+30.,topemlabel,alpha+'+NV',color='black',rotation=90,fontname='serif',verticalalignment='top')

    if OIem:
        plt.text(OI_emitted*(1+zem)-35.,topemlabel,'OI',color='black',rotation=90,fontname='serif',verticalalignment='top')

    if CIIem:
        plt.text(CII_emitted*(1+zem)-30.,topemlabel,'CII',color='black',rotation=90,fontname='serif',verticalalignment='top')

    if OVIem:
        plt.text(OVIll*(1+zem)-30.,topemlabel,'OVI',color='black',rotation=90,fontname='serif',verticalalignment='top')