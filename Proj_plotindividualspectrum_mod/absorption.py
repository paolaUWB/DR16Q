#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 09:03:57 2026

@author: jmaharaj
"""

import matplotlib.pyplot as plt


def plot_absorption_regions(zabs_min, zabs_max, vmin, ###function for plotting absorption regions on the spectrum
CIVll, NVll, NVllblue, NVllred,
OVIll, OVIllblue, OVIllred,
SiIVll, SiIVllblue, SiIVllred,
Lya, Lyb, CII_emitted, OI_emitted,
CIVabs, NVabs, OVIabs, SiIVabs, Lyaabs, Lybabs, CIIabs, OIabs,
colors):

    for k in range(0,len(vmin)): ###loops for each absorption system which is each pair of vmin/vmax
                                 ###converts velocity limits into absorption redshift vals
                                 ###uses redshifts to get observed wavelength ranges

        if CIVabs:
            plt.axvspan(CIVll*(1.+zabs_max[k]),CIVll*(1.+zabs_min[k]), alpha=0.2, color=colors['CIV']) ###shades CIV region in pink
            plt.text(CIVll*(1.+zabs_min[k])-30.,0.5-0.1*k,'CIV',color=colors['CIV'],fontname='serif',weight='bold') ###text for CIV

        if NVabs:
            plt.axvspan(NVllblue*(1.+zabs_max[k]),NVllred*(1.+zabs_min[k]), alpha=0.2, color=colors['NV'])
            plt.text(NVll*(1.+zabs_min[k]),1.3-0.1*k,'NV',color=colors['NV'],fontname='serif',weight='bold')

        if OVIabs:
            plt.axvspan(OVIllblue*(1.+zabs_max[k]),OVIllred*(1.+zabs_min[k]), alpha=0.2, color=colors['OVI'])
            plt.text(OVIll*(1.+zabs_min[k]),1.4-0.1*k,'OVI',color=colors['OVI'],fontname='serif',weight='bold')

        if SiIVabs:
            plt.axvspan(SiIVllblue*(1.+zabs_max[k]),SiIVllred*(1.+zabs_min[k]), alpha=0.2, color=colors['SiIV'])
            plt.text(SiIVll*(1.+zabs_min[k]),0.25*k+2,'SiIV',color=colors['SiIV'],fontname='serif',weight='bold')

        if Lyaabs:
            plt.axvspan(Lya*(1.+zabs_max[k]),Lya*(1.+zabs_min[k]), alpha=0.2, color=colors['Lya'])
            plt.text(Lya*(1.+zabs_min[k])-20.,1.45-0.1*k,'Lya',color=colors['Lya'],fontname='serif',weight='bold')
            
        if CIIabs:
            plt.axvspan(CII_emitted*(1.+zabs_max[k]),CII_emitted*(1.+zabs_min[k]), alpha=0.2, color=colors['CII'])
            plt.text(CII_emitted*(1.+zabs_min[k]),1.6-0.1*k,'CII',color=colors['CII'],fontname='serif',weight='bold')
            
        if OIabs:
            plt.axvspan(OI_emitted*(1.+zabs_max[k]),OI_emitted*(1.+zabs_min[k]), alpha=0.2, color=colors['OI'])
            plt.text(OI_emitted*(1.+zabs_min[k]),1.7-0.1*k,'OI',color=colors['OI'],fontname='serif',weight='bold')
            
        if Lybabs:
            plt.axvspan(Lyb*(1.+zabs_max[k]),Lyb*(1.+zabs_min[k]), alpha=0.2, color=colors['Lyb'])
            plt.text(Lyb*(1.+zabs_min[k]),1.8-0.1*k,'Lyb',color=colors['Lyb'],fontname='serif',weight='bold')