#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 12:03:36 2026

@author: lilianaflores
"""

import os
import sys
from os.path import exists
import csv
import pandas as pd
import scipy.stats as stat
import numpy as np
from pylab import*
from sympy import sympify
from matplotlib.backends.backend_pdf import PdfPages
from scipy import*
from astropy import*
from astropy.table import Table
from scipy.stats import ks_2samp
from astropy import constants as const
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import pearsonr
from astropy.io import fits
from astropy import stats
from scipy.optimize import curve_fit
import re
import math
import ndtest #pip installed from here: https://github.com/syrte/ndtest/blob/master/README.md
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'

# --------------- [Fig 1] CIV parameter space scatter plot -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def scatter_CivBlue_EW(x, y, ax, color, area, factor):
    ax.set_xlim([-2000, 7000])
    ax.set_ylim([0.7, 2.5])
    ax.scatter(x, y, s = area, color = color)
    ax.text(10.5,0.5,'')
    ax.set_xlabel(r'$\mathrm{C\,IV}$ blueshift [km s$^{-1}$]', fontsize=18)
    ax.set_ylabel(r'$\log_{10}(\mathrm{C\,IV}\ \mathrm{EW\ [\AA]})$', fontsize=18)
def figure_1(Blue_rankAll_good,EW_rankAll_good,Blue_rankAll_SNR,EW_rankAll_SNR,Blue_RH_parent,EW_RH_parent,rankAll_sample_name,rankAll_SNR_sample_name,RH_parent_sample_name):
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005
    rect_scatter = [left, bottom, width, height]
    fig = plt.figure(1)
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_axes(rect_scatter)
    scatter_CivBlue_EW(Blue_rankAll_good,EW_rankAll_good,ax, 'grey', 2, 0)
    scatter_CivBlue_EW(Blue_rankAll_SNR,EW_rankAll_SNR,ax,'cornflowerblue', 2, 0)
    scatter_CivBlue_EW(Blue_RH_parent,EW_RH_parent,ax, 'purple', 2, 0)
    ax.legend([rankAll_sample_name,rankAll_SNR_sample_name,RH_parent_sample_name],loc='upper right', markerscale=3, fontsize=12)
    plt.show()
    plt.close()

# -------------- [Fig 2] CIV parameter space scatter/hist plot -----------------------------------------------------------------------------------------------------------------------------------------------

def scatter_hist_CIV_outline(x, y, 
                             ax, 
                             ax_histx, ax_histy, 
                             color, area, mult, factor, 
                             x_outline=None, y_outline=None, 
                             cmap = None, cvals = None, 
                             outline_color='red',
                             x_size=1, y_size=1, 
                             scattermark = 'o', alpha=0.5):
    """
    This function creates a scatter plot along with its marginal histograms.

    Parameters:
    x (array): X data for scatter and histogram.
    y (array): Y data for scatter and histogram.
    ax (AxesSubplot): Axes object for the scatter plot.
    ax_histx (AxesSubplot): Axes object for the horizontal histogram.
    ax_histy (AxesSubplot): Axes object for the vertical histogram.
    color (str): Color for the scatter and histograms.
    area (int or float): Size of the scatter points.
    """
    ax_histx.tick_params(axis='x', labelbottom=False)
    ax_histy.tick_params(axis='y', labelleft=False)
    ax.set_xlim([-2000, 6999])
    
    if cmap != None:
        ax.scatter(x, y, s=area, c=cvals, cmap = cmap, marker = scattermark, edgecolor = 'k')
        lowy = np.nanmin(y) - 0.2
        topy = np.nanmax(y) + 0.2
    else:
        ax.scatter(x, y, s=area, color=color, marker = scattermark)
        lowy = 0.75
        topy = 2.5

    ax.set_ylim([lowy, topy])
    ax.text(10.5, 0.5, '')
    #ax.set_xlabel('C IV blueshift (km/s)', fontsize=18)
    #ax.set_ylabel('log C IV EW $(\AA)$', fontsize=18)
    ax.set_xlabel(r'$\mathrm{C\,IV}$ blueshift [km s$^{-1}$]', fontsize=18)
    ax.set_ylabel(r'$\log_{10}(\mathrm{C\,IV}\ \mathrm{EW\ [\AA]})$', fontsize=18)

    binwidth_x = 400.0
    binwidth_y = 0.05
    upperlimx = 6000
    lowerlimx = -2000
    upperlimy = 2.5
    lowerlimy = 0.7
    binsx = np.arange(lowerlimx, upperlimx + binwidth_x, binwidth_x)
    binsy = np.arange(lowerlimy, upperlimy + binwidth_y, binwidth_y)


    if mult == 'yes':
        for i in range(0, factor):
            x = np.append(x, x)
            y = np.append(y, y)
        x = np.hstack(x)
        y = np.hstack(y)
    weightsx = np.ones_like(x) / float(x_size)
    weightsy = np.ones_like(y) / float(y_size)
    ax_histx.hist(x, bins=binsx, weights=weightsx, color=color, alpha=alpha)
    ax_histy.hist(y, bins=binsy, weights=weightsy, orientation='horizontal', color=color, alpha=alpha)

    if x_outline is not None and y_outline is not None:
        outline_weights_x = np.ones_like(x_outline) / float(x_size)  # Normalized to the same total count as x
        outline_weights_y = np.ones_like(y_outline) / float(y_size)  # Normalized to the same total count as y
        ax_histx.hist(x_outline, bins=binsx, weights=outline_weights_x, histtype='step', color=outline_color)
        ax_histy.hist(y_outline, bins=binsy, weights=outline_weights_y, orientation='horizontal', histtype='step', color=outline_color)

def figure_2_outline(non_bal_Blue, non_bal_EW, bal_Blue, bal_EW, EHVO_Blue, EHVO_EW, non_bal_Sample_Name, bal_Sample_Name, EHVO_Sample_name, outline_Blue, outline_EW):
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]
    fig = plt.figure(1)
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_axes(rect_scatter)
    ax_histx = fig.add_axes(rect_histx, sharex=ax)
    ax_histy = fig.add_axes(rect_histy, sharey=ax)
    scatter_hist_CIV_outline(non_bal_Blue, non_bal_EW, ax, ax_histx, ax_histy, 'aqua', 2, 'no', 0, x_size=len(non_bal_Blue), y_size=len(non_bal_EW))
    scatter_hist_CIV_outline(bal_Blue, bal_EW, ax, ax_histx, ax_histy, 'cornflowerblue', 2, 'no', 0, x_size=len(bal_Blue), y_size=len(bal_EW))
    scatter_hist_CIV_outline(EHVO_Blue, EHVO_EW, ax, ax_histx, ax_histy, 'purple', 10, 'no', 0, x_outline=outline_Blue, y_outline=outline_EW, outline_color='purple', x_size=len(EHVO_Blue), y_size=len(EHVO_EW))
    #ax.scatter(801.759345146755, 2.02388722159978, marker='*', s=80, color='green') #FOR XRAY EHVO

    #Displaying the oulier differently
    outlier = np.where((EHVO_EW)>2.1) #finding outlier
    print()
    print(outlier)
    print(EHVO_Blue[71])
    ax.scatter(EHVO_Blue[71], EHVO_EW[71],facecolors='white', edgecolors='purple', s=12)
    outlier_name = 'Outlier'
    
    # Create proxy artists for the legend
    # Create proxy artists for the legend
    non_bal_proxy = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='aqua', markersize=6)
    bal_proxy = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='cornflowerblue', markersize=6)
    EHVO_proxy = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='purple', markersize=6)
    outline_proxy = plt.Line2D([0], [0], color='purple', linestyle='-')
    outlier_proxy = plt.Line2D([0], [0], marker='o', linestyle='none', color='purple', markerfacecolor='white', markersize=5)
    xray_EHVO_proxy = plt.Line2D([0], [0], marker='*', linestyle='none', color='green', markerfacecolor='green', markersize=8)


    ax.legend([non_bal_proxy, bal_proxy, EHVO_proxy, outlier_proxy, outline_proxy], [non_bal_Sample_Name, bal_Sample_Name, EHVO_Sample_name, outlier_name, 'EHVO (RH+2020)'], loc='upper right', fontsize=12)
    #ax.legend([non_bal_proxy, bal_proxy, EHVO_proxy, outlier_proxy, outline_proxy, xray_EHVO_proxy], [non_bal_Sample_Name, bal_Sample_Name, EHVO_Sample_name, outlier_name, 'EHVO (RH+2020)', 'X-ray selected EHVO'], loc='upper right')
    #FOR XRAY EHVO ^^^^
    
    plt.show()
    plt.close()

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- [Fig 4] CIV emission parameter space scatter/hist EHVOs only split by median speed w/velocity colormap -------------------------------------------------------------------------------------------------------------------------------------------------
def plot_EHVO_CIV_scatter_hist(
    x_data, y_data, v_data, small_speed_idx, large_speed_idx,
    scatter_label_low, scatter_label_high,
    v_label=r'$V [\mathrm{km\ s^{-1}}]$',
    cmap='Purples'):
    """
    Plot EHVO CIV emission parameter space with velocity subsamples.

    Parameters
    ----------
    x_data : array-like
        X-axis data (CivBlue values).
    y_data : array-like
        Y-axis data (log10 of CivEW values).
    v_data : array-like
        Velocity data (Vmin or Vmax).
    small_speed_idx : array-like
        Boolean/index array selecting "small" velocity objects.
    large_speed_idx : array-like
        Boolean/index array selecting "large" velocity objects.
    scatter_label_low : str
        Legend label for low velocity points.
    scatter_label_high : str
        Legend label for high velocity points.
    v_label : str, optional
        Colorbar label, by default r'$V [\mathrm{km\ s^{-1}}]$'
    cmap : str, optional
        Colormap for scatter points, by default 'Purples'.
    """
    
    # Define axes positions
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    # Create figure and axes
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_axes(rect_scatter)
    ax_histx = fig.add_axes(rect_histx, sharex=ax)
    ax_histy = fig.add_axes(rect_histy, sharey=ax)

    # Define histogram colors
    hist_color1 = 'grey'
    hist_color2 = 'k'

    # Scatter plots with histograms
    scatter_hist_CIV_outline(
        x_data[small_speed_idx],
        np.log10(y_data[small_speed_idx]),
        ax, ax_histx, ax_histy,
        hist_color1, 100, 'no', 0,
        cmap=cmap,
        cvals=np.abs(v_data[small_speed_idx]),
        scattermark='o'
    )

    scatter_hist_CIV_outline(
        x_data[large_speed_idx],
        np.log10(y_data[large_speed_idx]),
        ax, ax_histx, ax_histy,
        hist_color2, 100, 'no', 0,
        cmap=cmap,
        cvals=np.abs(v_data[large_speed_idx]),
        scattermark='v',
        x_outline=x_data[large_speed_idx],
        y_outline=np.log10(y_data[large_speed_idx]),
        outline_color='k',
        alpha=0
    )

    # Set color limits based on velocity
    plt.setp(ax.collections, clim=(np.min(np.abs(v_data)), np.max(np.abs(v_data))))

    # Colorbar
    divider = make_axes_locatable(ax_histy)
    cax = divider.append_axes("right", size="12%", pad=0.1)
    sc = ax.collections[0]
    cbar = plt.colorbar(sc, cax=cax)
    cbar.set_label(v_label, fontsize=16, rotation=270, labelpad=16)
    cbar.ax.yaxis.set_label_coords(6, 0.5)

    # Legend proxies
    low_proxy = plt.Line2D([0], [0], marker='o', color=hist_color1, markersize=8, linestyle='none')
    high_proxy = plt.Line2D([0], [0], marker='v', color=hist_color2, markersize=8, linestyle='-')
    ax.legend([low_proxy, high_proxy], [scatter_label_low, scatter_label_high], loc='upper right', fontsize=12)

    plt.show()
    plt.close()

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------- [Fig 5] HeII vs CIV distance scatter/hist plot  ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def scatter_hist_HeII_outline(x, y, ax, ax_histx, ax_histy, color, area, mult, factor, x_outline=None, y_outline=None, cmap = None, cvals = None, outline_color='red',x_size=1, y_size=1, scattermark = 'o', alpha=0.5, bisector='no'):
    """
    This function creates a scatter plot along with its marginal histograms.

    Parameters:
    x (array): X data for scatter and histogram.
    y (array): Y data for scatter and histogram.
    ax (AxesSubplot): Axes object for the scatter plot.
    ax_histx (AxesSubplot): Axes object for the horizontal histogram.
    ax_histy (AxesSubplot): Axes object for the vertical histogram.
    color (str): Color for the scatter and histograms.
    area (int or float): Size of the scatter points.
    """
    ax_histx.tick_params(axis='x', labelbottom=False)
    ax_histy.tick_params(axis='y', labelleft=False)
    ax.set_xlim([0, 1.4])
    
    if cmap != None:
        ax.scatter(x, y, s=area, c=cvals, cmap = cmap, marker = scattermark, edgecolor = 'k')
        lowy = np.nanmin(y) - 0.2
        topy = np.nanmax(y) + 0.2
    else:
        ax.scatter(x, y, s=area, color=color, marker = scattermark)
        lowy = -4
        topy = 1

    ax.set_ylim([lowy, topy])
    ax.text(10.5, 0.5, '')
    #ax.set_xlabel('C IV blueshift (km/s)', fontsize=15)
    #ax.set_ylabel('log C IV EW $(\AA)$', fontsize=15)
    ax.set_xlabel(r'$\mathrm{C\,IV}$ distance', fontsize=18)
    ax.set_ylabel(r'$\log_{10}(\mathrm{He\,II}\ \mathrm{EW\ [\AA]})$', fontsize=18)

    binwidth_x = 0.05
    binwidth_y = 0.25
    upperlimx = 1.4
    lowerlimx = 0
    upperlimy = 1
    lowerlimy = -4
    binsx = np.arange(lowerlimx, upperlimx + binwidth_x, binwidth_x)
    binsy = np.arange(lowerlimy, upperlimy + binwidth_y, binwidth_y)
    
    xticks = np.arange(0, 1.4 + 0.2, 0.2)  # x ticks every 0.2
    yticks = np.arange(-4, 1 + 1, 1)      # y ticks every 1
    
    ax.set_xticks(xticks[:-1])
    ax.set_yticks(yticks[:-1])
    
    if bisector == 'yes':

        def linear_model(x, m, b):
            return m * x + b
    
        mask = np.isfinite(x) & np.isfinite(y)
        x_fit = x[mask]
        y_fit = y[mask]
    
        (m, b), _ = curve_fit(linear_model, x_fit, y_fit)
    
        # ---- Black best-fit line ----
        x_line = np.linspace(min(x_fit), max(x_fit), 300)
        y_line = m * x_line + b
        ax.plot(x_line, y_line, color='k', linewidth=1.2, label='EHVO Fit Line')
    
        # ---- Adjustable red line ----
        slope_factor = 11    # change slope here
        m_red = (-1 / m) * slope_factor
    
        x0 = 0.8                # change horizontal position
        #y0 = -0.38               # vertical position option 1
        y0 = -0.15               # vertical position option 2
    
        redline_xvals = np.arange(0.4, 1.1, 0.1)
        
        x_vals = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 300)
        y_vals = m_red * (redline_xvals - x0) + y0
    
        ax.plot(redline_xvals, y_vals, 'r--', linewidth=1.5, label='EHVO Separation Line')
        
        indices_above = np.where(y > (m_red * (x - x0) + y0))[0]
        indices_below = np.where(y < (m_red * (x - x0) + y0))[0]
        
        print(indices_above)
        print(indices_below)
        
                

    if mult == 'yes':
        for i in range(0, factor):
            x = np.append(x, x)
            y = np.append(y, y)
        x = np.hstack(x)
        y = np.hstack(y)
    weightsx = np.ones_like(x) / float(x_size)
    weightsy = np.ones_like(y) / float(y_size)
    ax_histx.hist(x, bins=binsx, weights=weightsx, color=color, alpha=alpha)
    ax_histy.hist(y, bins=binsy, weights=weightsy, orientation='horizontal', color=color, alpha=alpha)

    if x_outline is not None and y_outline is not None:
        outline_weights_x = np.ones_like(x_outline) / float(x_size)  # Normalized to the same total count as x
        outline_weights_y = np.ones_like(y_outline) / float(y_size)  # Normalized to the same total count as y
        ax_histx.hist(x_outline, bins=binsx, weights=outline_weights_x, histtype='step', color=outline_color)
        ax_histy.hist(y_outline, bins=binsy, weights=outline_weights_y, orientation='horizontal', histtype='step', color=outline_color)

def figure_5_outline(non_bal_Blue, non_bal_EW, bal_Blue, bal_EW, EHVO_Blue, EHVO_EW, non_bal_Sample_Name, bal_Sample_Name, EHVO_Sample_name, outline_Blue, outline_EW):
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]
    fig = plt.figure(1)
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_axes(rect_scatter)
    ax_histx = fig.add_axes(rect_histx, sharex=ax)
    ax_histy = fig.add_axes(rect_histy, sharey=ax)
    
    scatter_hist_HeII_outline(non_bal_Blue, non_bal_EW, ax, ax_histx, ax_histy, 'aqua', 2, 'no', 0, x_size=len(non_bal_Blue), y_size=len(non_bal_EW))
    scatter_hist_HeII_outline(bal_Blue, bal_EW, ax, ax_histx, ax_histy, 'cornflowerblue', 2, 'no', 0, x_size=len(bal_Blue), y_size=len(bal_EW))
    scatter_hist_HeII_outline(EHVO_Blue, EHVO_EW, ax, ax_histx, ax_histy, 'purple', 10, 'no', 0, x_outline=outline_Blue, y_outline=outline_EW, outline_color='purple', x_size=len(EHVO_Blue), y_size=len(EHVO_EW), bisector='yes')
    #ax.scatter(801.759345146755, 2.02388722159978, marker='*', s=80, color='green') #FOR XRAY EHVO
    
    
    # Create proxy artists for the legend
    # Create proxy artists for the legend
    non_bal_proxy = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='aqua', markersize=5)
    bal_proxy = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='cornflowerblue', markersize=6)
    EHVO_proxy = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='purple', markersize=6)
    outline_proxy = plt.Line2D([0], [0], color='purple', linestyle='-')
    xray_EHVO_proxy = plt.Line2D([0], [0], marker='*', linestyle='none', color='green', markerfacecolor='green', markersize=8)
    fitline_proxy = plt.Line2D([0], [0], color='k', linestyle='-')
    sepline_proxy = plt.Line2D([0], [0], color='red', linestyle='--')


    ax.legend([non_bal_proxy, bal_proxy, EHVO_proxy, outline_proxy, fitline_proxy, sepline_proxy], [non_bal_Sample_Name, bal_Sample_Name, EHVO_Sample_name, 'EHVO (RH+2020)', 'EHVO Fit', 'EHVO Partition'], loc='lower left', fontsize=12)
    #ax.legend([non_bal_proxy, bal_proxy, EHVO_proxy, outlier_proxy, outline_proxy, xray_EHVO_proxy], [non_bal_Sample_Name, bal_Sample_Name, EHVO_Sample_name, outlier_name, 'EHVO (RH+2020)', 'X-ray selected EHVO'], loc='upper right')
    #FOR XRAY EHVO ^^^^
    
    plt.show()
    plt.close()

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------- [Figs 6 & 8] CIV Blushift/HeII EW colormap hex/scatter in physical property parameter spaces ----------------------------------------------------------------------------------------------------------------------------------------------------------
def hexfind(x,y,hexpos,parentval):
    #Loop hexpos
    pos = np.array([x,y])
    dist = []
    for i in hexpos:
        dist.append(math.dist(pos,i))
    minimum = np.min(dist)
    minpos = np.where(dist == minimum)

    value = parentval[minpos]
    return value

def plot_CIVBlue_hexbin_ehvo(x_parent, y_parent, c_parent,
                    x_ehvo, y_ehvo, c_ehvo,
                    ehvo_mask,
                    xlabel, ylabel, 
                    xlim= None, ylim = None,
                    gridsize=30,
                    show_left=True,
                    colormap='viridis',
                    fixed_vmin=500,
                    cmap_sample='CIVBlue',
                    hex_func=np.max,
                    left_label=None,
                    right_label=None):

    if show_left:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))
    else:
        fig, ax2 = plt.subplots(figsize=(6,5))
        
    if cmap_sample == 'CIVBlue': # want a set color scale min when plotting CIV Blueshift
        vmin_val = fixed_vmin
        vmax_val = max(np.concatenate([c_parent,c_ehvo]))
        
        cbar_label = r'$\mathrm{C\,IV}$ blueshift [km s$^{-1}$]'
        
        c_sample1 = c_ehvo[~ehvo_mask]
        c_sample2 = c_ehvo[ehvo_mask]
        
    elif cmap_sample == 'HeIIEW':
        combined = np.concatenate([c_parent,c_ehvo])
        combined = np.log10(combined[combined>0])
        
    
        vmin_val = np.nanmin(combined)
        vmax_val = np.nanmax(combined)
        
        cbar_label = r'$\log_{10}(\mathrm{He\,II}\ \mathrm{EW\ [\AA]})$'
        
        c_parent = np.log10(np.where(c_parent>0, c_parent, np.nan))
        c_sample1 = np.log10(np.where(c_ehvo[~ehvo_mask]>0, c_ehvo[~ehvo_mask], np.nan))
        c_sample2 = np.log10(np.where(c_ehvo[ehvo_mask]>0, c_ehvo[ehvo_mask], np.nan))


    # ---- LEFT PANEL (optional) ----
    if show_left:

        cm = ax1.hexbin(
            x_parent, y_parent,
            C=c_parent,
            gridsize=gridsize,
            cmap=colormap,
            reduce_C_function=hex_func,
            mincnt=4,
            vmin=vmin_val,
            vmax=vmax_val,
            edgecolor='white',
            linewidths=0.25,
            zorder=99
        )

        ax1.set_xlabel(xlabel, fontsize=18)
        ax1.set_ylabel(ylabel, fontsize=18)
        ax1.set_xlim(*xlim)
        ax1.set_ylim(*ylim)
        
        sc = ax1.scatter(
            x_ehvo[~ehvo_mask], # ~ gives the opposite mask condition
            y_ehvo[~ehvo_mask],
            c=c_sample1, 
            cmap=colormap,
            edgecolor='k',
            s=70,
            zorder=100,
            vmin=cm.norm.vmin,
            vmax=cm.norm.vmax,
            label=left_label
        )

    # ---- RIGHT PANEL ----
    cm = ax2.hexbin(
        x_parent, y_parent,
        C=c_parent,
        gridsize=gridsize,
        cmap=colormap,
        reduce_C_function=hex_func,
        mincnt=4,
        vmin=vmin_val,
        vmax=vmax_val,
        edgecolor='white',
        linewidths=0.25,
        zorder=99
    )

    sc = ax2.scatter(
        x_ehvo[ehvo_mask],
        y_ehvo[ehvo_mask],
        c=c_sample2,
        cmap=colormap,
        edgecolor='k',
        s=70,
        zorder=100,
        vmin=cm.norm.vmin,
        vmax=cm.norm.vmax,
        label=right_label
    )

    ax2.set_xlabel(xlabel, fontsize=18)
    ax2.set_xlim(*xlim)
    ax2.set_ylim(*ylim)

    if show_left:
        ax2.set_yticks([])
    else:
        ax2.set_ylabel(ylabel)

    ax1.legend(loc='upper left', fontsize=12)
    ax2.legend(loc='upper left', fontsize=12)

    # ---- COLORBAR ----
    cbarax = fig.add_axes([1, 0.125, 0.03, 0.845])
    cbar = plt.colorbar(cm, cax=cbarax)
    cbarax.set_ylabel(cbar_label, rotation=270,labelpad=25, fontsize=16)

    plt.tight_layout()
    plt.show()

    # ---- HEXBIN COMPARISON ----
    hexpositions = cm.get_offsets()
    hexvalues = cm.get_array()

    count = 0

    x_cut = x_ehvo[ehvo_mask]
    y_cut = y_ehvo[ehvo_mask]
    c_cut = c_ehvo[ehvo_mask]

    for i in range(len(c_cut)):

        ehvo_blue = c_cut[i]
        hex_blue = hexfind(x_cut[i], y_cut[i], hexpositions, hexvalues)

        if ehvo_blue >= hex_blue:
            count += 1

    print("EHVO above hex value:", count)

    return count

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

