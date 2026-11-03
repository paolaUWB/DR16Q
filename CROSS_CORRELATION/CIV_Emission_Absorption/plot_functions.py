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
from scipy.stats import pearsonr, t
from astropy.io import fits
from astropy import stats
from scipy.optimize import curve_fit
import re
import math
import ndtest #pip installed from here: https://github.com/syrte/ndtest/blob/master/README.md
import matplotlib as mpl
from CIVfunctions import project,CIV_distance
import matplotlib.patches as patches

mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'

# --------------- [Fig 1] CIV parameter space scatter plot -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bestfit = np.load("bestfit.npy")


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
    
    style = "Simple, tail_width=1.25, head_width=20, head_length=20"
    kw = dict(arrowstyle=style, color="k")
    a1 = patches.FancyArrowPatch((2000, 2.1), (3900, 1.7),
                                 connectionstyle="arc3,rad=.35", **kw)
    plt.gca().add_patch(a1)
    plt.text(x=2440,y=2.02,s="Increasing", fontsize=15, weight='extra bold')
    plt.text(x=2400,y=1.93,s="C$_\mathrm{IV}$ Distance", fontsize=15)
    
    plt.plot(bestfit[:,0], bestfit[:,1], 'k--')
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
                             scattermark = 'o', alpha=0.5, fixed_ylims=False, lowy = 0.6, topy = 2):
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
    
    if cmap != None and fixed_ylims == False:
        ax.scatter(x, y, s=area, c=cvals, cmap = cmap, marker = scattermark, edgecolor = 'k')
        lowy = np.nanmin(y) - 0.2
        topy = np.nanmax(y) + 0.2
        
    elif cmap != None and fixed_ylims == True:
        ax.scatter(x, y, s=area, c=cvals, cmap = cmap, marker = scattermark, edgecolor = 'k')
        lowy = 0.7 #fixing y limits so they are the same in both panels of Fig 4
        topy = 1.9
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
        scattermark='o', fixed_ylims=True
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
        alpha=0, fixed_ylims=True
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
    

    
    if cmap != None:
        ax.scatter(x, y, s=area, c=cvals, cmap = cmap, marker = scattermark, edgecolor = 'k')
        lowy = np.nanmin(y) - 0.2
        topy = np.nanmax(y) + 0.2
    else:
        ax.scatter(x, y, s=area, color=color, marker = scattermark)
        lowy = -4
        topy = 1

    axis_range = topy - lowy 
    x_mid = (0 + 1.4) / 2
    ax.set_xlim([0, 1.4])
    #ax.set_xlim(x_mid - axis_range/2, x_mid + axis_range/2) #for when aspect is equal
    
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
        '''
        
        x_fit = x
        y_fit = y
        print(x)
        print(y)
        
        print('EHVO sample before mask', len(x), len(y))
        print('EHVO sample being fit', len(x_fit), len(y_fit))
        '''
        (m, b), _ = curve_fit(linear_model, x_fit, y_fit)
    
        # ---- Black best-fit line ----
        x_line = np.linspace(min(x_fit), max(x_fit), 300)
        y_line = m * x_line + b
        
        
        # ---- Calculate 95% confidence band ----
        n = len(x_fit)
        y_fit_pred = m * x_fit + b
        residuals = y_fit - y_fit_pred
        
        dof = n - 2  # degrees of freedom for a linear fit
        
        if n > 2:
            residual_std = np.sqrt(np.sum(residuals**2) / dof) # Estimates the residual standard deviation, s
            x_mean = np.mean(x_fit)
            Sxx = np.sum((x_fit - x_mean)**2) # Measures how spread out the x-values are
        
            # Standard error of the fitted mean at each x
            se_fit = residual_std * np.sqrt(1/n + (x_line - x_mean)**2 / Sxx)
        
            t_crit = t.ppf(0.975, dof)  # 95% two-sided confidence
        
            lower = y_line - t_crit * se_fit
            upper = y_line + t_crit * se_fit
        
            # Shaded bow-tie confidence band
            ax.fill_between(x_line, lower, upper, color='black', alpha=0.18, linewidth=0, label='95% Confidence Band')
                
        # ------------------------------------------

        ax.plot(x_line, y_line, color='k', linewidth=1.2, label='EHVO Fit Line')

        # ---- Partition red line ----
        # Calculate the slope of a perpendicular partition line for an arbitrary plot aspect ratio.
        #
        # The standard negative reciprocal (m_red = -1/m) only gives a visually
        # perpendicular line when the x- and y-axes have equal physical scaling.
        #
        # Since our y-axis is log10(HeII EW) and our x-axis is linear (C IV distance),
        # we calculate the perpendicular slope in the plotted coordinate system,
        # accounting for the displayed axis ranges and physical dimensions.
        #
        # This ensures that the red partition line appears perpendicular to the
        # black EHVO fit line without requiring ax.set_aspect('equal').
        #
        # Note: Recalculate after figure layout/axis dimensions are finalized,
        # since resizing the figure can change the required perpendicular slope.
        #
        # Formula:
        # m_red = -(1/m) * (dy/dx)**2 * (W/H)**2
        #
        # where:
        # m  = slope of the black fit line
        # dx = x-axis data range
        # dy = y-axis data range
        # W  = physical width of plotting area
        # H  = physical height of plotting area
        
        # Get plotting dimensions and axis limits
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        
        bbox = ax.get_window_extent()
        
        W = bbox.width
        H = bbox.height
        
        dx = xlim[1] - xlim[0]
        dy = ylim[1] - ylim[0]
        
        # Calculate visually perpendicular slope
        m_red = -1 / m * (dy / dx)**2 * (W / H)**2
        
        '''
        slope_factor = 1    #[was 11] change slope here this makes partition line not perpendicluar to fit line
        m_red = (-1 / m) * slope_factor
        '''
        # OG method
        '''
        x0 = 0.8                # change horizontal position
        #y0 = -0.38               # vertical position option 1
        y0 = -0.15               # vertical position option 2
        '''
        
        # --- testing partition placement at first quartile of residuals plot ---------------
        # Desired intersection x-coordinate
        x_target = 0.7434 #CIV distance value at first quartile of the EHVO fit line residuals plot
        
        # Calculate y0 so the red line crosses the black line at x_target
        y0 = m * x_target + b
        
        #plot a line at the target point
        #ax.axvline(x=x_target, color='blue', linestyle=':', linewidth=1.5, label=f'x = {x_target}')
        x0 = x_target
        #---------------------------------------------------------------------
    
        #redline_xvals = np.arange(0.4, 1.1, 0.1)
        redline_xvals = np.array(ax.get_xlim()) # to plot over full plot limits
        #redline_xvals = np.sort(x) # to span EHVO x values
        y_vals = m_red * (redline_xvals - x0) + y0
    
        #ax.set_aspect('equal', adjustable='box') #aspect set equal
        ax.plot(redline_xvals, y_vals, 'r--', linewidth=1.5, label='EHVO Separation Line')
        
        # ---- Calculate intersection of red and black lines ----
        denom = m - m_red
        
        if not np.isclose(denom, 0):
            x_cross = (y0 - b - m_red * x0) / denom
            y_cross = m * x_cross + b
        
            print(f'EHVO lines intersect at: ({x_cross:.4f}, {y_cross:.4f})')
            '''
            # Mark the intersection
            ax.scatter(
                x_cross, y_cross,
                color='blue',
                marker='o',
                s=60,
                zorder=5,
                label='EHVO Line Intersection'
            )
        
            # Label the intersection coordinates
            ax.annotate(
                f'({x_cross:.3f}, {y_cross:.3f})',
                xy=(x_cross, y_cross),
                xytext=(8, 8),
                textcoords='offset points',
                color='blue',
                fontsize=9
            )
            '''
        else:
            print("The two lines are parallel or coincident.")
        
        indices_above = np.where(y > (m_red * (x - x0) + y0))[0]
        indices_below = np.where(y < (m_red * (x - x0) + y0))[0]
        
        print(indices_above)
        print(indices_below)
        print(f'EHVO Partition Line: ({m_red} * (x - {x0}) + {y0})')
    else:
        m_red = np.nan
        x0 = np.nan
        y0 = np.nan
        
                

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
        
    return m_red, x0, y0

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
    #m_red, x0, y0 = scatter_hist_HeII_outline(bal_Blue, bal_EW, ax, ax_histx, ax_histy, 'cornflowerblue', 2, 'no', 0, x_size=len(bal_Blue), y_size=len(bal_EW), bisector='yes')
    m_red, x0, y0 = scatter_hist_HeII_outline(EHVO_Blue, EHVO_EW, ax, ax_histx, ax_histy, 'purple', 10, 'no', 0, x_outline=outline_Blue, y_outline=outline_EW, outline_color='purple', x_size=len(EHVO_Blue), y_size=len(EHVO_EW), bisector='yes')
    #ax.scatter(801.759345146755, 2.02388722159978, marker='*', s=80, color='green') #FOR XRAY EHVO
    
    
    # Create proxy artists for the legend
    # Create proxy artists for the legend
    from matplotlib.patches import Patch

    confidence_proxy = Patch(facecolor='black', alpha=0.18, label='95% Confidence Band')
    non_bal_proxy = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='aqua', markersize=5)
    bal_proxy = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='cornflowerblue', markersize=6)
    EHVO_proxy = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='purple', markersize=6)
    outline_proxy = plt.Line2D([0], [0], color='purple', linestyle='-')
    xray_EHVO_proxy = plt.Line2D([0], [0], marker='*', linestyle='none', color='green', markerfacecolor='green', markersize=8)
    fitline_proxy = plt.Line2D([0], [0], color='k', linestyle='-')
    sepline_proxy = plt.Line2D([0], [0], color='red', linestyle='--')

    ax.legend([non_bal_proxy, bal_proxy, EHVO_proxy, outline_proxy, fitline_proxy, sepline_proxy, confidence_proxy], [non_bal_Sample_Name, bal_Sample_Name, EHVO_Sample_name, 'EHVO (RH+2020)', 'EHVO Fit', 'EHVO Partition', '95% Confidence Band'], loc='lower left', fontsize=12)
    #ax.legend([non_bal_proxy, bal_proxy, EHVO_proxy, outlier_proxy, outline_proxy, xray_EHVO_proxy], [non_bal_Sample_Name, bal_Sample_Name, EHVO_Sample_name, outlier_name, 'EHVO (RH+2020)', 'X-ray selected EHVO'], loc='upper right')
    #FOR XRAY EHVO ^^^^
    
    plt.show()
    plt.close()
    
    return m_red, x0, y0

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
                    parent_mask=None,
                    xlabel='add label', ylabel='add label', 
                    xlim= None, ylim = None,
                    gridsize=30,
                    show_left=True,
                    colormap='viridis',
                    fixed_vmin=500,
                    cmap_sample='CIVBlue',
                    hex_func=np.max,
                    left_label=None,
                    right_label=None, 
                    split_parent=False):

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
        
        if parent_mask is not None:
            c_parent1 = np.log10(np.where(c_parent[~parent_mask]>0, c_parent[~parent_mask], np.nan))
            c_parent2 = np.log10(np.where(c_parent[parent_mask]>0, c_parent[parent_mask], np.nan))
            
        else:
            c_parent = np.log10(np.where(c_parent>0, c_parent, np.nan))
        c_sample1 = np.log10(np.where(c_ehvo[~ehvo_mask]>0, c_ehvo[~ehvo_mask], np.nan))
        c_sample2 = np.log10(np.where(c_ehvo[ehvo_mask]>0, c_ehvo[ehvo_mask], np.nan))

    if parent_mask is None:
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
        
    else:
        
        extent = (*xlim, *ylim)

        # ---- LEFT PANEL (optional) ----
        if show_left:
    
            cm = ax1.hexbin(
                x_parent[~parent_mask], y_parent[~parent_mask],
                C=c_parent1,
                gridsize=gridsize,
                extent=extent,
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
        cm2 = ax2.hexbin(
            x_parent[parent_mask], y_parent[parent_mask],
            C=c_parent2,
            gridsize=gridsize,
            extent=extent,
            cmap=colormap,
            reduce_C_function=hex_func,
            mincnt=4,
            vmin=vmin_val,
            vmax=vmax_val,
            edgecolor='white',
            linewidths=0.25,
            zorder=99
        )
    
        sc2 = ax2.scatter(
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

        # LEFT PANEL: EHVOs above partition
        hexpositions1 = cm.get_offsets()
        hexvalues1 = cm.get_array()
        
        x_cut1 = x_ehvo[~ehvo_mask]
        y_cut1 = y_ehvo[~ehvo_mask]
        c_cut1 = c_ehvo[~ehvo_mask]
        
        count_left = 0
        
        print("Number of plotted left EHVOs:", np.sum(~ehvo_mask))
        print("Number of left EHVOs compared:", len(c_cut1))
        print("Number of left EHVOs above:", count_left)
        print("Number of left hexagons:", len(hexpositions1))

        
        for i in range(len(c_cut1)):
        
            ehvo_value = c_cut1[i]
        
            hex_value = hexfind(
                x_cut1[i],
                y_cut1[i],
                hexpositions1,
                hexvalues1
            )
        
            if ehvo_value >= hex_value:
                count_left += 1
        
        
        # RIGHT PANEL: EHVOs below partition
        hexpositions2 = cm2.get_offsets()
        hexvalues2 = cm2.get_array()
        
        x_cut2 = x_ehvo[ehvo_mask]
        y_cut2 = y_ehvo[ehvo_mask]
        c_cut2 = c_ehvo[ehvo_mask]
        
        count_right = 0
        
        for i in range(len(c_cut2)):
        
            ehvo_value = c_cut2[i]
        
            hex_value = hexfind(
                x_cut2[i],
                y_cut2[i],
                hexpositions2,
                hexvalues2
            )
        
            if ehvo_value >= hex_value:
                count_right += 1
        
        
        print("EHVOs above hex value - LEFT:", count_left)
        print("EHVOs above hex value - RIGHT:", count_right)
        
        count = [count_left, count_right]
    return count

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

