#!/usr/bin/env python3
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from astropy import stats
from astropy.io import fits
import os

infoRankineparent = os.getcwd() + "/DR16parent_DR14RankineInfo.csv"
infoRankineEHVO = os.getcwd() + "/DR16EHVO_DR14RankineInfo.csv"

SPACE_BETWEEN_GRAPHS = 0.01
# area w/ no points on graph
SPACE_X_AXIS = 0.5
SPACE_Y_AXIS = 0.1
HISTOGRAM_TICKS_ON = True

class GraphHandler:
    def __init__(self, left, bottom, width, height):
        rect_scatter = [left, bottom, width, height]
        # what does 0.2 mean in this context? TODO: remove magic number
        rect_histx = [left, bottom + height + SPACE_BETWEEN_GRAPHS, width, 0.2]
        rect_histy = [left + width + SPACE_BETWEEN_GRAPHS, bottom, 0.2, height]

        # this declaration seems redundant but it was in orignal code
        self.fig = plt.figure(1)
        self.fig = plt.figure(figsize=(8,8))

        # define graph boundaries based on dimensions
        self.ax = self.fig.add_axes(rect_scatter)
        self.ax_histx = self.fig.add_axes(rect_histx, sharex = self.ax)
        self.ax_histy = self.fig.add_axes(rect_histy, sharey = self.ax)

        self.x_vals     = []
        self.y_vals     = []
        self.color_vals = []
    
    # TODO: maybe throw an exception if any values out of range?
    def add_point(self, x, y, color, label, size=5, marker = None):
        self.ax.scatter(x, y, s=size, color=color, label=label, marker=marker)

        self.x_vals.append(x)
        self.y_vals.append(y)
        self.color_vals.append(color)
    
    def flatten(self, l):
        return [item for sublist in l for item in sublist]

    def set_axes(self, xlabel=None, ylabel=None):
       self.ax.set(xlabel=xlabel, ylabel=ylabel)

    def display_graph(self):
        # defining tick appearance
        self.ax_histx.tick_params(axis='x', labelbottom=False, bottom=HISTOGRAM_TICKS_ON)
        self.ax_histy.tick_params(axis='y', labelleft=False, left=HISTOGRAM_TICKS_ON)
    
        # defining x/y limits
        ymin, ymax = min(map(min, self.y_vals)), max(map(max, self.y_vals))
        xmin, xmax = min(map(min, self.x_vals)), max(map(max, self.x_vals))

        self.ax.set_ylim(ymin - SPACE_Y_AXIS, ymax + SPACE_Y_AXIS)
        self.ax.set_xlim(xmin - SPACE_X_AXIS, xmax + SPACE_X_AXIS)

        # find binwidth
        binwidth_x = stats.knuth_bin_width(self.flatten(self.x_vals))
        binwidth_y = stats.knuth_bin_width(self.flatten(self.y_vals))

        binsx = np.arange(xmin - SPACE_X_AXIS, xmax + SPACE_X_AXIS, binwidth_x)
        binsy = np.arange(ymin - SPACE_Y_AXIS, ymax + SPACE_Y_AXIS, binwidth_y)

        # iterate through all stored data to build histogram
        for x, y, color in zip(self.x_vals, self.y_vals, self.color_vals):
            weightsx = [np.ones_like(x)/float(len(x))]
            weightsy = [np.ones_like(y)/float(len(y))]

#            x = x + SPACE_X_AXIS
#            y = y + SPACE_Y_AXIS

            # TODO: remove magic number
            self.ax_histx.hist(x, bins=binsx, weights=weightsx, color=color, alpha=0.5)
            self.ax_histy.hist(y, bins=binsy, weights=weightsy, orientation='horizontal', color=color, alpha=0.5)

        # put finishing touches and display everything
        # TODO: set up legend config settings
        self.ax.legend()
        plt.show()

highZfile = pd.read_csv("highRedshiftQuasarCSV.csv")
highZmass = np.array(highZfile["bhm"])
highZedd =  np.array(highZfile["edd_ratio"])
highZbol_lum = np.array(highZfile["bol_lum"])

highZmass = np.log10(highZmass)
highZedd = np.log10(highZedd)
highZbol_lum = np.log10(highZbol_lum)

#Extracting values from Rankines info file for EHVO Parent
dfRPA = pd.read_csv(infoRankineparent, header=None)
dfRHV = pd.read_csv(infoRankineEHVO, header=None)

#for parent sample
Parentin14Rank_mbh=dfRPA[dfRPA.columns[16]].to_numpy()
Parentin14Rank_lbol=dfRPA[dfRPA.columns[17]].to_numpy()
Parentin14Rank_redd=dfRPA[dfRPA.columns[18]].to_numpy()

#for EHVOs
EHVOin14Rank_mbh=dfRHV[dfRHV.columns[16]].to_numpy()
EHVOin14Rank_lbol=dfRHV[dfRHV.columns[17]].to_numpy()
EHVOin14Rank_redd=dfRHV[dfRHV.columns[18]].to_numpy()


MBH_parentR = Parentin14Rank_mbh
Lbol_parentR = Parentin14Rank_lbol
Redd_parentR = Parentin14Rank_redd

MBH_EHVOR = EHVOin14Rank_mbh
Lbol_EHVOR = EHVOin14Rank_lbol
Redd_EHVOR = EHVOin14Rank_redd 


bhm_vs_eddratio = GraphHandler(0.1, 0.1, 0.65, 0.65)
bhm_vs_eddratio.set_axes(xlabel="X-axis", ylabel="Y-axis")

bhm_vs_eddratio.add_point(MBH_parentR, Redd_parentR, 'magenta', 'Renkine parent sample')
bhm_vs_eddratio.add_point(MBH_EHVOR, Redd_EHVOR, 'green', 'EHVO found by TODO')
bhm_vs_eddratio.add_point(highZmass, highZedd, 'red', 'QSO z>7.0', marker='*')

bhm_vs_eddratio.display_graph()
