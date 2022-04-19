import numpy as np
from matplotlib import pyplot as plt
import math

def plot_zem_histograms(x, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, plot, legend_location, FILE):
    '''
    Draws redshift histogram plots for DR16Q

    Parameters:
    -----------
    x: list
        List of x data to be plotted in histogram in following order: zem, zem_BAL, zem_EHVO
    bins:
    color: list
        List of colors to plot each x (input corresponding color in the same list index as in x)
    label: list
        List of labels for each x (input corresponding label in the same list index as in x)
    zem_BAL_size: int
        Number of zem_BAL to plot to scale to be visible on plot (i.e. zem_BAL_size = 2 --> zem_BAL_plot = [zem_BAL, zem_BAL])
    zem_EHVO_size:
        Number of zem_EHVO to plot to scale to be visible on plot (i.e. zem_EHVO_size = 2 --> zem_EHVO_plot = [zem_EHVO, zem_EHVO])
    x_label: string
        String to label x axis (i.e. xlabel = '$z_{em}$)
    y_label: string
        String to label y axis (i.e. ylabel = 'Number')
    plot: int
        plot = 'zem' (plot redshifts), plot = 'zem/class' (plot redshifts weighted by tot num in each class), plot = 'zem/parent' (plot redshift weighted by parent sample)
    FILE: file name?

    '''

    zem = x[0]
    zem_BAL = x[1]
    zem_EHVO = x[2]

    if plot == 'zem': 
        fig = plt.figure(1) 
        if zem_BAL_size > 1:
            zem_BAL = [zem_BAL] * zem_BAL_size

        if zem_EHVO_size > 1:
            zem_EHVO = [[zem_EHVO]] * zem_EHVO_size 
        
        x_data = np.array([zem,zem_BAL,zem_EHVO],dtype=object)


        plt.hist(x_data, bins, color=color, label=label, histtype='step')
        plt.legend(loc='upper ' + legend_location)
        plt.xlabel(x_label)
        plt.ylabel('Number')
        fig.savefig(FILE, dpi=200)

    elif plot == 'zem/class': 
        fig = plt.figure(2) 
        x = np.array([zem,zem_BAL,zem_EHVO],dtype=object)
        weights = [np.ones_like(zem)/float(len(zem)),np.ones_like(zem_BAL)/float(len(zem_BAL)),np.ones_like(zem_EHVO)/float(len(zem_EHVO))]

        plt.hist(x, bins, color=color, label=label, histtype='step', weights=weights)
        plt.legend(loc='upper ' + legend_location)
        plt.xlabel(x_label)
        plt.ylabel('$n/N_{tot}$')
        fig.savefig(FILE, dpi=200)
        
    elif plot == 'zem/parent': 
        hist_vals_PARENT, bins_PARENT, patches_PARENT = plt.hist(x[0], bins, color=color[0], label=label[0], histtype='step') 
        hist_vals_BALS, bin_edges_BALS, patchesBALS = plt.hist(x[1], bins, color=color[1], label=label[1], histtype='step')
        hist_vals_EHVO, bin_edges_EHVO, patchesEHVO = plt.hist(x[2], bins, color=color[2], label=label[2], histtype='step')

        hist_BALS = hist_vals_BALS 
        hist_EHVO = hist_vals_EHVO 

        # gives indeces of histogram where the parent sample is nonzero
        zem_index_nonzero, = np.where(hist_vals_PARENT != 0) 

        hist_BALS[zem_index_nonzero] = hist_vals_BALS[zem_index_nonzero]/hist_vals_PARENT[zem_index_nonzero] 
        hist_EHVO[zem_index_nonzero] = hist_vals_EHVO[zem_index_nonzero]/hist_vals_PARENT[zem_index_nonzero] 
        
        zem_index_zero, = np.where(hist_vals_PARENT == 0) 
        
        hist_BALS[zem_index_zero] = 0 
        hist_EHVO[zem_index_zero] = 0 

        hist_BALS = np.hstack(hist_BALS) 
        hist_EHVO = np.hstack(hist_EHVO) 


        fig=plt.figure(3)

        print(len(hist_EHVO))
        print(len(bin_edges_BALS))

        plt.step(bin_edges_BALS[1:len(bin_edges_EHVO)], hist_BALS, color=color[1], label=label[1])
        plt.step(bin_edges_EHVO[1:len(bin_edges_EHVO)], hist_EHVO, color=color[2], label=label[2])

        weights = [np.ones_like(zem_BAL)/float(len(zem)), np.ones_like(zem_EHVO)/float(len(zem))]

        plt.plot([bin_edges_BALS[0],bin_edges_BALS[0]],[0,hist_BALS[0]], color='blue')
        plt.plot([bin_edges_BALS[0],bin_edges_BALS[1]],[hist_BALS[0],hist_BALS[0]], color='blue')
        # plt.plot([bin_edges_EHVO[33],bin_edges_EHVO[33]],[0,hist_EHVO[32]], color='red')
        plt.plot([bin_edges_EHVO[len(hist_EHVO)-1],bin_edges_EHVO[len(hist_EHVO)-1]],[0,hist_EHVO[len(hist_EHVO)-2]], color='red')
        plt.plot([bin_edges_EHVO[len(hist_EHVO)],bin_edges_EHVO[len(hist_EHVO)]],[0,hist_EHVO[len(hist_EHVO)-1]], color='red') ## Comment out when old nhist used
        plt.plot([bin_edges_BALS[0],bin_edges_BALS[1]],[hist_EHVO[0],hist_EHVO[0]], color='red')

        plt.xlabel('$z_{em}$')
        plt.ylabel('$n/N_{parent}$')
        plt.legend(loc='upper left')

        fig.savefig(FILE, dpi=200)
    
    elif plot=='cdf': 
        fig = plt.figure(4)
        zem_sorted = np.sort(zem)
        zem_BAL_sorted = np.sort(zem_BAL)
        zem_EHVO_sorted = np.sort(zem_EHVO)
        print(zem_sorted)
        cdf_zem = np.arange(len(zem_sorted))/(len(zem_sorted)-1)
        cdf_zemBAL = np.arange(len(zem_BAL_sorted))/(len(zem_BAL_sorted)-1)
        cdf_zemEHVO = np.arange(len(zem_EHVO_sorted))/(len(zem_EHVO_sorted)-1)
        pdf_zem = [zem_sorted/sum(zem_sorted)]
        print(pdf_zem)
        plt.plot(zem_sorted, pdf_zem, color='green')
        plt.plot(zem_sorted, cdf_zem, color='black')
        plt.plot(zem_BAL_sorted, cdf_zemBAL, color='blue')
        plt.plot(zem_EHVO_sorted, cdf_zemEHVO, color='red')
        plt.xlabel('$z_{em}$')
        fig.savefig(FILE, dpi=200)