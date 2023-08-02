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
    ##------ data
    zem_EHVO = x[0] # EHVO redshift
    zem_EHVO = np.hstack(zem_EHVO)
    zem = x[1] # parent sample redshift
    zem = np.hstack(zem)
    zem_BAL = x[2] # BAL redshift
    zem_BAL = np.hstack(zem_BAL)


    ##------ plots redshift histograms w/ option of scaling BAL & EHVO to be able to see on plot
    if plot == 'zem': 
        fig = plt.figure(1) 

        if zem_BAL_size > 1:
            zem_BAL = [zem_BAL] * zem_BAL_size

        if zem_EHVO_size > 1:
            zem_EHVO = [zem_EHVO] * zem_EHVO_size 
        
        x_data = np.array([zem_EHVO, zem, zem_BAL], dtype=object) # check to see if there is a solution to not do this

        #--plots histograms w/ same bin size (uses parent sample bin size)
        plt.hist(x_data, bins=bins, color=color, label=label, histtype='step', linewidth=1.5)

        #-- plots histograms w/ different bin sizes for each class [bins input needs to be array of bin sizes]
        # plt.hist(x_data[0], bins=bins[0], color=color[0], label=label[0], histtype='step', linewidth=1.5) # EHVO 
        # plt.hist(x_data[1], bins=bins[1], color=color[1], label=label[1], histtype='step', linewidth=1.5) # parent sample
        # plt.hist(x_data[2], bins=bins[2], color=color[2], label=label[2], histtype='step', linewidth=1.5) # BAL

        plt.legend(loc='upper ' + legend_location)
        plt.xlabel(x_label)
        plt.ylabel('Number')
        
        fig.savefig(FILE, dpi=1200)


    ##------ plots redshift histograms --> scaled by dividing each class by itself
    elif plot == 'zem/class': 
        fig = plt.figure(2) 
        x = np.array([zem_EHVO, zem, zem_BAL],dtype=object)

        weights = [np.ones_like(zem_EHVO)/float(len(zem_EHVO)), np.ones_like(zem)/float(len(zem)), np.ones_like(zem_BAL)/float(len(zem_BAL))]

        #------ plots histograms w/ same bin size (uses parent sample bin size)
        plt.hist(x, bins, color=color, label=label, histtype='step', weights=weights, linewidth=1.5)
        
        #------ plots histograms w/ different bin sizes for each class [bins input needs to be array of bin sizes]
        # plt.hist(x[0], bins=bins[0], color=color[0], label=label[0], histtype='step', weights=weights[0], linewidth=1.5)
        # plt.hist(x[1], bins=bins[1], color=color[1], label=label[1], histtype='step', weights=weights[1], linewidth=1.5)
        # plt.hist(x[2], bins=bins[2], color=color[2], label=label[2], histtype='step', weights=weights[2], linewidth=1.5)

        plt.legend(loc='upper ' + legend_location)
        plt.xlabel(x_label)
        plt.ylabel('$n/N_{tot}$')

        handles, labels = plt.gca().get_legend_handles_labels()
        order = [1, 2, 0]
        plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order])
        
        plt.tight_layout()

        fig.savefig(FILE, dpi=1200)
        
    elif plot == 'zem/parent': 
        hist_vals_EHVO, bin_edges_EHVO, patchesEHVO = plt.hist(x[0], bins=bins, color=color[0], label=label[0], histtype='step', linewidth=1.5)
        hist_vals_PARENT, bins_PARENT, patches_PARENT = plt.hist(x[1], bins=bins, color=color[1], label=label[1], histtype='step', linewidth=1.5) 
        hist_vals_BALS, bin_edges_BALS, patchesBALS = plt.hist(x[2], bins=bins, color=color[2], label=label[2], histtype='step', linewidth=1.5)

        hist_EHVO = hist_vals_EHVO
        hist_BALS = hist_vals_BALS  

        #-- gives indeces of histogram where the parent sample is nonzero
        zem_index_nonzero, = np.where(hist_vals_PARENT != 0) 

        hist_EHVO[zem_index_nonzero] = hist_vals_EHVO[zem_index_nonzero]/hist_vals_PARENT[zem_index_nonzero] 
        hist_BALS[zem_index_nonzero] = hist_vals_BALS[zem_index_nonzero]/hist_vals_PARENT[zem_index_nonzero] 

        #-- gives indices of histogram where the parent sample is zero 
        zem_index_zero, = np.where(hist_vals_PARENT == 0) 
        
        hist_EHVO[zem_index_zero] = 0 
        hist_BALS[zem_index_zero] = 0 
 
        hist_EHVO = np.hstack(hist_EHVO) 
        hist_BALS = np.hstack(hist_BALS)

        fig=plt.figure(3)

        plt.step(bin_edges_BALS[1:len(bin_edges_EHVO)], hist_BALS, color=color[2], label=label[2], linewidth=1.5)
        plt.step(bin_edges_EHVO[1:len(bin_edges_EHVO)], hist_EHVO, color=color[0], label=label[0], linewidth=1.5)

        weights = [ np.ones_like(zem_EHVO)/float(len(zem)), np.ones_like(zem_BAL)/float(len(zem))]

        plt.plot([bin_edges_BALS[0],bin_edges_BALS[0]],[0,hist_BALS[0]], color=color[2], linewidth=1.5)
        plt.plot([bin_edges_BALS[0],bin_edges_BALS[1]],[hist_BALS[0],hist_BALS[0]], color=color[2], linewidth=1.5)
        plt.plot([bin_edges_EHVO[len(hist_EHVO)-1],bin_edges_EHVO[len(hist_EHVO)-1]],[0,hist_EHVO[len(hist_EHVO)-2]], color=color[0], linewidth=1.5)
        plt.plot([bin_edges_EHVO[len(hist_EHVO)],bin_edges_EHVO[len(hist_EHVO)]],[0,hist_EHVO[len(hist_EHVO)-1]], color=color[0], linewidth=1.5) ## Comment out when old nhist used
        plt.plot([bin_edges_BALS[0],bin_edges_BALS[1]],[hist_EHVO[0],hist_EHVO[0]], color=color[0], linewidth=1.5)
        plt.plot([bin_edges_BALS[0],bin_edges_BALS[0]],[0,hist_EHVO[0]], color=color[0], linewidth=1.5)

        plt.xlabel('$z_{em}$')
        plt.ylabel('$n/N_{parent}$')

        handles, labels = plt.gca().get_legend_handles_labels()
        order = [1, 0]
        plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc = 'upper left')
        
        plt.tight_layout()

        fig.savefig(FILE, dpi=1200)
    
    elif plot=='cdf': 
        fig = plt.figure(4)
        
        plt.hist(zem_EHVO, bins, density=True, histtype='step', cumulative=True, label='EHVO', color=color[0], linewidth=1.5)
        plt.hist(zem, bins, density=True, histtype='step', cumulative=True, label='parent', color=color[1], linewidth=1.5)
        plt.hist(zem_BAL, bins, density=True, histtype='step', cumulative=True, label='BAL', color=color[2], linewidth=1.5)
        
        plt.xlabel('$z_{em}$')

        handles, labels = plt.gca().get_legend_handles_labels()
        order = [1, 0, 2]
        plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc = 'upper left')
       
        plt.tight_layout()
       
        fig.savefig(FILE, dpi=1200)