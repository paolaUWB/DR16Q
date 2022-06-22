import os
import sys
import math
import numpy as np
from scipy.stats import ks_2samp, norm
from draw_histogram import plot_zem_histograms
sys.path.insert(0, os.getcwd() + '/../' + 'DR16Q')
from utility_functions import read_file


##------ number of spectra in sample
PARENT_SAMPLE = 18171  ## number of spectra in parent sample
EHVO_SAMPLE = 98 ## number of EHVO spectra
BAL_SAMPLE = 1913


##------ data files: DR16Q
info_DR16Q = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/DR16_parent_sample.csv"
info_EHVO_DR16Q = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/DR16Q_EHVO/DR16_EHVO_sorted_norm.csv"
info_BAL_DR16Q = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/DR16_BAL_parent_sample.csv"


##------ output file location
OUT_CONFIG = os.getcwd() + '/REDSHIFT/OUTPUT_FILES/'


##------ read DR16Q data: parent, EHVO, BAL
zem_DR16, snr_DR16, spectra_list_DR16 = read_file(info_DR16Q) # parent
zem_EHVO_DR16, snr_EHVO_DR16, spectra_list_EHVO_DR16 = read_file(info_EHVO_DR16Q) # EHVO
zem_BAL_DR16, bi_ratio_BAL_DR16, spectra_list_BAL_DR16 = read_file(info_BAL_DR16Q) # BAL


##------ plot zem
zem_DR16 = np.hstack(zem_DR16)
zem_EHVO_DR16 = np.hstack(zem_EHVO_DR16)
zem_BAL_DR16 = np.hstack(zem_BAL_DR16)

#-- Sturges - determines bin width based on parent sample
nhist = (1+(3.22*math.log(len(zem_DR16)))) # number of bins (float)
bins = np.linspace(min(zem_DR16), max(zem_DR16), int(nhist)) # bin edges

#-- inputs for plotting function
x_vals = [zem_EHVO_DR16, zem_DR16, zem_BAL_DR16]
color = ['xkcd:shocking pink', 'black', 'xkcd:purpleish blue']
label = ['EHVO', 'parent', 'BAL']

#-- scaling values for plotting 'zem'
zem_BAL_size = 1 # set to 1 if no BAL data
zem_EHVO_size = 1

x_label = '$z_{em}$'

plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'zem', 'right', OUT_CONFIG + 'zem_DR16_zem_test.png')
plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'zem/class', 'right', OUT_CONFIG + 'zem_DR16_zem-class_test.png')
plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'zem/parent', 'left', OUT_CONFIG + 'zem_DR16_zem-parent_test.png')
plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'cdf', 'left', OUT_CONFIG + 'zem_DR16_cdf_test.png')


##------ ks 2 sample statistical test: DR16Q
KS_EHVOall_zem=ks_2samp(zem_DR16,zem_EHVO_DR16)
KS_BALall_zem=ks_2samp(zem_DR16,zem_BAL_DR16)
KS_EHVOBAL_zem=ks_2samp(zem_BAL_DR16,zem_EHVO_DR16)

print('zem : p-value parent EHVO= '+str(KS_EHVOall_zem))
print('zem: p-value BAL parent= '+str(KS_BALall_zem))
print('zem: p-value BAL EHVO= '+str(KS_EHVOBAL_zem))


##------ to compare number of EHVOs in different zem ranges
'''
zem_EHVO_DR16_less = []
zem_DR16_less = []

zem_EHVO_DR16_more = []
zem_DR16_more = []

zem_EHVO_DR16_between = []
zem_DR16_between = []


for i in np.arange(2.0,5.0,0.5):
    for j in range(len(zem_EHVO_DR16)): 
        if (zem_EHVO_DR16[j] < i): 
            zem_EHVO_DR16_less.append(zem_EHVO_DR16[j])
    print('zem EHVO < ', i, '=', len(zem_EHVO_DR16_less))

    for k in range(len(zem_DR16)):
        if (zem_DR16[k] < i): 
            zem_DR16_less.append(zem_DR16[k])
    print('zem < ', i, '=', len(zem_DR16_less))

    print('% EHVO in parent sample for zem < ', i, '=', len(zem_EHVO_DR16_less)/len(zem_DR16_less))

    zem_EHVO_DR16_less = []
    zem_DR16_less=[]

    for l in range(len(zem_EHVO_DR16)): 
        if (zem_EHVO_DR16[l] >= i):
            zem_EHVO_DR16_more.append(zem_EHVO_DR16[l])
    print('zem EHVO >= ', i, '=', len(zem_EHVO_DR16_more))

    for m in range(len(zem_DR16)):
        if (zem_DR16[m] >= i):
            zem_DR16_more.append(zem_DR16[m])
    print('zem >= ', i, '=', len(zem_DR16_more))

    print('% EHVO in parent sample for zem >= ', i, '=', len(zem_EHVO_DR16_more)/len(zem_DR16_more))

    zem_EHVO_DR16_more = []
    zem_DR16_more = []

    if (i != 2.0):
        for n in range(len(zem_EHVO_DR16)): 
            if (zem_EHVO_DR16[n] < i) and (zem_EHVO_DR16[n] >= i-0.5):
                zem_EHVO_DR16_between.append(zem_EHVO_DR16[n])
        print(i-0.5, '=< zem EHVO < ', i, '=', len(zem_EHVO_DR16_between))

        for r in range(len(zem_DR16)):
            if (zem_DR16[r] < i) and (zem_DR16[r] >= i-0.5): 
                zem_DR16_between.append(zem_DR16[r])
        print(i-0.5, '=< zem < ', i, '=', len(zem_DR16_between))

        print('% EHVO in parent sample for ', i-0.5, ' =< zem < ', i, '=', len(zem_EHVO_DR16_between)/len(zem_DR16_between))

        zem_EHVO_DR16_between = []
        zem_DR16_between = []

print("number of BALs: ", len(zem_BAL_DR16))
print("zem BAL mean: ", np.mean(zem_BAL_DR16))
print("zem BAL median: ", np.median(zem_BAL_DR16))
print("zem BAL range: ", np.max(zem_BAL_DR16)-np.min(zem_BAL_DR16))
print("zem EHVO min: ", np.min(zem_EHVO_DR16))
print("zem EHVO max: ", np.max(zem_EHVO_DR16))
'''