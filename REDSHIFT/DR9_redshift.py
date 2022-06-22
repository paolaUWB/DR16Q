import os
import math
import numpy as np
from scipy.stats import ks_2samp
from draw_histogram import plot_zem_histograms

##---------- inputs/outputs to change
PARENT_SAMPLE_DR9 = 6743  ## number of spectra in parent sample (previously specnum1)
EHVO_SAMPLE_DR9 = 40 ## number of EHVO spectra (previously specnum3)


##------ data files

# DR9 ---
info_DR9Q = np.loadtxt(os.getcwd() + '/DR9Q_selection_minus17.dat',dtype=bytes, delimiter="\n").astype(str)
info_EHVO_DR9Q = np.loadtxt(os.getcwd() + '/EHVO_DR9.dat',dtype=bytes,delimiter="\n").astype(str)


##------ output locations
OUT_CONFIG = os.getcwd() + '/REDSHIFT/OUTPUT_FILES/'


##------ arrays: DR9 parent
zem_DR9 = np.zeros(PARENT_SAMPLE_DR9) - 1
PLATE_DR9 = np.zeros(PARENT_SAMPLE_DR9, dtype=int) - 1
MJD_DR9 = np.zeros(PARENT_SAMPLE_DR9, dtype=int) - 1
FIBERID_DR9 = np.zeros(PARENT_SAMPLE_DR9, dtype=int) - 1
BAL_DR9 = []
AI_DR9 = []
VMIN_DR9 = []
VMAX_DR9 = []
EW_DR9 = []


##------ arrays: DR9 EHVO
zem_EHVO_DR9 = []

#-- needed for BAL def?
QUASAR_NAME_EHVO_DR9 = []
EHVO_BALS_DR9 = [] 
PLATE_EHVO_DR9 = []
MJD_EHVO_DR9 = []
FIBERID_EHVO_DR9 = []
VMAX_EHVO_DR9 = []
VMIN_EHVO_DR9 = []
BI_EHVO_DR9 = [] # BI = balnicity index
EW_EHVO_DR9 = [] # EW = equivalent width
DEPTH_EHVO_DR9 = []


##------ arrays: DR9 BALS 
zem_BAL_DR9 = []

#-- needed for BAL def
VMAX_BAL_DR9 = []
VMIN_BAL_DR9 = []
BI_BAL_DR9 = []
EW_BAL_DR9 = []
DEPTH_BAL_DR9 = []


##------ read DR9 EHVO 
for m in range(0, EHVO_SAMPLE_DR9):
    DR9_EHVO_row = info_EHVO_DR9Q[m]
    DR9_EHVO_columns = DR9_EHVO_row.split()
    
    QUASAR_NAME_EHVO_DR9.append(str(DR9_EHVO_columns[0]))
    PLATE_EHVO_DR9.append(int(DR9_EHVO_columns[1]))
    MJD_EHVO_DR9.append(int(DR9_EHVO_columns[2]))
    FIBERID_EHVO_DR9.append(int(DR9_EHVO_columns[3]))
   
    VMAX_EHVO_DR9.append(int(DR9_EHVO_columns[4]))
    VMIN_EHVO_DR9.append(int(DR9_EHVO_columns[5]))
    BI_EHVO_DR9.append(int(DR9_EHVO_columns[6]))
    EW_EHVO_DR9.append(int(DR9_EHVO_columns[7]))
    DEPTH_EHVO_DR9.append(float(DR9_EHVO_columns[8]))

##------ read DR9 data
countBALinEHVO=0

for i in range(0, PARENT_SAMPLE_DR9):
    DR9_row = info_DR9Q[i]
    DR9_columns = DR9_row.split()
    
    zem_DR9[i] = float(DR9_columns[11])     
    BAL_DR9.append(int(DR9_columns[51]))
    AI_DR9.append(float(DR9_columns[54]))
    VMIN_DR9.append(float(DR9_columns[61]))  #vmin is vmax in Paris table
    VMAX_DR9.append(float(DR9_columns[60]))  #vmax is vmin in Paris table

    EW_DR9.append(float(DR9_columns[66]))
    
    PLATE_DR9[i]=int(DR9_columns[4])
    MJD_DR9[i]=int(DR9_columns[5])
    FIBERID_DR9[i]=int(DR9_columns[6])

    if (PLATE_DR9[i] != -1):
        nn,=np.where((PLATE_EHVO_DR9[:] == PLATE_DR9[i]) & (MJD_EHVO_DR9[:] == MJD_DR9[i]) & (FIBERID_EHVO_DR9[:] == FIBERID_DR9[i]))
    if (len(nn) != 0):
        zem_EHVO_DR9.append(zem_DR9[i])

    if ((BAL_DR9[i] == 1) and (AI_DR9[i] > 0) and (VMIN_DR9[i] < -5000)):
        zem_BAL_DR9.append(zem_DR9[i])
        VMAX_BAL_DR9.append(VMAX_DR9[i])
        VMAX_BAL_DR9.append(VMIN_DR9[i])
        BI_BAL_DR9.append(AI_DR9[i])
        EW_BAL_DR9.append(EW_DR9[i])
        
        qq,=np.where(np.array(QUASAR_NAME_EHVO_DR9) == DR9_columns[0])
        if (len(qq) != 0):
            countBALinEHVO = countBALinEHVO+1
            EHVO_BALS_DR9.append(DR9_columns[0])

##------ plot zem
zem_DR9 = np.hstack(zem_DR9)
zem_EHVO_DR9 = np.hstack(zem_EHVO_DR9)
zem_BAL_DR9 = np.hstack(zem_BAL_DR9)


#-- Sturges - determines bin width based on parent sample
nhist = (1+(3.22*math.log(len(zem_DR9)))) # number of bins (float)
bins = np.linspace(min(zem_DR9), max(zem_DR9), int(nhist)) # bin edges

#-- inputs for plotting function
x_vals = [zem_EHVO_DR9, zem_DR9, zem_BAL_DR9]
color = ['xkcd:shocking pink', 'black', 'xkcd:purpleish blue']

label = ['EHVO', 'parent', 'BAL']

#-- scaling values for plotting 'zem' 
zem_BAL_size = 1 # set to 1 if no BAL data
zem_EHVO_size = 1

x_label = '$z_{em}$'

plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'zem', 'right', OUT_CONFIG + '/zem_DR9_zem.png') # plots redshifts for each population (parent, EHVO, BAL)
plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'zem/class', 'right', OUT_CONFIG + '/zem_DR9_zem-class.png') # plots redshifts for each population, scaled per bin by number in each population/class per bin
plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'zem/parent', 'left', OUT_CONFIG + '/zem_DR9_zem-parent.png') # plots redshifts for each population, scaled per bin by number of parent sample per bin
plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'cdf', 'left', OUT_CONFIG + '/zem_DR9_cdf.png') # plots cumulative histogram of redshifts for each population (parent, EHVO, BAL)


#------ ks 2 sample statistical test - DR9
KS_EHVOall_zem=ks_2samp(zem_DR9, zem_EHVO_DR9)
KS_BALall_zem=ks_2samp(zem_DR9, zem_BAL_DR9)
KS_EHVOBAL_zem=ks_2samp(zem_BAL_DR9, zem_EHVO_DR9)

print('zem : p-value parent EHVO= '+str(KS_EHVOall_zem))
print('zem: p-value BAL parent= '+str(KS_BALall_zem))
print('zem: p-value BAL EHVO= '+str(KS_EHVOBAL_zem))


#------ to compare number of EHVOs in different zem ranges
'''
# this prints out spectra less than, greater than, or between certain redshift ranges and the percentage of EHVOs in the parent sample within that range

zem_EHVO_DR9_less = []
zem_DR9_less = []

zem_EHVO_DR9_more = []
zem_DR9_more = []

zem_EHVO_DR9_between = []
zem_DR9_between = []


for i in np.arange(2.0,5.0,0.5):
    for j in range(len(zem_EHVO_DR9)): 
        if (zem_EHVO_DR9[j] < i): 
            zem_EHVO_DR9_less.append(zem_EHVO_DR9[j])
    print('zem EHVO < ', i, '=', len(zem_EHVO_DR9_less))
    for k in range(len(zem_DR9)):
        if (zem_DR9[k] < i): 
            zem_DR9_less.append(zem_DR9[k])
    print('zem < ', i, '=', len(zem_DR9_less))
    print('% EHVO in parent sample for zem < ', i, '=', len(zem_EHVO_DR9_less)/len(zem_DR9_less))
    zem_EHVO_DR9_less = []
    zem_DR9_less=[]
    for l in range(len(zem_EHVO_DR9)): 
        if (zem_EHVO_DR9[l] >= i):
            zem_EHVO_DR9_more.append(zem_EHVO_DR9[l])
    print('zem EHVO >= ', i, '=', len(zem_EHVO_DR9_more))
    for m in range(len(zem_DR9)):
        if (zem_DR9[m] >= i):
            zem_DR9_more.append(zem_DR9[m])
    print('zem >= ', i, '=', len(zem_DR9_more))
    print('% EHVO in parent sample for zem >= ', i, '=', len(zem_EHVO_DR9_more)/len(zem_DR9_more))
    zem_EHVO_DR9_more = []
    zem_DR9_more = []
    if (i != 2.0):
        for n in range(len(zem_EHVO_DR9)): 
            if (zem_EHVO_DR9[n] < i) and (zem_EHVO_DR9[n] >= i-0.5):
                zem_EHVO_DR9_between.append(zem_EHVO_DR9[n])
        print(i-0.5, '=< zem EHVO < ', i, '=', len(zem_EHVO_DR9_between))
        for r in range(len(zem_DR9)):
            if (zem_DR9[r] < i) and (zem_DR9[r] >= i-0.5): 
                zem_DR9_between.append(zem_DR9[r])
        print(i-0.5, '=< zem < ', i, '=', len(zem_DR9_between))
        print('% EHVO in parent sample for ', i-0.5, ' =< zem < ', i, '=', len(zem_EHVO_DR9_between)/len(zem_DR9_between))
        zem_EHVO_DR9_between = []
        zem_DR9_between = []
    


print("number of BI: ", len(BI_EHVO_DR9))
print("zem BAL mean: ", np.mean(zem_BAL_DR9))
print("zem BAL median: ", np.median(zem_BAL_DR9))
print("zem BAL range: ", np.max(zem_BAL_DR9)-np.min(zem_BAL_DR9))
print("zem BI min: ", np.min(BI_EHVO_DR9))
print("zem BI max: ", np.max(BI_EHVO_DR9))
'''