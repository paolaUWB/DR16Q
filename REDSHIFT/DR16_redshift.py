import os
import sys
import math
import numpy as np
from scipy.stats import ks_2samp, norm
from draw_histogram import plot_zem_histograms
from utility_functions import read_file

# import seaborn as sns
# import matplotlib.pyplot as plt

##---------- inputs/outputs to change
PARENT_SAMPLE = 18171  ## number of spectra in parent sample
EHVO_SAMPLE = 98 ## number of EHVO spectra
BAL_SAMPLE = 1913

specnum1 = 6743
specnum3 = 40


##------ data files
# DR16 ---
info_DR16Q = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/DR16_parent_sample.csv"
info_EHVO_DR16Q = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/DR16_EHVO_sorted_norm.csv"
info_BAL_DR16Q = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/BAL_parent_sample.csv"

# DR9 ---
infoDR9 = np.loadtxt('/Users/mikelcharles/Documents/GitHub/DR16Q/DR9Q_selection_minus17.dat',dtype=bytes, delimiter="\n").astype(str)
infoEHVO = np.loadtxt('/Users/mikelcharles/Documents/GitHub/DR16Q/EHVO_DR9.dat',dtype=bytes,delimiter="\n").astype(str)


##------ arrays: DR16 parent
zem = np.zeros(PARENT_SAMPLE)-1
PLATE_DR16_in16 = np.zeros(PARENT_SAMPLE,dtype=int)-1
MJD_DR16_in16 = np.zeros(PARENT_SAMPLE,dtype=int)-1
FIBERID_DR16_in16 = np.zeros(PARENT_SAMPLE,dtype=int)-1
BAL_DR16_in16 = []
AI_DR16_in16 = []
vmin_DR16_in16 = []
vmax_DR16_in16 = []
EW_DR16_in16 = [] ## Do we need this for plotting redshift?? 


##------ arrays: DR16 EHVO
zem_EHVO = []

#-- needed for BAL def?
QUASAR_NAME_EHVO = []
EHVO_BALS = [] 
PLATE_inEHVO = []
MJD_inEHVO = []
FIBERID_inEHVO = []
vmax_inEHVO = []
vmin_inEHVO = []
BI_inEHVO = [] # BI = balnicity index
EW_inEHVO = [] # EW = equivalent width
depth_inEHVO = []


##------ arrays: DR16 BALS 
zem_BAL = []

#-- needed for BAL def
vmax_BAL = []
vmin_BAL = []
depth_BAL = []
bi_ratio_BAL = []


##------ read DR16 data (parent & EHVO)
zem, snr, spectra_list = read_file(info_DR16Q) # parent
zem_EHVO, snr_EHVO, spectra_list_EHVO = read_file(info_EHVO_DR16Q) # EHVO
zem_BAL, bi_ratio_BAL, spectra_list_BAL = read_file(info_BAL_DR16Q) # BAL


# with open(info_BAL_DR16Q) as f:  
#         for line in f:
#             each_row_in_file = line.split(" ")
#             spectra_list.append(each_row_in_file[0])
#             zem_BAL.append(np.float(each_row_in_file[1]))



##------ read BAL data
# for m in range(0, BAL_SAMPLE):
#     DR16_BAL_row = info_BAL_DR16Q[m]
#     DR16_BAL_columns = DR16_BAL_row.split()
    
#     zem_BAL.append(str(DR16_BAL_columns[1]))

# print('zem: ', len(zem_BAL))

##------ arrays: DR9 parent
# zem_dr9 = np.zeros(specnum1)-1
# PLATE_DR9_in9 = np.zeros(specnum1,dtype=int)-1
# MJD_DR9_in9 = np.zeros(specnum1,dtype=int)-1
# FIBERID_DR9_in9 = np.zeros(specnum1,dtype=int)-1
# BAL_DR9_in9 = []
# AI_DR9_in9 = []
# vmin_DR9_in9 = []
# vmax_DR9_in9 = []
# EW_DR9_in9 = []


# ##------ read DR9 data
# for i in range(0,specnum1):
#     row_DR9 = infoDR9[i]
#     columns_DR9 = row_DR9.split()
    
#     zem_dr9[i] = float(columns_DR9[11])     
#     BAL_DR9_in9.append(int(columns_DR9[51]))
#     AI_DR9_in9.append(float(columns_DR9[54]))
#     vmin_DR9_in9.append(float(columns_DR9[61]))  #vmin is vmax in Paris table
#     vmax_DR9_in9.append(float(columns_DR9[60]))  #vmax is vmin in Paris table

#     EW_DR9_in9.append(float(columns_DR9[66]))
    
#     PLATE_DR9_in9[i]=int(columns_DR9[4])
#     MJD_DR9_in9[i]=int(columns_DR9[5])
#     FIBERID_DR9_in9[i]=int(columns_DR9[6])



##------ arrays: DR9 EHVO
# zem_EHVO_DR9 = []

#-- needed for BAL def?
# QUASAR_NAME_EHVO_DR9 = []
# EHVO_BALS_DR9 = [] 
# PLATE_EHVO_DR9 = []
# MJD_EHVO_DR9 = []
# FIBERID_EHVO_DR9 = []
# VMAX_EHVO_DR9 = []
# VMIN_EHVO_DR9 = []
# BI_EHVO_DR9 = [] # BI = balnicity index
# EW_EHVO_DR9 = [] # EW = equivalent width
# DEPTH_EHVO_DR9 = []


##------ arrays: DR9 BALS 
# zem_BAL_DR9 = []

# #-- needed for BAL def
# VMAX_BAL_DR9 = []
# VMIN_BAL_DR9 = []
# BI_BAL_DR9 = []
# EW_BAL_DR9 = []
# DEPTH_BAL_DR9 = []


zem_DR9 = np.zeros(specnum1) - 1
PLATE_DR9 = np.zeros(specnum1, dtype=int) - 1
MJD_DR9 = np.zeros(specnum1, dtype=int) - 1
FIBERID_DR9 = np.zeros(specnum1, dtype=int) - 1
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
for m in range(0, specnum3):
    DR9_EHVO_row = infoEHVO[m]
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

for i in range(0, specnum1):
    DR9_row = infoDR9[i]
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
zem = np.hstack(zem)
zem_EHVO = np.hstack(zem_EHVO)
zem_BAL = np.hstack(zem_BAL)

#-- Sturges - determines bin width based on parent sample
nhist = (1+(3.22*math.log(len(zem)))) # number of bins (float)
bins = np.linspace(min(zem), max(zem), int(nhist)) # bin edges
print(bins)
# bins_parent = np.linspace(min(zem), max(zem), int(nhist)) # bin edges

#-- Sturges - determines bin width based on EHVO sample --> decide if we want to keep this or delete this 
# nhist_EHVO = (1+(3.22*math.log(len(zem_EHVO))))
# bins_EHVO = np.linspace(min(zem_EHVO), max(zem_EHVO), int(nhist)) ## number of bins or width of bins (change variable name)

#-- inputs for plotting function
# x_vals = [zem_EHVO, zem] # only plotting EHVO & parent --> delete once BAL def completed
x_vals = [zem_EHVO, zem, zem_BAL]
# x_vals_parent = [zem_DR9, zem, zem_BAL_DR9, zem_BAL, zem_EHVO_DR9, zem_EHVO] # DR9 & DR16 parent samples

# color = ['red', 'black']
color = ['xkcd:shocking pink', 'black', 'xkcd:purpleish blue']
# color_parent = ['black', 'red']

# label = ['EHVO', 'parent']
label = ['EHVO', 'parent', 'BAL']
# label_parent = ['DR9 Parent', 'DR16 Parent', 'DR9 BAL', 'DR16 BAL', 'DR9 EHVO', 'DR16 EHVO']
bins_parent = bins
# color_parent = ['xkcd:apple', 'black', 'xkcd:tangerine', 'red', 'xkcd:barney', 'blue']
# bins = [bins_EHVO, bins_parent]

zem_BAL_size = 1 # set to 1 if no BAL data
zem_EHVO_size = 1

x_label = '$z_{em}$'


# fig = plt.figure(4)

# zem_DR9_sorted = np.sort(zem_DR9)
# zem_DR16_sorted = np.sort(zem)

# zem_BAL_DR9_sorted = np.sort(zem_BAL_DR9)
# zem_BAL_DR16_sorted = np.sort(zem_BAL)

# zem_EHVO_DR9_sorted = np.sort(zem_EHVO_DR9)
# zem_EHVO_DR16_sorted = np.sort(zem_EHVO)


# cdf_DR9_zem = np.arange(len(zem_DR9_sorted))/(len(zem_DR9_sorted)-1)
# cdf_DR16_zem = np.arange(len(zem_DR16_sorted))/(len(zem_DR16_sorted)-1)

# cdf_DR9_zemBAL = np.arange(len(zem_BAL_DR9_sorted))/(len(zem_BAL_DR9_sorted)-1)
# cdf_DR16_zemBAL = np.arange(len(zem_BAL_DR16_sorted))/(len(zem_BAL_DR16_sorted)-1)

# cdf_DR9_zemEHVO = np.arange(len(zem_EHVO_DR9_sorted))/(len(zem_EHVO_DR9_sorted)-1)
# cdf_DR16_zemEHVO = np.arange(len(zem_EHVO_DR16_sorted))/(len(zem_EHVO_DR16_sorted)-1)


# plt.plot(zem_DR9_sorted, cdf_DR9_zem, color=color_parent[0], linewidth=1.5, label='DR9 Parent')#, label='dr9 parent')
# plt.plot(zem_DR16_sorted, cdf_DR16_zem, color=color_parent[1], linewidth=1.5, label='DR16 Parent')#, label='dr9 parent')

# plt.plot(zem_BAL_DR9_sorted, cdf_DR9_zemBAL, color=color_parent[2],linewidth=1.5, label='DR9 BAL')#, label='dr9 BAL')
# plt.plot(zem_BAL_DR16_sorted, cdf_DR16_zemBAL, color=color_parent[3],linewidth=1.5, label='DR16 BAL')#, label='dr9 BAL')

# plt.plot(zem_EHVO_DR9_sorted, cdf_DR9_zemEHVO, color=color_parent[4], linewidth=1.5, label='DR9 EHVO')#label='EHVO')
# plt.plot(zem_EHVO_DR16_sorted, cdf_DR16_zemEHVO, color=color_parent[5], linewidth=1.5, label='DR16 EHVO')#label='EHVO')



# plt.xlabel('$z_{em}$')

# handles, labels = plt.gca().get_legend_handles_labels()
# order = [0,1,2,3,4,5]
# plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc = 'upper left')
# plt.tight_layout()
# fig.savefig('zem_DR9_DR16_cdf.png', dpi=200)

# plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'zem', 'right', 'zem_DR16_zem.png')
# plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'zem/class', 'right', 'zem_DR16_zem-class.png')
# plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'zem/parent', 'left', 'zem_DR16_zem-parent.png')
plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'cdf', 'left', 'zem_DR16_cdf.png')

zem_compare = 2.0
zem_compare_2 = 4.5
zem_EHVO_highredshift = []
zem_highredshift = []
for i in range(len(zem_EHVO)): 
    if (zem_EHVO[i] <= zem_compare): # and (zem_EHVO[i] < zem_compare_2):
        zem_EHVO_highredshift.append(zem_EHVO[i])

for i in range(len(zem)):
    if (zem[i] <= zem_compare): # and (zem[i] < zem_compare_2): 
        zem_highredshift.append(zem[i])

print('zem EHVO < ', zem_compare, '=', len(zem_EHVO_highredshift))
print('zem < ', zem_compare, '=', len(zem_highredshift))

# zem_DR16_highredshift = []
# for i in range(len(zem)): 
#     if zem[i] >= zem_compare:
#         zem_DR16_highredshift.append(zem[i])

# print('min zem > ', str(zem_compare), '=', len(zem_DR16_highredshift))


# plt.figure(1)
# sns.kdeplot(zem, shade=True, bw_adjust=0.2, color='c', linewidth=2)

# # plt.boxplot(zem)
# # sns.boxplot(zem)

# # plt.figure(2)
# sns.kdeplot(zem_EHVO, shade=True, bw_adjust=0.2)
# plt.show()

# plt.figure(3)
# # sns.displot(zem, kde=True)
# sns.displot(zem, kind="ecdf")

# plt.figure(4)
# # sns.displot(zem_EHVO, kde=True)
# sns.displot(zem_EHVO, kind="ecdf")

print("number of BALs: ", len(zem_EHVO))
print("zem BAL mean: ", np.mean(zem_BAL))
print("zem BAL median: ", np.median(zem_BAL))
print("zem BAL range: ", np.max(zem_BAL)-np.min(zem_BAL))
print("zem EHVO min: ", np.min(zem_EHVO))
print("zem EHVO max: ", np.max(zem_EHVO))

KS_EHVOall_zem=ks_2samp(zem,zem_EHVO)
KS_BALall_zem=ks_2samp(zem,zem_BAL)
KS_EHVOBAL_zem=ks_2samp(zem_BAL,zem_EHVO)
print('zem : p-value parent EHVO= '+str(KS_EHVOall_zem))
print('zem: p-value BAL parent= '+str(KS_BALall_zem))
print('zem: p-value BAL EHVO= '+str(KS_EHVOBAL_zem))
## ------------------------------------------

'''
###### ZEM ######
# KS_EHVOall_zem=ks_2samp(zem,zem_EHVO)
# KS_BALall_zem=ks_2samp(zem,zem_BAL)
# KS_EHVOBAL_zem=ks_2samp(zem_BAL,zem_EHVO)
# print('zem : p-value parent EHVO= '+str(KS_EHVOall_zem))
# print('zem: p-value BAL parent= '+str(KS_BALall_zem))
# print('zem: p-value BAL EHVO= '+str(KS_EHVOBAL_zem))
'''