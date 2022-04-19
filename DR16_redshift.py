import math
import numpy as np
from scipy.stats import ks_2samp
from draw_histogram import plot_zem_histograms

######################### INPUTS/OUTPUTS TO CHANGE #########################

PARENT_SAMPLE = 18181 ## number of spectra in parent sample (previously specnum1)
EHVO_SAMPLE = 101 ## number of EHVO spectra (previously specnum3)

################################ DATA FILES ################################

## We need to create .dat files for the DR16 data 
info_DR16Q = np.loadtxt('/Users/mikelcharles/Documents/GitHub/DR16Q/DR9Q_selection_minus17.dat',dtype=bytes, delimiter="\n").astype(str)
info_EHVO_DR16Q = np.loadtxt('/Users/mikelcharles/Documents/GitHub/DR16Q/EHVO_DR9.dat',dtype=bytes,delimiter="\n").astype(str)

countBALinEHVO=0

############################## ARRAYS IN DR9 ###############################

zem = np.zeros(PARENT_SAMPLE)-1
PLATE_DR16_in16, MJD_DR16_in16, FIBERID_DR16_in16 = np.zeros(PARENT_SAMPLE,dtype=int)-1
BAL_DR16_in16, AI_DR16_in16 = []
vmin_DR16_in16, vmax_DR16_in16 = []
EW_DR16_in16 = [] ## Do we need this for plotting redshift?? 

############################## ARRAYS IN EHVO ##############################

QUASAR_NAME_EHVO, EHVO_BALS = [] 
PLATE_inEHVO, MJD_inEHVO, FIBERID_inEHVO = []
vmax_inEHVO, vmin_inEHVO = []
BI_inEHVO, EW_inEHVO, depth_inEHVO = []
zem_EHVO = []

############################## ARRAYS IN BALS ##############################

vmax_BAL, vmin_BAL, depth_BAL = []
zem_BAL=[]

################################# READ EHVO ################################

## Column numbers may need to change --> depends on .dat file for DR16 (maybe check paris+12 to compare??)
for m in range(0,EHVO_SAMPLE):
    ee=info_EHVO_DR16Q[m]
    columns3=ee.split()
    
    QUASAR_NAME_EHVO.append(str(columns3[0]))
    PLATE_inEHVO.append(int(columns3[1]))
    MJD_inEHVO.append(int(columns3[2]))
    FIBERID_inEHVO.append(int(columns3[3]))
   
    vmax_inEHVO.append(int(columns3[4]))
    vmin_inEHVO.append(int(columns3[5]))
    BI_inEHVO.append(int(columns3[6]))
    EW_inEHVO.append(int(columns3[7]))
    depth_inEHVO.append(float(columns3[8]))

################################ READ DR16Q ################################

for i in range(0,PARENT_SAMPLE):
    hh=info_DR16Q[i]
    columns=hh.split()
    
    zem[i]=float(columns[11])     
    BAL_DR16_in16.append(int(columns[51]))#
    AI_DR16_in16.append(float(columns[54]))
    vmin_DR16_in16.append(float(columns[61]))  #vmin is vmax --> may not be true for Lyke+2020
    vmax_DR16_in16.append(float(columns[60]))  #vmax is vmin --> may not be true for Lyke+2020

    EW_DR16_in16.append(float(columns[66]))
    
    PLATE_DR16_in16[i]=int(columns[4])
    MJD_DR16_in16[i]=int(columns[5])
    FIBERID_DR16_in16[i]=int(columns[6])
    
    ## Identifying BALs --> likely will need to be updated for Lyke+2020 table           
    if (PLATE_DR16_in16[i] != -1):
        nn,=np.where((PLATE_inEHVO[:] == PLATE_DR16_in16[i]) & (MJD_inEHVO[:] == MJD_DR16_in16[i]) & (FIBERID_inEHVO[:] == FIBERID_DR16_in16[i]))
        if (len(nn) != 0):
            zem_EHVO.append(zem[i])

    if ((BAL_DR16_in16[i] == 1) and (AI_DR16_in16[i] > 0) and (vmin_DR16_in16[i] < -5000)):
        zem_BAL.append(zem[i])
        vmax_BAL.append(vmax_DR16_in16[i]) ## are these necessary for redshift plots?
        vmin_BAL.append(vmin_DR16_in16[i]) ## are these necessary for redshift plots?
        
        qq,=np.where(np.array(QUASAR_NAME_EHVO) == columns[0])
        if (len(qq) != 0):
            countBALinEHVO = countBALinEHVO+1
            EHVO_BALS.append(columns[0])

################################# PLOT ZEM #################################

zem = np.hstack(zem)
zem_BAL = np.hstack(zem_BAL)
zem_EHVO = np.hstack(zem_EHVO)

x_vals = [zem, zem_BAL, zem_EHVO]

### Freedman & Diaconis
# iqr = np.subtract(*np.percentile(zem, [75, 25]))
# nhist = (max(zem)-min(zem))/(2*iqr*(len(zem)**(-1/3)))
# bins = np.linspace(min(zem), max(zem), int(nhist))

### Freedman & Diaconis (scaled by 0.75)
# bins = np.linspace(min(zem), max(zem), int(nhist * 0.75))

### Scott
# nhist = (max(zem)-min(zem))/((3.49 * np.std(zem)) * (len(zem)**(-1/3)))
# bins = np.linspace(min(zem), max(zem), int(nhist))

### Sturges
nhist = (1+(3.22*math.log(len(zem))))
bins = np.linspace(min(zem), max(zem), int(nhist))


color = ['black', 'blue', 'red']
label = ['parent', 'BAL', 'EHVO']
zem_BAL_size = 1
zem_EHVO_size = 1

x_label = '$z_{em}$'
y_label = 'Number'

# plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'zem', 'right', 'zem_DR9_test_1.png')
# plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'zem/class', 'right', 'zem_DR9_test_2_all_EHVO.png')
# plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'zem/parent', 'left', 'zem_DR9_test_3_sturge.png')
plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'cdf', 'left', 'cdf_DR9.png')


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