import numpy as np
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt
from draw_histogram import plot_zem_histograms
import math


specnum1=6743
specnum3=40

infoDR9 = np.loadtxt('/Users/mikelcharles/Documents/GitHub/DR16Q/DR9Q_selection_minus17.dat',dtype=bytes, delimiter="\n").astype(str)
infoEHVO = np.loadtxt('/Users/mikelcharles/Documents/GitHub/DR16Q/EHVO_DR9.dat',dtype=bytes,delimiter="\n").astype(str)

#----- Start arrays/lists:

countBALinEHVO=0

# --- Arrays in DR9

zem=np.zeros(specnum1)-1
PLATE_DR9_in9=np.zeros(specnum1,dtype=int)-1
MJD_DR9_in9=np.zeros(specnum1,dtype=int)-1
FIBERID_DR9_in9=np.zeros(specnum1,dtype=int)-1


BAL_DR9_in9=[]#
AI_DR9_in9=[]
vmin_DR9_in9=[]
vmax_DR9_in9=[]
EW_DR9_in9=[]



# -- Arrays in EHVO

QUASARNAME_EHVO=[]
EHVOBALS=[]

PLATE_inEHVO=[]
MJD_inEHVO=[]
FIBERID_inEHVO=[]

vmax_inEHVO=[]
vmin_inEHVO=[]
BI_inEHVO=[]
EW_inEHVO=[]
depth_inEHVO=[]


zem_EHVO=[]



# ---


zem_BAL=[]

vmax_BAL=[]
vmin_BAL=[]
BI_BAL=[]
EW_BAL=[]
depth_BAL=[]

    
# Read EHVO:
     
for m in range(0,specnum3):
    ee=infoEHVO[m]
    columns3=ee.split()
    
    QUASARNAME_EHVO.append(str(columns3[0]))
    PLATE_inEHVO.append(int(columns3[1]))
    MJD_inEHVO.append(int(columns3[2]))
    FIBERID_inEHVO.append(int(columns3[3]))
   
    vmax_inEHVO.append(int(columns3[4]))
    vmin_inEHVO.append(int(columns3[5]))
    BI_inEHVO.append(int(columns3[6]))
    EW_inEHVO.append(int(columns3[7]))
    depth_inEHVO.append(float(columns3[8]))

    
# Read DR9Q:
    
bbb=0
aaa=0

for i in range(0,specnum1):
    hh=infoDR9[i]
    columns=hh.split()
    
    zem[i]=float(columns[11])     
    BAL_DR9_in9.append(int(columns[51]))#
    AI_DR9_in9.append(float(columns[54]))
    vmin_DR9_in9.append(float(columns[61]))  #vmin is vmax
    vmax_DR9_in9.append(float(columns[60]))  #vmax is vmin

    EW_DR9_in9.append(float(columns[66]))
    
    PLATE_DR9_in9[i]=int(columns[4])
    MJD_DR9_in9[i]=int(columns[5])
    FIBERID_DR9_in9[i]=int(columns[6])
    
               
    if (PLATE_DR9_in9[i] != -1):
        nn,=np.where((PLATE_inEHVO[:] == PLATE_DR9_in9[i]) & (MJD_inEHVO[:] == MJD_DR9_in9[i]) & (FIBERID_inEHVO[:] == FIBERID_DR9_in9[i]))
        if (len(nn) != 0):
            zem_EHVO.append(zem[i])

    if ((BAL_DR9_in9[i] == 1) and (AI_DR9_in9[i] > 0) and (vmin_DR9_in9[i] < -5000)):
        zem_BAL.append(zem[i])
        vmax_BAL.append(vmax_DR9_in9[i])
        vmin_BAL.append(vmin_DR9_in9[i])
        BI_BAL.append(AI_DR9_in9[i])
        EW_BAL.append(EW_DR9_in9[i])
        
        qq,=np.where(np.array(QUASARNAME_EHVO) == columns[0])
        if (len(qq) != 0):
            countBALinEHVO = countBALinEHVO+1
            EHVOBALS.append(columns[0])

####################################################
zem = np.hstack(zem)
zem_BAL = np.hstack(zem_BAL)
zem_EHVO = np.hstack(zem_EHVO)
print(len(zem_BAL))
x_vals = [zem, zem_BAL, zem_EHVO]

# iqr = np.subtract(*np.percentile(zem, [75, 25]))
# nhist = (max(zem)-min(zem))/(2*iqr*(len(zem)**(-1/3)))
# bins = np.linspace(min(zem), max(zem), int(nhist))
# nhist = (max(zem)-min(zem))/((3.49 * np.std(zem)) * (len(zem)**(-1/3))) ## NEW (Scott)
nhist = (1+(3.22*math.log(len(zem))))
print(nhist)
bins = np.linspace(min(zem), max(zem), int(nhist)) ## NEW

color = ['black', 'blue', 'red']
label = ['parent', 'BAL', 'EHVO']
zem_BAL_size = 1
zem_EHVO_size = 1

x_label = '$z_{em}$'
y_label = 'Number'

# plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'zem', 'right', 'zem_DR9_test_1.png')
# plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'zem/class', 'right', 'zem_DR9_test_2_all_EHVO.png')
# plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'zem/parent', 'left', 'zem_DR9_test_3_sturge.png')
plot_zem_histograms(x_vals, bins, color, label, zem_BAL_size, zem_EHVO_size, x_label, 'cdf', 'keft', 'cdf_DR9.png')
