########################################################################################################
################## THIS CODE CREATES HISTOGRAM PLOTS FOR DR16 [FOR REDSHIFT ANALYSIS] ##################

import numpy as np
import pandas as pd
import os
import sys
from scipy.stats import ks_2samp
#from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from utility_functions import read_file, read_list_spectra

specnum1=18181 ## NUMBER OF SPECTRA IN PARENT SAMPLE  # data=6760
#specnum2=105783
specnum3=40  ## NUMBER OF EHVOS

#infoDR16 = np.loadtxt('/Users/Paola/QUASAR/Work_EHVO/DATA/DR9Q_selection_minus17.dat',dtype=bytes, delimiter="\n").astype(str)
#infoDR16 = np.loadtxt('/Users/mikelcharles/Documents/GitHub/DR16Q/DR16_sorted_norm.csv')

## FILE THAT CONTAINS DR16 DATA
# infoDR16_all = sys.argv[1] if len(sys.argv) > 1 else "DR16_sorted_norm.csv"
# infoDR16_snr = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/OUTPUT_FILES/NORMALIZATION/parent_sample.csv" ### ALL CASES WITH SNR > 10
infoDR16_parent = sys.argv[1] if len(sys.argv) > 1 else "DR16_parent_sample.csv"
#infoDR16_BAL = sys.argv[1] if len(sys.argv) > 1 else "DR16_BAL.csv"
#infoDR16_EHVO = sys.argv[1] if len(sys.argv) > 1 else "DR16_EHVO.csv"

#spectra_list, zem, calc_snr_list = read_list_spectra(infoDR16_snr, ["NORM SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR"]) 
# zem_orig, snr, spectra_name = read_file(infoDR16)
zem=np.zeros((specnum1)-1) ## Redshift
zem_EHVO=[] ## EHVO Redshift
zem_BAL=[] ## BAL Redshift

# zem, snr, spectra_list = read_file(infoDR16_snr)
# zem_all, snr_all, spectra_list_all = read_file(infoDR16_all)
zem, snr, spectra_list = read_file(infoDR16_parent)
#zem_BAL, snr_BAL, spectra_list_BAL = read_file(infoDR16_BAL)
#zem_EHVO, snr_EHVO, spectra_list_EHVO = read_file(infoDR16_EHVO)


## Read EHVO
# for m in range(0,specnum3):
#     ee=infoEHVO[m]
#     columns3=ee.split()
    
#     QUASARNAME_EHVO.append(str(columns3[0]))
#     PLATE_inEHVO.append(int(columns3[1]))
#     MJD_inEHVO.append(int(columns3[2]))
#     FIBERID_inEHVO.append(int(columns3[3]))
   
#     vmax_inEHVO.append(int(columns3[4]))
#     vmin_inEHVO.append(int(columns3[5]))
#     BI_inEHVO.append(int(columns3[6]))
#     EW_inEHVO.append(int(columns3[7]))
#     depth_inEHVO.append(float(columns3[8]))

## Read DR16

# for i in range(0, specnum1): 
#     index = infoDR16[i]
#     zem_orig = float(infoDR16[i][1])
#     zem.append(zem_orig)

# bbb=0
# aaa=0

# for i in range(0,specnum1):
#     hh=infoDR9[i]
#     columns=hh.split()
    
#     zem[i]=float(columns[11])   
#     orig_zem[i]=float(columns[11])   
#     errzem[i]=float(columns[12])
#     Mi[i]=float(columns[26])
#     BAL_DR9_in9.append(int(columns[51]))
#     AI_DR9_in9.append(float(columns[54]))
#     vmin_DR9_in9.append(float(columns[61]))  #vmin is vmax
#     vmax_DR9_in9.append(float(columns[60]))  #vmax is vmin
#     FIRST_fnu20cm[i]=float(columns[128])
#     SNFIRST_fnu20cm[i]=float(columns[129])

#     EW_DR9_in9.append(float(columns[66]))
    
#     SNR9.append(float(columns[30]))
#     SNR9b.append(float(columns[31]))
#     SNR9c.append(float(columns[32]))
    
#     PLATE_DR9_in9[i]=int(columns[4])
#     MJD_DR9_in9[i]=int(columns[5])
#     FIBERID_DR9_in9[i]=int(columns[6])
    
#     PLATE_DR7_in9[i]=int(columns[22])
#     MJD_DR7_in9[i]=int(columns[23])
#     FIBERID_DR7_in9[i]=int(columns[24])
    
    
#     if (PLATE_DR7_in9[i] != -1):
#         vv,=np.where((PLATE_DR7_inPBH[:] == PLATE_DR7_in9[i]) & (MJD_DR7_inPBH[:] == MJD_DR7_in9[i]) & (FIBERID_DR7_inPBH[:] == FIBERID_DR7_in9[i]))
#         if (len(vv) != 0):  
#             zem[i]=zem_inPBH[vv[0]]
#             errzem[i]=errzem_inPBH[vv[0]]
#             bbb=bbb+1
               
#     if (PLATE_DR9_in9[i] != -1):
#         nn,=np.where((PLATE_inEHVO[:] == PLATE_DR9_in9[i]) & (MJD_inEHVO[:] == MJD_DR9_in9[i]) & (FIBERID_inEHVO[:] == FIBERID_DR9_in9[i]))
#         if (len(nn) != 0):
               
#             MI_EHVO.append(Mi[i])
#             zem_EHVO.append(zem[i])
#             errzem_EHVO.append(errzem[i])
#             FIRST_inEHVO.append(FIRST_fnu20cm[i])
#             SNFIRST_inEHVO.append(SNFIRST_fnu20cm[i])
            
#             SNR9_EHVO.append(SNR9[i])

#     if ((BAL_DR9_in9[i] == 1) and (AI_DR9_in9[i] > 0) and (vmin_DR9_in9[i] < -5000)):
#         MI_BAL.append(Mi[i])
#         zem_BAL.append(zem[i])
#         errzem_BAL.append(errzem[i])
#         FIRST_BAL.append(FIRST_fnu20cm[i])
#         SNFIRST_BAL.append(SNFIRST_fnu20cm[i])
#         vmax_BAL.append(vmax_DR9_in9[i])
#         vmin_BAL.append(vmin_DR9_in9[i])
#         BI_BAL.append(AI_DR9_in9[i])
#         EW_BAL.append(EW_DR9_in9[i])
        
#         SNR9_BAL.append(SNR9[i])
        
#         qq,=np.where(np.array(QUASARNAME_EHVO) == columns[0])
#         if (len(qq) != 0):
#             countBALinEHVO = countBALinEHVO+1
#             EHVOBALS.append(columns[0])


# ----- plot zem ------------------------  4444444444444444

zem=np.hstack(zem)
#zem_BAL=np.hstack(zem_BAL)
#zem_EHVO=np.hstack(zem_EHVO)

#zem_BAL_toplot=[zem_BAL,zem_BAL] 
#zem_BAL_toplot=np.hstack(zem_BAL_toplot)

#zem_EHVO_toplot=[zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO]
#zem_EHVO_toplot=np.hstack(zem_EHVO_toplot)

fig=plt.figure(1) #### 4  

iqr = np.subtract(*np.percentile(zem, [75, 25]))
nhist=(max(zem)-min(zem))/(2*iqr*(len(zem)**(-1/3))) 

bins=np.linspace(min(zem),max(zem),int(nhist*0.75))#1.5))

# Original figure (Fig 12) without weighing, just multiplying by 2 and 10 the frequency of BALQSOs and EHVO
#plt.hist([zem,zem_BAL_toplot,zem_EHVO_toplot],bins,color=['black','blue','red'],label=['parent','BAL','EHVO'],histtype='step')

#plt.hist([zem,zem_all],bins,color=['black','blue'],label=['Parent w/o low SNR','Parent'], linewidth = 1, histtype='step')

plt.xlabel('$z_{em}$')
plt.ylabel('Number')
plt.legend(loc='upper right')
          
fig.savefig('zem_DR16.png', dpi=200) 

#KS_EHVOall_zem=ks_2samp(zem,zem_EHVO)
#KS_BALall_zem=ks_2samp(zem,zem_BAL)
#KS_EHVOBAL_zem=ks_2samp(zem_BAL,zem_EHVO)
#print('zem : p-value parent EHVO= '+str(KS_EHVOall_zem))
#print('zem: p-value BAL parent= '+str(KS_BALall_zem))
#print('zem: p-value BAL EHVO= '+str(KS_EHVOBAL_zem))


# ----- plot zem divided by total number in each class ----------  5555555

zem=np.hstack(zem)
zem_BAL=np.hstack(zem_BAL)
zem_EHVO=np.hstack(zem_EHVO)

zem_BAL_toplot=[zem_BAL,zem_BAL]
zem_BAL_toplot=np.hstack(zem_BAL_toplot)

zem_EHVO_toplot=[zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO]
zem_EHVO_toplot=np.hstack(zem_EHVO_toplot)

fig=plt.figure(2) 

iqr = np.subtract(*np.percentile(zem, [75, 25]))
nhist=(max(zem)-min(zem))/(2*iqr*(len(zem)**(-1/3)))

bins=np.linspace(min(zem),max(zem),int(nhist*0.75))  # before 1.5


# # If I want to weigh by the total number of cases in each category (Fig 13)
weights = [np.ones_like(zem)/float(len(zem)),np.ones_like(zem_BAL)/float(len(zem_BAL)),np.ones_like(zem_EHVO)/float(len(zem_EHVO))]
plt.hist([zem,zem_BAL,zem_EHVO],bins,color=['black','blue','red'],label=['parent','BAL','EHVO'],histtype='step',weights=weights)

plt.xlabel('$z_{em}$')
plt.ylabel('$n/N_{tot}$')
plt.legend(loc='upper right')
          
fig.savefig('zem_DR16_nNtot.png', dpi=200) 


#----- plot zem divided by parent sample in each bin ------ 6666666

zem=np.hstack(zem)
zem_EHVO=np.hstack(zem_EHVO)
zem_BAL=np.hstack(zem_BAL)

zem_BAL_toplot=[zem_BAL,zem_BAL]
zem_BAL_toplot=np.hstack(zem_BAL_toplot)

zem_EHVO_toplot=[zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO,zem_EHVO]
zem_EHVO_toplot=np.hstack(zem_EHVO_toplot)

iqr = np.subtract(*np.percentile(zem, [75, 25]))
nhist=(max(zem)-min(zem))/(2*iqr*(len(zem)**(-1/3)))

bins=np.linspace(min(zem),max(zem),int(nhist*0.75))  # before 1.5

# # If I want to weigh by total number of parent sample in that bin  (Fig 14)
nparent, binsparent, patchesparent = plt.hist([zem],bins,color=['black'],label=['parent'],histtype='step')
nBALS, binsBALS, patchesBALs = plt.hist([zem_BAL],bins,color=['blue'],label=['BAL'],histtype='step')
nEHVO, binsEHVO, patchesEHVO = plt.hist([zem_EHVO],bins,color=['red'],label=['EHVO'],histtype='step')

# ##### ??
# nNBALS=nBALS
# nNEHVO=nEHVO
# uu,=np.where(nparent != 0)
# nNBALS[uu]=nBALS[uu]/nparent[uu]
# nNEHVO[uu]=nEHVO[uu]/nparent[uu]
# uu0,=np.where(nparent == 0)
# nNBALS[uu0]=0
# nNEHVO[uu0]=0

# nNBALS=np.hstack(nNBALS)
# nNEHVO=np.hstack(nNEHVO)

# weights = [nNBALS,nNEHVO]
# weights = [np.ones_like(zem_BAL)/nparent,np.ones_like(zem_EHVO)/nparent]

fig=plt.figure(3)

weights = [np.ones_like(zem_BAL)/float(len(zem)),np.ones_like(zem_EHVO)/float(len(zem))]
plt.hist([zem_BAL,zem_EHVO],bins,color=['blue','red'],label=['BAL','EHVO'],histtype='step',weights=weights)

# plt.step(binsBALS[1:34],nNBALS[0:33],color='blue')
# plt.step(binsEHVO[1:34],nNEHVO[0:33],color='red')
# #plt.ylim(-0.1,0.4)

# plt.plot([binsBALS[0],binsBALS[0]],[0,nNBALS[0]], color='blue')
# plt.plot([binsBALS[0],binsBALS[1]],[nNBALS[0],nNBALS[0]], color='blue')
# plt.plot([binsEHVO[33],binsEHVO[33]],[0,nNEHVO[32]], color='red')
# plt.plot([binsBALS[0],binsBALS[1]],[nNEHVO[0],nNEHVO[0]], color='red')

# plt.xlabel('$z_{em}$')
# plt.ylabel('$n/N_{parent}$')
# #plt.text(2.0,0.35,'BALs',color='blue')
# plt.text(2.0,0.9,'BALs',color='blue')
# #plt.text(2.0,0.3,'EHVO',color='red')
# plt.text(2.0,0.8,'EHVO',color='red')
          
# fig.savefig('zem_DR9_nNparent2.png', dpi=200) 
fig.savefig('zem_DR16_weighted.png', dpi=200)