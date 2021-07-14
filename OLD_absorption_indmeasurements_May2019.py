#############################################################
#   absorption_code.py
#   Calculates absorption(BALnicity) for a list of spectra
#
#   Usage:
#      python absorption_code.py spectra_data_list.csv
#
#   Input file:
#        This program takes a CSV file with the format...
#
#        spectrum_name,z,snr
#
#   Input parameters:
#       Spectra base path(path to where spectra are stored on disk)
#          this needs to be changes when run on another computer:


# REWRITE! This code is a mess

#############################################################

import sys
import csv

from pylab import*
import numpy as np
from sympy import sympify
from scipy import signal,misc
from astropy import*
from matplotlib.pyplot import*
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages


# To run the code:

# %run absorption_indmeasurements_May2019.py temp_abs.csv

# where temp_abs.csv only has one spectrum

#__________________________________________________________
# Inputs & variables to change


# Location of normalized files:
spec_path = '/Users/Paola/QUASAR/Work_EHVO/ROUTINES/NORM/'

output_spec= '/Users/Paola/QUASAR/Work_EHVO/FIGURES/ABS_PLOTS/spec_'

# Name of output pdf with plots:
#pp=PdfPages('absorption_BI2000_test.pdf')

#pp=PdfPages('absorptiononly_BI1000_EHVOcases22Jan2019.pdf')
pp=PdfPages('absorptiononly_BI1000_EHVOcasescleaning_Dec2019.pdf')

#ffile='Absorption_EHVOcases22Jan2019.txt'
ffile='Absorption_cleaning_Dec2019.txt'

# Do you want to include all cases (if no, it only includes those with absorption)
plotall='yes'

countBI = 1000 # = lower limit of absorption width to be flagged 
  
maxvel = -30000.
minvel = -60000.

# Do you want to use smoothed norm flux/error instead of unsmoothed norm flux/error
sm='yes'
n=5  # Smooth boxcar size

# Necessary data from Verner table
wavelength_CIV_emit1=1550.7700
wavelength_CIV_emit2=1548.1950
avr_CIV_doublet = 1549.0524 #weighted average
avr_SiIV_doublet = 1396.747 # weighted average; individuals: 1402.770, 1393.755
CII_emitted = 1335.313 # (weighted average); individuals:
OI_emitted = 1303.4951 # weighted average; individuals pag 20 in Verner Table
avr_NV_doublet = 1240.15 # weighted average; individuals: 1242.80, 1238.82
avr_OVI_doublet=1033.8160 # weighted average; individuals: 1037.6167, 1031. 9261

#_______________________________________________________

fig=figure()

#...................................................
# Read list of spectra, zem, and snr 
config_file = sys.argv[1] #Set cfg file path (csv with spec_name,z,snr...)

spectra_list = list()
redshifts_list =  list() 
snr_list =  list() 

for line in open(config_file, 'r'):
	d = line.split(",")
	spectra_list.append(d[0])
	redshifts_list.append(np.float(d[1]))
	snr_list.append(np.float(d[2]))
    
    
#..................................................

# Variables .......................................

# XXX Check the variables are all necessary
brac_all=[]
deltav_all=[]

absspeccount=0
count=0
BI=0
vmins=[]
vmaxs=[]
vmins_all=[]
vmaxs_all=[]

final_depth_individual=[]
final_depth_all_individual=[]

BI_all=[]
BI_total=[]   # Total BI per spectrum
BI_ind_sum=[] # Each 
BI_individual=[]
BI_all_individual=[]
BI_ind=[]

EW_individual=[]
EW_ind=[]
EW_all_individual=[]
vlast=[]

#...................................................

def smooth(norm_flux, box_size):       
#    box = np.ones(box_size)/box_size
#    y_smooth = convolve(norm_flux, box, mode='same')
    y_smooth=signal.savgol_filter(norm_flux,box_size,2)  #linear
    return y_smooth



# Loops over the spectra
for i,j,k in zip (spectra_list, redshifts_list, snr_list):   
    
    count=count+1
    print(count)
    print(i)
    
    normalized_dr9 = loadtxt(spec_path + i) #Load in normalized spectrum
    wavelength = normalized_dr9[:,0] 
    norm_flux = normalized_dr9[:,1] 
    norm_error = normalized_dr9[:,2] 
    
   
    
    
    sm_flux=smooth(norm_flux,n) #Smooth the spectrum (3 point boxcar)
    sm_error=smooth(norm_error,n)/sqrt(n)
    
    norm_flux_used = norm_flux
    norm_error_used = norm_error
    
    if sm =='yes':
        norm_flux_used = sm_flux
        norm_error_used = sm_error
        
    # Error calculation -- Remove! XXXX
    #norm_flux_used=norm_flux_used-0.05
    
    # Masked values with another absorption present: 
    
    norm_flux_used[1272:1288] = 0.75
#    norm_flux_used[2673:2682] = 0.8
#    norm_flux_used[2686:2693] = 0.81
    
    
    zem=round (j,5) #Round the redshift
    
    snr = round(k, 5)

    #..................................
    # Initialize all variables for each spectrum
    vmins=[]
    vmaxs=[]
    BI_mid=[]    
    BI_individual=[]
    EW_individual=[]
    beta=[]
    index_depth_final=[]
    flux_depth=[]
    
    final_depth_individual = []
    non_trough_count = 100
  
    deltav = 0 #change in velocity
    part = 0
    bb = -1
                       
    count2=0   # variable initialization to get into vmin/vmax loop
    
    #..................................
    # Make beta array of velocities
    z_absC = (wavelength/avr_CIV_doublet)-1.
    RC=(1.+zem)/(1.+z_absC)
    betaC=((RC**2.)-1.)/((RC**2.)+1.)
    betaa = -betaC*(300000.)
    for ll in betaa:
        betas=round (ll,4)
        beta.append (betas)
    beta=array(beta)
        
    # Set the limits of beta based on minvel and maxvel.....................
    
    fst = 0

    if beta.any():
        try:
            fst = np.max(where(beta <= maxvel)) #index value of the starting point (on the very left) -- index value of minvel
        except:
            #fst = np.max(where(beta == maxvel))
            fst = 0
            
    try:
        lst = np.min(where(beta >= minvel)) #index value of the ending point (on the very right) -- index value of maxvel
    except:
        lst = where(beta == np.min(beta))
    
    jjj = arange(lst, fst)
    jjj = array(jjj)
    jjj = jjj[::-1]   # From right to left

    figure(count)
    
    for jjjs in jjj:
        
        # Initialize variables in each loop
        C = 0 
#        brac = (1. - (sm_flux[jjjs] / 0.9))  # [1 - f(v)/0.9] = brac > 0 when there is an absorption feature 
        brac = (1. - (norm_flux_used[jjjs] / 0.9))
        bracBAL= (1. - (norm_flux_used[jjjs] / 0.9))
#        print('brac'+str(brac))
        
        # Handle 3-point
        if brac > 0:
            non_trough_count = 0
        else:
            non_trough_count += 1
            brac = 0

        if((brac > 0) or (non_trough_count <= 3)):
            
            deltav = beta[jjjs] - beta[jjjs - 1]
            part = part + deltav
            brac_all.append(brac)
            deltav_all.append(deltav)
            
            EW = brac * deltav
            EW = round(EW, 4)
            EW_ind.append(EW)          
#           print('EW',EW)
        
            if part >= countBI:

                #print('I am in part > countBI')
                C = 1                
                BI = (brac * C) * (deltav) #Calculate BAL for this dv
                BI = round(BI, 4)
                BI_mid.append(BI) #Append to intermediate results

                BI_ind.append(BI)

                if non_trough_count == 0:
                    plot((beta[jjjs + 1],beta[jjjs]),(1.5,1.5),'k-')
#                axvspan(beta[jjjs+1],beta[jjjs], alpha=0.05, color='red')
                
                if (count2==0) and (non_trough_count==0):  
                    print('I am in vmin territory')

                    vmins_index=np.min(where(beta >= (beta[jjjs]+countBI)))  # vmins occurs current beta plus countBI
                    vvvmins=beta[vmins_index]
                    vvvmins=round (vvvmins,4)
                    vmins.append(vvvmins)
#                    print(vmins)
                    
                    plot ((beta[vmins_index], beta[vmins_index]), (-1,10),'r-')

                    # If the absorption is SiIV, this finds and plots where C, CII and OI would be
                    z_absSiIV = (wavelength[jjjs]/avr_SiIV_doublet)-1#
                    
                    obs_wavelength_C=(z_absSiIV+1)*(avr_CIV_doublet)#
                    obs_wavelength_C_index =np.min (where (wavelength>obs_wavelength_C))
                    obs_wavelength_C_vel=beta[obs_wavelength_C_index]+countBI
                    plot((obs_wavelength_C_vel, obs_wavelength_C_vel),(-1,10),'k-')

                    obs_wavelength_CII=(z_absSiIV+1)*(CII_emitted)#
                    obs_wavelength_CII_index =np.min (where (wavelength>obs_wavelength_CII))                  
                    obs_wavelength_CII_vel=beta[obs_wavelength_CII_index]+countBI
                    plot((obs_wavelength_CII_vel, obs_wavelength_CII_vel),(-1,10),'b-')

                    obs_wavelength_OI=(z_absSiIV+1)*(OI_emitted)#
                    obs_wavelength_OI_index =np.min (where (wavelength>obs_wavelength_OI))                  
                    obs_wavelength_OI_vel=beta[obs_wavelength_OI_index]+countBI
                    plot((obs_wavelength_OI_vel, obs_wavelength_OI_vel),(-1,10),'y-')

                    count2=1
                    
                nextbrac = (1. - (norm_flux_used[jjjs-1] / 0.9))
                nextnextbrac = (1. - (norm_flux_used[jjjs-2] / 0.9))
                nextnextnextbrac = (1. - (norm_flux_used[jjjs-3] / 0.9))
                nextnextnextnextbrac = (1. - (norm_flux_used[jjjs-4] / 0.9))
                
                if (((brac>0 and nextbrac<0 and nextnextbrac<0 and nextnextnextbrac<0 and nextnextnextnextbrac<0 and count2==1)) or (jjjs == lst)):  
                
                    print("I am vmax territory!")
                    vvmaxs = beta[jjjs]  
                    vmaxs_index = np.min (where (beta>= vvmaxs))
                    vvmaxs = round(vvmaxs,4)
                    vmaxs.append(vvmaxs)
#                    print(vvmaxs)
                    
                    axvspan(beta[vmins_index],beta[vmaxs_index], alpha=0.2, color='red')
                    print('vmins=',beta[vmins_index])
                    print('vmaxs=',beta[vmaxs_index])
                    z_absSiIV_final = (wavelength[vmaxs_index]/avr_SiIV_doublet)-1.
                    
                    obs_wavelength_Cfinal=(z_absSiIV_final+1.)*(avr_CIV_doublet)
                    obs_wavelength_Cfinal_index =np.min (where (wavelength>obs_wavelength_Cfinal))
                    obs_wavelength_C_final_vel=beta[obs_wavelength_Cfinal_index]
                    axvspan(obs_wavelength_C_vel,obs_wavelength_C_final_vel, alpha=0.2, color='grey')

                    obs_wavelength_CIIfinal=(z_absSiIV_final+1.)*(CII_emitted)
                    obs_wavelength_CIIfinal_index =np.min (where (wavelength>obs_wavelength_CIIfinal))
                    obs_wavelength_CII_final_vel=beta[obs_wavelength_CIIfinal_index]
                    axvspan(obs_wavelength_CII_vel,obs_wavelength_CII_final_vel, alpha=0.2, color='blue')

                    obs_wavelength_OIfinal=(z_absSiIV_final+1.)*(OI_emitted)
                    obs_wavelength_OIfinal_index =np.min (where (wavelength>obs_wavelength_OIfinal))
                    obs_wavelength_OI_final_vel=beta[obs_wavelength_OIfinal_index]
                    axvspan(obs_wavelength_OI_vel,obs_wavelength_OI_final_vel, alpha=0.2, color='yellow')

                    BI_ind_sum = round(sum(BI_ind),2)
                    BI_individual.append(BI_ind_sum)# this array contains one single BI value of each absortopn feature in a single spectrum
                    BI_ind = []
                    
                    EW_ind_sum = round(sum(EW_ind),2)
                    EW_individual.append(EW_ind_sum)
                    EW_ind = []
                                    
                    final_depth=round((1.-np.min(norm_flux_used[vmaxs_index:vmins_index])),2)
                    final_depth_individual.append(final_depth)
                    print('depth',final_depth_individual)
                    
                    count2=0                     
                                        
        else: #if the brac value is not more than zero (so if we don't have absorption feature)
            part=0 # this is so b/c we do not want to keep counting the width of the absorption feature if it is not wider than 600km/s
            count2=0# this is so b/c if the code encounters an other absorption feature which is wider than 600km/s, the code is going to go through the if statement on line 205
            EW_ind=[]
        
        if jjjs == lst:
            BI_total= round(sum(BI_mid),2)         
            BI_all.append(BI_total)    
            BI_all_individual.append(BI_individual)
            EW_all_individual.append(EW_individual)
   
    if (len(vmaxs) != 0) or (plotall == 'yes'): 
        absspeccount=absspeccount+1
        yes=(str(count)+';'+str(absspeccount)+'  name: ' + str(i) + '\n' + 'BI (-30000 > v > -60,000): ' + str(BI_total) + '\n' +  'vmins: ' + str(vmins) + '\n' + 'vmaxs: '+str(vmaxs) + '\n' + 'BI_individual: '+ str(BI_individual) + '\n' + 'EW_individual: '+ str(EW_individual) + '\n' + 'Depth: '+ str(final_depth_individual) +'\n'+'\n')
        vlast.append(yes)

    final_depth_all_individual.append(final_depth_individual)
 
    if (len(vmaxs) != 0) or (plotall == 'yes'):
           
        xlim(np.min(beta),0)# this is just seting how wide the graph should be (so we are setting the domain)
        title('Normalized flux vs velocity')
        xlabel('Velocity (km/s)')
        ylabel('Normalized flux')
        plot((np.min(beta),np.max(beta)),(1,1))
        plot((np.min(beta),np.max(beta)),(0.9,0.9),'r--')
        plot(beta, norm_flux_used, 'k-')
#        plot(beta, sm_flux, 'r-')
#        plot(beta, norm_flux, 'k-')
        plot (beta, norm_error_used,'k--')
        #plot (beta, norm_error,'k--')
        #plot (beta, sm_error,'k--')
        ylim(0,3)
        xlim(-70000,0)
        text(-60000, 2, str(i)+',     z='+str(zem)+' snr='+ str(snr), rotation=0, fontsize=9)
    
# FUTURE PLOTTING!    
#    if (len(vmins)) ne 0:
#        axvspan(beta[jjjs],beta[jjjs-1], alpha=0.05, color='red')
    
        pp.savefig()

        s = i.split('-')
        plateid = s[1]
        mjd = s[2]
        s = s[3].split('.')
        fiber = s[0]
        plt.savefig(output_spec + plateid + "-" + mjd + "-" + fiber + ".png")

    close(count)
    
    if (len(vmaxs) != 0) or (plotall == 'yes'):
        vmins_all.append(vmins)
        vmaxs_all.append(vmaxs)
          
BI_all= array(BI_all)

vmins = array(vmins)
vmaxs = array(vmaxs)
pp.close()
vmaxs_final=[]
vmins_final=[]

for loop in range (0, len (vmaxs_all)):
    vmaxs_final.append (str(vmaxs_all[loop])+ ',' )

for loop2 in range (0, len(vmins_all)):
    vmins_final.append (str(vmins_all[loop2])+ ',' )
                    

vmaxs_final = array(vmaxs_final)
vmins_final = array(vmins_final)    
savetxt(ffile,vlast,fmt='%s')

