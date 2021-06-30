"""
======================
absorption_code2021.py
======================

@author Paola Rodriguez Hidalgo, Wendy Garcia Naranjo, Mikel Charles, Nathnael Kahassi, 
Michael Parker, (not sure who else contributed before us but be add as needed)

Creates figure to visually inspect the absorption. Calculates absorption parameters 
(BALnicity Index BI, vmin and vmax) for a list of spectra.

Notes
-----
Usage:
    python absorption_code.py spectra_data_list.csv (the default is ...).

Input file:
    This program takes a CSV file with the format `spectrum_name`, `z`, `snr`.

Parameters
----------
    Spectra base path (path to where spectra are stored on disk).
"""

#############################################################################################
########################################## IMPORTS ##########################################

import os
import sys
import numpy as np 
from matplotlib import pyplot as plt
from numpy.lib.function_base import append
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
from utility_functions import print_to_file, clear_file, append_row_to_csv
from data_types import Range, RangesData, FigureData, FigureDataOriginal, FlaggedSNRData 
from useful_wavelength_flux_error_modules import wavelength_flux_error_for_points, wavelength_flux_error_in_range, calculate_snr
from file_reader import read_file
from scipy import signal

# imports currently copied from normalization.py
# add and delete imports as needed and as we go through the code.

#############################################################################################
############################## CHANGEABLE VARIABLES #########################################

# INPUT WHICH DATA RELEASE YOU ARE WORKING WITH [INPUT NUMBER ONLY i.e. '9']
DR = '16'

# read list of spectra, zem, and snr 
# ^^ set config_file ^^
config_file = sys.argv[1] #Set cfg file path (csv with spec_name,z,snr) waiting to hear ...
                          # ... what column they will be in good_normalization.csv from mikel
# IDEA: use file_reader.py to create own csv files with just spec_name,z,snr, as opposed ...
#       ... to certain columns from good_normalization.csv

# SETS THE DIRECTORY TO FIND THE DATA FILES (DR9, DR16)
SPEC_DIREC = os.getcwd() + "/OUTPUT_FILES/NORMALIZATION/good_normalization.csv"

# CREATES DIRECTORY FOR OUTPUT FILES
OUT_DIREC = os.getcwd() + "/OUTPUT_FILES/ABSORPTION/"

smooth ='yes' # do you want to use smoothed norm flux/error instead of unsmoothed norm flux/error
boxcar_size = 5  # boxcar_size must always be an odd integer.

# Set a variable to plot all cases or only those with absorption -- Do you want to include all ...
# ... cases (if no, it only includes those with absorption). It is called plotall in the old ...
# ... code, but I am not 100% sure of what it excludes. 

countBI = 2000 # = lower limit of absorption width to be flagged 
maxvel = -30000.
minvel = -60000.

#############################################################################################
######################################## OUTPUT FILES #######################################

# set name of output txt file with absorption values
LOG_FILE = OUT_DIREC + "/" + "log.txt"
ABSORPTION_VALUES = OUT_DIREC + "/" + "absorption_measurements.txt"

# set name of output pdf with plots 
ABSORPTION_OUTPUT_PLOT = PdfPages('absorption_BI2000_test.pdf') 

#############################################################################################
####################################### DO NOT CHANGE #######################################

# verner table data
wavelength_CIV_emit1=1550.7700
wavelength_CIV_emit2=1548.1950
avr_CIV_doublet = 1549.0524 #weighted average
avr_SiIV_doublet = 1396.747 # weighted average; individuals: 1402.770, 1393.755
CII_emitted = 1335.313 # (weighted average); individuals:
OI_emitted = 1303.4951 # weighted average; individuals pag 20 in Verner Table
avr_NV_doublet = 1240.15 # weighted average; individuals: 1242.80, 1238.82
avr_OVI_doublet=1033.8160 # weighted average; individuals: 1037.6167, 1031. 9261

#############################################################################################
######################################### FUNCTION(S) #########################################

def smooth(norm_flux, box_size):   
    """Function: 

    Parameters:
    -----------
    norm_flux : 

    box_size: int
        This is the number of points that are smoothed into one. Always be sure to use an odd 
        number, because we need the same amount of points on each side of the data point to be
        smoothed.

    Returns:
    --------
    y_smooth

    Examples:
    ---------

    """    
    y_smooth = signal.savgol_filter(norm_flux,box_size,2)  #linear
    return y_smooth

#############################################################################################
####################################### VARIABLES ###########################################

brac_all, deltav_all=[]
absspeccount=0
count=0
BI=0
vmins, vmaxs, vmins_all, vmaxs_all =[] # v = velocity
final_depth_individual, final_depth_all_individual =[]
BI_all, BI_total, BI_ind_sum, BI_individual, BI_all_individual, BI_ind=[]
EW_individual, EW_ind, EW_all_individual, vlast =[] #EW = equivalent width

#############################################################################################
######################################### MAIN CODE #########################################

# Clear files

# Loops over each spectra
# use and modify below from absorption_indmeasurements_May2019.py to do this.
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
    
    # Read the wavelength, norm_flux and norm_error, rounding the numbers.   
  
    # Include if statement for smoothing and smooth spectrum if so. Similar to normalization.py.   

    # Transform the wavelength array to velocity (called "beta" - we can change it) based on the CIV doublet: 
    z_absC = (wavelength/avr_CIV_doublet)-1.
    RC=(1.+zem)/(1.+z_absC)
    betaC=((RC**2.)-1.)/((RC**2.)+1.)
    betaa = -betaC*(300000.)
    beta=[]
    for ll in betaa:
        betas=round (ll,4)
        beta.append (betas)
    beta=array(beta)

    # Initialize all the variables 
      
    # Calculate BI, vmin and vmax by looping through the beta array in the velocity limits -- 
    ###################### 0 -> -
    # will call the module that does that. It in

	# It uses a loop: for jjjs in jjj:

    # Calculate depth of each individual absorption trough.  -- a module too?     ***** YES module***
    tmp = norm_flux_used[vmaxs_index:vmins_index]
    tmp_index = np.where(norm_flux_used[vmaxs_index:vmins_index] == np.min(norm_flux_used[vmaxs_index:vmins_index]))
                   
    # The depth is calculated around the minimum, instead of as the minimum point, which could be a spike.
    final_depth = round (np.median(1.-tmp[int(tmp_index[0])-10 : int(tmp_index[0])+10]),2)
    #final_depth=round((1.-np.min(norm_flux_used[vmaxs_index:vmins_index])),2) <-- kept it to show how it was done before with the minimum point. You can delete it later
    final_depth_individual.append(final_depth)
    print('depth',final_depth_individual)

    final_depth_all_individual.append(final_depth_individual)

    # Calculate where CIV, CII and OI would be for each pair of vmin and vmax *if* the EHVO absorption found were 
    # instead not EHVO and due to SiIV: 

    # If the absorption is SiIV, this finds and plots where CIV, CII and OI would be
                    z_absSiIV = (wavelength[jjjs]/avr_SiIV_doublet)-1    #<-- right now it does it in the loop value, ...
                    # ... that in the BI program is called jjjs ( for jjjs in jjj:). We want to rewrite this so i can be a module. 

			
                    axvspan(beta[vmins_index],beta[vmaxs_index], alpha=0.2, color='red')
                    	
                    z_absSiIV_final = (wavelength[vmaxs_index]/avr_SiIV_doublet)-1.
                    
                    obs_wavelength_Cfinal=(z_absSiIV_final+1.)*(avr_CIV_doublet)
                    obs_wavelength_Cfinal_index =np.min (where (wavelength>obs_wavelength_Cfinal))
                    obs_wavelength_C_final_vel=beta[obs_wavelength_Cfinal_index]
                    axvspan(obs_wavelength_C_vel,obs_wavelength_C_final_vel, alpha=0.2, color='grey')
		    plot((obs_wavelength_C_vel, obs_wavelength_C_vel),(-1,10),'k-')#  <-- plot a lie at the vmin and vmax of the same color

                    obs_wavelength_CIIfinal=(z_absSiIV_final+1.)*(CII_emitted)
                    obs_wavelength_CIIfinal_index =np.min (where (wavelength>obs_wavelength_CIIfinal))
                    obs_wavelength_CII_final_vel=beta[obs_wavelength_CIIfinal_index]
                    axvspan(obs_wavelength_CII_vel,obs_wavelength_CII_final_vel, alpha=0.2, color='blue')
		    plot((obs_wavelength_CII_vel, obs_wavelength_CII_vel),(-1,10),'b-')


                    obs_wavelength_OIfinal=(z_absSiIV_final+1.)*(OI_emitted)
                    obs_wavelength_OIfinal_index =np.min (where (wavelength>obs_wavelength_OIfinal))
                    obs_wavelength_OI_final_vel=beta[obs_wavelength_OIfinal_index]
                    axvspan(obs_wavelength_OI_vel,obs_wavelength_OI_final_vel, alpha=0.2, color='yellow')
		    plot((obs_wavelength_OI_vel, obs_wavelength_OI_vel),(-1,10),'y-')

     
    # Plot figure as if the absorption was SiIV, CII or OI (we will visually inspect this file). 

    if (len(vmaxs) != 0) or (plotall == 'yes'):

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
        text(-60000, 2, str(i)+',     z='+str(zem)+' snr='+ str(snr), rotation=0, fontsize=9) # <-- Different location?

	# the axvspan using beta need to be defined earlier to plot here f.e., axvspan(beta[jjjs],beta[jjjs-1], alpha=0.05, color='red')



    # Save absorption info in file. 
   
    if (len(vmaxs) != 0) or (plotall == 'yes'): # < -- I am confused about the "or" there
        absspeccount=absspeccount+1  # <-- I think this is a variable to show how many cases have absorption? 
        yes=(str(count)+';'+str(absspeccount)+'  name: ' + str(i) + '\n' + 'BI (-30000 > v > -60,000): ' + str(BI_total) + '\n' +  'vmins: ' + str(vmins) + '\n' + 'vmaxs: '+str(vmaxs) + '\n' + 'BI_individual: '+ str(BI_individual) + '\n' + 'EW_individual: '+ str(EW_individual) + '\n' + 'Depth: '+ str(final_depth_individual) +'\n'+'\n')
        vlast.append(yes)
      
    
    
        pp.savefig()

        s = i.split('-')
        plateid = s[1]
        mjd = s[2]
        s = s[3].split('.')
        fiber = s[0]
        plt.savefig(output_spec + plateid + "-" + mjd + "-" + fiber + ".png") # <-- I don't think we really want a million png... only if doing a single one this makes sense. 

    close(count)
    
    if (len(vmaxs) != 0) or (plotall == 'yes'):  #<-- Again, I am confused about what the if includes plotall. 
        vmins_all.append(vmins)
        vmaxs_all.append(vmaxs)
          
# Clean below this. Do we need it? Where is it saving it?

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

