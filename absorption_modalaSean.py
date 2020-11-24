#absorption_code.py
#Calculates absorption(BALnicity) for a list of spectra
import numpy as np 
from matplotlib import pyplot as plt
from scipy import signal
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys
from data_types import Range

VELOCITY_LIMIT = Range(-60000., -30000)

BI = 0 #change it to BI
BALNICITY_INDEX_LIMIT = 2000 #lower limit of absorption width

# Necessary data from Verner table
WAVELENGTH_CIV_EMIT_LIMIT = Range(1548.1950, 1550.7700) #never used?
AVERAGE_CIV_DOUBLET = 1549.0524 #weighted average
AVERAGE_SiIV_DOUBLET = 1396.747 # weighted average; individuals: 1402.770, 1393.755
AVERAGE_NV_DOUBLET = 1240.15 # weighted average; individuals: 1242.80, 1238.82
AVERAGE_OVI_DOUBLET=1033.8160 # weighted average; individuals: 1037.6167, 1031. 9261
CII_EMITTED = 1335.313 # (weighted average); individuals:
OI_EMITTED = 1303.4951 # weighted average; individuals pag 20 in Verner Table

non_trough_count = 100
boxcar_size = 5  
deltav = 0 #change in velocity
part = 0
count2 = 0   # variable initialization to get into vmin/vmax loop

SPECDIREC, OUTPUT_SPEC = os.getcwd() + "/files/", os.getcwd() + "/files/"
NORM_FILE_EXTENSION = "norm.dr9"
ABSORPTION_PDF = PdfPages('absorptiononly_BI1000_EHVOcasescleaning.pdf')
OUTPUT_CLEANING = 'Absorption_cleaning.txt'
NUM_OF_FILES = 10

# Do you want to include all cases (if no, it only includes those with absorption)
# then set plotall to 'yes'
plotall = 'yes'

# Do you want to use smoothed norm flux/error instead of unsmoothed norm flux/error,
# then set sm to 'yes'
sm = 'yes'

brac_all = []
deltav_all = []
vmins, vmaxs = [], []
vmins_all, vmaxs_all = [], []
final_depth_individual, final_depth_all_individual = [], []

#BI = balnicity index
#EW = Equivalent widths
BI_all , BI_total = [], []
BI_ind_sum = []
BI_individual = []
BI_all_individual = []
BI_ind = []
BI_mid = []  
EW_individual = []
EW_ind = []
EW_all_individual = []
vlast = []

def smooth(norm_flux, box_size):       
    y_smooth = signal.savgol_filter(norm_flux, box_size, 2)  #linear
    return y_smooth

fig= plt.figure()

CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else "sorted_norm.csv"

spectra_list, redshift_value_list, snr_value_list = [], [], []
for line in open(CONFIG_FILE, 'r'):
	each_row_in_file = line.split(",")
	spectra_list.append(each_row_in_file[0])
	redshift_value_list.append(np.float(each_row_in_file[1]))
	snr_value_list.append(np.float(each_row_in_file[2]))

file_count = 0
for spectrum_file_name, redshift_value, snr_value in zip(spectra_list, redshift_value_list, snr_value_list):   
    
    #########################################################
    ##If you dont have all files, point out stopping variable
    ##Change NUM_OF_FILES for different value
    if file_count == NUM_OF_FILES:
        break

    file_count += 1
    print(file_count, ":", spectrum_file_name)
      
    z = round (redshift_value, 5)
    snr = round(snr_value, 5)

    #read the norm file
    normalized_dr9 = np.loadtxt(SPECDIREC + spectrum_file_name[0:20] + NORM_FILE_EXTENSION) 
    wavelength = normalized_dr9[:, 0] 
    normalized_flux = normalized_dr9[:, 1] 
    error_normalized = normalized_dr9[:, 2]
    
    sm_flux = smooth(normalized_flux, boxcar_size) #Smooth the spectrum (3 point boxcar)
    sm_error = smooth(error_normalized, boxcar_size) / np.sqrt(boxcar_size)
    
    if sm == 'yes':
        normalized_flux = sm_flux
        norm_error_used = sm_error
        
    beta_list = []    # Make beta array of velocities
    z_absC = (wavelength/AVERAGE_CIV_DOUBLET) - 1.
    RC = (1. + z) / (1. + z_absC)
    betaC = ((RC**2.) - 1.) / ((RC**2.) + 1.)
    betaa = -betaC * (300000.)
    
    for values in betaa:
        betas = round (values,4)
        beta_list.append(betas)
        
    beta_list = np.array(beta_list)

    # Set the limits of beta based on minvel and maxvel.....................
    if beta_list.any():
        try:
            max_limit_value = np.max(np.where(beta_list <= VELOCITY_LIMIT.end)) #index value of the starting point (on the very left) -- index value of minvel
        except:
            max_limit_value = 0
            
    try:
        min_limit_value = np.min(np.where(beta_list >= VELOCITY_LIMIT.start)) #index value of the ending point (on the very right) -- index value of maxvel
    except:
        min_limit_value = np.where(beta_list == np.min(beta_list))
        

    beta_list_minmax_range = np.arange(min_limit_value, max_limit_value)
    beta_list_minmax_range = np.array(beta_list_minmax_range[::-1])  # From right to left [reversed list]

    plt.figure(file_count)

    for index in beta_list_minmax_range:     #change it to index   
        #BI formula = [1- (f(v)/0.9)]*C
        C = 0 #set the zero
        brac = (1. - (normalized_flux[index] / 0.9)) #Dividing f(v) by 0.9 avoids detections of shallow absorption;
        
        # Handle 3-point
        if brac > 0:
            non_trough_count = 0
        else:
            non_trough_count += 1
            brac = 0


        if brac > 0 or non_trough_count <= 3:
            
            deltav = beta_list[index] - beta_list[index - 1]
            part += deltav
            brac_all.append(brac)
            deltav_all.append(deltav)
            
            EW = brac * deltav
            EW = round(EW, 4)
            EW_ind.append(EW)          

            if part >= BALNICITY_INDEX_LIMIT:

                C = 1  #set to 1 only  square bracket is continuously positive over a velocity interval            
                BI = (brac * C) * (deltav) #Calculate BAL for this dv
                BI_mid.append(round(BI, 4)) #Append to intermediate results
                BI_ind.append(round(BI, 4))

                if non_trough_count == 0:
                    
                    plt.plot((beta_list[index + 1],beta_list[index]),(1.5,1.5),'k-')
                    #axvspan(beta[jjjs+1],beta[jjjs], alpha=0.05, color='red')

#vmin territory calculation and plotting                
                if count2 == 0 and non_trough_count == 0:  
                    print('I am in vmin territory')

                    vmins_index = np.min(np.where(beta_list >= (beta_list[index] + BALNICITY_INDEX_LIMIT)))  # vmins occurs current beta plus countBI
                    vmins.append(round (beta_list[vmins_index], 4))
                    
                    plt.plot((beta_list[vmins_index], beta_list[vmins_index]), (-1,10),'r-')

                    # If the absorption is SiIV, this finds and plots where C, CII and OI would be
                    z_absSiIV = (wavelength[index]/AVERAGE_SiIV_DOUBLET) - 1#
                    
                    observed_wavelength_C = (z_absSiIV + 1) * (AVERAGE_CIV_DOUBLET)#
                    obs_wavelength_C_index = np.min (np.where (wavelength > observed_wavelength_C))
                    obs_wavelength_C_vel = beta_list[obs_wavelength_C_index] + BALNICITY_INDEX_LIMIT
                    plt.plot((obs_wavelength_C_vel, obs_wavelength_C_vel), (-1,10),'k-')

                    obs_wavelength_CII = (z_absSiIV+1) * (CII_EMITTED)#
                    obs_wavelength_CII_index = np.min (np.where (wavelength>obs_wavelength_CII))                  
                    obs_wavelength_CII_vel = beta_list[obs_wavelength_CII_index] + BALNICITY_INDEX_LIMIT
                    plt.plot((obs_wavelength_CII_vel, obs_wavelength_CII_vel), (-1,10),'b-')

                    obs_wavelength_OI = (z_absSiIV+1) * (OI_EMITTED)#
                    obs_wavelength_OI_index = np.min (np.where (wavelength>obs_wavelength_OI))                  
                    obs_wavelength_OI_vel = beta_list[obs_wavelength_OI_index] + BALNICITY_INDEX_LIMIT
                    plt.plot((obs_wavelength_OI_vel, obs_wavelength_OI_vel), (-1,10), 'y-')

                    count2 = 1

                    
                nextbrac = (1. - (normalized_flux[index-1] / 0.9))
                nextnextbrac = (1. - (normalized_flux[index-2] / 0.9))
                nextnextnextbrac = (1. - (normalized_flux[index-3] / 0.9))
                nextnextnextnextbrac = (1. - (normalized_flux[index-4] / 0.9))

 #vmax territory calculation and plotting                               
                if (((brac > 0 and nextbrac < 0 and nextnextbrac < 0 and nextnextnextbrac < 0 and nextnextnextnextbrac < 0 and count2 == 1)) or (index == min_limit_value)):  
                
                    print("I am in vmax territory!")

                    vmaxs_index = np.min(np.where (beta_list >= beta_list[index]))
                    vmaxs.append(round (beta_list[index], 4))
                                       
                    plt.axvspan(beta_list[vmins_index], beta_list[vmaxs_index], alpha = 0.2, color = 'red')
                    print('vmins=', beta_list[vmins_index])
                    print('vmaxs=', beta_list[vmaxs_index])
                    z_absSiIV_final = (wavelength[vmaxs_index]/AVERAGE_SiIV_DOUBLET)-1.
                    
                    obs_wavelength_Cfinal = (z_absSiIV_final+1.) * (AVERAGE_CIV_DOUBLET)
                    obs_wavelength_Cfinal_index = np.min(np.where(wavelength>obs_wavelength_Cfinal))
                    obs_wavelength_C_final_vel = beta_list[obs_wavelength_Cfinal_index]
                    plt.axvspan(obs_wavelength_C_vel,obs_wavelength_C_final_vel, alpha = 0.2, color = 'grey')

                    obs_wavelength_CIIfinal = (z_absSiIV_final + 1.) * (CII_EMITTED)
                    obs_wavelength_CIIfinal_index = np.min (np.where (wavelength>obs_wavelength_CIIfinal))
                    obs_wavelength_CII_final_vel = beta_list[obs_wavelength_CIIfinal_index]
                    plt.axvspan(obs_wavelength_CII_vel,obs_wavelength_CII_final_vel, alpha = 0.2, color = 'blue')

                    obs_wavelength_OIfinal = (z_absSiIV_final + 1.) * (OI_EMITTED)
                    obs_wavelength_OIfinal_index = np.min (np.where (wavelength>obs_wavelength_OIfinal))
                    obs_wavelength_OI_final_vel = beta_list[obs_wavelength_OIfinal_index]
                    plt.axvspan(obs_wavelength_OI_vel,obs_wavelength_OI_final_vel, alpha = 0.2, color = 'yellow')

                    BI_ind_sum = round(sum(BI_ind), 2)
                    BI_individual.append(BI_ind_sum) # this array contains one single BI value of each absortopn feature in a single spectrum
                    BI_ind = []
                    
                    EW_ind_sum = round(sum(EW_ind), 2)
                    EW_individual.append(EW_ind_sum)
                    EW_ind = []
                                    
                    final_depth = round((1. - np.min(normalized_flux[vmaxs_index:vmins_index])), 2)
                    final_depth_individual.append(final_depth)
                    print('depth', final_depth_individual)
                    
                    count2 = 0                     
                                        
        else: #if the brac value is not more than zero (so if we don't have absorption feature)
            part = 0 # this is so b/c we do not want to keep counting the width of the absorption feature if it is not wider than 600km/s
            count2 = 0 # this is so b/c if the code encounters an other absorption feature which is wider than 600km/s, the code is going to go through the if statement on line 205
            EW_ind = []
        
        if index == min_limit_value:
            BI_total = round(sum(BI_mid),2)         
            BI_all.append(BI_total)    
            BI_all_individual.append(BI_individual)
            EW_all_individual.append(EW_individual)
   
    if (len(vmaxs) != 0) or (plotall == 'yes'):
        text = [f"{file_count};{spectrum_file_name}",
                f"BI ({VELOCITY_LIMIT.end} > v > {VELOCITY_LIMIT.start}): {BI_total}",
                f"vmins: {vmins}",
                f"vmaxs: {vmaxs}",
                f"BI_individual: {BI_individual}",
                f"EW_individual: {EW_individual}",
                f"Depth: {final_depth_individual}"]
        vlast.extend(['\n'.join(text), '\n'])

    final_depth_all_individual.append(final_depth_individual)
 
    if (len(vmaxs) != 0) or (plotall == 'yes'):
           
        plt.xlim(np.min(beta_list), 0) # this is just seting how wide the graph should be (so we are setting the domain)
        plt.title('Normalized Flux vs Velocity')
        plt.xlabel('Velocity (km/s)')
        plt.ylabel('Normalized Flux')
        plt.plot((np.min(beta_list), np.max(beta_list)), (1,1))
        plt.plot((np.min(beta_list), np.max(beta_list)), (0.9,0.9), 'r--')
        plt.plot(beta_list, normalized_flux, 'k-')
#       plot(beta, sm_flux, 'r-')
#       plot(beta, norm_flux, 'k-')
        plt.plot (beta_list, norm_error_used,'k--')
        #plot (beta, norm_error,'k--')
        #plot (beta, sm_error,'k--')
        plt.ylim(0, 3)
        plt.xlim(-70000, 0)
        plt.text(-60000, 2, str(spectrum_file_name)+',     z='+str(z)+' snr='+ str(snr), rotation = 0, fontsize = 9)
    

# FUTURE PLOTTING!    
#    if (len(vmins)) ne 0:
#        axvspan(beta[jjjs],beta[jjjs-1], alpha=0.05, color='red')
    
        ABSORPTION_PDF.savefig()
        s = spectrum_file_name.split('-')
        plateid = s[1]
        mjd = s[2]
        s = s[3].split('.')
        fiber = s[0]
        plt.savefig(OUTPUT_SPEC + plateid + "-" + mjd + "-" + fiber + ".png") #this saves a png file

    plt.close(file_count)
    
    if (len(vmaxs) != 0) or (plotall == 'yes'):
        vmins_all.append(vmins)
        vmaxs_all.append(vmaxs)
          
BI_all= np.array(BI_all)

vmins = np.array(vmins)
vmaxs = np.array(vmaxs)
ABSORPTION_PDF.close()
vmins_final, vmaxs_final = [], []

for loop in range(0, len(vmaxs_all)):
    vmaxs_final.append(str(vmaxs_all[loop]) + ',')

for loop2 in range(0, len(vmins_all)):
    vmins_final.append(str(vmins_all[loop2]) + ',')

vmaxs_final = np.array(vmaxs_final)
vmins_final = np.array(vmins_final)    
np.savetxt(OUTPUT_CLEANING,vlast,fmt='%s')