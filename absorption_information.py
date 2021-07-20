 
#Add keyword argument detailing the specific range over which BI is calculated. Outputs will need to be total BI_EHVO, BI per abs feature,
#vmax, and vmin per absorption feature, and depth of absorption feature
 
  # Initialize all variables for each spectrum (again, clean so it is not so many lines)
from absorption import BI_all


vmins, vmaxs=[]
BI_total, EW_individual=[] #EW is equivalent width = (vmax-vmin)
BI_individual=[] #BALnicity index will be single-valued per spectrum in terms of km/s
index_depth_final, flux_depth, final_depth_individual = []
non_trough_count = 100 ###what is this? Probably something just set as an arbitrary large number (like 99 or 999)

# verner table data
wavelength_CIV_emit1=  1550.7700
wavelength_CIV_emit2 = 1548.1950
avr_CIV_doublet = 1549.0524 #weighted average
avr_SiIV_doublet = 1396.747 # weighted average; individuals: 1402.770, 1393.755
CII_emitted = 1335.313 # (weighted average); individuals:
OI_emitted = 1303.4951 # weighted average; individuals pag 20 in Verner Table
avr_NV_doublet = 1240.15 # weighted average; individuals: 1242.80, 1238.82
avr_OVI_doublet = 1033.8160 # weighted average; individuals: 1037.6167, 1031. 9261


deltav = 0 #change in velocity
part = 0
bb = -1 #???
                    
count2=0   # variable initialization to get into vmin/vmax loop
        
# Set the limits of beta array based on minvel and maxvel.....................

first = 0

if beta.any(): #This "any" checks the truth value of any element, not what the actual element is
    try:
        first = np.max(np.where(beta <= maxvel)) #index value of the starting point (on the very left) 
#-- index value of minvel
    except:
        #first = np.max(where(beta == maxvel))
        first = 0
        
try:
    last = np.min(np.where(beta >= minvel)) #index value of the ending point (on the very right) -- index value of maxvel
except:
    last = np.where(beta == np.min(beta))

jjj = np.arange(last, first) #makes an array; inclusive on 'last', exclusive on 'first'
jjj = jjj[::-1]   #Reverses ordering from left to right; [a b c]->[c b a]

figure(count)

# Loop through the beta array in those limits:

for entry in jjj:
    # Initialize variables in each loop
    C = 0 
    """trough_cutoff has taken the place of brac"""
    # trough_cutoff = (1. - (sm_flux[entry] / 0.9))  # [1 - f(v)/0.9] = trough_cutoff > 0 when there is an absorption feature 
    trough_cutoff = (1. - (norm_flux_used[entry] / 0.9))
    bracBAL= (1. - (norm_flux_used[entry] / 0.9))   
    # print('trough_cutoff'+str(trough_cutoff))

# Handle 3-point spike
if trough_cutoff > 0:
    non_trough_count = 0
else:
    non_trough_count += 1
    trough_cutoff = 0

if((trough_cutoff > 0) or (non_trough_count <= 3)):
    
    deltav = beta[entry] - beta[entry - 1]
    part = part + deltav
    brac_all.append(trough_cutoff)
    deltav_all.append(deltav)
    
    EW = trough_cutoff * deltav
    EW = round(EW, 4)
    EW_ind.append(EW)          
    #print('EW', EW)
#***
    if part >= countBI:

        #print('I am in part > countBI')
        C = 1                
        BI = (trough_cutoff * C) * (deltav) #Calculate BAL for this dv
        BI = round(BI, 4)
        BI_total.append(BI) #Append to intermediate results

        BI_ind.append(BI)

        if non_trough_count == 0:
            plot((beta[entry + 1],beta[entry]),(1.5,1.5),'k-')
#                axvspan(beta[entry+1],beta[entry], alpha=0.05, color='red')
        
        if (count2==0) and (non_trough_count==0):  
            print('I am in vmin territory')

            vmins_index=np.min(where(beta >= (beta[entry]+countBI)))  # vmins occurs current beta plus countBI
            vvvmins=beta[vmins_index]
            vvvmins=round (vvvmins,4)
            vmins.append(vvvmins)
#                    print(vmins)

            plot ((beta[vmins_index], beta[vmins_index]) , (-1,10),'r-')

        # If the absorption is SiIV, this finds and plots where C, CII and OI would be [add proper if statement]
            z_absSiIV = (wavelength[entry]/avr_SiIV_doublet)-1#
            
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
            #This seems like the code was being tested in the short term to loop over each interval with "brac" (now trough cutoff)        
            nextbrac = (1. - (norm_flux_used[entry-1] / 0.9))
            nextnextbrac = (1. - (norm_flux_used[entry-2] / 0.9))
            nextnextnextbrac = (1. - (norm_flux_used[entry-3] / 0.9))
            nextnextnextnextbrac = (1. - (norm_flux_used[entry-4] / 0.9))
            
            if (((trough_cutoff>0 and nextbrac<0 and nextnextbrac<0 and nextnextnextbrac<0 and nextnextnextnextbrac<0 and count2==1)) or (entry == last)):  
            
                print("I am vmax territory!")
                vvmaxs = beta[entry]  
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
        BI_individual.append(BI_ind_sum)# this array contains one single BI value of each absorption feature in a single spectrum
                BI_ind = []
                
                EW_ind_sum = round(sum(EW_ind),2)
                EW_individual.append(EW_ind_sum)
                EW_ind = []
                                
                final_depth=round((1.-np.min(norm_flux_used[vmaxs_index:vmins_index])),2)
                final_depth_individual.append(final_depth)
                print('depth',final_depth_individual)
                
                count2=0                     
                                    
    else: #if the trough_cutoff value is not more than zero (so if we don't have absorption feature)
        part=0 # this is so b/c we do not want to keep counting the width of the absorption feature if it is not wider than 600km/s
        count2=0# this is so b/c if the code encounters an other absorption feature which is wider than 600km/s, the code is going to go through the if statement on line 205
        EW_ind=[]

    if entry == last:   
        BI_total= round(sum(BI_total),2)         
        BI_all.append(BI_total)    
        BI_all_individual.append(BI_individual)
        EW_all_individual.append(EW_individual)

"""Per spectrum, we can just take the individual BI for every trough and sum them for total BI;
I'm not sure what differentiates the BI_total and BI_all so I made BI_sum """
#This might be replacing the above work^
BI_sum = np.sum(BI_individual)

###################################### Working through #####################################
"""
test_beta = [1,2,3,4]

test_vel = np.linspace(-10,10,100)

test_vmins_index = [0,5,10,14]
test_vmaxs_index = [4,9,13,20]
b= np.linspace(0,4,5)
print(b)

#I need to generate a linspace of the indices so that I can find all the fluxes and get the minimum one.
for i in vmins_index:
It's actaully the elements of vmins_index and vmaxs_index that 
    for j in vmaxs_index:
        a=np.linspace(vmins_index[i],vmaxs_index[j],vmaxs_index[j]-vmins_index[i]+1)
        print(a)
#a= range(vmins_index[3])
"""
########################## Working through the above in spyder so it actually runs ###########################


def depth(vmins_index,vmaxs_index):
"""Returns the depth of each BAL from the given spectra. Must be evaluated on each spectrum.

Parameters
----------
vmins_index : float
    The starting point for the BAL trough
vmaxs_index : float
    The ending border of the BAL trough

Returns
-------
float
    "depth" will return the lowest flux for the region bordered by vmin and vmax. Be aware that broken pixels in the region may
    lead to erroneous depths.

"""
    #so the indices across each BAL need to be looped through:
    depths=[]
    for i in range(len(vmins_index)):
        BALregion=np.linspace(vmins_index[i],vmaxs_index[i],vmaxs_index[i]-vmins_index[i]+1)
        #print(BALregion)
        minimum = np.min(norm_flux_used(BALregion))
        #print(minimum)
        depths.append(minimum)
    return depths
