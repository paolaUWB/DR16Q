  # Initialize all variables for each spectrum (again, clean so it is not so many lines)
    vmins=[]
    vmaxs=[]
    BI_mid=[]    
    BI_individual=[]
    EW_individual=[]
    index_depth_final=[]
    flux_depth=[]
    
    final_depth_individual = []
    non_trough_count = 100
  
    deltav = 0 #change in velocity
    part = 0
    bb = -1
                       
    count2=0   # variable initialization to get into vmin/vmax loop
            
    # Set the limits of beta array based on minvel and maxvel.....................
    
    fst = 0

    if beta.any():
        try:
            fst = np.max(where(beta <= maxvel)) #index value of the starting point (on the very left) 
	#-- index value of minvel
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

# Loop through the beta array in those limits:

 for jjjs in jjj:
        
        # Initialize variables in each loop
        C = 0 
	# brac = (1. - (sm_flux[jjjs] / 0.9))  # [1 - f(v)/0.9] = brac > 0 when there is an absorption feature 
        brac = (1. - (norm_flux_used[jjjs] / 0.9))
        bracBAL= (1. - (norm_flux_used[jjjs] / 0.9))
	# print('brac'+str(brac))
        
        # Handle 3-point spike
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

  plot ((beta[vmins_index], beta[vmins_index]) , (-1,10),'r-')

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




