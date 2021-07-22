from matplotlib import pyplot as plt
import numpy as np

def draw_abs_figure(velocity, flux_normalized, error, savefile_name, spectra_name, redshift, snr, c0, c2, o1, vminz, vmaxz):
    """ Draws the normalized spectra graph.
    
    Parameters
    ----------
    flux_normalized: array
        The normalized flux to be graphed.
    velocity: array
        The value of the velocity calculated using the normalized flux.
    Returns
    -------
    None.
    
    Notes
    -----
    Creates a graph of the spectra and saves to the ``absorption_BI2000_test.pdf``
    """
    plt.plot(velocity, flux_normalized, color = 'k')
    plt.plot(velocity, error, color = 'grey')
    plt.title("Normalized Flux vs Velocity")
    plt.xlabel("Velocity (km/s)")
    plt.ylabel("Normalized Flux")
    plt.xlim(-70000, 0)
    max_peak = (np.mean(flux_normalized) * 1.9)
    min_peak = (np.min(error) - .5) 
    plt.text(-58000, (max_peak - .3), str(spectra_name) + ', z=' + str(redshift) + ', snr=' + str(snr), rotation = 0, fontsize = 8.5)
    plt.axhline(y=0.9, color='r', linestyle= '--')    
    plt.axhline(y=1.0)
    #max_peak = (np.max(flux_normalized[np.where(flux_normalized < 5)]))
    #for i in range(len(vminx)):
        #plt.axvspan(vminx[i], vmaxx[i], alpha=0.2, color='red')
    for i in range(len(vminz)):
        #plt.axvline(x=vminz[i], color='g', linestyle= '-')
        plt.axvspan(vminz[i],vmaxz[i], alpha = 0.2, color = 'yellow')

    plt.axvline(x=c0, color='k', linestyle= '-')
    plt.axvline(x=c2, color='b', linestyle= '-')
    plt.axvline(x=o1, color='y', linestyle= '-')
    plt.ylim((min_peak, max_peak))
    savefile_name.savefig()
    plt.close() 
 
 
 
"""
    # Calculate where CIV, CII and OI would be for each pair of vmin and vmax *if* the EHVO absorption found were 
    # instead not EHVO and due to SiIV: 
    # If the absorption is SiIV, this finds and plots where CIV, CII and OI would be
                    z_absSiIV = (wavelength[current_velocity_index]/avr_SiIV_doublet)-1    #<-- right now it does it in the loop value, ...
                    # ... that in the BI program is called current_velocity_index ( for jjjs in jjj:). We want to rewrite this so i can be a module. 
        
                    axvspan(beta[vmins_index],beta[vmaxs_index], alpha=0.2, color='red')
                        
                    z_absSiIV_final = (wavelength[vmaxs_index]/avr_SiIV_doublet)-1.
                    
                    obs_wavelength_Cfinal=(z_absSiIV_final+1.)*(avr_CIV_doublet)
                    obs_wavelength_Cfinal_index =np.min (where (wavelength>obs_wavelength_Cfinal))
                    obs_wavelength_C_final_vel=beta[obs_wavelength_Cfinal_index]
                    axvspan(obs_wavelength_C_vel,obs_wavelength_C_final_vel, alpha=0.2, color='grey')
            plot((obs_wavelength_C_vel, obs_wavelength_C_vel),(-1,10),'k-')#  <-- plot a line at the vmin and vmax of the same color
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
    # the axvspan using beta need to be defined earlier to plot here f.e., axvspan(beta[current_velocity_index],beta[current_velocity_index-1], alpha=0.05, color='red')
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
"""