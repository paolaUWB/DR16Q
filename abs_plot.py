########################################### IMPORTS ############################################################################
from matplotlib import pyplot as plt
import numpy as np
from data_types import Range

######################################## VERNER TABLE CONSTANTS ################################################################
# verner table data
WAVELENGTH_CIV_EMIT_LIMIT = Range(1548.1950, 1550.7700)                                    #never used?
AVERAGE_CIV_DOUBLET = 1549.0524 #weighted average
AVERAGE_SiIV_DOUBLET = 1396.747 # weighted average; individuals: 1402.770, 1393.755
AVERAGE_NV_DOUBLET = 1240.15 # weighted average; individuals: 1242.80, 1238.82             # not used
AVERAGE_OVI_DOUBLET=1033.8160 # weighted average; individuals: 1037.6167, 1031. 9261        # also not used
CII_EMITTED = 1335.313 # (weighted average); individuals:
OI_EMITTED = 1303.4951 # weighted average; individuals pag 20 in Verner Table
###############################################################################################################################


def draw_abs_figure(spectra_index, velocity, flux_normalized, error, savefile_name, spectra_name, redshift, snr):
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
    plt.xlabel("Velocity (km/s)")
    plt.ylabel("Normalized Flux")
    plt.xlim(-70000, 0)
    max_peak = (np.mean(flux_normalized) * 2)
    min_peak = (np.min(error) - .5) 
    plt.title(str(spectra_index) + ': ' + str(spectra_name) + ', z=' + str(redshift) + ' snr=' + str(snr))
    plt.axhline(y = 0.9, color='r', linestyle = '--')    
    plt.axhline(y = 1.0)
    plt.ylim(min_peak, max_peak)
    savefile_name.savefig()
    plt.close()

def vmin_plot(beta, wavelength, current_velocity_index, BALNICITY_INDEX_LIMIT):
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
    # Calculate where CIV, CII and OI would be for each pair of VMIN *if* the EHVO absorption found were 
    # instead not EHVO and due to SiIV: 
    z_absSiIV = (wavelength[current_velocity_index] / AVERAGE_SiIV_DOUBLET) - 1
        
    obs_wavelength_C = (z_absSiIV + 1) * (AVERAGE_CIV_DOUBLET)
    obs_wavelength_C_index = np.min(np.where(wavelength > obs_wavelength_C))
    obs_wavelength_C_vel = beta[obs_wavelength_C_index] + BALNICITY_INDEX_LIMIT
    plt.plot((obs_wavelength_C_vel, obs_wavelength_C_vel),(-1,10),'k-')

    obs_wavelength_CII = (z_absSiIV + 1) * (CII_EMITTED)
    obs_wavelength_CII_index = np.min(np.where(wavelength > obs_wavelength_CII))                  
    obs_wavelength_CII_vel = beta[obs_wavelength_CII_index] + BALNICITY_INDEX_LIMIT
    plt.plot((obs_wavelength_CII_vel, obs_wavelength_CII_vel),(-1,10),'b-')

    obs_wavelength_OI = (z_absSiIV + 1) * (OI_EMITTED)
    obs_wavelength_OI_index = np.min(np.where(wavelength > obs_wavelength_OI))                  
    obs_wavelength_OI_vel = beta[obs_wavelength_OI_index] + BALNICITY_INDEX_LIMIT
    plt.plot((obs_wavelength_OI_vel, obs_wavelength_OI_vel),(-1,10),'y-')

    return [obs_wavelength_C_vel, obs_wavelength_CII_vel, obs_wavelength_OI_vel]


def vmax_plot(beta, wavelength, vmaxs_index, obs_wavelength_C_vel, obs_wavelength_CII_vel, obs_wavelength_OI_vel):
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

    # Calculate where CIV, CII and OI would be for each pair of VMAX *if* the EHVO absorption found were 
    # instead not EHVO and due to SiIV: 
    # if the absorption is SiIV, this finds and plots where CIV, CII and OI would be ###########
    z_absSiIV_final = (wavelength[vmaxs_index] / AVERAGE_SiIV_DOUBLET) - 1.

    obs_wavelength_Cfinal = (z_absSiIV_final + 1.) * (AVERAGE_CIV_DOUBLET)
    obs_wavelength_Cfinal_index = np.min(np.where(wavelength > obs_wavelength_Cfinal))
    obs_wavelength_C_final_vel = beta[obs_wavelength_Cfinal_index]
    plt.axvspan(obs_wavelength_C_vel, obs_wavelength_C_final_vel, alpha = 0.2, color = 'grey')

    obs_wavelength_CIIfinal = (z_absSiIV_final + 1.) * (CII_EMITTED)
    obs_wavelength_CIIfinal_index = np.min (np.where (wavelength > obs_wavelength_CIIfinal))
    obs_wavelength_CII_final_vel = beta[obs_wavelength_CIIfinal_index]
    plt.axvspan(obs_wavelength_CII_vel,obs_wavelength_CII_final_vel, alpha = 0.2, color = 'blue')

    obs_wavelength_OIfinal = (z_absSiIV_final + 1.) * (OI_EMITTED)
    obs_wavelength_OIfinal_index = np.min(np.where (wavelength > obs_wavelength_OIfinal))
    obs_wavelength_OI_final_vel = beta[obs_wavelength_OIfinal_index]
    plt.axvspan(obs_wavelength_OI_vel,obs_wavelength_OI_final_vel, alpha = 0.2, color = 'yellow')

def vmin_line(beta, index):
    plt.plot((beta[index], beta[index]), (-1,10),'r-')

def span_vmin_vmax(beta, vmins, vmaxs):
    plt.axvspan(beta[vmins], beta[vmaxs], alpha = 0.2, color = 'red')

def black_line(beta, index):
    plt.plot((beta[index + 1], beta[index]), (1.5,1.5),'k-')
