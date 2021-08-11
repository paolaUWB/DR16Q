########################################### IMPORTS ############################################################################
from matplotlib import pyplot as plt
import numpy as np
from data_types import Range

######################################## VERNER TABLE CONSTANTS ################################################################
WAVELENGTH_CIV_EMIT_LIMIT = Range(1548.1950, 1550.7700)                                    #never used?
AVERAGE_CIV_DOUBLET = 1549.0524 #weighted average
AVERAGE_SiIV_DOUBLET = 1396.747 # weighted average; individuals: 1402.770, 1393.755
AVERAGE_NV_DOUBLET = 1240.15 # weighted average; individuals: 1242.80, 1238.82             # not used
AVERAGE_OVI_DOUBLET=1033.8160 # weighted average; individuals: 1037.6167, 1031. 9261        # also not used
CII_EMITTED = 1335.313 # (weighted average); individuals:
OI_EMITTED = 1303.4951 # weighted average; individuals pag 20 in Verner Table
###############################################################################################################################


def draw_abs_figure(spectra_count_abs, spectra_index, velocity, flux_normalized, error, savefile_name, spectra_name, redshift, snr, max_peak):
    """ Makes a flux vs velocity graph, that also has the error vs velocity on the same graph. Has text that identifies 
    what the graph number is, what the spectra name is, the signal to noise ratio, and the redshift value used.
    
    Parameters
    ----------
    spectra_count_abs: int or list
        Keeping record of how many absorption plots have been found.

    spectra_index: int or list
        The index number of the graph you are plotting.

    velocity: array or list
        The velocity values to plot on the x-axis.

    flux_normalized: array or list
        The flux values that will be plotted on the y-axis.

    error: array or list
        The error values to be plotted.

    savefile_name: string
        The name that you want the plot to be saved as.

    spectra_name: string
        The full name of the spectra. 

    redshift: int or list or array
        The redshift value used for calculations.

    snr: int or list or array
        The signal to noise ratio value.

    max_peak: int or array
        The max peak value which is used to scale the y-axis. 

    Returns
    -------
    None.
    """
    plt.plot(velocity, flux_normalized, color = 'k')
    plt.plot(velocity, error, color = 'grey')
    plt.xlabel("Velocity (km/s)")
    plt.ylabel("Normalized Flux")
    plt.xlim(-70000, 0)
    snr = round(snr, 2)
    min_peak = -0.1
    plt.title(str(spectra_count_abs) + ' abs |' + str(spectra_index) + ' tot: ' + str(spectra_name) + ', z=' + str(redshift) + ' snr=' + str(snr))
    plt.axhline(y = 0.9, color='r', linestyle = '--')    
    plt.axhline(y = 1.0)
    plt.ylim(min_peak, max_peak + (max_peak / 4))
    savefile_name.savefig()
    plt.close()

def vmin_plot_IF(beta, wavelength, current_velocity_index, BALNICITY_INDEX_LIMIT):
    """ Finds the minimum velocity index and value of where CIV, CII and OI would be *if* the EHVO absorption 
    found were instead not EHVO and due to SiIV, and plots a vertical line for that value.
    
    Parameters
    ----------
    beta: array or list
        The velocity value(s) to be used.

    wavelength: array or list or int
        The wavelength values to be used.

    current_velocity_index: list or int
        When looping through a velocity, the specific velocity index that you are on.

    BALNICITY_INDEX_LIMIT: int
        The minimum value of broad absorption width that we are looking for.

    Returns
    -------
    obs_wavelength_C_vel: array or list
        The minimum value of the velocity for where CIV would be for each pair of VMIN 
        *if* the EHVO absorption found were instead not EHVO and due to SiIV.

    obs_wavelength_CII_vel: array or list
        The minimum value of the velocity for where CII would be for each pair of VMIN 
        *if* the EHVO absorption found were instead not EHVO and due to SiIV.

    obs_wavelength_OI_vel: array or list
        The minimum value of the velocity for where OI would be for each pair of VMIN 
        *if* the EHVO absorption found were instead not EHVO and due to SiIV.
    """
    z_absSiIV = (wavelength[current_velocity_index] / AVERAGE_SiIV_DOUBLET) - 1
        
    obs_wavelength_C = (z_absSiIV + 1) * (AVERAGE_CIV_DOUBLET)
    obs_wavelength_C_index = np.min(np.where(wavelength > obs_wavelength_C))
    obs_wavelength_C_vel = beta[obs_wavelength_C_index] + BALNICITY_INDEX_LIMIT
    plt.plot((obs_wavelength_C_vel, obs_wavelength_C_vel), (-1,10), 'k-')

    obs_wavelength_CII = (z_absSiIV + 1) * (CII_EMITTED)
    obs_wavelength_CII_index = np.min(np.where(wavelength > obs_wavelength_CII))                  
    obs_wavelength_CII_vel = beta[obs_wavelength_CII_index] + BALNICITY_INDEX_LIMIT
    plt.plot((obs_wavelength_CII_vel, obs_wavelength_CII_vel), (-1,10), 'b-')

    obs_wavelength_OI = (z_absSiIV + 1) * (OI_EMITTED)
    obs_wavelength_OI_index = np.min(np.where(wavelength > obs_wavelength_OI))                  
    obs_wavelength_OI_vel = beta[obs_wavelength_OI_index] + BALNICITY_INDEX_LIMIT
    plt.plot((obs_wavelength_OI_vel, obs_wavelength_OI_vel), (-1,10), 'y-')

    return [obs_wavelength_C_vel, obs_wavelength_CII_vel, obs_wavelength_OI_vel]


def vmax_plot_span_IF(beta, wavelength, vmaxs_index, obs_wavelength_C_vel, obs_wavelength_CII_vel, obs_wavelength_OI_vel):
    """ Finds the maximum velocity index and value of where CIV, CII and OI would be *if* the EHVO absorption 
    found were instead not EHVO and due to SiIV, plots a vertical line for that value and creates a colored bar
    that spans from the minimum velocity previously found to this maximum velocity. 
    
    Parameters
    ----------
    beta: array or list
        The list velocity values to be used.

    wavelength: list
        A list of wavelength values to be used.

    vmaxs_index: int or list
        The index value of the maximum velocity found from the beta list to be used as a reference point to find
        the index value of the maximum velocity for where CIV, CII, OI would be *if* the EHVO absorption found 
        were instead not EHVO and due to SiIV.

    obs_wavelength_C_vel: array or list
        The minimum value of the velocity for where CIV would be for each pair of VMIN  *if* the EHVO absorption 
        found were instead not EHVO and due to SiIV.

    obs_wavelength_CII_vel: array or list
        The minimum value of the velocity for where CII would be for each pair of VMIN *if* the EHVO absorption 
        found were instead not EHVO and due to SiIV.

    obs_wavelength_OI_vel: array or list
        The minimum value of the velocity for where OI would be for each pair of VMIN *if* the EHVO absorption 
        found were instead not EHVO and due to SiIV. 

    Returns
    -------
    None.
    """
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
    """ Plots a red notable vertical line of the minimum velocity occurance if absorption found.

    Parameters
    ----------
    beta: array
        The velocity values you are looping through, should specifically be the minimum value of where
        absorption is found.
    index: list
        The value of the index of where the minimum value(s) of absorption is found.

    Returns
    -------
    None.
    """
    plt.plot((beta[index], beta[index]), (-1,10),'r-')

def span_vmin_vmax(beta, vmins, vmaxs):
    """ Plots a red bar that spans from the minimum velocity occurance and the maximum velocity 
    occurance if absorption found. 

    Parameters
    ----------
    beta: array
        The velocity values you are looping through.

    vmins: 
        The minimum velocity point you want to span.

    vmaxs: 
        The maximum velocity point you want to span.

    Returns
    -------
    None.
    """
    plt.axvspan(beta[vmins], beta[vmaxs], alpha = 0.2, color = 'red')

def black_line(beta, index):
    """ Plots a notable black horizontal line of when the minimum broad absorption width value has 
    passed if absorption found. This line is inside of span_vmin_vmax. 

    Parameters
    ----------
    beta: array
        The velocity values you are looping through.
    
    index: list or array
        The current index you are on.

    Returns
    -------
    None.
    """
    plt.plot((beta[index + 1], beta[index]), (1.5,1.5),'k-')
