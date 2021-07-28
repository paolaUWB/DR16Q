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

    max_peak = (np.mean(flux_normalized) * 1.9)
    min_peak = (np.min(error) - .5) 

    
    plt.plot(velocity, flux_normalized, color = 'k')
    plt.plot(velocity, error, color = 'grey')
    plt.xlabel("Velocity (km/s)")
    plt.ylabel("Normalized Flux")
    plt.xlim(-70000, 0)
    plt.title(str(spectra_name) + ', z=' + str(redshift) + ', snr=' + str(snr), rotation = 0, fontsize = 8.5)
    plt.axhline(y=0.9, color='r', linestyle= '--')    
    plt.axhline(y=1.0)

    for i in range(len(vminz)):
        #plt.axvline(x=vminz[i], color='g', linestyle= '-')
        plt.axvspan(vminz[i],vmaxz[i], alpha = 0.2, color = 'yellow')

    plt.axvline(x=c0, color='k', linestyle= '-')
    plt.axvline(x=c2, color='b', linestyle= '-')
    plt.axvline(x=o1, color='y', linestyle= '-')
    plt.ylim((min_peak, max_peak))
    savefile_name.savefig()
    plt.close() 
 