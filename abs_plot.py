from matplotlib import pyplot as plt
import numpy as np

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
    max_peak = (np.mean(flux_normalized) * 1.9)
    min_peak = (np.min(error) - .5) 
    plt.title(str(spectra_name) + ', z=' + str(redshift) + ' snr=' + str(snr))
    plt.axhline(y = 0.9, color='r', linestyle = '--')    
    plt.axhline(y = 1.0)
    plt.ylim((min_peak, max_peak))
    savefile_name.savefig()
    plt.close() 