## FIGURE OUT HOW OUTPUT FILES WILL WORK FOR USE IN BOTH NORMALIZATION AND ABSORPTION

import numpy as np 
from matplotlib import pyplot as plt
from data_types import RangesData, FigureData, FigureDataOriginal

def draw_dynamic_points(figure_index, wavelength, wavelength_observed_from, wavelength_observed_to, flux, test1: RangesData, test2: RangesData, number_of_anchor_points, anchor_pts, max_peak, bf, cf, z, snr, snr_mean_in_ehvo, spectrum_file_name, FILE):
    plt.figure(figure_index)
    plt.title(spectrum_file_name)
    plt.xlabel("Wavelength[A]")
    plt.ylabel("Flux[10^[-17]]cgs")
    subtitle_text = f"z={z} snr={snr} snr_mean_in_ehvo={snr_mean_in_ehvo}"
    plt.text(((wavelength_observed_from + wavelength_observed_to)/2.3), max_peak + 0.25, subtitle_text)
    plt.plot(wavelength, flux, color = "xkcd:ultramarine")
    plt.plot(test1.wavelength, test1.flux, color = "xkcd:green apple", linestyle = "-")
    plt.plot(test2.wavelength, test2.flux, color = "xkcd:bubblegum", linestyle = "-")
    for i in number_of_anchor_points:
        plt.plot(anchor_pts[i-1][0], anchor_pts[i-1][1], 'ro')
    plt.plot(wavelength, powerlaw(wavelength, bf, cf), color = "red", linestyle = "--")
    plt.xlim(wavelength_observed_from, wavelength_observed_to)
    plt.ylim(-2, max_peak)
    FILE.savefig()
    plt.show()

def draw_dynamic(wavelength, wavelength_observed_from, wavelength_observed_to, flux, test1: RangesData, test2: RangesData, max_peak):
    plt.plot(wavelength, flux, color = "xkcd:ultramarine")
    plt.plot(test1.wavelength, test1.flux, color = "xkcd:green apple", linestyle = "-")
    plt.plot(test2.wavelength, test2.flux, color = "xkcd:bubblegum", linestyle = "-")
    plt.xlim(wavelength_observed_from, wavelength_observed_to)
    plt.ylim(-2, max_peak)
    plt.show()

def powerlaw(wavelength, b, c):
    """ Calculates the power law. 

    Parameters:
    -----------
    wavelength: array
        Comes from RangesData().    
    b: int
        Initial parameter of powerlaw. 
    c: float
        Initial parameter of powerlaw.

    Returns:
    --------
    array
        Power law value in the form of an array.
    """
    return b * (np.power(wavelength, c))

def draw_original_figure(figure_index: int, original_ranges: RangesData, data: FigureDataOriginal, test1: RangesData, test2: RangesData, wavelength_observed_from, wavelength_observed_to, max_peak, FILE):
    """ Draws the original spectra graph.

    Parameters:
    -----------
    figure_index: int
        Makes a separate graph for each spectra. 
    original_ranges: RangesData
        Ranges of values for the original data.
    data: FigureDataOriginal
        Data from DR9Q (for now...).
    test1: RangesData
        Green highlighted area on graph. 
    test2: RangesData
        Pink highlighted area on graph.
    max_peak: any
        Max peak value of data per spectra.

    Returns:
    --------
    None.

    Note:
    -----
    Returns nothing, but draws the original spectra of the graph.
    """

    main_color = "xkcd:ultramarine"
    test_1_color, test_2_color = "xkcd:green apple", "xkcd:bubblegum"
    subtitle_text = f"z={data.FigureData.z} snr={data.FigureData.snr} snr_mean_in_ehvo={data.FigureData.snr_mean_in_ehvo}"
    plt.figure(figure_index)
    plt.title(data.FigureData.spectrum_file_name)
    plt.xlabel("Wavelength[A]")
    plt.ylabel("Flux[10^[-17]]cgs")
    plt.text(((data.FigureData.wavelength_from + data.FigureData.wavelength_to)/2.3), max_peak + 1, subtitle_text)
    plt.plot(original_ranges.wavelength, original_ranges.flux, color = main_color, linestyle = "-")
    plt.plot(data.power_law_data_x, data.power_law_data_y, 'ro')
    plt.plot(original_ranges.wavelength, original_ranges.error, color = "black", linestyle = "-")
    plt.plot(test1.wavelength, test1.flux, color = test_1_color, linestyle = "-")
    plt.plot(test2.wavelength, test2.flux, color = test_2_color, linestyle = "-")
    plt.plot(original_ranges.wavelength, powerlaw(original_ranges.wavelength, data.bf, data.cf), color = "red", linestyle = "--")
    plt.xlim(wavelength_observed_from, wavelength_observed_to)
    plt.ylim(-2, max_peak + 2)
    FILE.savefig()
    plt.close(figure_index)

def draw_normalized_figure(figure_index: int, original_ranges: RangesData, figure_data: FigureData, flux_normalized, error_normalized, #1 param more since I removed the tuple
                            test1: RangesData, test2: RangesData, normalized_flux_test_1, normalized_flux_test_2, wavelength_observed_from, wavelength_observed_to, max_peak, FILE):
    """ Draws the normalized spectra graph.

    Parameters:
    -----------
    figure_index: int
        Makes a separate graph for each spectra. 
    original_ranges: RangesData
        Ranges of values for the original data.
    figure_data: FigureData
        Data from DR9Q (for now...).
    flux_normalized: array
    error_normalized: array
    test1: RangesData
        Green highlighted area on graph. 
    test2: RangesData
        Pink highlighted area on graph.
    normalized_flux_test_1: any
    normalized_flux_test_2: any

    Returns:
    --------
    None.
    
    Notes:
    ------
    Creates a graph of the spectra and saves to the original_graphs.pdf
    """

    main_color = "xkcd:ultramarine"
    test_1_color, test_2_color = "xkcd:green apple", "xkcd:bubblegum"
    subtitle_text = f"z={figure_data.z} snr={figure_data.snr} snr_mean_in_ehvo={figure_data.snr_mean_in_ehvo}"
    plt.figure(figure_index) 
    plt.text(((figure_data.wavelength_from + figure_data.wavelength_to)/2.3), max_peak + 1, figure_data.spectrum_file_name)
    plt.text(((figure_data.wavelength_from + figure_data.wavelength_to)/2.3), max_peak + 0.5, subtitle_text)
    #plt.text(((figure_data.wavelength_from + figure_data.wavelength_to)/2.3), np.max(flux_normalized)/1.07, figure_data.spectrum_file_name)
    #plt.text(((figure_data.wavelength_from + figure_data.wavelength_to)/2.3), np.max(flux_normalized), subtitle_text)
    plt.title(figure_data.spectrum_file_name)
    plt.plot(original_ranges.wavelength, flux_normalized, color = main_color, linestyle = "-")
    plt.plot(original_ranges.wavelength, error_normalized, color = "black", linestyle = "-")
    plt.title("Normalized Data vs. Normalized Error")
    plt.xlabel("Wavelength [A]")
    plt.ylabel("Normalized Flux[10^[-17]]cgs")
    plt.plot(test1.wavelength, normalized_flux_test_1, color = test_1_color, linestyle = "-")
    plt.plot(test2.wavelength, normalized_flux_test_2, color = test_2_color, linestyle = "-")
    plt.plot((original_ranges.wavelength[0], original_ranges.wavelength[-1]), (1, 1), color = "red", linestyle = "-")
    plt.xlim(wavelength_observed_from, wavelength_observed_to)
    plt.ylim(-2, np.max(flux_normalized) + 1)
    FILE.savefig()
    plt.close(figure_index)
