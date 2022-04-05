## FIGURE OUT HOW OUTPUT FILES WILL WORK FOR USE IN BOTH NORMALIZATION AND ABSORPTION

import numpy as np 
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from data_types import RangesData, FigureData, FigureDataOriginal

def draw_dynamic_points(figure_index, wavelength, wavelength_observed_from, wavelength_observed_to, flux, test1: RangesData, test2: RangesData, number_of_anchor_points, anchor_pts, max_peak, bf, cf, z, snr, snr_mean_in_ehvo, spectrum_file_name, FILE, error):
    """ Draws a plot containing a dynamic number of anchor points 

    Parameters:
    -----------
    figure_index: int
        Makes a separate graph for each spectra 
    wavelength: list
        List of data that is being plotted on the x axis
    wavelength_observed_from: int
        The minimum x value to be plotted (xlim)
    wavelength_observed_to: int
        The maximum x value to be plotted (xlim)
    flux: list
        List of the data that is being plotted on the y axis
    test1: RangesData
        Green highlighted test region on the graph
    test2: RangesData
        Pink highlighted test region on the graph 
    number_of_anchor_points: int
        Variable containing the number of anchor points to be plotted
    anchor_pts: list
        A list of the wavelength, flux, and error for each anchor point
    max_peak: float
        Scaling value in the y direction
    bf: float
        Initial parameter of powerlaw curve fit
    cf: float
        Initial parameter of powerlaw curve fit
    z: float
        Redshift
    snr: float
        Signal to noise ratio provided by SDSS
    snr_mean_in_ehvo: float 
        Signal to noise ratio calculated in region we care about
    spectrum_file_name: string
        The name of the spectra
    FILE: path
        Path to file you would like these graphs to save to
    
    Returns
    -------
    None.

    Note:
    -----
    Returns nothing, but draws the spectra of the graph.

    """
    ## COMMENT OUT FOR PRESENTATION/PAPER FIGURES
    plt.title(spectrum_file_name)
    subtitle_text = f"z={z} snr={snr} snr_mean_in_ehvo={snr_mean_in_ehvo}"
    plt.text(((wavelength_observed_from + wavelength_observed_to)/2.7), max_peak - 0.25, subtitle_text)
    ########

    plt.figure(figure_index)
    plt.xlabel("Wavelength[$\AA$]")
    plt.ylabel("Flux[10$^{-17}$ erg/cm$^2$/$\AA$]")

    plt.plot(wavelength, flux, color = "black") ## PLOTS SPECTRA
    plt.plot(wavelength, error, color = "black") ## PLOTS ERROR
    plt.plot(test1.wavelength, test1.flux, color = "blue", linestyle = "-") ## PLOTS TEST REGION 1 (1315-1325) 
    plt.plot(test2.wavelength, test2.flux, color = "deeppink", linestyle = "-") ## PLOTS TEST REGION 2 (1350-1360)

    for i in number_of_anchor_points: ## PLOTS ANCHOR POINTS
        plt.plot(anchor_pts[i-1][0], anchor_pts[i-1][1], 'ro')
    plt.plot(wavelength, powerlaw(wavelength, bf, cf), color = "red", linestyle = "--")
    
    plt.xlim(wavelength_observed_from, wavelength_observed_to)
    plt.ylim(-2, max_peak + (max_peak / 1.5))
    
    if FILE == 'null': ## SAVES FIGURES
        plt.show()
    else:
        plt.tight_layout()
        FILE.savefig()
        plt.show()


def draw_dynamic(wavelength, wavelength_observed_from, wavelength_observed_to, flux, test1: RangesData, test2: RangesData, max_peak):
    ### ADD DOCUMENTATION - this plots the spectra to decide where to put anchor points -- not intended to be saved

    fig = plt.figure()
    ax = fig.gca()
    plt.plot(wavelength, flux, color = "black")
    plt.plot(test1.wavelength, test1.flux, color = "blue", linestyle = "-")
    plt.plot(test2.wavelength, test2.flux, color = "deeppink", linestyle = "-")
    plt.xlim(np.min(wavelength), np.max(wavelength)) #(wavelength_observed_from, wavelength_observed_to)
    np.arange(wavelength_observed_from, wavelength_observed_to, 100)
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.xaxis.set_major_formatter('{x:.0f}')
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.grid(which='minor', alpha=0.2)
    ax.grid(which='major', alpha=1)
    plt.ylim(-2, max_peak + (max_peak / 1.5))
    plt.show()

def powerlaw(wavelength, b, c):
    """ Calculates the power law. 

    Parameters
    ----------
    wavelength: array
        Comes from RangesData().    
    b: int
        Initial parameter of powerlaw. 
    c: float
        Initial parameter of powerlaw.

    Returns
    -------
    array
        Power law value in the form of an array.
    """
    return b * (np.power(wavelength, c))

def draw_original_figure(figure_index: int, original_ranges: RangesData, data: FigureDataOriginal, test1: RangesData, test2: RangesData, wavelength_observed_from, wavelength_observed_to, max_peak, FILE, flags, anchor_pts, z):
    """ Draws the original spectra graph.

    Parameters
    ----------
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

    Returns
    -------
    None.

    Note:
    -----
    Returns nothing, but draws the original spectra of the graph.
    """

    main_color = "black" # "midnightblue" 
    test_1_color, test_2_color = "blue", "deeppink"

    plt.figure(figure_index) 

    ## COMMENT OUT FOR PRESENTATION/PAPER FIGURES
    plt.title(data.FigureData.spectrum_file_name + flags)
    subtitle_text = f"z={data.FigureData.z} snr={data.FigureData.snr} snr_mean_in_ehvo={data.FigureData.snr_mean_in_ehvo}"
    plt.text(((data.FigureData.wavelength_from + data.FigureData.wavelength_to)/2.3), max_peak + 1, subtitle_text)
    ########

    plt.xlabel("Wavelength[$\AA$]")
    plt.ylabel("Flux[10$^{-17}$ erg/cm$^2$/$\AA$]")

    plt.plot(original_ranges.wavelength, original_ranges.flux, color = main_color, linestyle = "-") ## PLOTS SPECTRA
    plt.plot(original_ranges.wavelength, original_ranges.error, color = "black", linestyle = "-") ## PLOTS ERROR
    plt.plot(original_ranges.wavelength, powerlaw(original_ranges.wavelength, data.bf, data.cf), color = "red", linestyle = "--") ## PLOTS POWER LAW
    #plt.plot(original_ranges.wavelength, powerlaw(original_ranges.wavelength, data.bf, data.cf), color = "red", zorder = 3, linestyle = "--") ## PLOTS POWER LAW (DECIDE WHICH TO KEEP - THIS OR ABOVE LINE)

    plt.plot(data.power_law_data_x, data.power_law_data_y, 'ro') ## PLOTS ANCHOR POINTS ?? WHAT DOES THIS DO ??
    #plt.plot(anchor_pts[0][0], anchor_pts[0][1], 'ro')
    #plt.plot(anchor_pts[1][0], anchor_pts[1][1], 'ro')
    #plt.plot(anchor_pts[2][0], anchor_pts[2][1], 'ro')
    plt.plot(test1.wavelength, test1.flux, color = test_1_color, linestyle = "-") ## PLOTS TEST REGION 1 (1315-1325)
    plt.plot(test2.wavelength, test2.flux, color = test_2_color, linestyle = "-") ## PLOTS TEST REGION 2 (1350-1360)
    
    ##########################################################################################
    ############################### PLOTTING FOR PRESENTATIONS ###############################
    ##########################################################################################
    
    ### TO PLOT ANCHOR POINT *REGIONS*
    # ## LEFT POINT (1280-1290)
    # plt.axvline(1280 * (z + 1), 0, 1, color = "black")
    # plt.axvline(1290 * (z + 1), 0, 1, color = "black")
    # plt.axvspan(1280 * (z + 1), 1290 * (z + 1), 0, 1, alpha = 0.5, color = "lightskyblue")
    
    # ## MIDDLE POINT (1440-1450)
    # plt.axvline(1440 * (z + 1), 0, 1, color = "black")
    # plt.axvline(1450 * (z + 1), 0, 1, color = "black")
    # plt.axvspan(1440 * (z + 1), 1450 * (z + 1), 0, 1, alpha = 0.5, color = "lightskyblue")

    # ## TRY MIDDLE POINT - LEFT (1420-1430)
    # plt.axvline(1420 * (z + 1), 0, 1, color = "black")
    # plt.axvline(1430 * (z + 1), 0, 1, color = "black")
    # plt.axvspan(1420 * (z + 1), 1430 * (z + 1), 0, 1, alpha = 0.5, color = "lightskyblue")
   
    # ## TRY MIDDLE POINT - RIGHT (1460-1470)
    # plt.axvline(1460 * (z + 1), 0, 1, color = "black")
    # plt.axvline(1470 * (z + 1), 0, 1, color = "black")
    # plt.axvspan(1460 * (z + 1), 1470 * (z + 1), 0, 1, alpha = 0.5, color = "lightskyblue")    
    
    # ## RIGHT POINT (1690-1710)
    # plt.axvline(1690 * (z + 1), 0, 1, color = "black")
    # plt.axvline(1710 * (z + 1), 0, 1, color = "black")
    # plt.axvspan(1690 * (z + 1), 1710 * (z + 1), 0, 1, alpha = 0.5, color = "lightskyblue")

    ### PLOTTING REGION WE CARE ABOUT (1250-1400)
    # plt.axvline(1250 * (z + 1), 0, 1, color = "black")
    # plt.axvline(1400 * (z + 1), 0, 1, color = "black")
    # plt.axvspan(1250 * (z + 1), 1400 * (z + 1), 0, 1, alpha = 0.5, color = "lightskyblue")
    
    ##########################################################################################
    
    plt.xlim(wavelength_observed_from, wavelength_observed_to)
    plt.ylim(-2, max_peak + (max_peak / 1.5))
    FILE.savefig(bbox_inches='tight')
    plt.close(figure_index)

def draw_normalized_figure(figure_index: int, original_ranges: RangesData, figure_data: FigureData, flux_normalized, error_normalized,
                            test1: RangesData, test2: RangesData, normalized_flux_test_1, normalized_flux_test_2, wavelength_observed_from, wavelength_observed_to, max_peak, FILE, min_flux_green_region, min_flux_pink_region, max_flux_green_region, max_flux_pink_region, z, val):
    """ Draws the normalized spectra graph.

    Parameters
    ----------
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

    Returns
    -------
    None.
    
    Notes
    -----
    Creates a graph of the spectra and saves to the original_graphs.pdf
    """

    main_color = "black"
    test_1_color, test_2_color = "blue", "deeppink"
    
    plt.figure(figure_index) 

    ## COMMENT OUT FOR PRESENTATION/PAPER FIGURES
    plt.title("Normalized Data vs. Normalized Error")
    subtitle_text = f"z={figure_data.z} snr={figure_data.snr} snr_mean_in_ehvo={figure_data.snr_mean_in_ehvo}"
    plt.text(((figure_data.wavelength_from + figure_data.wavelength_to)/2.3), max_peak + 0.25, figure_data.spectrum_file_name)
    plt.text(((figure_data.wavelength_from + figure_data.wavelength_to)/2.3), np.max(flux_normalized)/1.07, figure_data.spectrum_file_name)
    # SHOULD WE USE THE LINES BELOW OR ABOVE TO PLACE THESE SUBTITLE_TEXTS ??
    # plt.text(((figure_data.wavelength_from + figure_data.wavelength_to)/2.3), max_peak, subtitle_text)
    # plt.text(((figure_data.wavelength_from + figure_data.wavelength_to)/2.3), np.max(flux_normalized), subtitle_text)
    ########

    plt.xlabel("Wavelength [$\AA$]")
    plt.ylabel("Normalized Flux")

    plt.plot(original_ranges.wavelength, flux_normalized, color = main_color, zorder = 1, linestyle = "-") ## PLOTS NORMALIZED SPECTRA
    plt.plot(original_ranges.wavelength, error_normalized, color = "black", linestyle = "-") ## PLOTS ERROR
    plt.plot((original_ranges.wavelength[0], original_ranges.wavelength[-1]), (1, 1), zorder = 2, color = "red", linestyle = "-", linewidth = 2) ## PLOTS LINE AT 1 ??

    plt.plot(test1.wavelength, normalized_flux_test_1, color = test_1_color, zorder = 1, linestyle = "-") ## PLOTS TEST REGION 1 (1315-1325)
    plt.plot(test2.wavelength, normalized_flux_test_2, color = test_2_color, zorder = 1, linestyle = "-") ## PLOTS TEST REGION 2 (1350-1360)

    ##########################################################################################
    ############################### PLOTTING FOR PRESENTATIONS ###############################
    ##########################################################################################
   
    ### TO PLOT RANGE ALLOWED FOR GOOD FIT (TEST#2)
    ## GREEN/BLUE
    # plt.axhspan(np.median(normalized_flux_test_1) - 0.05, np.median(normalized_flux_test_1) + 0.05, xmin = 0, xmax = 1, zorder = 1, alpha = 0.5, color = "grey")
    # plt.hlines(np.median(normalized_flux_test_1) + 0.05, wavelength_observed_from, wavelength_observed_to, zorder=2, color='black')
    # plt.hlines(np.median(normalized_flux_test_1) - 0.05, wavelength_observed_from, wavelength_observed_to, zorder=2, color='black')
    # plt.hlines(np.median(normalized_flux_test_1), wavelength_observed_from, wavelength_observed_to, zorder = 3, color = 'black', linestyle = 'dashed')
    # plt.axvline(1315 * (z + 1), 0, 1, color = "black", linestyle = "dashed")
    # plt.axvline(1325 * (z + 1), 0, 1, color = "black", linestyle = "dashed")

    ## PINK
    # plt.axhspan(np.median(normalized_flux_test_2) - 0.05, np.median(normalized_flux_test_2) + 0.05, xmin = 0, xmax = 1, zorder = 1, alpha = 0.5, color = "grey")
    # plt.hlines(np.median(normalized_flux_test_2) + 0.05, wavelength_observed_from, wavelength_observed_to, zorder=2, color='black')
    # plt.hlines(np.median(normalized_flux_test_2) - 0.05, wavelength_observed_from, wavelength_observed_to, zorder=2, color='black')
    # plt.hlines(np.median(normalized_flux_test_2), wavelength_observed_from, wavelength_observed_to, zorder = 3, color = 'black', linestyle = 'dashed')
    # plt.axvline(1350 * (z + 1), 0, 1, color = "black", linestyle = "dashed")
    # plt.axvline(1360 * (z + 1), 0, 1, color = "black", linestyle = "dashed")
    

    ### TO PLOT RANGE ALLOWED FOR GOOD FIT (TEST#3)
    # ## GREEN/BLUE
    # plt.hlines(min_flux_green_region, np.min(wavelength_observed_from), np.max(wavelength_observed_to), zorder=3, color='black')
    # plt.hlines(max_flux_green_region, np.min(wavelength_observed_from), np.max(wavelength_observed_to), zorder=3, color='black')
    # plt.axhspan(min_flux_green_region, max_flux_green_region, xmin = 0, xmax = 1, zorder = 2, alpha = 0.5, color = 'grey')
    # plt.hlines(np.median(normalized_flux_test_1), wavelength_observed_from, wavelength_observed_to, zorder = 3, color = 'black', linestyle = 'dashed', linewidth = 2)
    # plt.axvline(1315 * (z + 1), 0, 1, color = "black", linestyle = "dashed")
    # plt.axvline(1325 * (z + 1), 0, 1, color = "black",linestyle = "dashed")

    # ## PINK
    # plt.hlines(min_flux_pink_region, np.min(wavelength_observed_from), np.max(wavelength_observed_to), zorder=2, color='black')
    # plt.hlines(max_flux_pink_region, np.min(wavelength_observed_from), np.max(wavelength_observed_to), zorder=2, color='black')
    # plt.axhspan(min_flux_pink_region, max_flux_pink_region, xmin = 0, xmax = 1, zorder = 1, alpha = 0.5, color = "grey")
    # plt.hlines(np.median(normalized_flux_test_2), wavelength_observed_from, wavelength_observed_to, zorder = 3, color = 'black', linestyle = 'dashed', linewidth = 2)
    # plt.axvline(1350 * (z + 1), 0, 1, color = "black", linestyle = "dashed")
    # plt.axvline(1360 * (z + 1), 0, 1, color = "black", linestyle = "dashed")
    
    
    ##########################################################################################
    
    plt.xlim(1250 * (1 + z), 1400 * (1 + z))
    #plt.xlim(wavelength_observed_from, wavelength_observed_to)
    #plt.ylim(0, np.max(flux_normalized) + 1)
    #plt.ylim(0, max_peak + (max_peak / 4))
    plt.ylim(0,2)# max_peak)# + (max_peak / 1.5))
    FILE.savefig(bbox_inches='tight')
    plt.close(figure_index)

# def draw_normalized_figure(figure_index: int, original_ranges: RangesData, figure_data: FigureData, flux_normalized, error_normalized,
#                             test1: RangesData, test2: RangesData, normalized_flux_test_1, normalized_flux_test_2, wavelength_observed_from, wavelength_observed_to, max_peak, FILE):#, min_flux_green_region, min_flux_pink_region, max_flux_green_region, max_flux_pink_region):
#     """ Draws the normalized spectra graph.

#     Parameters
#     ----------
#     figure_index: int
#         Makes a separate graph for each spectra. 
#     original_ranges: RangesData
#         Ranges of values for the original data.
#     figure_data: FigureData
#         Data from DR9Q (for now...).
#     flux_normalized: array
#     error_normalized: array
#     test1: RangesData
#         Green highlighted area on graph. 
#     test2: RangesData
#         Pink highlighted area on graph.
#     normalized_flux_test_1: any
#     normalized_flux_test_2: any

#     Returns
#     -------
#     None.
    
#     Notes
#     -----
#     Creates a graph of the spectra and saves to the original_graphs.pdf
#     """

#     main_color = "midnightblue"
#     test_1_color, test_2_color = "xkcd:green apple", "xkcd:bubblegum"
#     subtitle_text = f"z={figure_data.z} snr={figure_data.snr} snr_mean_in_ehvo={figure_data.snr_mean_in_ehvo}"
#     plt.figure(figure_index) 
#     plt.text(((figure_data.wavelength_from + figure_data.wavelength_to)/2.3), max_peak + 0.25, figure_data.spectrum_file_name)
#     plt.text(((figure_data.wavelength_from + figure_data.wavelength_to)/2.3), max_peak, subtitle_text)
#     #plt.text(((figure_data.wavelength_from + figure_data.wavelength_to)/2.3), np.max(flux_normalized)/1.07, figure_data.spectrum_file_name)
#     #plt.text(((figure_data.wavelength_from + figure_data.wavelength_to)/2.3), np.max(flux_normalized), subtitle_text)
#     plt.title(figure_data.spectrum_file_name)
#     plt.plot(original_ranges.wavelength, flux_normalized, color = main_color, linestyle = "-")
#     plt.plot(original_ranges.wavelength, error_normalized, color = "black", linestyle = "-")
#     plt.title("Normalized Data vs. Normalized Error")
#     plt.xlabel("Wavelength [A]")
#     plt.ylabel("Normalized Flux[10^[-17]]cgs")
#     plt.plot(test1.wavelength, normalized_flux_test_1, color = test_1_color, linestyle = "-")
#     plt.plot(test2.wavelength, normalized_flux_test_2, color = test_2_color, linestyle = "-")
#     plt.plot((original_ranges.wavelength[0], original_ranges.wavelength[-1]), (1, 1), color = "red", linestyle = "-")

#     #plt.hlines(y=min_flux_green_region, xmin = np.min(wavelength_observed_from), xmax = np.max(wavelength_observed_to), zorder=1, linewidth = 0.5, color='xkcd:green apple')
#     #plt.hlines(y=max_flux_green_region, xmin = np.min(wavelength_observed_from), xmax = np.max(wavelength_observed_to), zorder=1, linewidth = 0.5, color='xkcd:green apple')
#     #plt.hlines(y=min_flux_pink_region, xmin = np.min(wavelength_observed_from), xmax = np.max(wavelength_observed_to), zorder=1, linewidth = 0.5, color='xkcd:bubblegum')
#     #plt.hlines(y=max_flux_pink_region, xmin = np.min(wavelength_observed_from), xmax = np.max(wavelength_observed_to), zorder=1, linewidth = 0.5, color='xkcd:bubblegum')
    
#     plt.xlim(wavelength_observed_from, wavelength_observed_to)
#     #plt.ylim(0, np.max(flux_normalized) + 1)
#     #plt.ylim(0, max_peak + (max_peak / 4))
#     plt.ylim(0, 6)
#     FILE.savefig()
#     plt.close(figure_index)
