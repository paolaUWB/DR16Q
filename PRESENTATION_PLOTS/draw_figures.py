import numpy as np 
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from data_types import RangesData, FigureData, FigureDataOriginal
from useful_wavelength_flux_error_modules import wavelength_flux_error_for_points


def draw_dynamic_points(figure_index, data,  wavelength_range, test1: RangesData, test2: RangesData, anchor_pts_data, max_peak, initial_parameters, user_input_wavelength, FILE, FILE2):
    """ Draws a plot containing a dynamic number of anchor points 

    Parameters:
    -----------
    figure_index: int
        Makes a separate graph for each spectra & index
    data: list
        data[0]     wavelength: list
                        List of data that is being plotted on the x axis
        data[1]     flux: list
                        List of the data that is being plotted on the y axis
        data[2]     error: list
                        List of the error for the plot - plotted at bottom of plot
        data[3]     z: float
                        Redshift
        data[4]     snr: float
                        Signal to noise ratio provided by SDSS
        data[5]     snr_mean_in_ehvo: float 
                        Signal to noise ratio calculated in region we care about
        data[6]     spectrum_file_name: string
                        The name of the spectra
    wavelength_range: list
        wavelength_range[0]     wavelength_observed_from: int
                                    The minimum x value to be plotted (xlim)
        wavelength_range[1]     wavelength_observed_to: int
                                    The maximum x value to be plotted (xlim)
    test1: RangesData
        Green highlighted test region on the graph
    test2: RangesData
        Pink highlighted test region on the graph 
    anchor_pts_data: list
        anchor_pts_data[0]      number_of_anchor_points: int
                                    Variable containing the number of anchor points to be plotted
        anchor_pts_data[1]      anchor_pts: list
                                    A list of the wavelength, flux, and error for each anchor point
    max_peak: float
        Scaling value in the y direction
    initial_parameters: list
        initial_parameters[0]    bf: float
                                    Initial parameter ("initial guess") of powerlaw curve fit (flux value when wavelength is 1)
        initial_parameters[1]    cf: float
                                    Initial parameter ("initial guess") for power law curve fit (slope)   
    FILE: path
        Path to file you would like these graphs to save to
    
    Returns
    -------
    None.

    Note:
    -----
    Returns nothing, but draws the spectra of the graph.

    """
    fig = plt.figure(figure_index)
    # plt.figure(figure_index)
    
    ## COMMENT OUT FOR PRESENTATION/PAPER FIGURES
    subtitle_text = f"z={data[3]} snr={data[4]} snr_mean_in_ehvo={data[5]}"
    title = f"{figure_index}: {data[6]}"
    plt.title(title + '\n' + subtitle_text)
    # plt.title(figure_index + data[6] + '\n' + subtitle_text)
    ## ---------------------------------------------------------------------

    ## PLOTS SPECTRUM W/ ANCHOR POINTS
    plt.title(title + '\n' + subtitle_text)
    plt.plot(data[0], data[1], color = "black") ## PLOTS SPECTRA
    plt.plot(data[0], data[2], color = "black") ## PLOTS ERROR
    plt.plot(test1.wavelength, test1.flux, color = "blue", linestyle = "-") ## PLOTS TEST REGION 1 (1315-1325) 
    plt.plot(test2.wavelength, test2.flux, color = "deeppink", linestyle = "-") ## PLOTS TEST REGION 2 (1350-1360)

    anchor_pts_wavelength = []
    anchor_pts_flux = []
    min_wavelength_pts = []
    max_wavelength_pts = []

    for i in anchor_pts_data[0]: ## PLOTS ANCHOR POINTS
        plt.plot(anchor_pts_data[1][i-1][0], anchor_pts_data[1][i-1][1], 'ro')

        anchor_pts_wavelength.append(anchor_pts_data[1][i-1][0])
        anchor_pts_flux.append(anchor_pts_data[1][i-1][1])

        min_wavelength = np.min(np.where(data[0] > anchor_pts_wavelength[i-1] - 100))
        max_wavelength = np.max(np.where(data[0] < anchor_pts_wavelength[i-1] + 100))

        min_wavelength_pts.append(min_wavelength)
        max_wavelength_pts.append(max_wavelength)


    plt.plot(data[0], powerlaw(data[0], initial_parameters[0], initial_parameters[1]), color = "red", linestyle = "--") ## PLOTS POWER LAW
    
    plt.xlabel("Wavelength[$\AA$]")
    plt.ylabel("Flux[10$^{-17}$ erg/cm$^2$/$\AA$]")

    plt.xlim(wavelength_range[0], wavelength_range[1])
    plt.ylim(-2, max_peak + (max_peak / 1.5))

    ## ---------------------------------------------------------------------
    graph_ylim = 10
    max_peak_pt1 = np.max(data[1][min_wavelength_pts[0] + graph_ylim : max_wavelength_pts[0] + graph_ylim])
    min_peak_pt1 = np.min(data[1][min_wavelength_pts[0] - graph_ylim : max_wavelength_pts[0] - graph_ylim])
    max_peak_pt2 = np.max(data[1][min_wavelength_pts[1] + graph_ylim : max_wavelength_pts[1] + graph_ylim])
    min_peak_pt2 = np.min(data[1][min_wavelength_pts[1] - graph_ylim : max_wavelength_pts[1] - graph_ylim])
    max_peak_pt3 = np.max(data[1][min_wavelength_pts[2] + graph_ylim : max_wavelength_pts[2] + graph_ylim])
    min_peak_pt3 = np.min(data[1][min_wavelength_pts[2] - graph_ylim : max_wavelength_pts[2] - graph_ylim])

    fig2 = plt.figure(figure_index + 1)

    ## FULL WAVELENGTH RANGE OF FIT SPECTRUM
    ax1 = fig2.add_subplot(2,3,(1,3))
    plt.title(title + '\n' + subtitle_text)
    plt.plot(data[0], data[1], color = "black")
    plt.plot(test1.wavelength, test1.flux, color = "blue", linestyle = "-") ## PLOTS TEST REGION 1 (1315-1325) 
    plt.plot(test2.wavelength, test2.flux, color = "deeppink", linestyle = "-") ## PLOTS TEST REGION 2 (1350-1360)
    plt.plot(anchor_pts_wavelength[0], anchor_pts_flux[0], 'ro')
    plt.plot(anchor_pts_wavelength[1], anchor_pts_flux[1], 'ro')
    plt.plot(anchor_pts_wavelength[2], anchor_pts_flux[2], 'ro')
    plt.plot(data[0], powerlaw(data[0], initial_parameters[0], initial_parameters[1]), color = "red", linestyle = "--") ## PLOTS POWER LAW
    ax1.xaxis.set_major_locator(MultipleLocator(500))
    ax1.xaxis.set_minor_locator(MultipleLocator(100))
    ax1.tick_params(labelsize=8)
    ax1.grid(which='minor', linewidth=0.3)
    ax1.grid(which='major', linewidth=1)
    plt.xlabel("Wavelength[$\AA$]", fontsize=8)
    plt.ylabel("Flux[10$^{-17}$ erg/cm$^2$/$\AA$]", fontsize=8)
    plt.xlim(np.min(data[0]), np.max(data[0]))
    plt.ylim(-2, max_peak + (max_peak / 3.5))
    
    ## ZOOMED IN PLOT OF FIRST ANCHOR POINT LOCATION
    ax2 = fig2.add_subplot(234) 
    plt.title("Left Anchor Point", fontsize=8)
    plt.plot(data[0], data[1], color = "black") ## PLOTS SPECTRA
    plt.plot(test1.wavelength, test1.flux, color = "blue", linestyle = "-") ## PLOTS TEST REGION 1 (1315-1325) 
    plt.plot(test2.wavelength, test2.flux, color = "deeppink", linestyle = "-") ## PLOTS TEST REGION 2 (1350-1360)
    plt.plot(anchor_pts_wavelength[0], anchor_pts_flux[0], 'ro')
    ax2.xaxis.set_major_locator(MultipleLocator(100))
    ax2.xaxis.set_minor_locator(MultipleLocator(20))
    ax2.tick_params(labelsize=8)
    ax2.grid(which='minor', linewidth=0.3)
    ax2.grid(which='major', linewidth=1)
    xlabel_2 = "input wavelength: " + str(user_input_wavelength[0])
    plt.xlabel(xlabel_2, fontsize=8)
    plt.ylabel("Flux", fontsize=8) 
    plt.xlim(anchor_pts_wavelength[0]-100, anchor_pts_wavelength[0]+100)
    plt.ylim(min_peak_pt1, max_peak_pt1) #max_peak + (max_peak / 1.5))
    
    ## ZOOMED IN PLOT OF SECOND ANCHOR POINT LOCATION
    ax3 = fig2.add_subplot(235) 
    plt.title("Middle Anchor Point", fontsize=8)
    plt.plot(data[0], data[1], color = "black") ## PLOTS SPECTRA
    plt.plot(test1.wavelength, test1.flux, color = "blue", linestyle = "-") ## PLOTS TEST REGION 1 (1315-1325) 
    plt.plot(test2.wavelength, test2.flux, color = "deeppink", linestyle = "-") ## PLOTS TEST REGION 2 (1350-1360)
    plt.plot(anchor_pts_wavelength[1], anchor_pts_flux[1], 'ro')
    ax3.xaxis.set_major_locator(MultipleLocator(100))
    ax3.xaxis.set_minor_locator(MultipleLocator(20))
    ax3.tick_params(labelsize=8)
    ax3.grid(which='minor', linewidth=0.3)
    ax3.grid(which='major', linewidth=1)
    xlabel_3 = "input wavelength: " + str(user_input_wavelength[1])
    plt.xlabel(xlabel_3, fontsize=8)
    plt.ylabel("Flux", fontsize=8) 
    plt.xlim(anchor_pts_wavelength[1]-100, anchor_pts_wavelength[1]+100)
    plt.ylim(min_peak_pt2, max_peak_pt2) 

    ## ZOOMED IN PLOT OF THIRD ANCHOR POINT LOCATION
    ax4 = fig2.add_subplot(236) 
    plt.title("Right Anchor Point", fontsize=8)
    plt.plot(data[0], data[1], color = "black") ## PLOTS SPECTRA
    plt.plot(test1.wavelength, test1.flux, color = "blue", linestyle = "-") ## PLOTS TEST REGION 1 (1315-1325) 
    plt.plot(test2.wavelength, test2.flux, color = "deeppink", linestyle = "-") ## PLOTS TEST REGION 2 (1350-1360)
    plt.plot(anchor_pts_wavelength[2], anchor_pts_flux[2], 'ro')
    ax4.xaxis.set_major_locator(MultipleLocator(100))
    ax4.xaxis.set_minor_locator(MultipleLocator(20))
    ax4.tick_params(labelsize=8)
    ax4.grid(which='minor', linewidth=0.3)
    ax4.grid(which='major', linewidth=1)   
    xlabel_4 = "input wavelength: " + str(user_input_wavelength[2])
    plt.xlabel(xlabel_4, fontsize=8)
    plt.ylabel("Flux", fontsize=8) 
    plt.xlim(anchor_pts_wavelength[2]-100, anchor_pts_wavelength[2]+100)
    plt.ylim(min_peak_pt3, max_peak_pt3) 

    plt.subplots_adjust(hspace=0.65)

    if FILE == 'null' and FILE2 == 'null': ## SHOWS FIGURE BUT DOES NOT SAVE
        plt.show()
    else: ## SAVES FIGURES
        FILE.savefig(fig)
        FILE2.savefig(fig2)
        plt.close(figure_index)
        plt.close(figure_index+1)


def draw_dynamic(wavelength, wavelength_observed_from, wavelength_observed_to, flux, test1: RangesData, test2: RangesData, max_peak):
    ### ADD DOCUMENTATION - this plots the spectra to decide where to put anchor points -- not intended to be saved
    ## combine this with draw_dynamic_points?
    fig = plt.figure()
    ax = fig.gca()
    plt.plot(wavelength, flux, color = "black")
    plt.plot(test1.wavelength, test1.flux, color = "blue", linestyle = "-")
    plt.plot(test2.wavelength, test2.flux, color = "deeppink", linestyle = "-")
    plt.xlim(wavelength_observed_from, wavelength_observed_to)
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

def draw_original_figure(figure_index: int, original_ranges: RangesData, data: FigureDataOriginal, test1: RangesData, test2: RangesData, wavelength_observed_from, wavelength_observed_to, max_peak, FILE, flags, anchor_pts, z, spectrum_file_name):
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
    subtitle_text = f"z={data.FigureData.z} snr={data.FigureData.snr} snr_mean_in_ehvo={data.FigureData.snr_mean_in_ehvo}"
    plt.title(spectrum_file_name + flags + '\n' + subtitle_text)

    ########

    plt.plot(original_ranges.wavelength, original_ranges.flux, color = main_color, linestyle = "-") ## PLOTS SPECTRA
    plt.plot(original_ranges.wavelength, original_ranges.error, color = "black", linestyle = "-") ## PLOTS ERROR
    plt.plot(original_ranges.wavelength, powerlaw(original_ranges.wavelength, data.bf, data.cf), color = "red", linestyle = "--") #, zorder = 3) ## PLOTS POWER LAW

    plt.plot(data.power_law_data_x, data.power_law_data_y, 'ro') ## PLOTS ANCHOR POINTS
    
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
    
    plt.xlabel("Wavelength[$\AA$]")
    plt.ylabel("Flux[10$^{-17}$ erg/cm$^2$/$\AA$]")
    plt.xlim(wavelength_observed_from, wavelength_observed_to)
    plt.ylim(-2, max_peak + (max_peak / 1.5))

    FILE.savefig(bbox_inches='tight')
    plt.close(figure_index)

def draw_normalized_figure(figure_index: int, original_ranges: RangesData, figure_data: FigureData, flux_normalized, error_normalized,
                            test1: RangesData, test2: RangesData, normalized_flux_test_1, normalized_flux_test_2, wavelength_observed_from, wavelength_observed_to, max_peak, FILE, z, spectrum_file_name):#, min_flux_green_region, min_flux_pink_region, max_flux_green_region, max_flux_pink_region, z, val):
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
    subtitle_text = f"z={figure_data.z} snr={figure_data.snr} snr_mean_in_ehvo={figure_data.snr_mean_in_ehvo}"
    plt.title("$\mathbf{Normalized ~Data ~vs. ~Normalized ~Error}$" + '\n' + spectrum_file_name + '\n' + subtitle_text)

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
    
    plt.xlim(wavelength_observed_from, wavelength_observed_to)
    plt.ylim(0, max_peak + (max_peak / 4))
    
    FILE.savefig(bbox_inches='tight')
    plt.close(figure_index)

def dynamic_find_anchor_points(spectra_data, z, user_anchors:list, user_delta:float, verbose=True):
    """
    Function based on 'define_three_anchor_points'. Defines a user-specified 
        number of anchor points. This function makes use of the function
        'wavelength_flux_error_for_points' to find the closest wavelength bin
        to the user-requested anchor point values.
    Parameters
    ----------
    spectra_data : tuple
        (wavelength, flux, error).
    z : float
        redshift.
    user_anchors : list
        User-defined desired anchor points.
    user_delta : float
        Wavelength range to search for anchor points.
    verbose : bool
        Print found wavelength bins for each user-defined anchor point.
    Returns
    -------
    anchor_pts : arr
        List of PointData objects.
    """

    anchor_pts = []
    
    
    
    for i, point in enumerate(user_anchors):
        
        llim = point - user_delta
        ulim = point + user_delta
        
        spec_point = wavelength_flux_error_for_points(llim, ulim, z, spectra_data)
        
        anchor_pts.append(spec_point)
        
        if verbose:
            print('Matched user requested point', np.round(point, 2))
            print('        with restframe point', np.round(spec_point[0]/(1+z), 2))

    return anchor_pts