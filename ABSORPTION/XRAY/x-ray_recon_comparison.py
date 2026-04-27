#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  6 15:25:15 2026

@author: lilianaflores and elijahfacklam
"""

import os
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
sys.path.insert(0, os.getcwd()+'/../')
from abs_function_module import wavelength_to_velocity
sys.path.insert(0, os.getcwd()+'/../../')
from utility_functions import read_list_spectra


######################################## PATH FINDING ########################################

#defining the config file
CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()+"/spec_lists/xray_list_SNR10_z1.9.csv"
#SDSS_CSV = sys.argv[2] if len(sys.argv) > 2 else os.getcwd()+"/ordered_sdss_names.csv"
SDSS_CSV = CONFIG_FILE

# range of spectra you are working with from the good_fit.csv file
STARTS_FROM, ENDS_AT = 1, 250 # eventually 1 -> 250

norm_spectra_list, redshift_list, calc_snr_list = read_list_spectra(CONFIG_FILE, ["NORM SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR"]) 
sdss_spectra_list, sdss_redshift_list, sdss_calc_snr_list = read_list_spectra(SDSS_CSV, ["SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR"]) 

# Defines an output folder and makes it if it does not exist
output_folder = os.path.join(os.getcwd(), "OUTPUT_FILES_PDF")
os.makedirs(output_folder, exist_ok=True)

########################################## CHANGEABLE VARIABLES ##########################################

output_morphed_pdf = os.path.join(output_folder, "spectra_recon_morphed.pdf")
output_og_pdf = os.path.join(output_folder, "spectra_og.pdf")
output_sdss_pdf = os.path.join(output_folder, "spectra_full_og.pdf")

# Defining x-limits
xlims = -70000, 0

# Modes can be, "all", "morphed", or "og"
MODE = "morphed"
MODE = "og"
MODE = "all" # includes SDSS downloaded all epoch and the Hiremath recon and og spectra

# Do you want to produce plots and pdfs for all three modes ^^
run_all_modes = True


# Do you want to display the plots as they are made in the pdf?
show_plot = True
show_plot = False

# Do you want to scale the SDSS spectra to match the Hiremath spectra? - Still need to determine why they have different scales. Maybe Hiremath has corrected for galactic extinction?
scale_SDSS = True
#scale_SDSS = False

########################################## FUNCTIONS ##########################################

def plot_func(files, file_types, xlims=None, ylims=None, redshift=0, scale_SDSS=True, mode='all'):

    # establish plot
    fig, ax = plt.subplots()
    
    for i in range(len(files)):
        file = files[i]
        file_type = file_types[i]
        
        if file_type == 'Hiremath':
            
            # -------- OG DATA --------
            data = fits.open(file)
            data_tab = data[1].data
        
            wavelength = data_tab['wave']
            flux = data_tab['flux']
            morph = data_tab['morph']
            recon = data_tab['recon']
            error = data_tab['noise']
        
            data.close()
            
            # -------- VELOCITY --------
            beta  = wavelength_to_velocity(0, wavelength)
            
        
        elif file_type == 'SDSS_full':
            
            # -------- SDSS DATA --------
            sdss = fits.open(file)
            sdss_tab = sdss[1].data
            
            """
            SDSS Column Names
            ColDefs(
                name = 'FLUX'; format = 'E'; unit = '10^-17 ergs/s/cm^2/Angs'
                name = 'LOGLAM'; format = 'E'; unit = 'log10(Angs)'
                name = 'IVAR'; format = 'E'
                name = 'AND_MASK'; format = 'J'
                name = 'OR_MASK'; format = 'J'
                name = 'WDISP'; format = 'E'; unit = 'Pixels'
                name = 'SKY'; format = 'E'; unit = '10^-17 ergs/s/cm^2/Angs'
                name = 'MODEL'; format = 'E'
                name = 'WRESL'; format = 'E'; unit = 'Angs'
            """
        
            sdss_wavelength = 10**sdss_tab['LOGLAM']
            sdss_flux = sdss_tab['FLUX']
        
            sdss.close()
    
            # -------- VELOCITY --------
            beta2 = wavelength_to_velocity(redshift, sdss_wavelength)
    
    
    if scale_SDSS == True:
        try: 
            og = flux * morph
            # --- SCALE SDSS TO MATCH OG ---
            #np.interp(x_new, x_known, y_known)
            og_i = np.interp(beta2, beta, og)
            
            scale = np.nanmedian(og) / np.nanmedian(sdss_flux)
            sdss_flux = sdss_flux * scale
            
            
            
            xmin, xmax = -45000, -2000
           
            mask_sdss = (beta2 >= xmin) & (beta2 <= xmax)
            mask_og   = (beta  >= xmin) & (beta  <= xmax)
    
            y_values_sdss = []
            y_values_og = []
    
    
            if np.any(mask_sdss):
                y_values_sdss.append(sdss_flux[mask_sdss])
    
            if np.any(mask_og):
                y_values_og.append((og_i)[mask_sdss])
    
            avg_diff = np.nanmean((np.array(y_values_sdss)-np.array(y_values_og)))
            
            sdss_flux = sdss_flux-avg_diff
        except:
            pass
    else:
        pass
        

    # -------- PLOT --------
    if mode.lower() == 'all':
        ax.plot(beta2, sdss_flux, color='red', label='SDSS Full allepoch')
        ax.plot(beta, flux*morph, color='blue', label='Original (Hiremath+2025)')
        ax.plot(beta, recon*morph, color='green', label='Recon (Hiremath+2025)')
        
        # -------- Y-LIMITS --------
        #play around with the y-limits, possibly take the median value within a range
        if xlims is not None:
            xmin, xmax = -63000, 0

            mask_sdss = (beta2 >= xmin) & (beta2 <= xmax)
            mask_og   = (beta  >= xmin) & (beta  <= xmax)

            y_values = []

            if np.any(mask_sdss):
                y_values.append(sdss_flux[mask_sdss])

            if np.any(mask_og):
                y_values.append((flux*morph)[mask_og])

            if len(y_values) > 0:
                y_all = np.concatenate(y_values)
                ymin = np.min(y_all)
                ymax = np.max(y_all)

                ax.set_ylim(ymin * 0.95, ymax * 1.05)
        
    elif mode.lower() == 'og':
        ax.plot(beta, flux*morph, color='blue', label='Original (Hiremath+2025)')
        ax.plot(beta, recon*morph, color='red', label='Recon (Hiremath+2025)')
        ax.plot(beta, error*morph, color='grey', label='Error (Hiremath+2025)')

        # -------- Y-LIMITS --------      
        xmin, xmax = -63000, 0

        mask_og   = (beta  >= xmin) & (beta  <= xmax)

        y_values = []

        y_values.append((flux*morph)[mask_og])
        y_values.append((error*morph)[mask_og])


        ymin = np.min(y_values)
        ymax = np.max(y_values)

        ax.set_ylim(ymin * 1, ymax * 1.05)
        
    elif mode.lower() == 'morphed':
        ax.plot(beta, flux/recon, color='blue', label='Original (Hiremath+2025)')
        ax.plot(beta, error/recon, color='grey', label='Error (Hiremath+2025)')
        ax.plot(beta, np.ones_like(beta), color='k', linestyle='--')
        
        # -------- Y-LIMITS --------      
        xmin, xmax = -63000, 0

        mask_og   = (beta  >= xmin) & (beta  <= xmax)

        y_values = []

        y_values.append((flux/recon)[mask_og])
        y_values.append((error/recon)[mask_og])

        ymin = np.min(y_values)
        ymax = np.max(y_values)

        ax.set_ylim(ymin * 0.95, ymax * 1.5)
        
    ax.set_ylabel('Flux')
    ax.set_xlabel('Velocity km/s')
    ax.legend(loc='upper right')
    ax.set_title(os.path.basename(file))


    if xlims is not None:
        ax.set_xlim(xlims)

    return fig



def Spectra_Comparison_Generate_PDF(output_path, xlims, scale_SDSS, show_plot = False, mode='all'):
    with PdfPages(output_path) as pdf:
        for spectra_index in range(STARTS_FROM, ENDS_AT + 1):
            #closes the global storage of figures
            plt.close('all')
    
            plt.figure()
            
            #keeping track of how many spectra do not have data points within the velocity limits
            norm_spectrum_file_name = norm_spectra_list[spectra_index - 1]
            
            #gets the cwd and and names each files corectly
            recon_file = os.getcwd() + "/Recons_Hiremath2025/" + str(norm_spectrum_file_name)
            
            sdss_spectrum_file_name = sdss_spectra_list[spectra_index - 1]
            sdss_file = os.path.join(os.getcwd(), "Downloading_SDSS_specs/downloaded_SDSS_spectra", str(sdss_spectrum_file_name))
                
                
            redshift = redshift_list[spectra_index - 1]
            fig = plot_func([recon_file, sdss_file], ['Hiremath', 'SDSS_full'], xlims, redshift=redshift, scale_SDSS=scale_SDSS, mode=mode)
            
            #Saves the current figure into a page of the pdf
            pdf.savefig(fig)
            if show_plot == True:
                plt.show(fig)
            plt.close(fig)
    
if run_all_modes==False:
    Spectra_Comparison_Generate_PDF(output_sdss_pdf, xlims, scale_SDSS ,show_plot, MODE)

elif run_all_modes==True:
    modes = ["all", "morphed", "og"]
    outnames = [output_sdss_pdf, output_morphed_pdf, output_og_pdf]
    for i in range(len(modes)):
        Spectra_Comparison_Generate_PDF(outnames[i], xlims, scale_SDSS ,show_plot, modes[i])
