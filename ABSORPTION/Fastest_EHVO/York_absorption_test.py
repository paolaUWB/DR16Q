"""
absorption
==========

@author Wendy Garcia Naranjo, Mikel Charles, Nathnael Kahassai, Michael Parker 
based on code prepared by Abdul Khatri and Paola Rodriguez Hidalgo

Short description:
    Creates text file and plot for BI.

Extended description:
    Loops through a list of spectra and uses absorption_parameters_with_plot from basic_absorption_parameters.py to receive: 
    BI_total, vmins, vmaxs, BI_individual, EW_individual, final_depth_individual values and a plot is created to show where 
    CIV, CII, and OI would be *if* the EHVO absorption found was due to SiIV. From those values a plot and text file 
    are created and saved.

Input file:
    This program takes a CSV file with the format ``spectrum_name``, ``z``, ``snr``. In this paticular case it is good_fit.csv. 
    From those spectra names it reads in the wavelength, flux, and error from NORM_DR16Q.
"""

###############################################################################################################################
########################################## IMPORTS ############################################################################
import os
import sys
import numpy as np 
import math
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
sys.path.insert(0, os.path.dirname(os.getcwd()+'/../../')) # changes the directory to the DR16Q --> all paths after this will need to be written as if this was in the top level of the DR16Q
from utility_functions import clear_file, read_list_spectra, read_spectra, append_row_to_csv
from data_types import Range
from abs_function_module import smooth, abs_parameters_plot_optional
from abs_plot_module import draw_abs_figure

###############################################################################################################################
############################## CHANGEABLE VARIABLES ###########################################################################

# data from Verner table
wavelength_CIV_emit1 = 1548.1950
wavelength_CIV_emit2 = 1550.7700
wavelength_SiIV_emit1 = 1393.755
wavelength_SiIV_emit2 = 1402.770
avr_CIV_doublet = 1549.0524 # weighted average
avr_SiIV_doublet = 1396.747 # weighted average

#choose your reference wavelength for wavelength to velocity conversion
ref = 'blue'
#ref = 'red'
#ref = 'avg'

#PREVIOUS trial settings not in use
#files = 'og_data'
#files = 'alt'
#files = 'norm_add_altdiff'
#files = 'norm_sub_altdiff'
#files = 'divided' #not currently usable, finalize divided fits before. Right now divided is not usable revisit idea with Pat.
#files = 'sigma_shiftup'
#files = 'sigma_shiftdown'

#CIV settings
files = 'og_data_all' # running speftra fittings with no alterations to flux
#files = 'CIV_add_0.5err' # using .txt files produced from shifting spectra up/down by 50% of the error and then refitting - [THIS WAS THE METHOD PUBLISHED FOR ERROR CONSTRAINING!]
#files = 'CIV_sub_0.5err'

#files = 'og_data_all_shift_down' # shiffting fitting up - [NOT USED IN PAPER]
#files = 'og_data_all_shift_up' # shifting fitting - [NOT USED IN PAPER]

#SiIV settings
#files = 'SiIV_fit'
#files = 'SiIV_shift_spec_up' # using .txt files produced from shifting spectra up/down by 50% of the error and then refitting - [THIS WAS THE METHOD PUBLISHED FOR ERROR CONSTRAINING!]
#files = 'SiIV_shift_spec_down'
#files = 'SiIV57_upperlim' # only done for 57328 upper limit using .txt file produced from shifting spectra down by 50% of the error and then refitting - [THIS WAS THE METHOD PUBLISHED FOR ERROR CONSTRAINING!]

#files = 'SiIV_fit_up' # shiffting fitting up - [NOT USED IN PAPER]
#files = 'SiIV_fit_down' # shiffting fitting up - [NOT USED IN PAPER]



shift_ogfit = False



#defining the config file
CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/fit_spectra_list_ogdata_V3.csv"
CONFIG_FILE2 = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/fit_spectra_list_alternates_V3.csv" 
CONFIG_FILE3 = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/fit_spectra_list_SiIV.csv" 
CONFIG_FILE4 = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/fit_spectra_list_SiIV57_upperlim.csv" 

# directory of where normalized data files are
# data NOT on github but local computer

NORM_DIREC = os.getcwd() + '/output_files/'

NORM_DIREC2 = os.getcwd() + "/output_files/alternate_normalized_fits/"

NORM_DIREC3 = "/Users/lilianaflores/GitHub/DR16Q/ABSORPTION/York_paper/Finding_Absorption.py_Errors/output_files/CIV_shift_spec/"



# creates directory for output files


OUT_DIREC = os.getcwd() + "/output_files/"



# do you want to use smoothed norm flux/error
# boxcar_size must always be an odd integer
want_to_smooth = 'no' 
boxcar_size = 11

# plot all cases or only those with absorption
# and provide text file for all cases or only those with absorption 
# yes for everything, no for only absorption
all_plot_and_text = 'yes'

# lower limit of absorption width to be flagged 
BALNICITY_INDEX_LIMIT = 2000 #450 # changed to 450 to test for the narrow abs 

# limits on velocity

xlow = -100000
xhigh = -70000

VELOCITY_LIMIT = Range(xhigh, xlow)

'''

# limits on velocity     min,   max
VELOCITY_LIMIT = Range(-70000, -100000.)
'''

# what percentage value you want to go below the continuum
percent = 0.9

# whether you want to output a csv table of your run
want_csv = 'yes'
want_csv = 'no'


want_LaTex = 'no'
errors = True #do you want error values in the LaTex table
# what kind of errors
shift = 'spec'

ref_title = ''

if files == 'og_data':
    CONFIG_FILE = CONFIG_FILE
    NORM_DIREC = NORM_DIREC
    
    # range of spectra you are working with spectra list 'fit_spectra_list_ogdata_V3.csv'
    STARTS_FROM, ENDS_AT = 1, 8 #for norm york spec to compare to alt
    
    out_name = 'og_comp'
    
    ref_wave = wavelength_CIV_emit1
    ref_title = ''

    
    
elif files == 'alt':
    CONFIG_FILE = CONFIG_FILE2
    NORM_DIREC = NORM_DIREC2
    
    # with spectra list 'fit_spectra_list_alternates_V3.csv':
    STARTS_FROM, ENDS_AT = 1, 8 #for alternate york spec (norm2)
    
    out_name = 'alt'
    
    ref_wave = wavelength_CIV_emit1
    ref_title = ''


    
elif files == 'norm_add_altdiff':
    CONFIG_FILE = CONFIG_FILE2
    NORM_DIREC = NORM_DIREC2
    
    # with spectra list 'fit_spectra_list_alternates_V3.csv':
    STARTS_FROM, ENDS_AT = 9, 16 #norm og data + flux difference between og data and alternate normalized data (norm1 + (norm1-norm2))

    out_name = 'added_altdiff'
    
    ref_wave = wavelength_CIV_emit1
    ref_title = ''



elif files == 'norm_sub_altdiff':
    CONFIG_FILE = CONFIG_FILE2
    NORM_DIREC = NORM_DIREC2
    
    # with spectra list 'fit_spectra_list_alternates_V3.csv':
    STARTS_FROM, ENDS_AT = 17, 24 #norm og data + flux difference between og data and alternate normalized data (norm1 + (norm1-norm2))

    out_name = 'sub_altdiff'

    file = 'diffs_df_09_sub_diff.csv'
    
    ref_wave = wavelength_CIV_emit1
    ref_title = ''



elif files == 'divided': #not actively in use
    CONFIG_FILE = CONFIG_FILE2
    NORM_DIREC = NORM_DIREC2
    
    # with spectra list 'fit_spectra_list_alternates_V3.csv':
    STARTS_FROM, ENDS_AT = 25, 32 #norm og data/alternate normalized data (norm1/norm2)
    
    out_name = 'divided'
    
    ref_wave = wavelength_CIV_emit1
    ref_title = ''


    
elif files == 'sigma_shiftup':
    CONFIG_FILE = CONFIG_FILE2
    NORM_DIREC = NORM_DIREC2
    
    # with spectra list 'fit_spectra_list_alternates_V3.csv':
    STARTS_FROM, ENDS_AT = 33, 40 #norm og data plus xsigma to shift up if noise around continuum is 3sigma
    
    out_name = 'sigma_shiftup'
    
    ref_wave = wavelength_CIV_emit1
    ref_title = ''


    
elif files == 'sigma_shiftdown':
    CONFIG_FILE = CONFIG_FILE2
    NORM_DIREC = NORM_DIREC2
    
    # with spectra list 'fit_spectra_list_alternates_V3.csv':
    STARTS_FROM, ENDS_AT = 41, 48 #norm og data plus xsigma to shift up if noise around continuum is 3sigma
    
    out_name = 'sigma_shiftdown'
    
    ref_wave = wavelength_CIV_emit1
    ref_title = ''


    
elif files == 'og_data_all':
    CONFIG_FILE = CONFIG_FILE
    NORM_DIREC = NORM_DIREC

    # range of spectra you are working with spectra list 'fit_spectra_list_ogdata_V3.csv'
    STARTS_FROM, ENDS_AT = 1, 12 #for normal york spec
    
    out_name = 'all_og'
    
    if ref == 'red':
        ref_wave = wavelength_CIV_emit2
        OUT_DIREC = OUT_DIREC + '/other_ref_wavelength/'
        ref_title = '_red'
        
    elif ref == 'blue':
        ref_wave = wavelength_CIV_emit1

        
    elif ref == 'avg':
        ref_wave = avr_CIV_doublet
        OUT_DIREC = OUT_DIREC + '/other_ref_wavelength/'
        ref_title = '_avg'
    
elif files == 'CIV_add_0.5err':
    CONFIG_FILE = CONFIG_FILE2
    NORM_DIREC = NORM_DIREC3
    
    # with spectra list 'fit_spectra_list_alternates_V3.csv':
    STARTS_FROM, ENDS_AT = 49,60  #norm og data plus 0.5 err
    
    out_name = 'CIV_add_0.5err'
    
    ref_wave = wavelength_CIV_emit1
    ref_title = ''
    


    
elif files == 'CIV_sub_0.5err':
    CONFIG_FILE = CONFIG_FILE2
    NORM_DIREC = NORM_DIREC3
    
    # with spectra list 'fit_spectra_list_alternates_V3.csv':
    STARTS_FROM, ENDS_AT = 61,72  #norm og data plus 0.5 err
    
    out_name = 'CIV_sub_0.5err'
    
    ref_wave = wavelength_CIV_emit1
    ref_title = ''


    
elif files == 'og_data_all_shift_down':
    CONFIG_FILE = CONFIG_FILE
    NORM_DIREC = NORM_DIREC

    # range of spectra you are working with spectra list 'fit_spectra_list_ogdata_V3.csv'
    STARTS_FROM, ENDS_AT = 1, 12 #for normal york spec
    
    out_name = 'all_og_shift_down'
    
    ref_wave = wavelength_CIV_emit1
    ref_title = ''
    shift_ogfit = 'down'

    
elif files == 'og_data_all_shift_up':
    CONFIG_FILE = CONFIG_FILE
    NORM_DIREC = NORM_DIREC

    # range of spectra you are working with spectra list 'fit_spectra_list_ogdata_V3.csv'
    STARTS_FROM, ENDS_AT = 1, 12 #for normal york spec
    
    out_name = 'all_og_shift_up'
    
    ref_wave = wavelength_CIV_emit1
    ref_title = ''
    shift_ogfit = 'up'
   

    
elif files == 'SiIV_fit':
    CONFIG_FILE = CONFIG_FILE3
    NORM_DIREC = NORM_DIREC

    # range of spectra you are working with spectra list 'fit_spectra_list_ogdata_V3.csv'
    STARTS_FROM, ENDS_AT = 1, 8 #for SiIV fittings
    
    out_name = 'all_SiIV'
        
    if ref == 'red':
        ref_wave = wavelength_SiIV_emit2
        OUT_DIREC = OUT_DIREC + '/other_ref_wavelength/'
        ref_title = '_red'
        
    elif ref == 'blue':
        ref_wave = wavelength_SiIV_emit1
        
    elif ref == 'avg':
        ref_wave = avr_SiIV_doublet
        OUT_DIREC = OUT_DIREC + '/other_ref_wavelength/'
        ref_title = '_avg'


elif files == 'SiIV_shift_spec_up':
    CONFIG_FILE = CONFIG_FILE3
    NORM_DIREC = NORM_DIREC

    # range of spectra you are working with spectra list 'fit_spectra_list_ogdata_V3.csv'
    STARTS_FROM, ENDS_AT = 9, 16 #for SiIV fittings
    
    out_name = 'SiIV_shift_spec_up'
    
    if ref == 'red':
        ref_wave = wavelength_SiIV_emit2
        OUT_DIREC = OUT_DIREC + '/other_ref_wavelength/'
        ref_title = '_red'
        
    elif ref == 'blue':
        ref_wave = wavelength_SiIV_emit1
        
    elif ref == 'avg':
        ref_wave = avr_SiIV_doublet
        OUT_DIREC = OUT_DIREC + '/other_ref_wavelength/'
        ref_title = '_avg'


elif files == 'SiIV_shift_spec_down':
    CONFIG_FILE = CONFIG_FILE3
    NORM_DIREC = NORM_DIREC

    # range of spectra you are working with spectra list 'fit_spectra_list_ogdata_V3.csv'
    STARTS_FROM, ENDS_AT = 17, 24 #for SiIV fittings
    
    out_name = 'SiIV_shift_spec_down'
    
    ref_wave = wavelength_SiIV_emit1
    ref_title = ''

    
elif files == 'SiIV_fit_up':
    CONFIG_FILE = CONFIG_FILE3
    NORM_DIREC = NORM_DIREC

    # range of spectra you are working with spectra list 'fit_spectra_list_ogdata_V3.csv'
    STARTS_FROM, ENDS_AT = 1, 8 #for SiIV fittings
    
    out_name = 'SiIV_fit_up'
    
    ref_wave = wavelength_SiIV_emit1
    ref_title = ''
    shift_ogfit = 'up'
    


elif files == 'SiIV_fit_down':
    CONFIG_FILE = CONFIG_FILE3
    NORM_DIREC = NORM_DIREC

    # range of spectra you are working with spectra list 'fit_spectra_list_ogdata_V3.csv'
    STARTS_FROM, ENDS_AT = 1, 8 #for SiIV fittings
    
    out_name = 'SiIV_fit_down'
    
    ref_wave = wavelength_SiIV_emit1

    ref_title = ''
    shift_ogfit = 'down'

elif files == 'SiIV57_upperlim':
    CONFIG_FILE = CONFIG_FILE4
    NORM_DIREC = NORM_DIREC

    # range of spectra you are working with spectra list 'fit_spectra_list_ogdata_V3.csv'
    STARTS_FROM, ENDS_AT = 1, 4 #special fitting shifted down for upper limit constraint
    
    out_name = 'SiIV_fit57_upperlim'
    
    
    
    if ref == 'red':
        ref_wave = wavelength_SiIV_emit2
        OUT_DIREC = OUT_DIREC + '/other_ref_wavelength/'
        ref_title = '_red'
        
    elif ref == 'blue':
        ref_wave = wavelength_SiIV_emit1
        
    elif ref == 'avg':
        ref_wave = avr_SiIV_doublet
        OUT_DIREC = OUT_DIREC + '/other_ref_wavelength/'
        ref_title = '_avg'



###############################################################################################################################
######################################## OUTPUT FILES #########################################################################

# set name of output .txt file with absorption values
ABSORPTION_VALUES = OUT_DIREC + "/" + 'BI' + str(BALNICITY_INDEX_LIMIT) + '_P' + str(percent) + out_name  + str(ref_title) + '.txt'

# set name of output pdf with plots 
ABSORPTION_OUTPUT_PLOT_PDF = PdfPages(OUT_DIREC + 'BI' + str(BALNICITY_INDEX_LIMIT) + '_P' + str(percent) + out_name + str(ref_title) + '.pdf') 

ABSORPTION_TABLE = OUT_DIREC + 'absorption_table' + '_P' + str(percent) + out_name + 'BI' + str(BALNICITY_INDEX_LIMIT)  + str(ref_title) +  '.csv'

###############################################################################################################################
######################################### MAIN CODE ###########################################################################

# clear files
if __name__ == "__main__":
    clear_file(ABSORPTION_VALUES)
    if (want_csv == 'yes'):
        clear_file(ABSORPTION_TABLE)

# read list of normalized spectra, zem, and calculated snr from csv file (in this case good_normalization.csv)
# and set variable name to each value
norm_spectra_list, redshift_list, calc_snr_list = read_list_spectra(CONFIG_FILE, ["NORM SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR"]) 
vlast = []
latex_table_rows = [] #Added for LaTex table 
SiIV_latex_table_rows = []
paper_tab_rows = [] #added for shorter paper version of Latex table
SiIV_paper_tab_rows = [] 
# whether abs_count or all_count is used is based on the value of all_plot_and_text
abs_count = 0 # counter for amount of spectra that have absorption when all_plot_and_text = no and for text files
all_count = 0 # counter for all spectra ran when all_plot_and_text = yes

if (want_csv == 'yes'):
    field = ['NORM SPECTRA FILE NAME','BI TOTAL','BI INDIVIDUAL','VMINS', 'VMAXS', 'EW INDIVIDUAL', 'DEPTH']
    append_row_to_csv(ABSORPTION_TABLE, field)

# loops over each spectra from a specified starting and ending point
for spectra_index in range(STARTS_FROM, ENDS_AT + 1):

    # rounding the numbers of the redshift, calculated snr and setting the norm file name to the current file name from the csv
    z = round(redshift_list[spectra_index - 1], 5)
    calc_snr = round(calc_snr_list[spectra_index - 1], 5)
    norm_spectrum_file_name = norm_spectra_list[spectra_index - 1]

    # from the norm spectra name retrieving it's wavelength, normalized flux, and normalized error (in this case from NORM_DRXQ)
    print(str(spectra_index), "current spectra file name:", norm_spectrum_file_name)
    norm_spectra_data = np.loadtxt(NORM_DIREC + norm_spectrum_file_name)

    # setting a variable for each of those values from the spectra
    wavelength, normalized_flux, normalized_error = read_spectra(norm_spectra_data)

    
    # test of shifting just fitting up 
    if shift_ogfit == 'up':
        normalized_flux = normalized_flux + (0.5*normalized_error)
    
    elif shift_ogfit == 'down':
        normalized_flux = normalized_flux - (0.5*normalized_error)
        
    else:
        print('Flux not shifted')
        #continue
        pass
    


    # smoothing the flux and error based on what the user wants (yes or no)
    if want_to_smooth == 'yes':
        normalized_flux = smooth(normalized_flux, boxcar_size)
        normalized_error = smooth(normalized_error, boxcar_size) / math.sqrt(boxcar_size)

    # getting various BI-related values from the absorption_parameters_with_plot function
    BI_total, BI_individual, BI_all, vmins, vmaxs, EW_individual, final_depth_individual, final_depth_all_individual, beta, vminindex_for_range, vmaxindex_for_range, masked_regions_all = abs_parameters_plot_optional(
        z, wavelength, normalized_flux, BALNICITY_INDEX_LIMIT, VELOCITY_LIMIT, ref_wavelength = ref_wave, percent=percent)
                                                                #velocity_limits, ref_wavelength = avr_CIV_doublet, percent = 0.9, plots = 'yes', flag = 'N', norm_spectrum_file_name = 'Manual Mask Regions', manual_depth_masking = False, masks=[]
    
    max_peak = np.max(normalized_flux[vmaxindex_for_range + 1 : vminindex_for_range + 1])

    ############################# putting things into a text file or plot #######################################
    fields = [norm_spectrum_file_name, BI_total, BI_individual, vmins, vmaxs, EW_individual, final_depth_individual]
    
    if (all_plot_and_text == 'yes'): # plot all is yes, graph everything but only text file for when abs is found
        all_count += 1
        if (len(vmaxs) != 0): # text file created only when absorption is found
            abs_count += 1
            text = [f"{str(abs_count)} abs | {str(all_count)} tot",
                    f"{norm_spectrum_file_name}",
                    f"BI ({VELOCITY_LIMIT.start} > v > {VELOCITY_LIMIT.end}): {BI_total}",
                    f"vmins: {vmins}",
                    f"vmaxs: {vmaxs}",
                    f"BI_individual: {BI_individual}",
                    f"EW_individual: {EW_individual}",
                    f"Depth: {final_depth_individual}"]
            vlast.extend(['\n'.join(text), '\n'])
            abs = abs_count
        else: # create graph no matter what
            abs = 'no'
        draw_abs_figure(
            abs, all_count, beta, normalized_flux, normalized_error, ABSORPTION_OUTPUT_PLOT_PDF, norm_spectrum_file_name, z, calc_snr, max_peak, VELOCITY_LIMIT, percent, xlow, xhigh) 
        # whether you want to create a master csv table or not
        if (want_csv == 'yes'):
            append_row_to_csv(ABSORPTION_TABLE, fields)  
    else: # plot all is no and only create text file and graph of cases where absorption is found
        all_count += 1
        if (len(vmaxs) != 0):
            abs_count += 1
            text = [f"{str(abs_count)} abs | {str(all_count)} tot",
                    f"{norm_spectrum_file_name}",
                    f"BI ({VELOCITY_LIMIT.start} > v > {VELOCITY_LIMIT.end}): {BI_total}",
                    f"vmins: {vmins}",
                    f"vmaxs: {vmaxs}",
                    f"BI_individual: {BI_individual}",
                    f"EW_individual: {EW_individual}",
                    f"Depth: {final_depth_individual}"]
            vlast.extend(['\n'.join(text), '\n'])
            abs = abs_count
            draw_abs_figure(
                abs_count, all_count, beta, normalized_flux, normalized_error, ABSORPTION_OUTPUT_PLOT_PDF, norm_spectrum_file_name, z, calc_snr, max_peak, VELOCITY_LIMIT, percent, xlow, xhigh)
        
        # whether you want to create a master csv table or not
        if (want_csv == 'yes'):
            append_row_to_csv(ABSORPTION_TABLE, fields)  
            
            
    #####################################################################################################################
    
    final_depth_all_individual.append(final_depth_individual)

#added to write LaTex table .txt file ---------------------------------------------------------------------------------------------------------------------------------
    #file = 'diffs_df_09_alt.csv'
    #file = 'diffs_df_09_sigma_shift_up.csv'
    #file = 'diffs_df_09_sigma_shift_down.csv'
    #file = 'diffs_df_09_sub_diff.csv'
    #file = 'diffs_df_09_add_diff.csv' 
    #file = 'diffs_df_09_add_err.csv'
    
    if want_LaTex == 'yes':
        if shift == 'spec':
            file = 'diffs_CIV_add_0.5err_rounded_P0.9_shiftspec.csv' # - error
            file2 = 'diffs_CIV_sub_0.5err_rounded_P0.9_shiftspec.csv' # + error
            
            file_P099 = 'diffs_CIV_add_0.5err_rounded_P0.99_shiftspec.csv' # - error
            file2_P099 = 'diffs_CIV_sub_0.5err_rounded_P0.99_shiftspec.csv' # + error
            
            SiIV_file = 'diffs_SiIV_add_0.5err_rounded_P0.9_shiftspec.csv' # - error
            SiIV_file2 = 'diffs_SiIV_sub_0.5err_rounded_P0.9_shiftspec.csv' # + error
            
            SiIV_file_P099 = 'diffs_SiIV_add_0.5err_rounded_P0.99_shiftspec.csv' # - error
            SiIV_file2_P099 = 'diffs_SiIV_sub_0.5err_rounded_P0.99_shiftspec.csv' # + error
            
        elif shift == 'fit': #not in use rn
            file = 'diffs_df_09_add_0.5err_rounded_P0.9_shiftfit.csv' # - error
            file2 = 'diffs_df_09_sub_0.5err_rounded_P0.9_shiftfit.csv' # + error
            
            file_P099 = 'diffs_df_09_add_0.5err_rounded_P0.99_shiftfit.csv' # - error
            file2_P099 = 'diffs_df_09_sub_0.5err_rounded_P0.99_shiftfit.csv' # + error
            
            SiIV_file = 'diffs_df_09_add_0.5err_rounded_P0.9_shiftfit.csv' # - error
            SiIV_file2 = 'diffs_df_09_sub_0.5err_rounded_P0.9_shiftfit.csv' # + error
            
            SiIV_file_P099 = 'diffs_df_09_add_0.5err_rounded_P0.99_shiftfit.csv' # - error
            SiIV_file2_P099 = 'diffs_df_09_sub_0.5err_rounded_P0.99_shiftfit.csv' # + error
            
        else:
            print('Choose error version')
            
            
   
        df = pd.read_csv(OUT_DIREC+ 'diffs/' +file)
        df2 = pd.read_csv(OUT_DIREC+ 'diffs/' +file2)
        
        df_P099 = pd.read_csv(OUT_DIREC+ 'diffs/' +file_P099)
        df2_P099 = pd.read_csv(OUT_DIREC+ 'diffs/' +file2_P099)
        
        #SiIV
        SiIV_df = pd.read_csv(OUT_DIREC+ 'diffs/' + SiIV_file)
        SiIV_df2 = pd.read_csv(OUT_DIREC+ 'diffs/' + SiIV_file2)
        
        SiIV_df_P099 = pd.read_csv(OUT_DIREC+ 'diffs/' + SiIV_file_P099)
        SiIV_df2_P099 = pd.read_csv(OUT_DIREC+ 'diffs/' + SiIV_file2_P099)
        
        
    


        def round_by_error(value, error):
            if error == 0:
                print('Check errors. Error is 0')
                return value
            
            if error >= 100000:
                if str(error)[0] == '1' and str(error)[1] != '0' :
                    whole = 10000
                    
                else:
                    whole = 100000
                    
                    
            elif error >= 10000:
                if str(error)[0] == '1' and str(error)[1] != '0' :
                    whole = 1000
                    
                else:
                    whole = 10000
                    
    
            elif error >= 1000:
                if str(error)[0] == '1' and str(error)[1] != '0' :
                    whole = 100
                    
                else:
                    whole = 1000
                    
    
            elif error >= 100:
                if str(error)[0] == '1' and str(error)[1] != '0' :
                    whole = 10
                    
                else:
                    whole = 100
                    
    
            elif error >= 10:
                if str(error)[0] == '1' and str(error)[1] != '0' :
                    whole = 1
                    
                else:
                    whole = 10
            
    
            else:
                print('Your value is 1-9?')
                print('WARNING IS YOUR ERROR ZERO? CHECK VALIDITY!!')
             
                
            try:
                rounded_value = round(value/whole)*whole #rounds to nearest multiple of another number
                return int(rounded_value)
            except:
                return np.nan
        
        
        
        def get_values_str(string):
            string = string.strip('[]').replace(',', ' ')
            values = string.split()
            return values
        
        def get_max(str_values):
            values = get_values_str(str_values)
                
            numbers = []
            for x in values:
                numbers.append(float(x))
            
            max_value = max(np.abs(numbers))
            return max_value
    
        def get_min(str_values):
            values = get_values_str(str_values)
                
            numbers = []
            for x in values:
                numbers.append(float(x))
            
            min_value = min(np.abs(numbers))
            return min_value
    
        def get_sum(str_values):
            values = get_values_str(str_values)
                
            numbers = []
            for x in values:
                numbers.append(float(x))
            
            sum_values = sum(np.abs(numbers))
            return sum_values
    
        
    
        #CIV errors
        BI_Err_sub = int(df2["BI_TOT_DIFF_09"].iloc[spectra_index-1])
        vmax_Err_sub = int(df2["VMAX_DIFF_09"].iloc[spectra_index-1])
        vmin_Err_sub = int(df2["VMIN_DIFF_09"].iloc[spectra_index-1])
        EW_Err_sub = int(df2["EW_IND_DIFF_09"].iloc[spectra_index-1])
                            
        BI_Err_add = int(df["BI_TOT_DIFF_09"].iloc[spectra_index-1])
        vmax_Err_add = int(df["VMAX_DIFF_09"].iloc[spectra_index-1])
        vmin_Err_add = int(df["VMIN_DIFF_09"].iloc[spectra_index-1])
        EW_Err_add = int(df["EW_IND_DIFF_09"].iloc[spectra_index-1])
        
        
        
        vmax_Err_sub_P099 = int(df2_P099["VMAX_DIFF_09"].iloc[spectra_index-1])
        vmin_Err_sub_P099 = int(df2_P099["VMIN_DIFF_09"].iloc[spectra_index-1])
        EW_Err_sub_P099 = int(df2_P099["EW_IND_DIFF_09"].iloc[spectra_index-1])
                            
        
        vmax_Err_add_P099 = int(df_P099["VMAX_DIFF_09"].iloc[spectra_index-1])
        vmin_Err_add_P099 = int(df_P099["VMIN_DIFF_09"].iloc[spectra_index-1])
        EW_Err_add_P099 = int(df_P099["EW_IND_DIFF_09"].iloc[spectra_index-1])
        
        
        #SiIV errors
        
        '''
        some values for the errors (diffs) are np.nan as absorption was not there for some SiIV
        so this function helps to make values nan and not try to convert to int
        
        '''
        '''
        def nan_or_int(val):
            if math.isnan(val):
                return np.nan
            else:
                return int(val)
            '''
        '''
        SiIV_BI_Err_sub = SiIV_df2["BI_TOT_DIFF_09"].iloc[spectra_index-1]
        SiIV_vmax_Err_sub = int(SiIV_df2["VMAX_DIFF_09"].iloc[spectra_index-1])
        SiIV_vmin_Err_sub = int(SiIV_df2["VMIN_DIFF_09"].iloc[spectra_index-1])
        SiIV_EW_Err_sub = int(SiIV_df2["EW_IND_DIFF_09"].iloc[spectra_index-1])
                            
        SiIV_BI_Err_add = int(SiIV_df["BI_TOT_DIFF_09"].iloc[spectra_index-1])
        SiIV_vmax_Err_add = int(SiIV_df["VMAX_DIFF_09"].iloc[spectra_index-1])
        SiIV_vmin_Err_add = int(SiIV_df["VMIN_DIFF_09"].iloc[spectra_index-1])
        SiIV_EW_Err_add = int(SiIV_df["EW_IND_DIFF_09"].iloc[spectra_index-1])
        
        
        
        SiIV_vmax_Err_sub_P099 = int(SiIV_df2_P099["VMAX_DIFF_09"].iloc[spectra_index-1])
        SiIV_vmin_Err_sub_P099 = int(SiIV_df2_P099["VMIN_DIFF_09"].iloc[spectra_index-1])
        SiIV_EW_Err_sub_P099 = int(SiIV_df2_P099["EW_IND_DIFF_09"].iloc[spectra_index-1])
                            
        
        SiIV_vmax_Err_add_P099 = int(SiIV_df_P099["VMAX_DIFF_09"].iloc[spectra_index-1])
        SiIV_vmin_Err_add_P099 = int(SiIV_df_P099["VMIN_DIFF_09"].iloc[spectra_index-1])
        SiIV_EW_Err_add_P099 = int(SiIV_df_P099["EW_IND_DIFF_09"].iloc[spectra_index-1])
        '''
        
        SiIV_BI_Err_sub = SiIV_df2["BI_TOT_DIFF_09"].iloc[spectra_index-1]
        SiIV_vmax_Err_sub = SiIV_df2["VMAX_DIFF_09"].iloc[spectra_index-1]
        SiIV_vmin_Err_sub = SiIV_df2["VMIN_DIFF_09"].iloc[spectra_index-1]
        SiIV_EW_Err_sub = SiIV_df2["EW_IND_DIFF_09"].iloc[spectra_index-1]
                                
        SiIV_BI_Err_add = SiIV_df["BI_TOT_DIFF_09"].iloc[spectra_index-1]
        SiIV_vmax_Err_add = SiIV_df["VMAX_DIFF_09"].iloc[spectra_index-1]
        SiIV_vmin_Err_add = SiIV_df["VMIN_DIFF_09"].iloc[spectra_index-1]
        SiIV_EW_Err_add = SiIV_df["EW_IND_DIFF_09"].iloc[spectra_index-1]
        
        SiIV_vmax_Err_sub_P099 = SiIV_df2_P099["VMAX_DIFF_09"].iloc[spectra_index-1]
        SiIV_vmin_Err_sub_P099 = SiIV_df2_P099["VMIN_DIFF_09"].iloc[spectra_index-1]
        SiIV_EW_Err_sub_P099 = SiIV_df2_P099["EW_IND_DIFF_09"].iloc[spectra_index-1]
        
        SiIV_vmax_Err_add_P099 = SiIV_df_P099["VMAX_DIFF_09"].iloc[spectra_index-1]
        SiIV_vmin_Err_add_P099 = SiIV_df_P099["VMIN_DIFF_09"].iloc[spectra_index-1]
        SiIV_EW_Err_add_P099 = SiIV_df_P099["EW_IND_DIFF_09"].iloc[spectra_index-1]

        
        obs = norm_spectrum_file_name[6:11]
    
    
        
        
        BI = round_by_error(BI_total, max(BI_Err_add, BI_Err_sub)) 
        vmax = round_by_error(np.min(vmaxs), max(vmax_Err_add, vmax_Err_sub)) 
        vmin = round_by_error(np.max(vmins), max(vmin_Err_add, vmin_Err_sub)) 
        EW = round_by_error(np.sum(EW_individual), max(EW_Err_add, EW_Err_sub)) 
        
        #need to read in csv for these values measured at 0.99 for the CIV
        P099_abs_values = 'absorption_table_P0.99all_ogBI2000.csv' 
    
        df_P099_abs_values = pd.read_csv(OUT_DIREC+P099_abs_values)
        
        
        vmax_P099 = get_max(df_P099_abs_values["VMAXS"].iloc[spectra_index-1])
        vmin_P099 = get_min(df_P099_abs_values["VMINS"].iloc[spectra_index-1])
        EW_P099 = get_sum(df_P099_abs_values["EW INDIVIDUAL"].iloc[spectra_index-1])
        
        if shift == 'fit':
            vmax_Err_sub_P099 = np.nan
            vmin_Err_sub_P099 = np.nan
            EW_Err_sub_P099 = np.nan
        
        vmax_P099 = round_by_error(vmax_P099, max(vmax_Err_add_P099, vmax_Err_sub_P099)) 
        vmin_P099 = round_by_error(vmin_P099, max(vmin_Err_add_P099, vmin_Err_sub_P099)) 
        EW_P099 = round_by_error(EW_P099, max(EW_Err_add_P099, EW_Err_sub_P099)) 
        
        
        if norm_spectrum_file_name[12:14] == 'UZ':
            smoothing = 'N'
            
        else:
            smoothing = 'Y'
            
            
        if norm_spectrum_file_name[16:19] == 'RP2':
            power = 'CS'
            
        else:
            power = 'SMC'
            
        
        if shift == 'fit':
            vmax_Err_sub_P099 = ''
            vmin_Err_sub_P099 = ''
            EW_Err_sub_P099 = ''
        
        #need to read in csv for these values measured at 0.9 and 0.99 for the SiIV
        
        P09_SiIVabs_values = 'absorption_table_P0.9all_SiIVBI2000.csv' 

        
        P099_SiIVabs_values = 'absorption_table_P0.99all_SiIVBI2000.csv' 
        
        P099_SiIVabs_57upperlim_values = 'absorption_table_P0.99SiIV_fit57_upperlimBI2000.csv'
    
        df_P09_SiIVabs_values = pd.read_csv(OUT_DIREC+P09_SiIVabs_values)
        df_P099_SiIVabs_values = pd.read_csv(OUT_DIREC+P099_SiIVabs_values)
        df_P099_SiIVabs_57upperlim_values = pd.read_csv(OUT_DIREC+P099_SiIVabs_57upperlim_values)


        df_P099_SiIVabs_values.iloc[:4] = df_P099_SiIVabs_57upperlim_values.iloc[:4].values




        BI_SiIV = int(df_P09_SiIVabs_values["BI TOTAL"].iloc[spectra_index-1])

        SiIV_vmax_P09 = get_max(df_P09_SiIVabs_values["VMAXS"].iloc[spectra_index-1])
        SiIV_vmin_P09 = get_min(df_P09_SiIVabs_values["VMINS"].iloc[spectra_index-1])
        SiIV_EW_P09 = get_sum(df_P09_SiIVabs_values["EW INDIVIDUAL"].iloc[spectra_index-1])
        
        SiIV_vmax_P099 = get_max(df_P099_SiIVabs_values["VMAXS"].iloc[spectra_index-1])
        SiIV_vmin_P099 = get_min(df_P099_SiIVabs_values["VMINS"].iloc[spectra_index-1])
        SiIV_EW_P099 = get_sum(df_P099_SiIVabs_values["EW INDIVIDUAL"].iloc[spectra_index-1])
        
        if shift == 'fit':
            vmax_Err_sub_P099 = np.nan
            vmin_Err_sub_P099 = np.nan
            EW_Err_sub_P099 = np.nan
            print('The SiIV values with errors from fitting shifted are not formated for this table!!!!!!!!')
        #stop
        
        
        '''
        if BI_SiIV == SiIV_BI_Err_sub:
            BI_SiIV = BI_SiIV
            SiIV_BI_Err_sub = 0
        elif BI_SiIV == SiIV_BI_Err_add:
            BI_SiIV = BI_SiIV
            SiIV_BI_Err_add = 0
        else:
            BI_SiIV = round_by_error(BI_SiIV, max(SiIV_BI_Err_add, SiIV_BI_Err_sub)) 
        
        print()
        print(BI_SiIV)
        print(SiIV_BI_Err_add)
        print(SiIV_BI_Err_sub)
        print()
        SiIV_BI_Err_add = np.where(BI_SiIV == SiIV_BI_Err_add, np.nan, SiIV_BI_Err_add)#not working as intended as SiIV BI err are rounded before being uploaded into this program
        SiIV_BI_Err_sub = np.where(BI_SiIV == SiIV_BI_Err_sub, np.nan, SiIV_BI_Err_sub)
        '''
        BI_SiIV = round_by_error(BI_SiIV, np.nanmax([SiIV_BI_Err_add, SiIV_BI_Err_sub])) 
        
        SiIV_vmax_P09 = round_by_error(SiIV_vmax_P09, np.nanmax([SiIV_vmax_Err_add, SiIV_vmax_Err_sub])) 
        SiIV_vmin_P09 = round_by_error(SiIV_vmin_P09, np.nanmax([SiIV_vmin_Err_add, SiIV_vmin_Err_sub])) 
        SiIV_EW_P09 = round_by_error(SiIV_EW_P09, np.nanmax([SiIV_EW_Err_add, SiIV_EW_Err_sub])) 
        
        SiIV_vmax_P099 = round_by_error(SiIV_vmax_P099, np.nanmax([SiIV_vmax_Err_add_P099, SiIV_vmax_Err_sub_P099])) 
        SiIV_vmin_P099 = round_by_error(SiIV_vmin_P099, np.nanmax([SiIV_vmin_Err_add_P099, SiIV_vmin_Err_sub_P099])) 
        SiIV_EW_P099 = round_by_error(SiIV_EW_P099, np.nanmax([SiIV_EW_Err_add_P099, SiIV_EW_Err_sub_P099])) 
        
        
        if norm_spectrum_file_name[12:14] == 'UZ':
            smoothing = 'N'
            
        else:
            smoothing = 'Y'
            
            
        if norm_spectrum_file_name[16:19] == 'RP2':
            power = 'CS'
            
        else:
            power = 'SMC'
            
        
        if shift == 'fit':
            vmax_Err_sub_P099 = ''
            vmin_Err_sub_P099 = ''
            EW_Err_sub_P099 = ''
        
        #latex_table_row = f"{obs} & {power} & {smoothing} & \\ensuremath{BI^{BI_Err_add}_{BI_Err_sub}} & ${vmax} \\pm {vmax_Err}$ & ${vmin} \\pm {vmin_Err}$ & ${EW} \\pm {EW_Err}$ \\\\ "
                                                            #\ensuremath{12.3^{+0.5}_{-0.2}}
        latex_table_row = (
            f"{obs} & CIV & {smoothing} & {power} &"
            f"\\ensuremath{{{BI}^{{+{BI_Err_sub}}}_{{-{BI_Err_add}}}}} & "
            f"\\ensuremath{{{vmax}^{{-{vmax_Err_sub}}}_{{+{vmax_Err_add}}}}} & "
            f"\\ensuremath{{{vmin}^{{+{vmin_Err_sub}}}_{{-{vmin_Err_add}}}}} & "
            f"\\ensuremath{{{EW}^{{+{EW_Err_sub}}}_{{-{EW_Err_add}}}}} & "
            f"\\ensuremath{{-{vmax_P099}^{{-{vmax_Err_sub_P099}}}_{{+{vmax_Err_add_P099}}}}} & "
            f"\\ensuremath{{-{vmin_P099}^{{+{vmin_Err_sub_P099}}}_{{-{vmin_Err_add_P099}}}}} & "
            f"\\ensuremath{{{EW_P099}^{{+{EW_Err_sub_P099}}}_{{-{EW_Err_add_P099}}}}} \\\\"
            )
        
        SiIV_latex_table_row = (
            f"{obs} & SiIV &  {smoothing} & {power} &"
            f"\\ensuremath{{{BI_SiIV}^{{+{SiIV_BI_Err_sub}}}_{{-{SiIV_BI_Err_add}}}}} & "
            f"\\ensuremath{{{SiIV_vmax_P09}^{{-{SiIV_vmax_Err_sub}}}_{{+{SiIV_vmax_Err_add}}}}} & "
            f"\\ensuremath{{{SiIV_vmin_P09}^{{+{SiIV_vmin_Err_sub}}}_{{-{SiIV_vmin_Err_add}}}}} & "
            f"\\ensuremath{{{SiIV_EW_P09}^{{+{SiIV_EW_Err_sub}}}_{{-{SiIV_EW_Err_add}}}}} & "
            f"\\ensuremath{{-{SiIV_vmax_P099}^{{-{SiIV_vmax_Err_sub_P099}}}_{{+{SiIV_vmax_Err_add_P099}}}}} & "
            f"\\ensuremath{{-{SiIV_vmin_P099}^{{+{SiIV_vmin_Err_sub_P099}}}_{{-{SiIV_vmin_Err_add_P099}}}}} & "
            f"\\ensuremath{{{SiIV_EW_P099}^{{+{SiIV_EW_Err_sub_P099}}}_{{-{SiIV_EW_Err_add_P099}}}}} \\\\"
            )
        
        latex_table_rows.append(latex_table_row)

        SiIV_latex_table_rows.append(SiIV_latex_table_row)
        
        
        
        #df version to filter latexx table for paper info-------
        
        

        #CIV paper table rows
        paper_tab_rows.append({
                "MJD": obs,
                "Trough": "CIV",
                "Smooth": smoothing,
                "Power": power,
                "BI": f"\\ensuremath{{{BI}^{{+{BI_Err_sub}}}_{{-{BI_Err_add}}}}}",
                "vmax_09": f"\\ensuremath{{{vmax}^{{-{vmax_Err_sub}}}_{{+{vmax_Err_add}}}}}",
                "vmin_09": f"\\ensuremath{{{vmin}^{{+{vmin_Err_sub}}}_{{-{vmin_Err_add}}}}}",
                "EW_09": f"\\ensuremath{{{EW}^{{+{EW_Err_sub}}}_{{-{EW_Err_add}}}}}",
                "vmax_099": f"\\ensuremath{{-{vmax_P099}^{{-{vmax_Err_sub_P099}}}_{{+{vmax_Err_add_P099}}}}}",
                "vmin_099": f"\\ensuremath{{-{vmin_P099}^{{+{vmin_Err_sub_P099}}}_{{-{vmin_Err_add_P099}}}}}",
                "EW_099": f"\\ensuremath{{{EW_P099}^{{+{EW_Err_sub_P099}}}_{{-{EW_Err_add_P099}}}}}",
                    })

         #SiIV paper table rows
        SiIV_paper_tab_rows.append({
                 "MJD": obs,
                 "Trough": "SiIV",
                 "Smooth": smoothing,
                 "Power": power,
                 "BI": f"\\ensuremath{{{BI_SiIV}^{{+{SiIV_BI_Err_sub}}}_{{-{SiIV_BI_Err_add}}}}}",
                 "vmax_09": f"\\ensuremath{{-{SiIV_vmax_P09}^{{-{SiIV_vmax_Err_sub}}}_{{+{SiIV_vmax_Err_add}}}}}",
                 "vmin_09": f"\\ensuremath{{-{SiIV_vmin_P09}^{{+{SiIV_vmin_Err_sub}}}_{{-{SiIV_vmin_Err_add}}}}}" ,
                 "EW_09": f"\\ensuremath{{{SiIV_EW_P09}^{{+{SiIV_EW_Err_sub}}}_{{-{SiIV_EW_Err_add}}}}}",
                 "vmax_099": f"\\ensuremath{{-{SiIV_vmax_P099}^{{-{SiIV_vmax_Err_sub_P099}}}_{{+{SiIV_vmax_Err_add_P099}}}}}",
                 "vmin_099": f"\\ensuremath{{-{SiIV_vmin_P099}^{{+{SiIV_vmin_Err_sub_P099}}}_{{-{SiIV_vmin_Err_add_P099}}}}}",
                 "EW_099": f"\\ensuremath{{{SiIV_EW_P099}^{{+{SiIV_EW_Err_sub_P099}}}_{{-{SiIV_EW_Err_add_P099}}}}}",
                     })

if want_LaTex == 'yes':

    paper_tab_df = pd.DataFrame(paper_tab_rows+SiIV_paper_tab_rows)
            
    paper_tab_df_filtered = paper_tab_df[(paper_tab_df['Smooth'] == 'N') & (paper_tab_df['Power'] == 'CS')]
    #print(paper_tab_df_filtered)
    
    paper_tab_df_filtered = paper_tab_df_filtered.drop(columns=['Smooth', 'Power'])
    #print(paper_tab_df_filtered)
            
    paper_latex_rows = []
    for _, row in paper_tab_df_filtered.iterrows(): #iterate through each row. index as _ since not using the info and then row data
        paper_latex_row = " & ".join(row.values) + r" \\" #joining the values of the cells in the rows, in this case strings, with & and ending each line with \\ to build table
        paper_latex_rows.append(paper_latex_row) #append to table
            
        
    
    paper_table_header = '''
    \\begin{tabular}{ccccccccc}
    \\hline
    \\text{\\parbox[t]{1cm}{\\centering MJD}} & \\text{\\parbox[t]{1cm}{\\centering Trough}}  & \\text{\\parbox[t]{1.3cm}{\\centering $BI_{EHVO}$ \\\ (km s$^{-1}$)}} & \\text{\\parbox[t]{1.3cm}{\\centering $v_{max, 0.9}$ \\\ (km s$^{-1}$)}} & \\text{\\parbox[t]{1.3cm}{\\centering $v_{min, 0.9}$ \\\ (km s$^{-1}$)}} & \\text{\\parbox[t]{1.3cm}{\\centering $EW, 0.9$ \\\ (km s$^{-1}$)}} \\\
        & \\text{\\parbox[t]{1.3cm}{\\centering $v_{max, 0.99}$ \\\ (km s$^{-1}$)}} & \\text{\\parbox[t]{1.3cm}{\\centering $v_{min, 0.99}$ \\\ (km s$^{-1}$)}} & \\text{\\parbox[t]{1.3cm}{\\centering $EW, 0.99$ \\\ (km s$^{-1}$)}} \\\\\hline
    '''
    
    paper_table_footer = r"\hline \end{tabular}"
    
    # Combine everything
    paper_latex_table = paper_table_header + "\n" + "\n".join(paper_latex_rows) + "\n" + paper_table_footer
    #print(paper_latex_table)
            
            
            #-------------------------------------------------------
            
    
    
    
    latex_table = '''
    \\begin{tabular}{ccccccccccc}
    \\hline
    \\text{\\parbox[t]{1cm}{\\centering MJD}} & \\text{\\parbox[t]{1cm}{\\centering Trough}} & \\text{\\parbox[t]{1cm}{\\centering Smooth}} & \\text{\\parbox[t]{1cm}{\\centering Power}} & \\text{\\parbox[t]{1.3cm}{\\centering $BI_{EHVO}$ \\\ (km s$^{-1}$)}} & \\text{\\parbox[t]{1.3cm}{\\centering $v_{max, 0.9}$ \\\ (km s$^{-1}$)}} & \\text{\\parbox[t]{1.3cm}{\\centering $v_{min, 0.9}$ \\\ (km s$^{-1}$)}} & \\text{\\parbox[t]{1.3cm}{\\centering $EW, 0.9$ \\\ (km s$^{-1}$)}} \\\
        & \\text{\\parbox[t]{1.3cm}{\\centering $v_{max, 0.99}$ \\\ (km s$^{-1}$)}} & \\text{\\parbox[t]{1.3cm}{\\centering $v_{min, 0.99}$ \\\ (km s$^{-1}$)}} & \\text{\\parbox[t]{1.3cm}{\\centering $EW, 0.99$ \\\ (km s$^{-1}$)}} \\\\\hline
    '''
    
    
    latex_table_end ='''
    \\end{tabular}
    '''
    
    
    full_latex_table = latex_table + ''.join(latex_table_rows) + ''.join(SiIV_latex_table_rows) + latex_table_end

#if want_LaTex == 'yes':

    if shift == 'spec':
        with open(os.path.join(OUT_DIREC, 'abs_LaTex_table_P' + str(percent) + str() + 'shiftspec.txt'), 'w') as f:f.write(full_latex_table)
        with open(os.path.join(OUT_DIREC, 'abs_LaTex_table_P' + str(percent) + str() + 'paper.txt'), 'w') as f:f.write(paper_latex_table)

    elif shift == 'fit':
        with open(os.path.join(OUT_DIREC, 'abs_LaTex_table_P' + str(percent) + str(files) + 'shiftfit.txt'), 'w') as f:f.write(full_latex_table)

else:
    print('No Latex table produced')


#---------------------------------------------------------------------------------------------------------------------------------
        



BI_all= np.array(BI_all)

vmins = np.array(vmins)
vmaxs = np.array(vmaxs)

ABSORPTION_OUTPUT_PLOT_PDF.close()

vmins_final, vmaxs_final = [], []

'''
# creating list of all vmaxs
for loop in range(0, len(vmaxs_all)):
    vmaxs_final.append(str(vmaxs_all[loop])+ ',' )

# creating list of all vins
for loop2 in range(0, len(vmins_all)):
    vmins_final.append(str(vmins_all[loop2])+ ',' ) 
                    
vmaxs_final = np.array(vmaxs_final)
vmins_final = np.array(vmins_final)
'''

np.savetxt(ABSORPTION_VALUES, vlast, fmt='%s')



