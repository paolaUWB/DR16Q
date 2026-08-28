# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 11:33:16 2017

This code plots pretty spectra for posters/presentations/papers/etc.

@author: Paola
"""

import yaml
import numpy as np
import matplotlib.pyplot as plt 
import os
from lines import CIVll, SiIVll, NVll, OVIll, CIVllred, CIVllblue, SiIVllred, SiIVllblue, CII_emitted, OI_emitted, NVllred, NVllblue, Lya, OVIllred, OVIllblue, Lyb
from absorption import plot_absorption_regions
from emission import plot_emission_labels
from matplotlib.backends.backend_pdf import PdfPages

##------ Inputs/Outputs to change
script_dir = os.path.dirname(os.path.abspath(__file__)) ###pathing

with open(os.path.join(script_dir, "config.yaml"), "r") as file: ###yaml code
    config = yaml.safe_load(file)

with open(os.path.join(script_dir, "config.yaml"), "r") as file:
    config_text = file.read()

save_format = 'pdf' # 'pdf' to save as pdf file, 'png' to save as png file

#-- TO CHANGE EVERY TIME:
save_file_name = config["save_file_name"] 
norm_spectra = config["norm_spectra"] ###Spectra file with all the data (3 columns: Observed wavelength, normalized flux, error)

zem = config["zem"] # redshift of norm_spectra

topylim = config["topylim"] ###Max Norm flux on the graph
topemlabel = topylim - 0.03 # where you want to place the ion labels ###Like OI, CII, and CIV

zem_label_x = config["zem_label_x"] # x-coordinate for zem label (in restframe)
zem_label_y = config["zem_label_y"] # y-coordinate for zem label

wavelength_emit1_initial = config["restframe_min"]  # left xlim in restframe
wavelength_emit2_initial = config["restframe_max"]  # right xlim in restframe

vmin = config["vmin"] # make smaller to move right line right (bigger to move right line left)
vmax = config["vmax"] # make bigger to move left line left (smaller to move left line right)

#-- absorption shading: 'yes' to include
CIVabs = config['CIVabs']
NVabs = config['NVabs']
OVIabs = config['OVIabs']
SiIVabs = config['SiIVabs']
Lyaabs = config['Lyaabs']
Lybabs = config['Lybabs']
CIIabs = config['CIIabs']
OIabs = config['OIabs']

CIVem = config['CIVem']
SiIVem = config['SiIVem']
LyaNVem = config['LyaNVem']
OIem = config['OIem']
CIIem = config['CIIem']
OVIem = config["OVIem"] # OVI emission label: 'true' to include

n = config["smooth_box"] # smooth box car

#------ saving files
#-- save png
if save_format == 'png': 
    pp2 = save_file_name + '.png'

#-- save pdf
if save_format == 'pdf':
    pp2 = save_file_name + '.pdf'


#------ location (wavelength) for doublets
c = 300000. # speed of light

##----- functions
def smooth(norm_flux, box_pts):   #<--- reduces noise ###norm flux = og flux val, box_pts = smooth box car -> n=3
    box = np.ones(box_pts)/box_pts ###np.ones makes array with ones, box_pts = 3 so np.ones(box_pts) is [1, 1, 1] then [1, 1, 1]/3 = [0.33, 0.33, 0.33]
    y_smooth = np.convolve(norm_flux, box, mode='same') ###np.convolve does moving average smoothing, np.convolve(data, smoothing filter, mode), mode='same' includes data in box while mode='full' does all data included outside of box bounds  
    return y_smooth

def zabs(v): ###calculates absorption redshift of absorp line corresponding to the outflow velocity v
    beta = -v/c ###v=-37190.96072 so -v=37190.96072 and c=300000 so beta=0.123969 which mean gas is moving 0.123969 times speed of light
    Rc = np.sqrt((beta+1.)/(1.-beta)) ###Rc is relativistic Doppler factor which is 1.132706
    za = ((1.+zem)/Rc)-1. ###za is absorption redshift of 1.57537, Quasar emission lines are at z=zem, Absorbing gas produces lines at z=a
    return za

spectra_count = 1 ###Number of spectra being processed, right now just 1 graph so = 1
print('spec_name' + ": " + norm_spectra)
print(spectra_count)
print('zem=',zem) ###All 3 will print spec_name: spec-10431-58137-0135norm.dr16, 1, and zem= 1.91713968767695 <- prints in terminal when code is ran
    
data='' ###sets data as an empty string but gets overwritten in next line
data = np.loadtxt(os.path.join(script_dir, norm_spectra)) ###specdirec = /home/jmaharaj/, norm_spectra = "spec-10431-58137-0135norm.dr16", so full path becomes /home/jmaharaj/spec-10431-58137-0135norm.dr16
                                          ###np.loadtxt() reads txt files numbers and stores it as numpy array -> 1, 2, 3 becomes [1, 2, 3]
                                          ###data becomes array so data = [1, 2, 3]

wavelength_observe1 = (zem+1.)*wavelength_emit1_initial #Shift start wavelength into frame |<--This makes our wavelength range for
wavelength_observe2 = (zem+1.)*wavelength_emit2_initial #Shift end wavelength into frame   |     the region we want to look at
                                                        ###We only want to look at the spectra from wave_obs1 to wave_obs2 like 2900 to 4600
                                                        ###2900 to 4600 is 1000 to 1600 in restframe like seen on the graph

wavelength_lower_limit = np.where(data[:,0] > wavelength_observe1) ###in data theres column 0, 1, and 2 -> column 0 = wavelength vals
                                                                   ###data[:,0] pulls all wavelength vals
                                                                   ###wave_obs1=2917 so wave_low_lim is positions where wavelength is > 2917
wavelength_upper_limit = np.where(data[:,0] < wavelength_observe2) ###same thing but with less than
                                                                   ###so wave_up_lim is positions where wavelength is < 4667

minwave= np.min(wavelength_lower_limit[0]) ###finds the smallest index from wave_low_lim which is in row 0 -> wavelength is 3594
maxwave= np.max(wavelength_upper_limit[0]) ###finds the biggest index from wave_up_lim which is in row 1134 -> wavelength is 4666
                                          
wavelength = data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit [0]),0] #Get wavelengths in our data set that fall into our region of study
                              ###sets wavelength equal to the range of data we need which is from row 0 to 1134 or wavelength vals 3594 to 4666
actual_wavelength= wavelength ###renames variable

flux = data[np.min(wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),1] #Get flux values in our region
###does same as last two lines but for flux over the range we are looking at, uses 1 for the 2nd column since 2 is flux (0 is wave in 1st column)
actual_flux = flux ###renames

error= data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2] #Get error values in our region
###gives error in same format over range of what we are observing and uses 2 for column 3 which is error
messed_up_error = np.where ( data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2]  >3) #Get inexes of points with error > 3
###messed up error makes any val in error array ignored if it is greater than 3, this seems to not be used as mentioned later
plerror=error ###renames
   
wavelength_emit = wavelength/(zem+1) #Unshift(?) the wavelength, back to a rest frame ###observed to restframe wavelength

normflux=flux ###renamed
error_normflux=error ###renamed

aa= np.where(error_normflux > 2) ###Finds where in the array error is greater than 2
if len(aa) > 0: ###if there are any error values
    error_normflux[aa]=0 ###replace the error vals found with val of 0

zabs_min=vmin ###renamed, vmin=outlfow vel
zabs_max=vmax ###renamed

for i in range(0,len(vmin)): ###saying to repeat code for all vals of vmin (only 1 val)
    zabs_min[i]=zabs(vmin[i]) ###velocity limits converted to absorption redshifts, vmin to zabs_min, so absorption features are drawn on graph
    zabs_max[i]=zabs(vmax[i])
    

# ------ 
    #SOMETIMES, THERE ARE PIXEL PROBLEMS, AND WE MIGHT GET AN ERROR OF 30 IN FLUX. TO AVOID THAT, WE HAVE DONE THIS. MESSED UP ERROR IS 
    ####       DEFINED ABOVE.
    #if len (messed_up_error[0]) > 0:######################################################original
        #plerror[messed_up_error[0]]=0####################################################
	#flux[messed_up_error[0]]=0

fig, ay1 = plt.subplots() ###makes figure and axes 

# ay1 = fig.add_subplot(1, 1, 1)

# plt.title(spectrum)
ay1.set_xlabel(r"Observed Wavelength [$\rm \AA$]") ###x axis label
ay1.set_ylabel(r"Normalized Flux") ###y axis label
     
ay1.plot (wavelength, smooth(normflux,n),'k-') ###plots (x,y) which is (wave,flux) and in smooth(normflux,n), n=3 like before, style is k- (normal)
ay1.plot (wavelength, error_normflux,'k--') ###plots the error line, style is k-- which is dashed

plt.plot([wavelength_observe1,wavelength_observe2],[1,1],'r--') ###draws normalized flux line at y=1
colors = config['colors']



plot_absorption_regions(zabs_min, zabs_max, vmin,
                        CIVll, NVll, NVllblue, NVllred,
                        OVIll, OVIllblue, OVIllred,
                        SiIVll, SiIVllblue, SiIVllred,
                        Lya, Lyb, CII_emitted, OI_emitted,
                        CIVabs, NVabs, OVIabs, SiIVabs,
                        Lyaabs, Lybabs, CIIabs, OIabs,
                        colors)



#matplotlib.rcParams['font.sans-serif'] = ['Source Han Sans TW', 'sans-serif']

plot_emission_labels(
    zem, topemlabel,
    CIVem, SiIVem, LyaNVem, OIem, CIIem, OVIem,
    CIVll, SiIVll, NVll, OVIll,
    CIVllred, CIVllblue,
    SiIVllred, SiIVllblue,
    CII_emitted, OI_emitted,
    NVllred, NVllblue,
    Lya,
    OVIllred, OVIllblue,
    Lyb
)

coem=zem+1. ###stores redshift + 1 for later
plt.xlim(wavelength_observe1,wavelength_observe2) ###sets the x axis limits from 2917 to 4667

ay2 = ay1.twiny() ###creates 2nd x axis for the restframe wavelength at the top

ay2.plot(wavelength/coem,1000*smooth(flux,n)+1000.) ###converts observed to restframe for the top x axis
ay2.set_xlabel(r"Restframe Wavelength [$\rm \AA$]") ###adds label for restframe wavelegth in top x axis

ay1.xaxis.set_label_coords(0.48, -0.08) ###location of bottom x axis label
ay2.xaxis.set_label_coords(0.48, 1.11) ###location of top x axis label
ay2.xaxis.set_major_locator(plt.MaxNLocator(5)) ###limits top x axis to 5 tick marks

zem_plot = "z = " + str(zem) ###converts z val into text for redshift label
plt.text(zem_label_x, zem_label_y, zem_plot, bbox=dict(facecolor='none', edgecolor='black', pad=7.0)) ###places box around redshoft label
plt.xlim(wavelength_observe1/coem,wavelength_observe2/coem) ###sets x axis limits in range of the restframe from 1000 to 1600
plt.ylim(0,topylim) ###sets y axis range of normalized flux from 0 to 2.7

fig.tight_layout() ###adjusts spacing and fit to ensure everything is properly visible

outdir = os.path.join(script_dir, "OUTPUT_FILES") ###creates path for where output files will be saved
os.makedirs(outdir, exist_ok=True) ###creates output folder if not already exisitng

pdf_metadata = {                          ### meta data code
    'Title': save_file_name,
    'Subject': 'Normalized Spectrum Plot',
    'Keywords': config_text,
    'Creator': 'plotindividualspectrum.py'
}

pdf_path = os.path.join(outdir, pp2) ### Save the spectrum plot and configuration as a multi-page PDF

with PdfPages(pdf_path) as pdf:

    pdf.savefig(fig, dpi=100) ### Page 1: spectrum plot

    config_fig = plt.figure(figsize=(8.5, 11)) ### Page 2: configuration parameters
    
    config_fig.text(
        0.05,
        0.95,
        'Configuration Parameters',
        fontsize=16,
        fontweight='bold',
        verticalalignment='top'
    )

    config_fig.text(
        0.05,
        0.91,
        config_text,
        fontsize=10,
        fontfamily='monospace',
        verticalalignment='top'
    )

    pdf.savefig(config_fig, dpi=100)

    plt.close(config_fig)