# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from utility_functions import read_spectra, read_list_spectra

#print ( os.getcwd()) 

CONFIG_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/overplot.csv"

# directory where normalized data files are
# Directory to download files from: J2318
#Files to download:
NORM_DIREC = os.getcwd() + '/../spec/' #+ "/../" + "/DR16Q_EHVO/NORM_DR16Q_EHVO/"

# creates directory for output files
OUT_DIREC = os.getcwd() + "/output/"

#specfiles = os.getcwd()
specfiles = os.getcwd() + '/../spec/'


#-- saving files
#-- save png
saveformat = 'png'
if saveformat == 'png': 
    pp2 = '.png'

#-- save pdf
if saveformat == 'pdf':
    pp2 = '.pdf'

quasarspecname = 'test'
quasarspec = specfiles
file_list = os.listdir(quasarspec)

# #DR16 spec files, normalized, 2/108

# spec_1 = specfiles + 'spec-0628-52083-0304norm.dr16'
# spec_2 = specfiles + 'spec-1272-52989-0064norm.dr16'


# #Get the data from the files

# data_1 = np.loadtxt(spec_1)
# data_2 = np.loadtxt(spec_2)


# #Define the data columns [wavelength, flux, error]
# #Get the wavelength, flux, and error from each spec file

# wavelength_1 = data_1[:,0]
# wavelength_2 = data_2[:,0]

# flux_1 = data_1[:,1]
# flux_2 = data_2[:,1]

# error_1 = data_1[:,2]
# error_2 = data_2[:,2]


# #plot the spectra on the same plot

# plt.plot(wavelength_1, flux_1, color = 'purple', alpha = 0.8)
# plt.plot(wavelength_2, flux_2, color = 'green', alpha = 0.2)

#inside the for loop, one line of plotting
#plt.plot( wavelength, flux )

norm_spectra_list, redshift_list, calc_snr_list = read_list_spectra(CONFIG_FILE, ["normspec", "redshift", "snr"]) 

for i in range ( len(norm_spectra_list)): 
    norm_spectrum_file_name = norm_spectra_list [i]
    norm_spectra_data = np.loadtxt(NORM_DIREC + norm_spectrum_file_name)
    wavelength, normalized_flux, normalized_error = read_spectra(norm_spectra_data)
    plt.plot( wavelength, normalized_flux )

#outside the for loop
#Label the plot

plt.xlabel ('Wavelength')
plt.ylabel ('Normalized Flux')
plt.title ('Overplot for Average of EHVO')

#See the plot and close the plot after the loop is over

plt.show ()
plt.close ()

#Save the plot on the current working directory, still working this part 

plt.savefig( os.getcwd() + pp2 )


#put file list in a csv file w/ column of list of names that is spec files being pulled
#are these in restframe wavelength
# not in restframe, need redshift included for correction
#line 115, set up like absorption.py
#copy line 126 
#copy line 131 exactly from the for loop
#somthing like 135 as well
#copy 138, uses from like 135 to get the wavelength flux and error using the utility functions program 

