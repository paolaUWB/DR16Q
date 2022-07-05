# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 18:45:19 2022

@author: bunge
"""

import numpy as np
import os

parent_data_raw = np.loadtxt(os.getcwd() + "/DR9Q_EHVO/DR9Q_selection_minus17.dat", dtype='str', delimiter='\n').T
directory = os.getcwd() + '/VARIABILITY/DATA_VARIABILITY/'

###Choose which data release the data is from

DR9 = True
DR16 = False

NORM = True

parent_data = []
EHVO_data = []
plate = ''
mjd = ''
fiber = ''
parent_names = []
EHVO_names = []
previous_name = []
objid = []

z = ''
snr = ''


for i in range(len(parent_data_raw)):
    parent_data.append(parent_data_raw[i].split())
    plate = parent_data[i][4].zfill(4)
    mjd = parent_data[i][5]
    fiber = parent_data[i][6].zfill(4)
    parent_names.append('spec-' + plate + '-' + mjd + '-' + fiber)
    objid.append('J' + parent_data[i][0])
    
if DR9:
    cutoff = 4
    if NORM:
        cutoff += 4
    for f in os.listdir(directory):
        print(f)
        for file in os.listdir(directory + f):
            print(file)
            if file[:-cutoff] in parent_names:
                previous_name = file[:-cutoff]
                mask = np.where(previous_name == np.array(parent_names))[0][0]
                z = parent_data[mask][7]
                snr = parent_data[mask][29]
                
            header = ['#name', 'z_pipe', 'snr']
            save = np.array([header, [previous_name, z, snr]])
            np.savetxt(directory + f + '/info.txt', save, fmt = '%s, ')

