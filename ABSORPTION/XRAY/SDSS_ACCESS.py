#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 13:28:27 2026

@author: lilianaflores


NOTE: This program only works if you have followed steps to run SDSS Access.
You should not need to run this program as it has already produced 'xray_sample_field.csv'
which is in spec_lists folder in this directory.
"""

import sdss_access
import numpy as np
import os
import sys
import astropy.io.fits
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, os.getcwd()+'/../')
from abs_function_module import wavelength_to_velocity

sdss_path = sdss_access.Path(release='dr19', verbose=True)
path = sdss_access.Path(release='dr19')
access = sdss_access.Access(release='dr19', verbose=True)


allspec_file = sdss_path.full('allspec', vers='1.0.1', release='dr19')

if not sdss_path.exists('',full=allspec_file):
    # if the file does not exist locally, this code will download the data.
    access.remote()
    access.add('allspec', vers='1.0.1', release='dr19')
    access.set_stream()
    access.commit()



allspec_hdus = astropy.io.fits.open(allspec_file)
data = allspec_hdus[1].data

mjd_SDSS = data["mjd"]
run2d_SDSS = data["run2d"]
catalogid_SDSS = data["catalogid"]
sdss_id_SDSS = data["sdss_id"]

dic_df = {'mjd':mjd_SDSS,
          'run2d':run2d_SDSS,
          'catalogid':catalogid_SDSS,
          'sdss_id':sdss_id_SDSS}

df_sdss = pd.DataFrame(dic_df, index=None)
allspec_hdus.close()


def find_field(run2d, mjd, catalogid, sdss_id):
    mask = (
        (df_sdss['run2d'] == run2d) &
        (df_sdss['mjd'] == mjd) &
        (df_sdss['catalogid'] == catalogid) &
        (df_sdss['sdss_id'] == sdss_id)
    )


    matched = data[mask]
    field = matched['plate_or_fps_field']
    return field

input_file = os.getcwd()+'/spec_lists/xray_SDSSCross_SNR10_z1.9.csv'
data_in = pd.read_csv(input_file)

run2d_col = data_in['RUN2D']
mjd_col = data_in['MJD']
catalogid_col = data_in['CATALOGID']
sdssid_col = data_in['SDSS_ID']


field_list = []

for i in range(len(data_in)):
    
    run2d = run2d_col[i]
    mjd = mjd_col[i]
    catalogid = catalogid_col[i]
    sdss_id = sdssid_col[i]
    
    field = find_field(run2d, mjd, catalogid, sdss_id)
    
    if len(field) == 0:
        field_list.append(None)
        print('Warning. No plate or field number.')
    else:
        field_list.append(field[0])
    
    print(field)  

data_in['plate_or_fps_field'] = field_list

output_dir = os.path.join(os.getcwd(), 'spec_lists')

# create the directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

output_file = os.path.join(output_dir, 'xray_sample_field.csv')

data_in.to_csv(output_file, index=False)














