#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 16:18:26 2026

@author: lilianaflores
"""
import sdss_access
import numpy as np
import os
import sys
import astropy.io.fits
import pandas as pd
import matplotlib.pyplot as plt
import shutil
sys.path.insert(0, os.getcwd() + '/../')
from abs_function_module import wavelength_to_velocity


download_dir = os.path.join(os.getcwd(), "downloaded_SDSS_spectra")

os.makedirs(download_dir, exist_ok=True)

from sdss_access import HttpAccess
http_access = HttpAccess(release='DR19', verbose=True)

# set to use remote
http_access.remote()

input_file = os.getcwd()+ '/spec_lists/xray_sample_field.csv'
data_in = pd.read_csv(input_file)

df = pd.DataFrame(data_in)
#print(df.columns)
#Index(['MJD', 'RUN2D', 'CATALOGID', 'SDSS_ID', 'plate_or_fps_field'], dtype='object')

run2d_col = df['RUN2D']
mjd_col = df['MJD']
catalogid_col = df['CATALOGID']
sdssid_col = df['SDSS_ID']
field_col = df['plate_or_fps_field']
z_col = df['Z']

for i in range(len(df)):
#for i in range(0,5):

    run2d = run2d_col[i]
    mjd = mjd_col[i]
    catalogid = catalogid_col[i]
    sdss_id = sdssid_col[i]
    field = field_col[i]
    z = z_col[i]

    path = sdss_access.Path(release='DR19', force_modules=True)
    spec_path = path.full('specFull', catalogid=catalogid, run2d=run2d, mjd=mjd, fieldid=field)
    
    
    filename = os.path.basename(spec_path)
    new_path = os.path.join(download_dir, filename)

    if os.path.exists(new_path):
        print(f"Already exists in flat dir: {filename}")
        
        spec_file = astropy.io.fits.open(new_path)
        spec_dat = np.array(spec_file[1].data)

        #print(spec_file[1].columns)#run to print column names
     
        flux = spec_dat['FLUX']
        loglam = spec_dat['LOGLAM']
        sky = spec_dat['SKY'] #is this the error?

        dumb = np.zeros_like(flux)
        output_data = np.column_stack((10**loglam, flux, dumb))
        
        new_filename = filename[:-5]
        
        np.savetxt(os.getcwd()+'/modified_data/'+f"{new_filename}.txt", output_data, header="#wavelengths flux error", comments='')
        
        '''
        beta = wavelength_to_velocity(z, 10**loglam)
        
        plt.plot(beta,flux, color='k')
        plt.plot(beta,sky, color='grey')
        plt.title(f'mjd={mjd} | catalogid={catalogid} | z = {z}')
        
        plt.ylim(-5,40)
        plt.xlim(-70000,-20000)
        plt.ylabel('Flux 10^-17 ergs/s/cm^2/Angs')
        plt.xlabel('Velocity km/s')
        plt.show()
        plt.close()
        '''
        
        continue
    else:
        # otherwise download
        try:
            http_access.get('specFull',
                            catalogid=catalogid,
                            run2d=run2d,
                            mjd=mjd,
                            fieldid=field)
        
            # move into flat directory
            shutil.move(spec_path, new_path)
        except:
            print('!!!!! Could not find. Download FAILED!!!!!!')

    print(i)


    
    












