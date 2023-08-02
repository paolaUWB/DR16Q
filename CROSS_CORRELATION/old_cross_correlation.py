# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 12:13:58 2022

@author: Dakota Bunger
"""

from astropy.io import fits

file = "C:/Users/bunge/downloads/dr14q_spec_prop.fits"

hdu = fits.open(file)
data = hdu[1].data

sdss_name = data['SDSS_NAME']
plate = data['PLATE   ']
mjd = data['MJD     ']
fiber = data['FIBERID ']
log_mbh = data['LOG_MBH ']
log_mbh_err = data['LOG_MBH_ERR']
q_mbh = data['QUALITY_MBH']
log_lbol = data['LOG_LBOL']
q_lbol = data['QUALITY_LBOL']
log_redd = data['LOG_REDD']
q_redd = data['QUALITY_REDD']
bi_civ = data['BI_CIV  ']
bi_civ_err = data['ERR_BI_CIV']
bal_flag = data['BAL_FLAG']
hdu.close()

mbh = 10**log_mbh
mbh_err = 10**log_mbh_err
lbol = 10**log_lbol
redd = 10**log_redd