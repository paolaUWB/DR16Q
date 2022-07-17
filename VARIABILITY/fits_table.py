import os
from astropy.io import fits
from astropy.table import Table

FITS_DIREC = os.getcwd() + '/../DR16Q_v4.fits'

t = Table.read(FITS_DIREC)

with fits.open(FITS_DIREC) as hdu_list:
    hdu_list.info()
    table_data = hdu_list[1].data
    print('Column names: \n', table_data.names)
    print('\nRow 1: \n', table_data[1])
    print('\nNumber of rows: \n', len(table_data))