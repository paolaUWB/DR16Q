
import os
import sys
from matplotlib.pyplot import table
import numpy as np
from astropy.io import fits
from astropy.table import Table
import shutil
# from VARIABILITY.EHVO_files import FITS_DUPLICATES_FILE
sys.path.insert(0, os.getcwd() + '/../' + 'DR16Q') # changes the directory to the DR16Q --> all paths after this will need to be written as if this was in the top level of the DR16Q
from utility_functions import print_to_file, read_file, append_row_to_csv

##------ fits table directory
FITS_DIREC = os.getcwd() + '/../../DR16Q_v4.fits'
SPEC_DIREC = os.getcwd() + '/DR16Q_EHVO/DR16Q_EHVO_FILES/'
EHVO_DR16_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + "/../DR16Q_EHVO/DR16_EHVO_sorted_norm.csv"

FITS_DUPLICATES_FILE = os.getcwd() + '/fits_duplicates_no-1.txt'

zem_DR16, snr_DR16, EHVO_spectra_list_DR16 = read_file(EHVO_DR16_FILE)

#-- opens fits table and reads data
hdu = fits.open(FITS_DIREC)
data = hdu[1].data

#-- get column data from fits table
fits_PLATE = data['PLATE'].astype(str)
fits_MJD = data['MJD'].astype(str)
fits_FIBER = data['FIBERID'].astype(str)
fits_OBJID = data['OBJID']
fits_duplicate_PLATE = data['PLATE_DUPLICATE'].astype(str)
fits_duplicate_MJD = data['MJD_DUPLICATE'].astype(str)
fits_duplicate_FIBER = data['FIBERID_DUPLICATE'].astype(str)
hdu.close()

#%%
#-- initializing variables to store spectra names & duplicates
fits_names = []
fits_duplicates = []
# duplicates = [] # is this necessary? (since it is also initialized below in for loop)

# count = 0
# count2 = 0
##------
print(len(fits_PLATE))
for i in range(len(fits_PLATE)): #750414 spectra total in DR16Q fits file
    duplicates = []
    fits_names.append('spec-' + fits_PLATE[i].zfill(4) + '-' + fits_MJD[i] + '-' + fits_FIBER[i].zfill(4))
    
    if (i % 50000) == 0:
        print('You have loaded ', i, ' spectra!')

    if np.any(np.array(fits_duplicate_PLATE[i]) != '-1'): # is this necessary? 
        for j in range(74): # 74 *possible* observations per object
            if fits_duplicate_PLATE[i][j] != '-1':
                duplicates.append('spec-' + fits_duplicate_PLATE[i][j].zfill(4) + '-' + fits_duplicate_MJD[i][j].zfill(4) + '-' + fits_duplicate_FIBER[i][j].zfill(4))
    
    if duplicates != []:
        fits_duplicates.append(duplicates)

    else:
        fits_duplicates.append('-1')
        
    # fields = [fits_names[i], fits_duplicates[i]]
    print(fits_duplicates[i])
    if fits_duplicate_PLATE[i][0] != '-1':
        fields = [fits_names[i]]
    
        for k in range(len(duplicates)): 
            fields.append(duplicates[k])
    
    for m in range(len(fields)): 
        
        print(fields)
        append_row_to_csv(os.getcwd() + '/fits_duplicates.csv', fields)
    #     count += 1
    #     # print_to_file(fields, os.getcwd() + '/fits_duplicates5.txt')
    
    # if fits_duplicate_PLATE[i][0] == '-1':
    #     count2 += 1
    
    
    # print_to_file(fields, os.getcwd() + '/fits_duplicates4.txt')
    # print_to_file(fits_names[i], os.getcwd() + '/fits_names2.txt')

filenames = []
foldernames = []
index = []

# for j in range(len(EHVO_spectra_list_DR16)): # i starts at 0
#     EHVO_DR16_minus_dered = EHVO_spectra_list_DR16[j][:-11] # -11 for dr16, -10 for dr9
#     try: 
#         loc = np.where(EHVO_DR16_minus_dered == )


for file in os.listdir(SPEC_DIREC):
    loc = 0
    if file.startswith('spec'):
        filenames.append(file)
    if file.endswith('dered.dr16'):
        try:
            loc = np.where(file[:-10] == np.array(fits_names))[0][0]########### not getting here ever
            print('loc 1: ', loc)
        except IndexError:
            pass
        if loc == 0:
            for i in range(len(fits_duplicates)):
                if type(fits_duplicates[i]) == list:
                    for j in range(len(fits_duplicates[i])):
                        if file[:-10] == fits_duplicates[i][j]:
                            loc = i###### not getting here ever
                            print('loc 2: ', loc)
    if loc != 0:
        index.append(loc)

    if file.startswith('spec'):
        foldernames.append('J' + fits_OBJID[loc])


og_names = np.array(fits_names)[index]
data_dir = os.getcwd() + '/DR16Q_EHVO/DR16Q_EHVO_FILES/' #possibly: /../ + /DATA/DR9Q_SNR10 #'C:/Users/Dakota/Documents/GitHub/DR16Q/DATA/DR9Q_SNR10'

for i in range(len(filenames)): 
    try: 
        os.mkdir(os.getcwd() + '/VARIABILITY/DATA_VARIABILITY/' + foldernames[i])
    except FileExistsError:
        pass
    # os.replace(os.getcwd() + '/' + filenames[i], os.getcwd() + '/VARIABILITY/DATA_VARIABILITY/' + foldernames[i] + '/' + filenames[i])
    shutil.copyfile(data_dir + filenames[i], os.getcwd() + '/VARIABILITY/DATA_VARIABILITY/' + foldernames[i] + '/' + filenames[i])

for i in range(len(og_names)):
    try:
        loc = np.where(og_names[i] + '-dered.dr16' == np.array(os.listdir(data_dir)))[0][0]
        print('yes')
        
    except IndexError:
        pass

# t = Table.read(FITS_DIREC)

# with fits.open(FITS_DIREC) as hdu_list:
#     hdu_list.info()
#     table_data = hdu_list[1].data
#     print('Column names: \n', table_data.names)
#     print('\nRow 1: \n', table_data[1])
#     print('\nNumber of rows: \n', len(table_data))
# %%
