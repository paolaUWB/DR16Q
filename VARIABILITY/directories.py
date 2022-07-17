'''  '''
#%%
import os
import sys
import numpy as np
from astropy.io import fits
import shutil
# from VARIABILITY.EHVO_files import FITS_DUPLICATES_FILE
sys.path.insert(0, os.getcwd() + '/../' + 'DR16Q') # changes the directory to the DR16Q --> all paths after this will need to be written as if this was in the top level of the DR16Q
from utility_functions import print_to_file, read_file, clear_file

#read EHVOs
#open fits table
#create directories

DR = '16' # which data release are you working with? DRX (i.e. '9' or '16') 

save_txt_files = 'no' # do you want to save new output files? 'yes'/'no'
save_ehvo_txt_files = 'yes'
##------ fits table directory
FITS_DIREC = os.getcwd() + '/../DR16Q_v4.fits'

#-- new directory path
MAKE_DIREC = os.getcwd() + '/VARIABILITY/DATA_VARIABILITY/'

##------ output files
FITS_NAMES_FILE = os.getcwd() + '/VARIABILITY/fits_names.txt'
FITS_NAMES_DUPLICATES_FILE = os.getcwd() + '/VARIABILITY/fits_names_duplicates.txt'
FITS_DUPLICATES_FILE = os.getcwd() + '/VARIABILITY/fits_duplicates.txt'
EHVO_DUPLICATES_FILE = os.getcwd() + '/VARIABILITY/DR' + DR + '_ehvo_duplicates.txt'

#-- DR9 paths
if DR == '9':
    SPEC_DIREC = os.getcwd() + '/DR9Q_EHVO/DR9Q_EHVO_FILES/'
    EHVO_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + '/DR9Q_EHVO/DR9_EHVO_sorted_norm.csv'
    zem, snr, EHVO_spectra_list = read_file(EHVO_FILE)

#-- DR16 paths
if DR == '16':
    SPEC_DIREC = os.getcwd() + '/DR16Q_EHVO/DR16Q_EHVO_FILES/'
    NORM_SPEC_DIREC = os.getcwd() + '/DR16Q_EHVO/NORM_DR16Q_EHVO/'
    EHVO_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + '/DR16Q_EHVO/DR16_EHVO_sorted_norm.csv'
    zem, snr, EHVO_spectra_list = read_file(EHVO_FILE)
    # EHVO_NORM_FILE = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + '/DR16Q_EHVO/good_fit_EHVO.csv'
    # EHVO_spectra_list, EHVO_norm_spectra_list, idx = read_file(EHVO_NORM_FILE)

#%%
#-- opens fits table and reads data
hdu = fits.open(FITS_DIREC)
data = hdu[1].data

#-- get column data from fits table
fits_SDSS_NAME = data['SDSS_NAME'].astype(str)
fits_PLATE = data['PLATE'].astype(str)
fits_MJD = data['MJD'].astype(str)
fits_FIBER = data['FIBERID'].astype(str)
fits_OBJID = data['OBJID']
fits_nspecBOSS = data['NSPEC_BOSS']
fits_duplicate_PLATE = data['PLATE_DUPLICATE'].astype(str)
fits_duplicate_MJD = data['MJD_DUPLICATE'].astype(str)
fits_duplicate_FIBER = data['FIBERID_DUPLICATE'].astype(str)
hdu.close()

#%%
#-- clear files
if save_txt_files == 'yes':
    clear_file(FITS_NAMES_FILE)
    clear_file(FITS_NAMES_DUPLICATES_FILE)
    clear_file(FITS_DUPLICATES_FILE)

if save_ehvo_txt_files == 'yes':
    clear_file(EHVO_DUPLICATES_FILE)

#%%
#-- initializing variables to store spectra names & duplicates
fits_names = []
fits_duplicates = [] # any spectra with a duplicate

fits_duplicates_ehvo = [] # duplicate spectra of EHVO
ehvo_duplicates_index = [] # index of EHVO in original list
ehvo_duplicates = [] # original EHVO names (that have a duplicate)
foldernames = []
filenames = []
full_filenames = []

EHVO_norm_spectra = []

index = []

count = 0
#%%
# for file in os.listdir(os.getcwd() + '/VARIABILITY/'):
#     if file.startswith('spec-'):
#         full_filenames.append(file)
#         if file.endswith('-dered.dr16'):
#             file_name = file[:-11]
#             filenames.append(file_name)
#         if file.endswith('-dered.txt'):
#             file_name = file[:-10]
#             filenames.append(file_name)

# for i in range(len(EHVO_norm_spectra_list)):
#     EHVO_norm_spectra.append(EHVO_norm_spectra_list[i][:-9])
    
# for file in os.listdir(os.getcwd() + '/DR16Q_EHVO/NORM_DR16Q_EHVO/'):
#     if file.startswith('spec-'):
#         full_filenames.append(file)
        
#         if file.endswith('norm.dr16'):
#             file_name = file[:-9]
#             filenames.append(file_name)
            
#%%

##------
for i in range(len(fits_PLATE)): #750414 spectra total in DR16Q fits file
    duplicates = []
    fields = []
    fits_names.append('spec-' + fits_PLATE[i].zfill(4) + '-' + fits_MJD[i] + '-' + fits_FIBER[i].zfill(4))
    foldernames.append('J' + fits_SDSS_NAME[i])
    
    if save_txt_files == 'yes':
        print_to_file(fits_names[i], FITS_NAMES_FILE)
    
    if (i % 50000) == 0:
        print('You have loaded ', i, ' spectra!')

    for j in range(74): # 74 *possible* observations per object
        if fits_duplicate_PLATE[i][j] != '-1':
            duplicates.append('spec-' + fits_duplicate_PLATE[i][j].zfill(4) + '-' + fits_duplicate_MJD[i][j].zfill(4) + '-' + fits_duplicate_FIBER[i][j].zfill(4))

    if duplicates != []:
        fits_duplicates.append(duplicates)
        fields = [i, fits_names[i], duplicates]
        
        for k in range(len(EHVO_spectra_list)):
            if DR == '16':
                EHVO_spec = EHVO_spectra_list[k][:-11]
            if DR == '9': 
                EHVO_spec = EHVO_spectra_list[k][:-4]
            
            for p in range(len(duplicates)):
                if (duplicates[p] == EHVO_spec): 
                    ehvo_duplicates_index.append(i)
                    ehvo_duplicates.append(fits_names[i])
                    fits_duplicates_ehvo.append(duplicates)
                    count+=1
                    
                    if save_ehvo_txt_files == 'yes':
                        print_to_file(fields, EHVO_DUPLICATES_FILE)
                        
                    try: 
                        os.mkdir(MAKE_DIREC + foldernames[i])
                    except FileExistsError:
                        pass
                
                    if DR == '16':
                        shutil.copyfile(SPEC_DIREC + EHVO_spectra_list[k], MAKE_DIREC + foldernames[i] + '/' + EHVO_spectra_list[k])
                        # shutil.copyfile(NORM_SPEC_DIREC + EHVO_norm_spectra_list[k], MAKE_DIREC + foldernames[i] + '/' + EHVO_norm_spectra_list[k])

                    # if DR == '9': 
                        # EHVO_spectra_list_norm = EHVO_spectra_list[k][:-4] + 'norm.dr9'
                        # shutil.copyfile(SPEC_DIREC + '/' + EHVO_spectra_list_norm, os.getcwd() + '/VARIABILITY/DATA_VARIABILITY/' + foldernames[i] + '/' + EHVO_spectra_list_norm)
            
                # for x in range(len(full_filenames)):
                    # if (duplicates[p] == filenames[x]):
                        # try:
                        #     shutil.move(os.getcwd() + '/VARIABILITY/' + full_filenames[x], MAKE_DIREC + foldernames[i] + '/' + full_filenames[x])
                        # except:
                        #     pass
                        
            if (fits_names[i] == EHVO_spec):
                ehvo_duplicates_index.append(i)
                ehvo_duplicates.append(fits_names[i])
                fits_duplicates_ehvo.append(duplicates)
                
                if save_ehvo_txt_files == 'yes':
                    print_to_file(fields, EHVO_DUPLICATES_FILE)
                
                try: 
                    os.mkdir(MAKE_DIREC + foldernames[i])
                except FileExistsError:
                    pass
                
                if DR == '16':
                    shutil.copyfile(SPEC_DIREC + EHVO_spectra_list[k], MAKE_DIREC + foldernames[i] + '/' + EHVO_spectra_list[k])
                    # shutil.copyfile(NORM_SPEC_DIREC + EHVO_norm_spectra_list[k], MAKE_DIREC + foldernames[i] + '/' + EHVO_norm_spectra_list[k])
                # if DR == '9': 
                #     EHVO_spectra_list_norm = EHVO_spectra_list[k][:-4] + 'norm.dr9'
                #     shutil.copyfile(SPEC_DIREC + '/' + EHVO_spectra_list_norm, os.getcwd() + '/VARIABILITY/DATA_VARIABILITY/' + foldernames[i] + '/' + EHVO_spectra_list_norm)
                
    else:
        fits_duplicates.append('-1')
    
    if fits_duplicate_PLATE[i][0] != '-1':   
        index.append(i) # check that this matches the index file (but +1 in file)
        
    if save_txt_files == 'yes':
        if duplicates != []:
            print_to_file(fields, FITS_NAMES_DUPLICATES_FILE) # saves index in fits_names file, name of the original spectra, and any duplicates of that spectra
        print_to_file(fits_duplicates[i], FITS_DUPLICATES_FILE) # saves ONLY the duplicates
    

# for i in range(len(filenames)):
#     os.replace(os.getcwd() + '/VARIABILITY/' + filenames[i], MAKE_DIREC + foldernames[i] + '/' + filenames[i])

###################################

# filenames = []
# foldernames = []
# index = []

# original = []
# duplicates_only = []
# fits_duplicate_index = []
'''
EHVO_duplicates = []
for q in range(len(original)): # original contains 119924 spectra (only spectra with duplicate observations)
# for q in range(len(original)): 
    for j in range(len(EHVO_spectra_list_DR16)): # i starts at 0
        # EHVO_DR16_minus_dered = EHVO_spectra_list_DR16[j][:-11] # -11 for dr16, -10 for dr9
        fields2 = [q, j, EHVO_spectra_list_DR16[j][:-11], original[q]]
        if original[q] == EHVO_spectra_list_DR16[j][:-11]:
            # append_row_to_csv(os.getcwd() + '/EHVO_duplicates.csv', fields2)
            foldernames.append('J' + fits_OBJID[fits_duplicate_index[q]]) # possibly fits_OBJID[q][:4] instead (so just first 4 of RA)
            # foldernames = 'J' + fits_OBJID[q]
            data_dir = os.getcwd() + '/DR16Q_EHVO/DR16Q_EHVO_FILES/'
            try: 
                os.mkdir(os.getcwd() + '/VARIABILITY/DATA_VARIABILITY/' + foldernames[j])
            except FileExistsError:
                pass
            shutil.copyfile(data_dir + '/' + EHVO_spectra_list_DR16[j], os.getcwd() + '/VARIABILITY/DATA_VARIABILITY/' + foldernames[i] + '/' + EHVO_spectra_list_DR16[j])
'''
            
    # try: 
    #     loc = np.where(EHVO_DR16_minus_dered == original[q])
    #     # print(loc)
    # except:
    #     pass
'''
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
data_dir = os.getcwd() + '/../DR16Q_EHVO/DR16Q_EHVO_FILES/' #possibly: /../ + /DATA/DR9Q_SNR10 #'C:/Users/Dakota/Documents/GitHub/DR16Q/DATA/DR9Q_SNR10'

for i in range(len(filenames)): 
    try: 
        os.mkdir(os.getcwd() + '/DATA_VARIABILITY/' + foldernames[i])
    except FileExistsError:
        pass
    # os.replace(os.getcwd() + '/' + filenames[i], os.getcwd() + '/VARIABILITY/DATA_VARIABILITY/' + foldernames[i] + '/' + filenames[i])
    shutil.copyfile(data_dir + filenames[i], os.getcwd() + '/DATA_VARIABILITY/' + foldernames[i] + '/' + filenames[i])

for i in range(len(og_names)):
    try:
        loc = np.where(og_names[i] + '-dered.dr16' == np.array(os.listdir(data_dir)))[0][0]
        print('yes')
        
    except IndexError:
        pass
'''

# t = Table.read(FITS_DIREC)

# with fits.open(FITS_DIREC) as hdu_list:
#     hdu_list.info()
#     table_data = hdu_list[1].data
#     print('Column names: \n', table_data.names)
#     print('\nRow 1: \n', table_data[1])
#     print('\nNumber of rows: \n', len(table_data))
#%%
