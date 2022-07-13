'''  '''
#%%
import os
import sys
import numpy as np
from astropy.io import fits
import shutil
# from VARIABILITY.EHVO_files import FITS_DUPLICATES_FILE
sys.path.insert(0, os.getcwd() + '/../' + 'DR16Q/') # changes the directory to the DR16Q --> all paths after this will need to be written as if this was in the top level of the DR16Q
from utility_functions import print_to_file, read_file, clear_file

#read EHVOs  -
#open fits table -
#create directories


DR = '16' # data release you are working with (i.e. '9' or '16')

fits_direc = os.getcwd() + '/../DR16Q_v4.fits'
make_direc = os.getcwd() + '/VARIABILITY/DATA_VARIABILITY/'

data_req_file = os.getcwd() + '/VARIABILITY/data_request.txt'


##------ read EHVOs
#-- dr9
if DR == '9':
    ehvo_file = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + '/DR9Q_EHVO/DR9_EHVO_sorted_norm.csv'
    zem, snr, ehvo_spectra_list = read_file(ehvo_file)
    
#-- dr16
if DR == '16':
    ehvo_file = sys.argv[1] if len(sys.argv) > 1 else os.getcwd() + '/DR16Q_EHVO/DR16_EHVO_sorted_norm.csv'
    zem, snr, ehvo_spectra_list = read_file(ehvo_file)


#%%
##------ open fits table
hdu = fits.open(fits_direc)
data = hdu[1].data

fits_SDSS_NAME = data['SDSS_NAME'].astype(str)
fits_PLATE = data['PLATE'].astype(str)
fits_MJD = data['MJD'].astype(str)
fits_FIBER = data['FIBERID'].astype(str)
fits_OBJID = data['OBJID']
fits_nspecSDSS = data['NSPEC_SDSS'] # The number of additional spectra for an object from SDSS-I/II
fits_nspecBOSS = data['NSPEC_BOSS'] # The number of additional spectra for an object from BOSS/eBOSS
fits_duplicate_PLATE = data['PLATE_DUPLICATE'].astype(str)
fits_duplicate_MJD = data['MJD_DUPLICATE'].astype(str)
fits_duplicate_FIBER = data['FIBERID_DUPLICATE'].astype(str)
hdu.close()


#%%
##------ create directories
fits_names = []
fits_duplicates = []

foldernames = []

duplicates = []

count1 = 0
count2 = 0 

for i in range(len(fits_PLATE)):
    # duplicates = []
    
    fits_names.append('spec-' + fits_PLATE[i].zfill(4) + '-' + fits_MJD[i] + '-' + fits_FIBER[i].zfill(4))
    foldernames.append('J' + fits_SDSS_NAME[i])
    
    if (i % 50000) == 0:
        print('You have loaded ', i, ' spectra!')

    for j in range(74): # 74 *possible* observations per object
        if fits_duplicate_PLATE[i][j] != '-1':
            duplicates.append('spec-' + fits_duplicate_PLATE[i][j].zfill(4) + '-' + fits_duplicate_MJD[i][j].zfill(4) + '-' + fits_duplicate_FIBER[i][j].zfill(4))
    
    if duplicates != []:
        fits_duplicates.append(duplicates)
        
        for k in range(len(ehvo_spectra_list)):
            if DR == '16':
                ehvo_spec = ehvo_spectra_list[k][:-11]
            if DR == '9': 
                ehvo_spec = ehvo_spectra_list[k][:-4]
        
        try: 
            os.mkdir(MAKE_DIREC + foldernames[i])
        except FileExistsError:
            pass
    

#%%



save_txt_files = 'no' # do you want to save new output files? 'yes'/'no'
save_ehvo_txt_files = 'yes'

##------ output files
FITS_NAMES_FILE = os.getcwd() + '/VARIABILITY/fits_names.txt'
FITS_NAMES_DUPLICATES_FILE = os.getcwd() + '/VARIABILITY/fits_names_duplicates.txt'
FITS_DUPLICATES_FILE = os.getcwd() + '/VARIABILITY/fits_duplicates.txt'
EHVO_DUPLICATES_FILE = os.getcwd() + '/VARIABILITY/DR' + DR + '_ehvo_duplicates.txt'

#%%
#-- clear files
if save_txt_files == 'yes':
    clear_file(FITS_NAMES_FILE)
    clear_file(FITS_NAMES_DUPLICATES_FILE)
    clear_file(FITS_DUPLICATES_FILE)

if save_ehvo_txt_files == 'yes':
    clear_file(EHVO_DUPLICATES_FILE)

#%%
fits_duplicates_ehvo = [] # duplicate spectra of EHVO
ehvo_duplicates_index = [] # index of EHVO in original list
ehvo_duplicates = [] # original EHVO names (that have a duplicate)
filenames = []
full_filenames = []

EHVO_norm_spectra = []

index = []

count = 0
            
#%%

##------
for i in range(len(fits_PLATE)): #750414 spectra total in DR16Q fits file
    fields = []
    
    if save_txt_files == 'yes':
        print_to_file(fits_names[i], FITS_NAMES_FILE)
    
    if (i % 50000) == 0:
        print('You have loaded ', i, ' spectra!')


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
 