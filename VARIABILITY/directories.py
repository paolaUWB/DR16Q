'''  '''
#%%
import os
import sys
from astropy.io import fits
sys.path.insert(0, os.getcwd() + '/../' + 'DR16Q/') # changes the directory to the DR16Q --> all paths after this will need to be written as if this was in the top level of the DR16Q
from utility_functions import read_file, clear_file, append_row_to_csv#, print_to_file


DR = '16' # data release you are working with (i.e. '9' or '16')

fits_direc = os.getcwd() + '/../DR16Q_DATA/DR16Q_v4.fits'
make_direc = os.getcwd() + '/VARIABILITY/DATA_VARIABILITY/'

data_req_file = os.getcwd() + '/VARIABILITY/data_request.csv'
ehvo_dup_file = os.getcwd() + '/VARIABILITY/ehvo_duplicate_file.csv'


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

fits_nspec = data['NSPEC'].astype(int)

fits_duplicate_PLATE = data['PLATE_DUPLICATE'].astype(str)
fits_duplicate_MJD = data['MJD_DUPLICATE'].astype(str)
fits_duplicate_FIBER = data['FIBERID_DUPLICATE'].astype(str)

hdu.close()


#%%
##------ clear csv files
clear_file(ehvo_dup_file)
clear_file(data_req_file)

##------ create directories
fits_names = []
fits_duplicates = []

ehvo_names = []
ehvo_duplicates = []

foldernames = []

ehvo_fields = ['sdssName', 'ehvoPlate', 'ehvoMJD', 'ehvoFiber', 'duplicate_spectra']
append_row_to_csv(ehvo_dup_file, ehvo_fields)

data_req_fields = ['sdssName', 'Plate', 'MJD', 'Fiber']
append_row_to_csv(data_req_file, data_req_fields)


#%%
for i in range(len(fits_duplicate_PLATE)):
    duplicates = []
    test_duplicates = []
    
    fits_names.append('spec-' + fits_PLATE[i].zfill(4) + '-' + fits_MJD[i] + '-' + fits_FIBER[i].zfill(4))
    foldernames.append('J' + fits_SDSS_NAME[i]) ## move this lower (only needs to be created for ehvo duplicates)
    
    if (i % 50000) == 0:
        print('You have loaded ', i, ' spectra!')

    if fits_duplicate_PLATE[i][0] != '-1':
        
        for ehvo in range(len(ehvo_spectra_list)):
            if DR == '16':
                ehvo_spec = ehvo_spectra_list[ehvo][:-11]
            if DR == '9': 
                ehvo_spec = ehvo_spectra_list[ehvo][:-4]


            if ehvo_spec == fits_names[i]:
                ehvo_PLATE = fits_PLATE[i].zfill(4)
                ehvo_MJD = fits_MJD[i].zfill(4)
                ehvo_FIBER = fits_FIBER[i].zfill(4)
                
                for j in range(fits_nspec[i]):
                    dup_spec_name = 'spec-' + fits_duplicate_PLATE[i][j].zfill(4) + '-' + fits_duplicate_MJD[i][j].zfill(4) + '-' + fits_duplicate_FIBER[i][j].zfill(4)
                    duplicates.append(dup_spec_name)
                    data_req_fields = [fits_SDSS_NAME[i], fits_duplicate_PLATE[i][j].zfill(4), fits_duplicate_MJD[i][j].zfill(4), fits_duplicate_FIBER[i][j].zfill(4)]
                    append_row_to_csv(data_req_file, data_req_fields)
                
                ehvo_fields = [fits_SDSS_NAME[i], ehvo_PLATE, ehvo_MJD, ehvo_FIBER, duplicates]
                append_row_to_csv(ehvo_dup_file, ehvo_fields)
                
