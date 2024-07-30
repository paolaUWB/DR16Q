# Overview of files

## find_weird_variability.py
This is the 'main' program. As far as I am aware it must be ran from the command line. The first command line argument should be to absorption_table.csv (this is produced by the code in ABSORPTION.) The second command line argument should be to ehvo_duplicate_file.csv, produced by make_directories.py in VARIABILITY



This should produce a CSV and PDF file (the latter can be disabled by setting GRAPH_RESULTS to false.) 

## graph_utilities.py
This file contains all the methods used for displaying spectra with MatPlotLib.

## strange_utilities.py
This file contains miscellaneous utilities that other projects may want to use. This includes things like getSpectrum

# How to run

## Building sample
This program uses outputs from VARIABILITY/make_directories and ABSORPTION/absorption.py . It is recommended that you start with absorption

### Absorption
You need to use a modified version of the ABSORPTION code stored on the Google Drive (email/DM me if you don't have access to it)

In ABSORPTION run
`
python absorption_modification_variability.py good_fit.csv
`
Where good_fit.csv is a list of normalized spectra which have a 'good fit' for normalization. There should be 16,512 elements in good_fit.csv . If you don't have access to it on the Google Drive, please dm or email me

### Variability
In VARIABILITY use make_directories.py to generate a list of repeated observations of the same objects
`
python make_directories.py ../ABSORPTION/OUTPUT_FILES/absorption_table.csv
`
## Executing find_weird_variability from command line
In the STRANGE_VARIABILITY directory, run
`
python find_weird_variability.py ../ABSORPTION/OUTPUT_FILES/absorption_table.csv ../VARIABILITY/ehvo_duplicate_file.csv
`

## Executing find_weird_variability without command line
If you can't run Python from the command line for whatever reason, hardcode the values of absorptionCsvFilename and duplicateCsvFilename. An example of this would be
`
absorptionCsvFilename = "absorption.csv"
duplicateCsvFilename  = "../../examples/duplicates.csv"
`
## Notes
If you want graphs to be generated, ensure that GRAPH_RESULTS is set to true

# Common issues
If you're not using version v0.4.7 of AstroQuery the program may be unable to download the spectra
