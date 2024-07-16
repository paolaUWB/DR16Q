##find_weird_variability.py
This is the 'main' program. As far as I am aware it must be ran from the command line. The first command line argument should be to absorption_table.csv (this is produced by the code in ABSORPTION.) The second command line argument should be to ehvo_duplicate_file.csv, produced by make_directories.py in VARIABILITY

`
python find_weird_variability.py ../ABSORPTION/OUTPUT_FILES/absorption_table.csv ../VARIABILITY/ehvo_duplicate_file.csv
`

This should produce a CSV and PDF file (the latter can be disabled by setting GRAPH_RESULTS to false.) 

##graph_utilities.py
This file contains all the methods used for displaying spectra with MatPlotLib.

##strange_utilities.py
This file contains miscellaneous utilities that other projects may want to use. This includes things like getSpectrum
