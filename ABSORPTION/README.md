# DR16Q --> ABSORPTION

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# put some kind of overview/intro to what the code does & its "importance" here

Searches for absorption. 

-Default is it searches for non-shallow absorption, 10% below the continuum. 
    -We have the contiuum set to 1, and 0.9 to show where it flags.  
-Has a 3-spike limit, so after the spectrum exits 0.9 value three times it 
stops counting.
-BI value is set to start calculating at 2,000 km/s.
-Searches between 30,000 and 60,000 km/s


-These defaults can be changed to other values.

-outputs a text file with :BI value, max v, min v, depth, and equivalent width

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### RUNNING ABSORPTION FILE [absorption.py]

-Files needed to run absorption.py:
    -good_fit.csv [from NORMALIZATION]

# add "instructions" and/or important things to know about running the code


# add any variables to change and describe them 

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### FILES IN ABSORPTION DIRECTORY
absorption.py
-contains a way to run multiple spectra through the functions necessary 
to calculate BI, vmin, vmax, EW, and depth.

abs_plot_module.py
-modules that contriute to plotting spectra in terms of velocities and 
explicitly plotting various absorption features

abs_function_module.py
-contains various functions to run absorption code

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### OUTPUT FILES 

# possibly split between OUTPUT FILES (csv/txt) and OUTPUT FILES (pdf) if there 
# are many files. 
--> include a brief description of what each file contains. If txt or csv file, 
include what each column contains as well. 

text file with :BI value, max v, min v, depth, and equivalent width

pdf file of graph depending on user input:
    -graph can be smoothed or not
    -graph can include every spectrum analyzed or only those containig absorption