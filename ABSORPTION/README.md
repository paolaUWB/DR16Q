# DR16Q --> ABSORPTION
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## QUICK OVERVIEW OF ABSORPTION FOLDER

### `absorption.py`
- Contains a way to run multiple spectra through the functions necessary to calculate BI, vmin, vmax, EW, and depth. Also, outputs graphs of the spectra when absorption is found.

### `abs_plot_module.py`
- Modules that contriute to plotting spectra in terms of velocities and explicitly plotting various absorption features. 

### `abs_function_module.py`
- Contains various functions to run absorption code.

### `OUTPUT_FILES`
- What in this folder will vary but from absorption.py the output files will be placed here, that includes the text file and pdf-graphs.
- If the user changes the name of the files multiple graphs and text files will be created. 

### `README.md`
- this file.
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## RUNNING ABSORPTION FILE

- Files needed to run `absorption.py`:
    - `good_fit.csv` from **NORMALIZATION** folder

<!--- add "instructions" and/or important things to know about running the code pending, might delete this seciton in general. --->

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## FILES IN ABSORPTION DIRECTORY

### absorption.py
- Searches for absorption based on the parameters listed below. It outputs pdf-graphs of the spectra and also with various absorption characteristics. It also outputs a text file with the values of various absorption characteristics which are listed below.   
- Default is it searches for non-shallow absorption, we have this as being 10% below the continuum. 
    - We have the contiuum set to a value of 1 from the normalization, so 0.9 would be 10% below the continuum, this is where it flags.  
    - Has a 3-spike limit, so after the spectrum crosses the 0.9 threshold three times it stops counting.
- Looking for broad absorption which we characterized at 2,000 km/s.
- BI value will start being counted after that braod absorption is met.
- Searches between 30,000 km/s and 60,000 km/s.
- These defaults can be changed to other values.
- When the pdf-graphs of spectra are produced the user can state if they would like to plot all of the spectra even if there are no absorption features or only spectra that contain absorption features based on the prerequistes entered above (i.e. the range of speed, non-shallow value, and broad absorption value).
- The user can also state whether to 'smooth' the spectra and can set a value for the amount of smoothness.
      
### OUTPUT_FILES --> textFILES
- Text file with :BI, v mins, v maxs, BI individual, equivalent width individual and depth.
    - If there are multiple troughs found the text file will have multiple values for everything except the total BI value. 
        
### OUTPUT_FILES --> pdfFILES
- pdf graph of spectra which shows normalized flux vs velociity.
- There is a horizontal black line that plots at a value of 1.37 to indicate where the broad absorption value is met.
- Visually shows where where CIV, CII and OI would be _if_ the EHVO absorption found were instead not EHVO and due to SiIV, plots a vertical line for that value and creates a colored bar that spans from the minimum velocity previously found to the maximum velocity. 
    - red is used for SiIV
    - gray is used for CIV
    - blue is used for CII
    - yellow is used for OI
        
### `abs_function_module.py`
- This is the meat of `absorption.py`, most of the code is in here. `absorption.py` uses this module to run it over many spectra.
- The `wavelength_to_velocity function` is in here.
- The `smooth` function is in here.

### `abs_plot_module.py`
- Contains all of the code necessary to plot spectra and various absorption characteristics. 
- Also contains a function for plots we utilize only for presentation purposes.










    