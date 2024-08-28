# README File for Taylor Gibbons' Mathematical Study of Strange Variability

## Table of Contents:
* [Project Information](#project-information)
* [Code](#code)
* [Setup and Launch](#setup-and-launch)
* [Uncommited Files](#currently-uncommited-files)


# Project Information

Quasars are the most luminous type of Active Galactic Nuclei. Both consist of a supermassive blackhole at the center with a surrounding luminous accretion disk. Quasars can expel clouds of gas outwards, which are called outflows. Historically quasar outflows have been detected by the presence of absorption lines in quasar spectra, with the most typical ion detected being CIV (Carbon ionized 3 times, or C3+). Some of these outflows reach relativistic speeds beyond 10% of the speed of light. Our team has built the largest database of these outflows, called extremely high velocity outflows. All outflows have been found to be variable between observations, but the exact cause of this variability is unknown. The most prominent theories are variability due to motion in and out of our line-of-sight, and the recombination/ionization of the outflow gas. While carrying out our study of quasar outflow variability, we found that several cases have a common minimum flux value but vary almost everywhere else within their absorption trough; we call this ‘strange’ variability. Its importance is that it can help us constrain the cause of variability since it sets restrictions on several absorption parameters: coverage fraction, how much the source of light is covered, and the optical depth, how easily light is able to travel through the gas.

# Code
This is the list of code that is included in the folder (in order of use):
* ## StrangeVariabilityCoverageOptical.py
  - This set of code is front end of the Strange Variability Mathematical Suite (SVMS). It contains no functions but is instead the code that calls the functions that have been made, shown below.
* ## StrangeVariabilityCalculationFunctions.py
  - This set of code is the backend of calculating the values for the minimum coverage fraction and optical depth. This code contains:
    - ### get_spectral_values
      - Calculates and returns the minimum coverage fraction for both the first and second observation (Spectral_value_1, Spectral_value_2), as well as a ratio between them (Spectral_value_Cf2Cf1)
    - ### min_value_finder
      - Calculates the list of coverage fraction values based on the observed spectra data
    - ### tau_finder
      - Calculates the list of optical depth values that will explain the values we see from the observed spectra data.
* ## StrangeVariabilityPlottingFunctions.py
  - This set of code is the backend of all the plotting functions for the SVMS. This code contains:
    - ### visualize
      - A very early testing of visualizing the observed and calculated data. Still in progress
    - ### Cf_tau_grapher
      - The main plotting function of the SVMS. It takes in a data frame and will result in graphs for each of the detections of Strange Variability. There is also an option of grouping the plots by the alpha values within an interval
    - ### alpha_group_grapher
      - The first function to be called when trying to group values within an interval. This function returns a new dataframe that is based on a premade template of I01,I02, I1,I2 and alpha values where the calculated data fills another column. This dataframe is then sent into Cf_tau_grapher with alpha_group = 'yes' when calling. (See below for more details)
    - ### alpha_testing
      - A testing function to find out if we need to be super strict on the numbering of our observations, or if it is okay to just know that, in general, there is one observation with a larger I0 than the other. This code results in graphs interchanging the equations and plots that show the observations are just inversed. 
* ## AlphaChangeTest.py
  - This very simple code calls the alpha_testing function to produce the graphs. This is not a needed program to get results

# Setup and Launch
This is how you would set up your system to be able to use the code that has been made first by Dr. Paola Rodriguez Hidalgo, myself: Taylor Gibbons, with help from Alex Salley and Imad Morsli.
1. In order to start running the code, make sure to have the output csv file from Alex Salley's portion: Detecting Strange Variability. This will be a csv file containing the name of the object that the detection is from, the spectra names, two flux values from around the trough (I01,I02) and two values at the bottom of the trough (I1,I2). After downloading the folder STRANGEVARIABILITY-MATHEMATICAL STUDY, the out.csv file will need to be inserted into the CSV Files folder.
2. At default, the process of making graphs, grouped graphs and saving figures is set to no. This means the code will only do the calculations and create a dataframe/xlsx file.
3. Make sure the `filename` filepath is correctly set to find the out.csv or there is also an option for the out.xlsx if the file format permits. At this point, the code should be ready to run by using all the folders within the parent SV folder.
4. Running for the first time will only dave the dataframe into an `.xlsx` file and into a pickle (`.pkl`). The excel files are the ones to be read by people, however if read through a code, it will end up giving some of the columns in the form of strings, instead of the desired list of floats. The pickle files are python read only and can be handled similarly to how excel files are read and saved.
5. To run the graphing function, change the `graph = 'no'` to `graph = 'yes'`
6. If you would like to create plots based on an interval of the alpha value; default is set to $`interval = \left(\alpha \pm 0.05 | \alpha \in [0.25,0.35,0.45, ... ,0.95]\right)`$

# Currently Uncommited Files:
IN PROGRESS
`git status`
