# PLOT INDIVIDUAL SPECTRUM
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## QUICK OVERVIEW OF PROJECT FOLDER

This project contains a modularized version of the individual spectrum plotting code. The program reads a normalized spectrum, processes the wavelength/flux/error data, and produces a plot showing the normalized spectrum, error spectrum, absorption regions, and emission-line labels.

The project is designed so that most parameters can be changed through `config.yaml` without modifying the main Python code.

### `plotindividualspectrum.py`
- Main Python script used to process and plot an individual normalized spectrum.
- Loads the configuration from `config.yaml`.
- Loads the normalized spectrum.
- Selects the desired wavelength range.
- Converts between observed-frame and rest-frame wavelengths.
- Calculates absorption redshifts from the specified outflow velocities.
- Smooths the normalized spectrum.
- Plots the normalized flux and error.
- Calls the absorption plotting function from `absorption.py`.
- Calls the emission-line labeling function from `emission.py`.
- Creates the rest-frame wavelength axis.
- Adds the quasar redshift label.
- Saves the final figure into `OUTPUT_FILES`.

### `absorption.py`
- Contains the `plot_absorption_regions()` function.
- Used to plot shaded absorption regions and their labels.
- Absorption lines can be individually turned on or off using Boolean values in `config.yaml`.
- Currently supports:
    - CIV
    - NV
    - OVI
    - SiIV
    - Lyα
    - Lyβ
    - CII
    - OI
- Absorption-region colors are controlled through the `colors` section of `config.yaml`.
- This allows the default colors to be changed without modifying `absorption.py`.

### `emission.py`
- Contains the `plot_emission_labels()` function.
- Used to place emission-line labels on the spectrum.
- Emission lines can be individually turned on or off using Boolean values in `config.yaml`.
- Currently supports:
    - CIV
    - SiIV + OIV]
    - Lyα + NV
    - OI
    - CII
    - OVI

### `lines.py`
- Contains the wavelength values used for the spectral lines.
- Keeps the wavelength constants separate from the plotting code.
- Contains weighted-average wavelengths and individual doublet wavelengths.
- Currently contains values for:
    - CIV
    - SiIV
    - NV
    - OVI
    - CII
    - OI
    - Lyα
    - Lyβ
- Also contains the red and blue wavelengths of the relevant doublets.

### `config.yaml`
- Configuration file used to control the spectrum plot.
- Allows the user to change input/output settings and plotting parameters without modifying the main Python script.
- Contains:
    - Spectrum name
    - Quasar redshift
    - Rest-frame wavelength limits
    - Plot y-axis limit
    - Redshift label location
    - Absorption velocity limits
    - Absorption line Boolean settings
    - Emission line Boolean settings
    - Smoothing value
    - Absorption colors

### `OUTPUT_FILES`
- Contains the figures produced by `plotindividualspectrum.py`.
- The output directory is automatically created if it does not already exist.
- The current code is configured to save the spectrum as a PDF.
- PNG output can also be enabled in the main script.

### `README.md`
- This file.
- Contains documentation for the project and instructions for running and modifying the code.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## RUNNING THE CODE

### Files needed to run `plotindividualspectrum.py`

The project directory should contain:

- `plotindividualspectrum.py`
- `absorption.py`
- `emission.py`
- `lines.py`
- `config.yaml`
- The normalized spectrum file being plotted
- `OUTPUT_FILES` (this can be created automatically)
- `README.md` is not required for the code to run

Example:

```text
Proj_plotindividualspectrum_mod/
│
├── plotindividualspectrum.py
├── absorption.py
├── emission.py
├── lines.py
├── config.yaml
├── spec-10431-58137-0135norm.dr16
├── README.md
│
└── OUTPUT_FILES/