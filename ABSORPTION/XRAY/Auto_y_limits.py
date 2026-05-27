# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 11:41:04 2026

@author: elijahf
"""
import numpy as np
import matplotlib.pyplot as plt

def auto_ylim(ax, x_values=None, y_values=None, x_limits=None, Normalization=False, padding=0.05):
    """
    TO DO: Update the Normalization section by either simplifying it to a single
    function or testing spikes and creating logic for noise
    
    Automatically sets the y-axis limits for a matplotlib plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis object where the y-limits will be applied.

    x_values : list or array, optional
        The x-values corresponding to the plotted data. Can be a single array
        or a list of arrays. If not provided, values will be extracted from ax.

    y_values : list or array, optional
        The y-values corresponding to the plotted data. Can be a single array
        or a list of arrays. If not provided, values will be extracted from ax.

    x_limits : tuple, optional
        A tuple of (xmin, xmax). If provided, only y-values within this x-range
        will be used to determine limits.

    padding : float, optional
        padding is added to the min and max y-values (default is 5%).

    Returns
    -------
    None
        Sets the y-limits directly on the provided axis.
    """
    
        
     
    # Extract data from axis if not provided
    if x_values is None or y_values is None:
        lines = ax.get_lines()

        x_values = []
        y_values = []

        for line in lines:
            x_values.append(line.get_xdata())
            y_values.append(line.get_ydata())

     
    # Ensure inputs are lists for consistent handling
    if not isinstance(x_values, (list, tuple)):
        x_values = [x_values]

    if not isinstance(y_values, (list, tuple)):
        y_values = [y_values]

    # Collect relevant y-values based on x-range
    y_collection = []
    
    if Normalization == False:
        for x, y in zip(x_values, y_values):
            
            # Ensures Array
            x = np.asarray(x)
            y = np.asarray(y)
            
            if x_limits is not None:
                xmin, xmax = x_limits

                mask = (x >= xmin) & (x <= xmax)

                if np.any(mask):
                    y_collection.append(y[mask])

            else:
                y_collection.append(y)

         
        # Handle case where no valid data is found
        if len(y_collection) == 0:
            print("Warning: No valid y-values found for given x-range.")
            return
        
    if Normalization == True:
        for x, y in zip(x_values, y_values):
            
            # Ensures Array
            x = np.asarray(x)
            y = np.asarray(y)
            
            if x_limits is not None:
                xmin, xmax = x_limits

                mask = (x >= xmin) & (x <= xmax)

                if np.any(mask):
                    y_collection.append(y[mask])

            else:
                y_collection.append(y)

         
        # Handle case where no valid data is found
        if len(y_collection) == 0:
            print("Warning: No valid y-values found for given x-range.")
            return
            
    
    # Compute min and max y-values
    y_all = np.concatenate(y_collection)

    ymin = np.min(y_all)
    ymax = np.max(y_all)

     
    # Handle flat-line case
    
    if ymin == ymax:
        ymin = ymin - 1
        ymax = ymax + 1

     
    # Apply padding and set limits
     
    ymin_final = ymin * (1.0 - padding)
    ymax_final = ymax * (1.0 + padding)

    ax.set_ylim(ymin_final, ymax_final)

    return

#######################################################################
"""
Example use
"""
if __name__ == "__main__":

    # Create synthetic spectral-style data
    x = np.linspace(-70000, 0, 5000)

    # Simulated absorption feature
    signal = 1 - 0.4 * np.exp(-((x + 25000)**2) / (2 * (4000**2)))

    # Add noise
    noise = np.random.normal(0, 0.03, len(x))

    # Simulated observed spectrum
    spectrum = signal + noise

    # Simulated error array
    error = 0.05 + 0.01 * np.sin(x / 5000)

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 5))

    ax.plot(x, spectrum, label='Spectrum')
    ax.plot(x, error, label='Error')
    ax.axhline(1.0, linestyle='--', color='k', label='Continuum')

    # Visible x-range
    ax.set_xlim(-70000, 0)

    # Auto-scale y-limits using only the absorption region
    auto_ylim(ax,
              x_limits=(-45000, -5000),
              padding=0.05)

    plt.title("Auto Y-Limit Example")
    plt.xlabel("Velocity (km/s)")
    plt.ylabel("Normalized Flux")
    plt.legend()

    plt.show()
