# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 11:41:04 2026

@author: elijahfacklam
"""
import numpy as np
import matplotlib.pyplot as plt

def auto_ylim(ax, x_values=None, y_values=None, x_limits=None,
              Normalization=False, padding=0.05):
    """
    Automatically sets the y-axis limits for a matplotlib plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The plot where the y-limits will be applied.

    x_values, y_values : list or array, optional
        Data used to calculate the limits.
        Can be one dataset or a list of datasets.
        If neither is provided, use the lines already on the plot.
        If providing data, both x_values and y_values are required.

    x_limits : tuple, optional
        Only use data between (xmin, xmax) to calculate the y-limits.
        This does not change the displayed x-axis range.

    Normalization : bool, optional
        Kept so older code using this option still works.
        Currently does not change the calculation.

    padding : float, optional
        Extra space above and below the data.
        A value of 0.05 adds 5% of the y-range to each end.

    Returns
    -------
    tuple or None
        The new (lower, upper) y-limits.
        Returns None if no usable data is found.
    """

    # Check that padding is a valid number.
    padding = float(padding)

    if not np.isfinite(padding) or padding < 0:
        raise ValueError("padding must be a finite number of 0 or greater.")

    # If providing data, both x and y are needed.
    if x_values is None and y_values is not None:
        raise ValueError("Please provide x_values as well as y_values.")

    if y_values is None and x_values is not None:
        raise ValueError("Please provide y_values as well as x_values.")

    # Check the optional x-range.
    if x_limits is not None:
        x_limits = np.asarray(x_limits, dtype=float)

        if x_limits.shape != (2,):
            raise ValueError("x_limits must contain two numbers: (xmin, xmax).")

        if not np.all(np.isfinite(x_limits)):
            raise ValueError("x_limits cannot contain NaN or infinity.")

        xmin, xmax = x_limits

        if xmin > xmax:
            raise ValueError("xmin must be less than or equal to xmax.")

    # ---------------------------------------------------------------
    # Get data from the plot if no data was provided.
    # ---------------------------------------------------------------
    if x_values is None and y_values is None:
        lines = ax.get_lines()

        x_values = []
        y_values = []

        for line in lines:
            x_values.append(line.get_xdata())
            y_values.append(line.get_ydata())

    # ---------------------------------------------------------------
    # Put the data into a consistent format:
    # a list containing one array for each dataset.
    # ---------------------------------------------------------------

    # Handle x-values.
    if isinstance(x_values, np.ndarray):
        if x_values.ndim == 1:
            x_values = [x_values]
        elif x_values.ndim == 2:
            x_values = list(x_values)
        else:
            raise ValueError("x_values must contain one-dimensional datasets.")

    elif isinstance(x_values, (list, tuple)):
        # A plain list such as [1, 2, 3] is one dataset.
        if len(x_values) > 0:
            if all(np.isscalar(value) for value in x_values):
                x_values = [x_values]

    else:
        x_values = [x_values]

    # Handle y-values in the same way.
    if isinstance(y_values, np.ndarray):
        if y_values.ndim == 1:
            y_values = [y_values]
        elif y_values.ndim == 2:
            y_values = list(y_values)
        else:
            raise ValueError("y_values must contain one-dimensional datasets.")

    elif isinstance(y_values, (list, tuple)):
        if len(y_values) > 0:
            if all(np.isscalar(value) for value in y_values):
                y_values = [y_values]

    else:
        y_values = [y_values]

    # Each x dataset must have a matching y dataset.
    if len(x_values) != len(y_values):
        raise ValueError(
            "x_values and y_values must have the same number of datasets."
        )

    # ---------------------------------------------------------------
    # Collect valid y-values within the chosen x-range.
    # ---------------------------------------------------------------
    y_collection = []

    for x, y in zip(x_values, y_values):

        # Convert each dataset into a numerical array.
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)

        if x.ndim != 1 or y.ndim != 1:
            raise ValueError("Each x and y dataset must be one-dimensional.")

        if len(x) != len(y):
            raise ValueError("Each x and y pair must have the same length.")

        # Keep only points where both x and y are finite.
        # This removes NaN and infinity.
        mask = np.isfinite(x) & np.isfinite(y)

        # If an x-range was provided, keep only points in that range.
        if x_limits is not None:
            mask = mask & (x >= xmin) & (x <= xmax)

        if np.any(mask):
            y_collection.append(y[mask])

    # Leave the plot unchanged if no usable points were found.
    if len(y_collection) == 0:
        print("Warning: No valid data found for the selected x-range.")
        return None

    # ---------------------------------------------------------------
    # Calculate the minimum and maximum y-values.
    # ---------------------------------------------------------------
    y_all = np.concatenate(y_collection)

    ymin = np.min(y_all)
    ymax = np.max(y_all)

    # Handle a flat line so the limits are not identical.
    if ymin == ymax:
        # Use 5% of the constant value, with a minimum of 0.05.
        extra_space = max(abs(ymin) * 0.05, 0.05)

        ymin = ymin - extra_space
        ymax = ymax + extra_space

    # ---------------------------------------------------------------
    # Add space above and below the data.
    # ---------------------------------------------------------------
    y_range = ymax - ymin
    extra_space = padding * y_range

    ymin_final = ymin - extra_space
    ymax_final = ymax + extra_space

    ax.set_ylim(ymin_final, ymax_final)

    return ymin_final, ymax_final

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
