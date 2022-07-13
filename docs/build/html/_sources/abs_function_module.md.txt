# `abs_function_module`

```{eval-rst}
.. function:: abs_function_module.wavelength_to_velocity

    Reads in a list of wavelength values to be converted to velocity.
    
    :param kind: redshift
    :type kind: list[int]
    :param kind: wavelength
    :type kind: list[int]
    :returns: velocity
    :rtype: array[int]
    
```

```{eval-rst}
.. function:: abs_function_module.smooth 

    Smooths out values entered.   
     
    :param kind: smooth_this
    :type kind: list[int]
    :param kind: box_size
    :type kind: int
    :returns: now_smooth
    :rtype: list[int]
    
```

```{eval-rst}
.. function:: abs_function_module.abs_parameters_plot_optional 

    Based off and does what find_absorption_parameters does, but also includes plotting.

    Reads in a list of redshift, wavelength, velocity limit (your integral bounds), broad absorption width, and percentage value 
    you want to go below the continuum to calculate BI. Returns BI for each indivdual trough, the total BI from all troughs, 
    vmin/vmaxs of the trough found if there are more than one, the equivalent width, the depth of the trough and the beta values
    that were converted from wavelength. Plots where CIV, CII, and OI would be *if* the EHVO absorption found was due to SiIV and 
    plots the spectra as normalzied flux vs velocity, and error vs velocity. 
    
    :param kind: z
    :type kind: list[int]
    
    :param kind: wavelength
    :type kind: list[int]
    
    :param kind: normalized_flux
    :type kind: list[int]
    
    :param kind: BALNICITY_INDEX_LIMIT
    :type kind: int
    
    :param kind: velocity_limits
    :type kind: namedtuple
    
    :param kind: percent
    :type kind: float
    
    :param kind: plots
    :type kind: string
    
    :returns: BI_total
    :rtype: list[int]
    
    :returns: BI_individual
    :rtype: list[int]
    
    :returns: BI_all
    :rtype: list[int]
    
    :returns: vmins
    :rtype: array[int]
    
    :returns: vmaxs
    :rtype: array[int]
    
    :returns: EW_individual
    :rtype: list[int]
    
    :returns: final_depth_individual
    :rtype: list[int]
    
    :returns: beta
    :rtype: array[int]
    
    :returns: vminindex_for_range
    :rtype: int
    
    :returns: vmaxindex_for_range
    :rtype: int
    
```














