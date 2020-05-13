# OMEGA-Py documentation - v1.2

## `omegapy.omega_plots`

Display of OMEGAdata cubes.

`show_omega(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, title='auto', xlim=(None, None), ylim=(None, None), Nfig=None)`

`show_omega_v2(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, alpha=None, title='auto', lonlim=(None, None), latlim=(None, None), Nfig=None, polar=False, cbar=True, grid=True)`

`show_omega_interactif(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, title='auto', autoyscale=True, xlim=(None, None), ylim=(None, None))`

`show_omega_interactif_v2(omega, lam=1.085, refl=True, lam_unit='m', data=None, cmap='Greys_r', cb_title='IBD', title='auto', vmin=None, vmax=None, autoyscale=True, ylim_sp=(None, None), alpha=None, lonlim=(None, None), latlim=(None, None), polar=False, cbar=True, grid=True)`

`show_ibd_v2(omega, ibd, cmap='viridis', vmin=None, vmax=None, alpha=None, title='auto', cb_title = 'IBD', lonlim=(None, None), latlim=(None, None), Nfig=None, polar=False, cbar=True, grid=True)`

`show_omega_list_v2(omega_list, lam=1.085, lat_min=-90, lat_max=90, lon_min=0, lon_max=360, pas_lat=0.1, pas_lon=0.1, cmap='Greys_r', vmin=None, vmax=None, title='auto', Nfig=None, polar=False, cbar=True, cb_title='auto', data_list=None, mask_list=None, plot=True, grid=True, out=False, negatives_longitudes=False, **kwargs)`

`save_map_omega_list(omega_list, lat_min=-90, lat_max=90, lon_min=0, lon_max=360, pas_lat=0.1, pas_lon=0.1, lam=1.085, data_list=None, data_desc='', mask_list=None, sav_filename='auto', base_folder='../data/OMEGA/sav_map_list_v2/', ext='')`

`load_map_omega_list(filename)`

`show_omega_list_v2_man(data, grid_lat, grid_lon, infos, cmap='Greys_r', vmin=None, vmax=None, title='auto', Nfig=None, polar=False, cbar=True, cb_title='auto', grid=True, negatives_longitudes=False, **kwargs)`


### Display cube

~~~python
omegapy.omega_plots.show_omega(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None,
               title='auto', xlim=(None, None), ylim=(None, None), Nfig=None):
    Display an OMEGA/MEx observation in a rectangular pixel grid.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA/MEx observation
    lam : float
        The selected wavelength.
    refl : bool, optional (default True)
        True -> The reflectance is display.
        False -> The radiance is display.
    lam_unit : str, optional (default 'm')
        The unit of the `lam` parameter:
        | 'm' -> `lam` is the wavelength value (in µm).
        | else -> `lam` is the index of the wavelength in the omega.lam array (must be int).
    cmap : str, optional (default 'Greys_r')
        The matplotlib colormap.
    vmin : float or None, optional (default None)
        The lower bound of the coloscale.
    vmax : float or None, optional (default None)
        The upper bound of the colorscale.
    title : str, optional (default 'auto')
        The title of the figure.
    xlim : tuple of int or None, optional (default (None, None))
        The bounds of the x-axis of the figure.
    ylim : tuple of int or None, optional (default (None, None))
        The bounds of the y-axis of the figure.
    Nfig : int or str or None, optional (default None)
        The target figure ID.
~~~

~~~python
omegapy.omega_plots.show_omega_v2(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None,
                  alpha=None, title='auto', lonlim=(None, None), latlim=(None, None), Nfig=None,
                  polar=False, cbar=True, grid=True):
    Display an OMEGA/MEx observation with respect of the lat/lon coordinates of the pixels,
    and allows to use a polar projection if desired.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA/MEx observation
    lam : float
        The selected wavelength.
    refl : bool, optional (default True)
        True -> The reflectance is display.
        False -> The radiance is display.
    lam_unit : str, optional (default 'm')
        The unit of the `lam` parameter:
        | 'm' -> `lam` is the wavelength value (in µm).
        | else -> `lam` is the index of the wavelength in the omega.lam array (must be int).
    cmap : str, optional (default 'Greys_r')
        The matplotlib colormap.
    vmin : float or None, optional (default None)
        The lower bound of the coloscale.
    vmax : float or None, optional (default None)
        The upper bound of the colorscale.
    alpha : float or None, optional (default None)
        Opacity of the plot.
    title : str, optional (default 'auto')
        The title of the figure.
    lonlim : tuple of int or None, optional (default (None, None))
        The longitude bounds of the figure.
    latlim : tuple of int or None, optional (default (None, None))
        The latitude bounds of the y-axis of the figure.
    Nfig : int or str or None, optional, default None)
        The target figure ID.
    polar : bool, optional (default False)
        If True -> Use a polar projection for the plot.
    cbar : bool, optional (default True)
        If True -> Diplay the colorbar.
    grid : bool, optional (default True)
        Enable the display of the lat/lon grid.
~~~

### Display cube interactive version

~~~python
omegapy.omega_plots.show_omega_interactif(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', 
                          vmin=None, vmax=None, title='auto', autoyscale=True,
                          xlim=(None, None), ylim=(None, None)):
    Affichage interactif d'un cube de données.
    Possibilité d'afficher le spectre associé à un pixel en cliquant dessus
    (maintenir Ctrl pour supperposer plusieurs spectres), ou en se déplaçant avec les flèches.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA/MEx observation
    lam : float
        The selected wavelength.
    refl : bool, optional (default True)
        True -> The reflectance is display.
        False -> The radiance is display.
    lam_unit : str, optional (default 'm')
        The unit of the `lam` parameter:
        | 'm' -> `lam` is the wavelength value (in µm).
        | else -> `lam` is the index of the wavelength in the omega.lam array (must be int).
    cmap : str, optional (default 'Greys_r')
        The matplotlib colormap.
    vmin : float or None, optional (default None)
        The lower bound of the coloscale.
    vmax : float or None, optional (default None)
        The upper bound of the colorscale.
    title : str, optional (default 'auto')
        The title of the figure.
    xlim : tuple of int or None, optional (default (None, None))
        The bounds of the x-axis of the figure.
    ylim : tuple of int or None, optional (default (None, None))
        The bounds of the y-axis of the figure.
~~~

~~~python
omegapy.omega_plots.show_omega_interactif_v2(omega, lam=1.085, refl=True, lam_unit='m', data=None, 
                             cmap='Greys_r', cb_title='IBD', title='auto',
                             vmin=None, vmax=None, autoyscale=True, ylim_sp=(None, None),
                             alpha=None, lonlim=(None, None), latlim=(None, None),
                             polar=False, cbar=True, grid=True):
    Affichage interactif d'un cube de données.
    Possibilité d'afficher le spectre associé à un pixel en cliquant dessus
    (maintenir Ctrl pour supperposer plusieurs spectres), ou en se déplaçant avec les flèches.

    Display an OMEGA/MEx observation with respect of the lat/lon coordinates of the pixels,
    and allows to use a polar projection if desired.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA/MEx observation
    lam : float, optional (default 1.085)
        The selected wavelength.
    refl : bool, optional (default True)
        True -> The reflectance is display.
        False -> The radiance is display.
    lam_unit : str, optional (default 'm')
        The unit of the `lam` parameter:
        | 'm' -> `lam` is the wavelength value (in µm).
        | else -> `lam` is the index of the wavelength in the omega.lam array (must be int).
    data : 2D array or None, optional (default None)
        Array of high-level data (e.g. IBD map) computed from the omega observation.
    cmap : str, optional (default 'Greys_r')
        The matplotlib colormap.
    cb_title : str,  optional (default 'IBD')
        The title of the colorbar.
        Note : Only for the `data` plots.
    title : str, optional (default 'auto')
        The title of the figure.
    vmin : float or None, optional (default None)
        The lower bound of the coloscale.
    vmax : float or None, optional (default None)
        The upper bound of the colorscale.
    autoyscale : bool, optional (default True)
        | True -> Enable the auto-scaling of the spectra y-axis.
        | False -> Force use of the (vmin, vmax) bounds for the spectra plots.
    ylim_sp : tuble of float or None, optional (default (None, None))
        If autoyscale is False, can specify the bound values for the spectrum y-axis,
        other that (vmin, vmax).
    alpha : float or None, optional (default None)
        Opacity of the plot.
    lonlim : tuple of int or None, optional (default (None, None))
        The longitude bounds of the figure.
    latlim : tuple of int or None, optional (default (None, None))
        The latitude bounds of the y-axis of the figure.
    polar : bool, optional (default False)
        If True -> Use a polar projection for the plot.
    cbar : bool, optional (default True)
        If True -> Diplay the colorbar.
    grid : bool, optional (default True)
        Enable the display of the lat/lon grid.
~~~

### Display derived data map from OMEGA observation (ex IBD)

~~~python
omegapy.omega_plots.show_ibd_v2(omega, ibd, cmap='viridis', vmin=None, vmax=None, alpha=None, title='auto', 
                cb_title = 'IBD', lonlim=(None, None), latlim=(None, None), Nfig=None, 
                polar=False, cbar=True, grid=True):
    Affichage IBD avec pcolormesh.
    Display an OMEGA/MEx observation with respect of the lat/lon coordinates of the pixels,
    and allows to use a polar projection if desired.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA/MEx observation
    ibd : 2D array
        The array of the computed IBD values from the omega observation
    cmap : str, optional (default 'Greys_r')
        The matplotlib colormap.
    vmin : float or None, optional (default None)
        The lower bound of the coloscale.
    vmax : float or None, optional (default None)
        The upper bound of the colorscale.
    alpha : float or None, optional (default None)
        Opacity of the plot.
    title : str, optional (default 'auto')
        The title of the figure.
    cb_title : str,  optional (default 'IBD')
        The title of the colorbar.
    lonlim : tuple of int or None, optional (default (None, None))
        The longitude bounds of the figure.
    latlim : tuple of int or None, optional (default (None, None))
        The latitude bounds of the y-axis of the figure.
    Nfig : int or str or None, optional (default None)
        The target figure ID.
    polar : bool, optional (default False)
        If True -> Use a polar projection for the plot.
    cbar : bool, optional (default True)
        If True -> Display the colorbar.
    grid : bool, optional (default True)
        Enable the display of the lat/lon grid.
~~~

### Display composite map of several OMEGA observations, sample on a lat/lon grid

~~~python
omegapy.omega_plots.show_omega_list_v2(omega_list, lam=1.085, lat_min=-90, lat_max=90, lon_min=0, lon_max=360,
                       pas_lat=0.1, pas_lon=0.1, cmap='Greys_r', vmin=None, vmax=None, 
                       title='auto', Nfig=None, polar=False, cbar=True, cb_title='auto',
                       data_list=None, mask_list=None, plot=True, grid=True, out=False, 
                       negatives_longitudes=False, **kwargs):
    Display an composite map from a list OMEGA/MEx observations, sampled on a new lat/lon grid.

    Parameters
    ==========
    omega_list : array of OMEGAdata
        The list of OMEGA/MEx observations.
    lam : float, optional (default 1.085)
        The selected wavelength (in µm).
    lat_min : float, optional (default -90)
        The minimal latitude of the grid.
    lat_max : float, optional (default 90)
        The maximum latitude of the grid.
    lon_min : float, optional (default 0)
        The minimal longitude of the grid.
    lon_max : float, optional (default 360)
        The maximal longitude of the grid.
    pas_lat : float, optional (default 0.1)
        The latitude intervals of the grid.
    pas_lon : float, optional (default 0.1)
        The longitude intervals of the grid.
    cmap : str, optional (default 'Greys_r')
        The matplotlib colormap.
    vmin : float or None, optional (default None)
        The lower bound of the coloscale.
    vmax : float or None, optional (default None)
        The upper bound of the colorscale.
    title : str, optional (default 'auto')
        The title of the figure.
    Nfig : int or str or None, optional (default None)
        The target figure ID.
    polar : bool, optional (default False)
        If True -> Use a polar projection for the plot.
    cbar : bool, optional (default True)
        If True -> Diplay the colorbar.
    cb_title : str, optional (default 'auto')
        The title of the colorbar.
    data_list : 3D array or None, optional (default None)
        1D array of the same dimension of `omega_list` containing alternative maps (2D arrays),
        in the **same order** than the observations of `omega_list`.
    mask_list : 3D array
        1D array of the same dimension of `omega_list` containing the masks to remove the
        corrupted pixels of each observaiton, in the **same order** than the observations of 
        `omega_list`.
        Each mask is a 2D array, filled with 1 for good pixels and NaN for bad ones.
    plot : bool, optional (default True)
        If True -> Diplay the final figure.
    grid : bool, optional (default True)
        Enable the display of the lat/lon grid.
    out : bool, optional (default False)
        If True -> Return output.
    negatives_longitudes : bool, optional (default False)
        Argument for non-polar plots.
        | True -> longitudes between 0° and 360°.
        | False -> longitudus between -180° and 180°.
    **kwargs:
        Optional arguments for the plt.pcolormesh() function.

    Returns (if out=True)
    =======
    data : 2D array (dim : Nlon x Nlat)
        The omega reflectance at lam, sampled on the new lat/lon grid.
    mask : 2D array
        The array indicating where the new grid has been filled by the OMEGA data.
    grid lat : 2D array
        The new latitude grid.
    grid lon : 2D array
        The new longitude grid.
    mask_obs : 2D array of str
        The array indicating which observations have been used to fill each grid position.
~~~

### Save & restore composite map of several OMEGA observations, sample on a lat/lon grid
~~~python
omegapy.omega_plots.save_map_omega_list(omega_list, lat_min=-90, lat_max=90, lon_min=0, lon_max=360,
                        pas_lat=0.1, pas_lon=0.1, lam=1.085, data_list=None, data_desc='', 
                        mask_list=None, sav_filename='auto', 
                        base_folder='../data/OMEGA/sav_map_list_v2/', ext=''):
    Save the output of the omega_plots.show_omega_list_v2() function with the requested
    parameters as a dictionary.

    Parameters
    ==========
    omega_list : array of OMEGAdata
        The list of OMEGA/MEx observations.
    lat_min : float, optional (default -90)
        The minimal latitude of the grid.
    lat_max : float, optional (default 90)
        The maximum latitude of the grid.
    lon_min : float, optional (default 0)
        The minimal longitude of the grid.
    lon_max : float, optional (default 360)
        The maximal longitude of the grid.
    pas_lat : float, optional (default 0.1)
        The latitude intervals of the grid.
    pas_lon : float, optional (default 0.1)
        The longitude intervals of the grid.
    lam : float, optional (default 1.085)
        The selected wavelength (in µm).
    data_list : 3D array or None, optional (default None)
        1D array of the same dimension of `omega_list` containing alternative maps (2D arrays),
        in the **same order** than the observations of `omega_list`.
    data_desc : str, optional (default '')
        Description of the data contained in data_list (if used).
    mask_list : 3D array
        1D array of the same dimension of `omega_list` containing the masks to remove the
        corrupted pixels of each observaiton, in the **same order** than the observations of 
        `omega_list`.
        Each mask is a 2D array, filled with 1 for good pixels and NaN for bad ones.
    sav_filename : str, optional (default 'auto')
        The saving file name.
        | If 'auto' -> Automatically generated.
    base_folder : str, optional (default '../data/OMEGA/sav_map_list_v2/')
        The base folder to save the data.
    ext : str, optional (default '')
        Extension to add at the end of the filename (useful in case of automatic generation).
~~~

~~~python
omegapy.omega_plots.load_map_omega_list(filename):
    Load and return the result of omega_plots.show_omega_list_v2() previously saved
    with save_map_omega_list().

    Parameters
    ==========
    filename : str
        The file path.

    Returns
    =======
    data : 2D array
        The omega reflectance at lam, sampled on the new lat/lon grid.
    mask : 2D array
        The array indicating where the new grid has been filled by the OMEGA data.
    grid lat : 2D array
        The new latitude grid.
    grid lon : 2D array
        The new longitude grid.
    mask_obs : 2D array of str
        The array indicating which observations have been used to fill each grid position.
    infos : dict
        The informations about the computation of the data.
~~~

~~~python
omegapy.omega_plots.show_omega_list_v2_man(data, grid_lat, grid_lon, infos, cmap='Greys_r', vmin=None, vmax=None, 
                           title='auto', Nfig=None, polar=False, cbar=True, cb_title='auto',
                           grid=True, negatives_longitudes=False, **kwargs):
    Display an composite map from a list OMEGA/MEx observations, previously sampled on 
    a new lat/lon grid with show_omega_list_v2() and saved with save_map_omega_list().

    Parameters
    ==========
    data : 2D array
        The omega reflectance at lam, sampled on the new lat/lon grid.
    grid lat : 2D array
        The new latitude grid.
    grid lon : 2D array
        The new longitude grid.
    infos : dict
        The informations about the computation of the data.
    cmap : str, optional (default 'Greys_r')
        The matplotlib colormap.
    vmin : float or None, optional (default None)
        The lower bound of the coloscale.
    vmax : float or None, optional (default None)
        The upper bound of the colorscale.
    title : str, optional (default 'auto')
        The title of the figure.
    Nfig : int or str or None, optional (default None)
        The target figure ID.
    polar : bool, optional (default False)
        If True -> Use a polar projection for the plot.
    cbar : bool, optional (default True)
        If True -> Diplay the colorbar.
    cb_title : str, optional (default 'auto')
        The title of the colorbar.
    grid : bool, optional (default True)
        Enable the display of the lat/lon grid.
    negatives_longitudes : bool, optional (default False)
        Argument for non-polar plots.
        | True -> longitudes between 0° and 360°.
        | False -> longitudus between -180° and 180°.
    **kwargs:
        Optional arguments for the plt.pcolormesh() function.
~~~
