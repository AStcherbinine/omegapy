# OMEGA-Py documentation - v1.2

## `omegapy.omega_plots`

Display of OMEGAdata cubes.

`show_omega(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, title='auto', xlim=(None, None), ylim=(None, None), Nfig=None)`

`show_omega_v2(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, alpha=None, title='auto', lonlim=(None, None), latlim=(None, None), Nfig=None, polar=False, cbar=True)`

`show_omega_interactif(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, title='auto', autoyscale=True, xlim=(None, None), ylim=(None, None))`

`show_omega_interactif_v2(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, title='auto', autoyscale=True, alpha=None, lonlim=(None, None), latlim=(None, None), polar=False)`

`show_ibd_v2(omega, ibd, cmap='viridis', vmin=None, vmax=None, alpha=None, title='auto', cb_title = 'IBD', lonlim=(None, None), latlim=(None, None), Nfig=None, polar=False, cbar=True)`


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
    refl : bool (default True)
        True -> The reflectance is display.
        False -> The radiance is display.
    lam_unit : str (optional, default 'm')
        The unit of the `lam` parameter:
        | 'm' -> `lam` is the wavelength value (in µm).
        | else -> `lam` is the index of the wavelength in the omega.lam array (must be int).
    cmap : str (optional, default 'Greys_r')
        The matplotlib colormap.
    vmin : float or None (optional, default None)
        The lower bound of the coloscale.
    vmax : float or None (optional, default None)
        The upper bound of the colorscale.
    title : str (optional, default 'auto')
        The title of the figure.
    xlim : tuple of int or None (optional, default (None, None))
        The bounds of the x-axis of the figure.
    ylim : tuple of int or None (optional, default (None, None))
        The bounds of the y-axis of the figure.
    Nfig : int or str or None (optional, default None)
        The target figure ID.
~~~

~~~python
omegapy.omega_plots.show_omega_v2(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None,
                  alpha=None, title='auto', lonlim=(None, None), latlim=(None, None), Nfig=None,
                  polar=False, cbar=True):
    Display an OMEGA/MEx observation with respect of the lat/lon coordinates of the pixels,
    and allows to use a polar projection if desired.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA/MEx observation
    lam : float
        The selected wavelength.
    refl : bool (default True)
        True -> The reflectance is display.
        False -> The radiance is display.
    lam_unit : str (optional, default 'm')
        The unit of the `lam` parameter:
        | 'm' -> `lam` is the wavelength value (in µm).
        | else -> `lam` is the index of the wavelength in the omega.lam array (must be int).
    cmap : str (optional, default 'Greys_r')
        The matplotlib colormap.
    vmin : float or None (optional, default None)
        The lower bound of the coloscale.
    vmax : float or None (optional, default None)
        The upper bound of the colorscale.
    alpha : float or None (optional, default None)
        Opacity of the plot.
    title : str (optional, default 'auto')
        The title of the figure.
    lonlim : tuple of int or None (optional, default (None, None))
        The longitude bounds of the figure.
    latlim : tuple of int or None (optional, default (None, None))
        The latitude bounds of the y-axis of the figure.
    Nfig : int or str or None (optional, default None)
        The target figure ID.
    polar : bool (optional, default False)
        If True -> Use a polar projection for the plot.
    cbar : bool (optional, default True)
        If True -> Diplay the colorbar.
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
    refl : bool (default True)
        True -> The reflectance is display.
        False -> The radiance is display.
    lam_unit : str (optional, default 'm')
        The unit of the `lam` parameter:
        | 'm' -> `lam` is the wavelength value (in µm).
        | else -> `lam` is the index of the wavelength in the omega.lam array (must be int).
    cmap : str (optional, default 'Greys_r')
        The matplotlib colormap.
    vmin : float or None (optional, default None)
        The lower bound of the coloscale.
    vmax : float or None (optional, default None)
        The upper bound of the colorscale.
    title : str (optional, default 'auto')
        The title of the figure.
    xlim : tuple of int or None (optional, default (None, None))
        The bounds of the x-axis of the figure.
    ylim : tuple of int or None (optional, default (None, None))
        The bounds of the y-axis of the figure.
~~~

~~~python
omegapy.omega_plots.show_omega_interactif_v2(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', 
                          vmin=None, vmax=None, title='auto', autoyscale=True,
                          alpha=None, lonlim=(None, None), latlim=(None, None),
                          polar=False):
    Affichage interactif d'un cube de données.
    Possibilité d'afficher le spectre associé à un pixel en cliquant dessus
    (maintenir Ctrl pour supperposer plusieurs spectres), ou en se déplaçant avec les flèches.

    Display an OMEGA/MEx observation with respect of the lat/lon coordinates of the pixels,
    and allows to use a polar projection if desired.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA/MEx observation
    lam : float
        The selected wavelength.
    refl : bool (default True)
        True -> The reflectance is display.
        False -> The radiance is display.
    lam_unit : str (optional, default 'm')
        The unit of the `lam` parameter:
        | 'm' -> `lam` is the wavelength value (in µm).
        | else -> `lam` is the index of the wavelength in the omega.lam array (must be int).
    cmap : str (optional, default 'Greys_r')
        The matplotlib colormap.
    vmin : float or None (optional, default None)
        The lower bound of the coloscale.
    vmax : float or None (optional, default None)
        The upper bound of the colorscale.
    title : str (optional, default 'auto')
        The title of the figure.
    autoyscale : bool (optional, default True)
        | True -> Enable the auto-scaling of the spectra y-axis.
        | False -> Force use of the (vmin, vmax) bounds for the spectra plots.
    alpha : float or None (optional, default None)
        Opacity of the plot.
    lonlim : tuple of int or None (optional, default (None, None))
        The longitude bounds of the figure.
    latlim : tuple of int or None (optional, default (None, None))
        The latitude bounds of the y-axis of the figure.
    polar : bool (optional, default False)
        If True -> Use a polar projection for the plot.
~~~

### Display derived data map from OMEGA observation (ex IBD)

~~~python
omegapy.omega_plots.show_ibd_v2(omega, ibd, cmap='viridis', vmin=None, vmax=None, alpha=None, title='auto', 
                cb_title = 'IBD', lonlim=(None, None), latlim=(None, None), Nfig=None, 
                polar=False, cbar=True):
    Affichage IBD avec pcolormesh.
    Display an OMEGA/MEx observation with respect of the lat/lon coordinates of the pixels,
    and allows to use a polar projection if desired.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA/MEx observation
    ibd : 2D array
        The array of the computed IBD values from the omega observation
    cmap : str (optional, default 'Greys_r')
        The matplotlib colormap.
    vmin : float or None (optional, default None)
        The lower bound of the coloscale.
    vmax : float or None (optional, default None)
        The upper bound of the colorscale.
    alpha : float or None (optional, default None)
        Opacity of the plot.
    title : str (optional, default 'auto')
        The title of the figure.
    cb_title : str (optional, default 'IBD')
        The title of the colorbar.
    lonlim : tuple of int or None (optional, default (None, None))
        The longitude bounds of the figure.
    latlim : tuple of int or None (optional, default (None, None))
        The latitude bounds of the y-axis of the figure.
    Nfig : int or str or None (optional, default None)
        The target figure ID.
    polar : bool (optional, default False)
        If True -> Use a polar projection for the plot.
    cbar : bool (optional, default True)
        If True -> Display the colorbar.
~~~

### Display composite map of several OMEGA observations, sample on a lat/lon grid

~~~python
omegapy.omega_plots.show_omega_list_v2(omega_list, lam, lat_min=-90, lat_max=90, lon_min=0, lon_max=360,
                       pas_lat=0.1, pas_lon=0.1, cmap='Greys_r', vmin=None, vmax=None, 
                       title='auto', Nfig=None, polar=False, cbar=True, cb_title='auto',
                       data_list=None, plot=True, out=False, **kwargs):
    Display an composite map from a list OMEGA/MEx observations, sampled on a new lat/lon grid.

    Parameters
    ==========
    omega_list : array of OMEGAdata
        The list of OMEGA/MEx observations.
    lam : float
        The selected wavelength (in µm).
    lat_min : float (optional, default -90)
        The minimal latitude of the grid.
    lat_max : float (optional, default 90)
        The maximum latitude of the grid.
    lon_min : float (optional, default 0)
        The minimal longitude of the grid.
    lon_max : float (optional, default 360)
        The maximal longitude of the grid.
    pas_lat : float (optional, default 0.1)
        The latitude intervals of the grid.
    pas_lon : float (optional, default 0.1)
        The longitude intervals of the grid.
    cmap : str (optional, default 'Greys_r')
        The matplotlib colormap.
    vmin : float or None (optional, default None)
        The lower bound of the coloscale.
    vmax : float or None (optional, default None)
        The upper bound of the colorscale.
    title : str (optional, default 'auto')
        The title of the figure.
    Nfig : int or str or None (optional, default None)
        The target figure ID.
    polar : bool (optional, default False)
        If True -> Use a polar projection for the plot.
    cbar : bool (optional, default True)
        If True -> Diplay the colorbar.
    cb_title : str (optional, default 'auto')
        The title of the colorbar
    data_list : 3D array or None (optional, default None)
        1D array of the same dimension of `omega_list` containing alternative maps (2D arrays),
        in the **same order** than the observations of `omega_list`.
    plot : bool (optional, default True)
        If True -> Diplay the final figure.
    out : bool (optional, default False)
        If True -> Return output.
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
~~~

