#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## omega_plots.py
## Created by Aurélien STCHERBININE
## Last modified by Aurélien STCHERBININE : 20/04/2022

##----------------------------------------------------------------------------------------
"""Display of OMEGAdata cubes.
"""
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
## Packages
# Global
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from copy import deepcopy
from tqdm import tqdm
import datetime
import os
# Local
from . import useful_functions as uf
from . import omega_data as od
from .omega_data import OMEGAdata

# if __name__ == '__main__':
# Activation of the interactive mode
plt.ion()
# De-activation of default matplotlib keybindings for keyboard arrows
if 'left' in mpl.rcParams['keymap.back']:
    mpl.rcParams['keymap.back'].remove('left')
if 'right' in mpl.rcParams['keymap.forward']:
    mpl.rcParams['keymap.forward'].remove('right')

# Name of the current file
_py_file = 'omega_plots.py'

##----------------------------------------------------------------------------------------
## Initialisation variables globales
picked_spectra = {}

##----------------------------------------------------------------------------------------
## Affichage cube
def show_cube(cube, i_lam, cmap='Greys_r', vmin=None, vmax=None, cb_title='', Nfig=None):
    """Display the cube from an OMEGA/MEx observation.

    Parameters
    ==========
    cube : 3D array
        The data cube (X,Y,wvl).
    i_lam : int
        The index of the selected wavelength.
    cmap : str, optional (default 'Greys_r')
        The matplotlib colormap.
    vmin : float or None, optional (default None)
        The lower bound of the coloscale.
    vmax : float or None, optional (default None)
        The upper bound of the colorscale.
    cb_title : str, optional (default '')
        The title of the colorbar.
    Nfig : int or str or None, optional (default None)
        The target figure ID.
    """
    fig = plt.figure(Nfig)
    plt.imshow(cube[:,:,i_lam], cmap=cmap, vmin=vmin, vmax=vmax,
               aspect='equal', origin='lower', interpolation=None)
    cb = plt.colorbar()
    cb.set_label(cb_title)
    plt.tight_layout()

def show_omega(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None,
               title='auto', xlim=(None, None), ylim=(None, None), Nfig=None, mask=None):
    """Display an OMEGA/MEx observation in a rectangular pixel grid.

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
    mask : 2D array or None, optional (default None)
        The array that identify the bad/corrupted pixels to remove.
        If None, all the pixels are conserved.
        | 1 -> Good pixel
        | NaN -> Bad pixel
    """
    if ((lam_unit == 'm') or isinstance(lam, float)) and (lam < 10):
        i_lam = uf.where_closer(lam, omega.lam)
    else:
        i_lam = deepcopy(lam)
    lam = omega.lam[i_lam]
    if refl:
        cube = deepcopy(omega.cube_rf)
        cb_title = r'Reflectance @ $\lambda$' + ' = {0:.2f} µm'.format(lam)
    else:
        cube = deepcopy(omega.cube_i)
        cb_title = (r'Radiance [W.m$^{-2}$.sr$^{-1}$.µm$^{-1}$] at $\lambda$' + 
                    ' = {0:.2f} µm'.format(lam))
    if not (mask is None):
        cube = (cube.T * mask.T).T      # apply mask to remove bad pixels (turned to NaN)
    if title == 'auto':
        title = 'OMEGA/MEx observation {0}\n'.format(omega.name) 
    show_cube(cube, i_lam, cmap, vmin, vmax, cb_title, Nfig)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.title(title)
    plt.tight_layout()

def show_omega_v2(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None,
                  alpha=None, title='auto', lonlim=(None, None), latlim=(None, None), Nfig=None,
                  polar=False, cbar=True, grid=True, mask=None, negatives_longitudes='auto',
                  **kwargs):
    """Display an OMEGA/MEx observation with respect of the lat/lon coordinates of the pixels,
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
    mask : 2D array or None, optional (default None)
        The array that identify the bad/corrupted pixels to remove.
        If None, all the pixels are conserved.
        | 1 -> Good pixel
        | NaN -> Bad pixel
    negatives_longitudes : str or bool, optional (default 'auto')
        Argument for non-polar plots.
        | True -> longitudes between 0° and 360°.
        | False -> longitudus between -180° and 180°.
        | 'auto' -> automatic detection of the best case.
    **kwargs:
        Optional arguments for the plt.pcolormesh() function.
    """
    if ((lam_unit == 'm') or isinstance(lam, float)) and (lam < 10):
        i_lam = uf.where_closer(lam, omega.lam)
    else:
        i_lam = deepcopy(lam)
    lam = omega.lam[i_lam]
    if refl:
        cube = deepcopy(omega.cube_rf)
        cb_title = r'Reflectance @ $\lambda$' + ' = {0:.2f} µm'.format(lam)
    else:
        cube = deepcopy(omega.cube_i)
        cb_title = (r'Radiance [W.m$^{-2}$.sr$^{-1}$.µm$^{-1}$] at $\lambda$' + 
                    ' = {0:.2f} µm'.format(lam))
    cube_map = cube[:, :, i_lam]    # extracted map to display
    if not (mask is None):
        cube_map *= mask    # apply mask to remove bad pixels (turned to NaN)
    if title == 'auto':
        title = 'OMEGA/MEx observation {0}'.format(omega.name) 
    if isinstance(negatives_longitudes, str):
        mask_lat = (np.abs(omega.lat) < 85)
        if (omega.lon[mask_lat] < 10).any() and (omega.lon[mask_lat] > 350).any():
            negatives_longitudes = True
    fig = plt.figure(Nfig)
    Nfig = fig.number   # get the actual figure number if Nfig=None
    if len(fig.get_axes()) != 0:    # If presence of axes
        ax0 = fig.get_axes()[0]
        is_ax0_polar = hasattr(ax0, 'set_theta_offset') # Test if ax has polar projection
        if not polar == is_ax0_polar:
            raise ValueError("Can not mix polar and non-polar projections in the same plot")
    if polar:
        if len(fig.get_axes()) == 0:    # Test presence of axes in the figure
            ax = plt.axes(polar=True)
        else:
            ax = fig.get_axes()[0]  # Do not create new axes instance
        plt.pcolormesh(omega.lon_grid*np.pi/180, omega.lat_grid, cube_map, cmap=cmap, 
                       alpha=alpha, vmin=vmin, vmax=vmax, **kwargs)
        ax.set_yticklabels([])  # remove the latitude values in the plot
        if latlim[0] is None:
            if np.max(omega.lat) > 0:
                latlim = (90, np.min(omega.lat_grid)-1)
            else:
                latlim = (-90, np.max(omega.lat_grid)+1)
        if latlim[0] > 0:   # Northern hemisphere
            ax.set_theta_offset(-np.pi/2)   # longitude origin at the bottom
        else:               # Southern hemisphere
            ax.set_theta_offset(np.pi/2)    # longitude origin at the top
            ax.set_theta_direction(-1)      # clockwise theta
        plt.xlim(lonlim)
        plt.ylim(latlim)
    else:
        lon_grid2 = deepcopy(omega.lon_grid)
        if negatives_longitudes:
            lon_grid2[lon_grid2 > 180] -= 360
        plt.pcolormesh(lon_grid2, omega.lat_grid, cube_map, cmap=cmap, alpha=alpha,
                       vmin=vmin, vmax=vmax, **kwargs)
        plt.gca().axis('equal')
        plt.xlim(lonlim)
        plt.ylim(latlim)
        plt.gca().set_adjustable('box')
        plt.xlabel('Longitude [°]')
        plt.ylabel('Latitude [°]')
    if cbar:
        cb = plt.colorbar()
        cb.set_label(cb_title)
    if grid:
        ax = plt.figure(Nfig).get_axes()[0]
        lonlim = ax.get_xlim()
        latlim = ax.get_ylim()
        lon_sgn = np.sign(lonlim[1] - lonlim[0])
        lat_sgn = np.sign(latlim[1] - latlim[0])
        lon_grid = np.arange(np.round(lonlim[0]/10)*10, np.round(lonlim[1]/10)*10+lon_sgn, 
                    10 * lon_sgn)   # 10° grid in longitude
        lat_grid = np.arange(np.round(latlim[0]/10)*10, np.round(latlim[1]/10)*10+lat_sgn, 
                    10 * lat_sgn)   # 10° grid in latitude
        plt.grid()
        if polar:
            ax.set_rticks(lat_grid)
        else:
            ax.set_xticks(lon_grid)
            ax.set_yticks(lat_grid)
    plt.title(title)
    plt.tight_layout()

##----------------------------------------------------------------------------------------
## Affichage cube interactif
def show_omega_interactif(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', 
                          vmin=None, vmax=None, title='auto', autoyscale=True,
                          xlim=(None, None), ylim=(None, None)):
    """Affichage interactif d'un cube de données.
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
    """
    # Initialisation
    if refl:
        yaxis = 'Reflectance'
        cube = omega.cube_rf
    else:
        yaxis = r'Radiance [W.m$^{-2}$.sr$^{-1}$.µm$^{-1}$]'
        cube = omega.cube_i
    ny, nx, nlam = cube.shape
    xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
    fig1, ax1 = plt.subplots(1,1)
    nfig = fig1.number
    ax1.scatter(xx, yy, marker='s', s=1, picker=True, alpha=0)
    show_omega(omega, lam, refl, lam_unit, cmap, vmin, vmax, title, 
               xlim, ylim, nfig)
    sc_pos = []
    if xlim[0] is None:
        xcoord = 0
    else:
        xcoord = deepcopy(xlim[0])
    if ylim[0] is None:
        ycoord = 0
    else:
        ycoord = deepcopy(ylim[0])
    
    #---------------------------------
    # Plot spectra fig2 function
    def plot_sp(xcoord, ycoord, clear=True):
        nonlocal sc_pos
        fig2 = plt.figure(-nfig)
        if clear:
            fig2.clf()
            for sc in sc_pos:
                sc.remove()
            sc_pos = []
        line = plt.plot(omega.lam, cube[ycoord, xcoord], 
                        label='lat = {0:.2f}° | lon= {1:.2f}°'.format(omega.lat[ycoord, xcoord], 
                                                                      omega.lon[ycoord, xcoord]))
        plt.xlabel(r'$\lambda$ [µm]')
        plt.ylabel(yaxis)
        plt.title('OMEGA/MEx observation {0}'.format(omega.name))
        plt.legend(loc='best')
        ymin, ymax = fig2.get_axes()[0].get_ylim()
        # Rescale ordonnées
        if autoyscale:
            if (vmin!=None) and (ymin > vmin):
                ymin = vmin
            if (vmax!=None) and (ymax < vmax):
                ymax = vmax
        else:
            ymin, ymax = vmin, vmax
        plt.ylim(ymin, ymax)
        fig2.canvas.draw()
        fig2.canvas.flush_events()
        fig2.tight_layout()
        last_plot_color = line[0].get_color()
        sc_pos.append(ax1.scatter(xcoord, ycoord, marker='s', s=20, color=last_plot_color))
        fig1.canvas.draw()
        fig1.canvas.flush_events()
    
    #---------------------------------
    # Picking function clic souris
    def pick_pos(event):
        nonlocal xcoord, ycoord
        ctrl = event.mouseevent.key == 'control'    # test si la touche ctrl est enfoncée
        artist = event.artist
        ind = event.ind[0]
        xcoord, ycoord = artist.get_offsets()[ind]
        xcoord, ycoord = int(xcoord), int(ycoord)
        # if not ctrl:        # Si Ctrl enfoncée, pas le plot précédent est conservé
            # plt.clf()
        plot_sp(xcoord, ycoord, not ctrl)

    #---------------------------------
    # Picking function keyboard
    def change_pos(event):
        nonlocal xcoord, ycoord
        key = event.key
        if (0 < xcoord) and (key=='left'):
            xcoord -= 1
            plot_sp(xcoord, ycoord, clear=True)
        elif (xcoord < nx-1) and (key=='right'):
            xcoord += 1
            plot_sp(xcoord, ycoord, clear=True)
        if (0 < ycoord) and (key=='down'):
            ycoord -= 1
            plot_sp(xcoord, ycoord, clear=True)
        elif (ycoord < ny-1) and (key=='up'):
            ycoord += 1
            plot_sp(xcoord, ycoord, clear=True)

    #---------------------------------
    # Lien avec la figure
    cid = fig1.canvas.mpl_connect('pick_event', pick_pos)
    cid2 = fig1.canvas.mpl_connect('key_press_event', change_pos)

def show_omega_interactif2(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', 
                          vmin=None, vmax=None, title='auto', 
                          xlim=(None, None), ylim=(None, None)):
    """Affichage interactif d'un cube de données.
    Possibilité d'afficher le spectre associé à un pixel en cliquant dessus
    (maintenir Ctrl pour supperposer plusieurs spectres), ou en se déplaçant avec les flèches.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA/MEx observation
    lam : float
        The selected wavelength.
    refl : bool optional (default True)
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
    """
    # Initialisation
    if refl:
        yaxis = 'Reflectance'
        cube = omega.cube_rf
    else:
        yaxis = r'Radiance [W.m$^{-2}$.sr$^{-1}$.µm$^{-1}$]'
        cube = omega.cube_i
    Lam = omega.lam
    ny, nx, nlam = cube.shape
    xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
    fig1, ax1 = plt.subplots(1,1)
    nfig = fig1.number
    ax1.scatter(xx, yy, marker='s', s=1, picker=True, alpha=0)
    show_omega(omega, lam, refl, lam_unit, cmap, vmin, vmax, title, 
               xlim, ylim, nfig)
    sc_pos = []
    if xlim[0] is None:
        xcoord = 0
    else:
        xcoord = deepcopy(xlim[0])
    if ylim[0] is None:
        ycoord = 0
    else:
        ycoord = deepcopy(ylim[0])
    # Sliders
    axcolor = 'lightgoldenrodyellow'
    axlam = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    slam = Slider(axlam, r'$\lambda$ [µm]', Lam[0], Lam[-1], valinit=lam, valfmt='%1.2f')

    def update_img(lam):
        # ax1.collections[0].remove()
        # print(lam)
        lam = slam.val
        if ((lam_unit == 'm') or isinstance(lam, float)) and (lam < 10):
            i_lam = uf.where_closer(lam, omega.lam)
        else:
            i_lam = deepcopy(lam)
        lam = omega.lam[i_lam]
        ax1.images[0].set_array(cube[:,:,i_lam])
        title = r'$\lambda$' + ' = {0:.2f} µm'.format(Lam[i_lam])
        ax1.set_title(title)
        fig1.canvas.draw_idle()
    slam.on_changed(update_img)
    
    #---------------------------------
    # Plot spectra fig2 function
    def plot_sp(xcoord, ycoord, clear=True):
        nonlocal sc_pos
        fig2 = plt.figure(-nfig)
        if clear:
            fig2.clf()
            for sc in sc_pos:
                sc.remove()
            sc_pos = []
        line = plt.plot(omega.lam, cube[ycoord, xcoord], 
                        label='lat = {0:.2f}° | lon= {1:.2f}°'.format(omega.lat[ycoord, xcoord], 
                                                                      omega.lon[ycoord, xcoord]))
        plt.xlabel(r'$\lambda$ [µm]')
        plt.ylabel(yaxis)
        plt.title('OMEGA/MEx observation {0}'.format(omega.name))
        plt.legend(loc='best')
        ymin, ymax = fig2.get_axes()[0].get_ylim()
        # Rescale ordonnées
        if (vmin!=None) and (ymin > vmin):
            ymin = vmin
        if (vmax!=None) and (ymax < vmax):
            ymax = vmax
        plt.ylim(ymin, ymax)
        fig2.canvas.draw()
        fig2.canvas.flush_events()
        fig2.tight_layout()
        last_plot_color = line[0].get_color()
        sc_pos.append(ax1.scatter(xcoord, ycoord, marker='s', s=20, color=last_plot_color))
        fig1.canvas.draw()
        fig1.canvas.flush_events()
    
    #---------------------------------
    # Picking function clic souris
    def pick_pos(event):
        nonlocal xcoord, ycoord
        ctrl = event.mouseevent.key == 'control'    # test si la touche ctrl est enfoncée
        artist = event.artist
        ind = event.ind[0]
        xcoord, ycoord = artist.get_offsets()[ind]
        xcoord, ycoord = int(xcoord), int(ycoord)
        # if not ctrl:        # Si Ctrl enfoncée, pas le plot précédent est conservé
            # plt.clf()
        plot_sp(xcoord, ycoord, not ctrl)

    #---------------------------------
    # Picking function keyboard
    def change_pos(event):
        nonlocal xcoord, ycoord
        key = event.key
        if (0 < xcoord) and (key=='left'):
            xcoord -= 1
            plot_sp(xcoord, ycoord, clear=True)
        elif (xcoord < nx-1) and (key=='right'):
            xcoord += 1
            plot_sp(xcoord, ycoord, clear=True)
        if (0 < ycoord) and (key=='down'):
            ycoord -= 1
            plot_sp(xcoord, ycoord, clear=True)
        elif (ycoord < ny-1) and (key=='up'):
            ycoord += 1
            plot_sp(xcoord, ycoord, clear=True)
    
    #---------------------------------
    # Lien avec la figure
    cid = fig1.canvas.mpl_connect('pick_event', pick_pos)
    cid2 = fig1.canvas.mpl_connect('key_press_event', change_pos)

def show_omega_interactif_v2(omega, lam=1.085, refl=True, lam_unit='m', data=None, 
                             cmap='Greys_r', cb_title='data', title='auto',
                             vmin=None, vmax=None, autoyscale=True, ylim_sp=(None, None),
                             alpha=None, lonlim=(None, None), latlim=(None, None),
                             polar=False, cbar=True, grid=True, mask=None, lam_mask=None,
                             negatives_longitudes='auto', **kwargs):
    """Affichage interactif d'un cube de données.
    Possibilité d'afficher le spectre associé à un pixel en cliquant dessus
    (maintenir Ctrl pour supperposer plusieurs spectres), ou en se déplaçant avec les flèches.
    Les spectres affichés sont stockés dans le dictionnaire `picked_spectra[nfig]`.

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
    cb_title : str,  optional (default 'data')
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
    mask : 2D array or None, optional (default None)
        The array that identify the bad/corrupted pixels to remove.
        If None, all the pixels are conserved.
        | 1 -> Good pixel
        | NaN -> Bad pixel
    lam_mask : 1D array or None, optional (default None)
        The array that identify the bad/corrupted spectels to remove.
        If None, all the spectels are conserved.
        | True -> Good spectel
        | False -> Bad spectel
    negatives_longitudes : str or bool, optional (default 'auto')
        Argument for non-polar plots.
        | True -> longitudes between 0° and 360°.
        | False -> longitudus between -180° and 180°.
        | 'auto' -> automatic detection of the best case.
    **kwargs:
        Optional arguments for the plt.pcolormesh() function.
    """
    # Initialisation
    if refl:
        yaxis = 'Reflectance'
        cube = omega.cube_rf
    else:
        yaxis = r'Radiance [W.m$^{-2}$.sr$^{-1}$.µm$^{-1}$]'
        cube = omega.cube_i
    ny, nx, nlam = cube.shape
    if isinstance(negatives_longitudes, str):
        mask_lat = (np.abs(omega.lat) < 85)
        if (omega.lon[mask_lat] < 10).any() and (omega.lon[mask_lat] > 350).any():
            negatives_longitudes = True
    if polar:
        lon, lat = omega.lon*np.pi/180, omega.lat
    else:
        lon = deepcopy(omega.lon)
        if negatives_longitudes:
            lon[lon > 180] -= 360
        # lon, lat = omega.lon, omega.lat
        lat = omega.lat
    bij = np.zeros((ny, nx), dtype=int)
    for j in range(ny):
        for i in range(nx):
            bij[j,i] = 10000*j + i
    if lam_mask is None:
        lam_mask = np.ones(len(omega.lam), dtype=bool)
    elif len(lam_mask) != len(omega.lam):
        raise ValueError('omega.lam and lam_mask must have the same dimension')
    lam2 = deepcopy(omega.lam)[lam_mask]
    #---------------------------------
    # Display map
    fig1 = plt.figure()
    nfig = fig1.number
    if data is None:
        show_omega_v2(omega, lam, refl, lam_unit, cmap, vmin, vmax, alpha, title, 
                      lonlim, latlim, nfig, polar, cbar, grid, mask, negatives_longitudes, 
                      **kwargs)
    else:
        show_data_v2(omega, data, cmap, vmin, vmax, alpha, title, cb_title, 
                     lonlim, latlim, nfig, polar, cbar, grid, mask, negatives_longitudes,
                     **kwargs)
    ax1 = fig1.gca()
    ax1.scatter(lon, lat, c=bij, marker='s', s=1, picker=True, alpha=0)
    sc_pos = []
    if lonlim[0] is None:
        xcoord = 0
    else:
        xcoord = deepcopy(lonlim[0])
    if latlim[0] is None:
        ycoord = 0
    else:
        ycoord = deepcopy(latlim[0])
    ylim_sp = np.array(ylim_sp)
    if ylim_sp[0] is None:
        ylim_sp[0] = vmin
    if ylim_sp[1] is None:
        ylim_sp[1] = vmax

    picked_spectra[nfig] = [lam2]
    
    #---------------------------------
    # Plot spectra fig2 function
    def plot_sp(xcoord, ycoord, clear=True):
        nonlocal sc_pos
        global picked_spectra
        fig2 = plt.figure(-nfig)
        if clear:
            fig2.clf()
            picked_spectra[nfig] = [lam2]
            for sc in sc_pos:
                sc.remove()
            sc_pos = []
        sp_i = cube[ycoord, xcoord, lam_mask]
        lati = omega.lat[ycoord, xcoord]
        longi = omega.lon[ycoord, xcoord]
        picked_spectra[nfig].append(sp_i)   # Stockage spectres dans variable globale
        line = plt.plot(lam2, sp_i, 
                label='lat = {0:.2f}°N | lon = {1:.2f}°E | pixel coord = ({2:d}, {3:d})'.format(
                                                lati, longi, ycoord, xcoord))
        plt.xlabel(r'$\lambda$ [µm]')
        plt.ylabel(yaxis)
        plt.title('OMEGA/MEx observation {0}'.format(omega.name))
        plt.legend(loc='best')
        ymin, ymax = fig2.get_axes()[0].get_ylim()
        # Rescale ordonnées
        if autoyscale:
            if (vmin!=None) and (ymin > vmin):
                ymin = vmin
            if (vmax!=None) and (ymax < vmax):
                ymax = vmax
        else:
            ymin, ymax = ylim_sp[0], ylim_sp[1]
        plt.ylim(ymin, ymax)
        fig2.canvas.draw()
        fig2.canvas.flush_events()
        fig2.tight_layout()
        longi2 = lon[ycoord, xcoord]    # longitude adaptée si projection polaire
        last_plot_color = line[0].get_color()
        sc_pos.append(ax1.scatter(longi2, lati, marker='s', s=20, color=last_plot_color))
        fig1.canvas.draw()
        fig1.canvas.flush_events()
    
    #---------------------------------
    # Picking function clic souris
    def pick_pos(event):
        nonlocal xcoord, ycoord
        ctrl = event.mouseevent.key == 'control'    # test si la touche ctrl est enfoncée
        artist = event.artist
        ind = event.ind[0]
        bij_value = artist.get_array()[ind]
        xcoord = int(bij_value % 10000)
        ycoord = int(bij_value // 10000)
        plot_sp(xcoord, ycoord, not ctrl)

    #---------------------------------
    # Picking function keyboard
    def change_pos(event):
        nonlocal xcoord, ycoord
        key = event.key
        if (0 < xcoord) and (key=='left'):
            xcoord -= 1
            plot_sp(xcoord, ycoord, clear=True)
        elif (xcoord < nx-1) and (key=='right'):
            xcoord += 1
            plot_sp(xcoord, ycoord, clear=True)
        if (0 < ycoord) and (key=='down'):
            ycoord -= 1
            plot_sp(xcoord, ycoord, clear=True)
        elif (ycoord < ny-1) and (key=='up'):
            ycoord += 1
            plot_sp(xcoord, ycoord, clear=True)

    #---------------------------------
    # Lien avec la figure
    cid = fig1.canvas.mpl_connect('pick_event', pick_pos)
    cid2 = fig1.canvas.mpl_connect('key_press_event', change_pos)

##----------------------------------------------------------------------------------------
## Affichage données haut-niveau
def show_data_v2(omega, data, cmap='viridis', vmin=None, vmax=None, alpha=None, title='auto', 
                cb_title = 'data', lonlim=(None, None), latlim=(None, None), Nfig=None, 
                polar=False, cbar=True, grid=True, mask=None, negatives_longitudes='auto',
                **kwargs):
    """Affichage données haut-niveau avec pcolormesh.
    Display an OMEGA/MEx observation with respect of the lat/lon coordinates of the pixels,
    and allows to use a polar projection if desired.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA/MEx observation
    data : 2D array
        The array of the computed data values from the omega observation
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
    cb_title : str,  optional (default 'data')
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
    mask : 2D array or None, optional (default None)
        The array that identify the bad/corrupted pixels to remove.
        If None, all the pixels are conserved.
        | 1 -> Good pixel
        | NaN -> Bad pixel
    negatives_longitudes : str or bool, optional (default 'auto')
        Argument for non-polar plots.
        | True -> longitudes between 0° and 360°.
        | False -> longitudus between -180° and 180°.
        | 'auto' -> automatic detection of the best case.
    **kwargs:
        Optional arguments for the plt.pcolormesh() function.
    """
    if isinstance(negatives_longitudes, str):
        mask_lat = (np.abs(omega.lat) < 85)
        if (omega.lon[mask_lat] < 10).any() and (omega.lon[mask_lat] > 350).any():
            negatives_longitudes = True
    if title == 'auto':
        title = ('OMEGA/MEx observation {0}'.format(omega.name))
    fig = plt.figure(Nfig)
    Nfig = fig.number   # get the actual figure number if Nfig=None
    if not (mask is None):
        data = deepcopy(data) * mask     # apply mask to remove bad pixels (turned to NaN)
    if polar:
        ax = plt.axes(polar=True)
        plt.pcolormesh(omega.lon_grid*np.pi/180, omega.lat_grid, data, cmap=cmap, 
                       alpha=alpha, vmin=vmin, vmax=vmax, **kwargs)
        ax.set_yticklabels([])  # remove the latitude values in the plot
        if latlim[0] is None:
            if np.max(omega.lat) > 0:
                latlim = (90, np.min(omega.lat_grid)-1)
            else:
                latlim = (-90, np.max(omega.lat_grid)+1)
        if latlim[0] > 0:   # Northern hemisphere
            ax.set_theta_offset(-np.pi/2)   # longitude origin at the bottom
        else:               # Southern hemisphere
            ax.set_theta_offset(np.pi/2)    # longitude origin at the top
            ax.set_theta_direction(-1)      # clockwise theta
        plt.xlim(lonlim)
        plt.ylim(latlim)
    else:
        lon_grid2 = deepcopy(omega.lon_grid)
        if negatives_longitudes:
            lon_grid2[lon_grid2 > 180] -= 360
        plt.pcolormesh(lon_grid2, omega.lat_grid, data, cmap=cmap, alpha=alpha,
                       vmin=vmin, vmax=vmax, **kwargs)
        plt.gca().axis('equal')
        plt.xlim(lonlim)
        plt.ylim(latlim)
        plt.gca().set_adjustable('box')
        plt.xlabel('Longitude [°]')
        plt.ylabel('Latitude [°]')
    if cbar:
        cb = plt.colorbar()
        cb.set_label(cb_title)
    if grid:
        ax = plt.figure(Nfig).get_axes()[0]
        lonlim = ax.get_xlim()
        latlim = ax.get_ylim()
        lon_sgn = np.sign(lonlim[1] - lonlim[0])
        lat_sgn = np.sign(latlim[1] - latlim[0])
        lon_grid = np.arange(np.round(lonlim[0]/10)*10, np.round(lonlim[1]/10)*10+lon_sgn, 
                    10 * lon_sgn)   # 10° grid in longitude
        lat_grid = np.arange(np.round(latlim[0]/10)*10, np.round(latlim[1]/10)*10+lat_sgn, 
                    10 * lat_sgn)   # 10° grid in latitude
        plt.grid()
        if polar:
            ax.set_rticks(lat_grid)
        else:
            ax.set_xticks(lon_grid)
            ax.set_yticks(lat_grid)
    plt.title(title)
    plt.tight_layout()

##----------------------------------------------------------------------------------------
## Projection grille
def proj_grid(omega, data, lat_min=-90, lat_max=90, lon_min=0, lon_max=360,
              pas_lat=0.1, pas_lon=0.1, negative_values=False):
    """Sample the data from the input OMEGA/MEx observation on a given lat/lon grid.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA/MEx observation
    data : 2D array
        The initial array of values associated to the OMEGAdata observation.
        e.g.: Refelectance at selected wvl, spectra, or derived data such as IBD map.
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
    negative_values : bool, optional (default False)
        Set if the negative values are considered as relevant data or not.

    Returns
    =======
    grid_data : 2D array (dim : Nlon x Nlat)
        The data values, sampled on the new lat/lon grid.
    mask : 2D array
        The array indicating where the new grid has been filled by the OMEGA data.
    grid lat : 2D array
        The new latitude grid.
    grid lon : 2D array
        The new longitude grid.
    """
    # Initialisation
    lat_array = np.arange(lat_min, lat_max+pas_lat, pas_lat)
    lon_array = np.arange(lon_min, lon_max+pas_lon, pas_lon)
    Nlon, Nlat = len(lon_array)-1, len(lat_array)-1
    grid_lat, grid_lon = np.meshgrid(lat_array, lon_array)
    grid_data = np.zeros((Nlon, Nlat))
    mask = np.zeros((Nlon, Nlat))
    lat2 = np.floor(np.clip(omega.lat, lat_min, lat_max-0.1*pas_lat) /pas_lat) * pas_lat
    lon2 = np.floor(np.clip(omega.lon, lon_min, lon_max-0.1*pas_lon) /pas_lon) * pas_lon
    nx, ny = lat2.shape
    # Sampling on the new grid
    for j in range(ny):
        for i in range(nx):
            longi, lati = lon2[i,j], lat2[i,j]
            i_lon = int(longi/pas_lon - lon_min/pas_lon)
            j_lat = int(lati/pas_lat - lat_min/pas_lat)
            data_tmp = data[i,j]
            # if (not np.isnan(data_tmp)) & (data_tmp > 0):   # Filtrage régions sans données
            if negative_values:
                if (not np.isnan(data_tmp)) & (not np.isinf(data_tmp)):   # Filtrage régions sans données
                    grid_data[i_lon,j_lat] += data_tmp
                    mask[i_lon,j_lat] += 1
            else:   # negative values = No data
                if (not np.isnan(data_tmp)) & (data_tmp > 0) & (not np.isinf(data_tmp)):   # Filtrage régions sans données
                    grid_data[i_lon,j_lat] += data_tmp
                    mask[i_lon,j_lat] += 1
    grid_data[grid_data==0] = np.nan
    grid_data2 = grid_data / mask       # Normalisation
    mask2 = np.clip(mask, 0, 1)
    return grid_data2, mask2, grid_lat, grid_lon

def check_list_data_omega(omega_list, data_list, disp=True):
    """Check the compatibility between data_list and the list of OMEGA/MEx observations.
    Raise ValueError if uncompatibility.

    Parameters
    ==========
    omega_list : array of OMEGAdata
        List of OMEGA/MEx observations.
    data_list : 3D array
        List of high-level map associated to the observations of omega_list.
    disp : bool, optional (default True)
        Enable the display of the result of the test.
    """
    if len(omega_list) != len(data_list):
        raise ValueError("omega_list and data_list must have the same size")
    else:
        for i in range(len(omega_list)):
            if omega_list[i].lat.shape != data_list[i].shape:
                raise ValueError("The shapes of items {0} of omega_list and data_list does not match.".format(i))
    if disp:
        print("\033[01;32mCompatibility between omega_list and data_list OK\033[0m")

def check_list_mask_omega(omega_list, mask_list, disp=True):
    """Check the compatibility between mask_list and the list of OMEGA/MEx observations.
    Raise ValueError if uncompatibility.

    Parameters
    ==========
    omega_list : array of OMEGAdata
        List of OMEGA/MEx observations.
    mask_list : 3D array
        List of masks to remove the corrupted pixels of each OMEGA/MEx observation.
    disp : bool, optional (default True)
        Enable the display of the result of the test.
    """
    if len(omega_list) != len(mask_list):
        raise ValueError("omega_list and mask_list must have the same size")
    else:
        for i in range(len(omega_list)):
            if omega_list[i].lat.shape != mask_list[i].shape:
                raise ValueError("The shapes of items {0} of omega_list and mask_list does not match.".format(i))
    if disp:
        print("\033[01;32mCompatibility between omega_list and mask_list OK\033[0m")

def show_omega_list_v2(omega_list, lam=1.085, lat_min=-90, lat_max=90, lon_min=0, lon_max=360,
                       pas_lat=0.1, pas_lon=0.1, cmap='Greys_r', vmin=None, vmax=None, 
                       title='auto', Nfig=None, polar=False, cbar=True, cb_title='auto',
                       data_list=None, mask_list=None, negative_values=False, plot=True, 
                       grid=True, out=False, negatives_longitudes=False, 
                       edgecolor='face', lw=0.1, **kwargs):
    """Display an composite map from a list OMEGA/MEx observations, sampled on a new lat/lon grid.

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
    negative_values : bool, optional (default False)
        Set if the negative values are considered as relevant data or not.
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
    edgecolor : {'none', None, 'face', color', color sequence}, optional (default 'face')
        The color of the edges, see documentation of plt.pcolormesh for more details.
        > Added in version 2.2.8 to fix display due for new version of matplotlib.
    lw : float, optional (default 0.1)
        The line width of the edges (if diplayed).
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
    """
    # Sampling on same grid
    lat_array = np.arange(lat_min, lat_max+pas_lat, pas_lat)
    lon_array = np.arange(lon_min, lon_max+pas_lon, pas_lon)
    Nlon, Nlat = len(lon_array)-1, len(lat_array)-1
    grid_lat, grid_lon = np.meshgrid(lat_array, lon_array)
    data, mask = np.zeros((Nlon, Nlat)), np.zeros((Nlon, Nlat))
    mask_obs = np.ndarray((Nlon, Nlat), dtype=object)
    mask_obs.fill('')
    if not (mask_list is None):
        check_list_mask_omega(omega_list, mask_list, disp=True)
    if data_list is None:
        for i, omega in enumerate(tqdm(omega_list)):
            i_lam = uf.where_closer(lam, omega.lam)
            if mask_list is None:
                data_tmp = omega.cube_rf[:,:,i_lam]     # Reflectance without mask
            else:
                data_tmp = omega.cube_rf[:,:,i_lam] * mask_list[i]  # Reflectance with mask
            data0, mask0 = proj_grid(omega, data_tmp, lat_min, lat_max,
                                     lon_min, lon_max, pas_lat, pas_lon, negative_values)[:2]
            data += np.nan_to_num(data0)    # Conversion NaN -> 0 pour somme des images
            mask += mask0
            mask_obs[mask0 == 1] += (omega.name + ',')
    else:
        check_list_data_omega(omega_list, data_list, disp=True)
        for i, omega in enumerate(tqdm(omega_list)):
            if mask_list is None:
                data_tmp = data_list[i]     # Data without mask
            else:
                data_tmp = data_list[i] * mask_list[i]  # Data with mask
            data0, mask0 = proj_grid(omega, data_tmp, lat_min, lat_max,
                                     lon_min, lon_max, pas_lat, pas_lon, negative_values)[:2]
            data += np.nan_to_num(data0)    # Conversion NaN -> 0 pour somme des images
            mask += mask0
            mask_obs[mask0 == 1] += (omega.name + ',')
    data[mask == 0] = np.nan
    data2 = data/mask   # Normalisation
    # Affichage figure
    if plot:
        if title == 'auto':
            title = 'Composite map from OMEGA/MEx observations' 
        fig = plt.figure(Nfig)
        Nfig = fig.number   # get the actual figure number if Nfig=None
        if polar:
            ax = plt.axes(polar=True)
            plt.pcolormesh(grid_lon*np.pi/180, grid_lat, data2, cmap=cmap, 
                        vmin=vmin, vmax=vmax, edgecolor=edgecolor, lw=lw, **kwargs)
            ax.set_yticklabels([])  # remove the latitude values in the plot
            plt.xlim(0, 2*np.pi)
            if np.abs(lat_max) >= np.abs(lat_min):
                latlim = (lat_max, lat_min)
            else:
                latlim = (lat_min, lat_max)
            if latlim[0] > 0:   # Northern hemisphere
                ax.set_theta_offset(-np.pi/2)   # longitude origin at the bottom
            else:               # Southern hemisphere
                ax.set_theta_offset(np.pi/2)    # longitude origin at the top
                ax.set_theta_direction(-1)      # clockwise theta
            plt.ylim(latlim)
        else:
            if negatives_longitudes and (lon_max > 180):
                n_neg_lon = np.sum(grid_lon[:,0] > 180) # nb of negative longitudes (>180°)
                i_lon180 = np.where(grid_lon[:,0] > 180)[0][0] # 1st index of lon > 180°
                grid_lon_nl = deepcopy(grid_lon)        # new longitude grid [-180°, 180°]
                grid_lon_nl[:n_neg_lon] = grid_lon[i_lon180-1:-1] - 360
                grid_lon_nl[n_neg_lon:] = grid_lon[:i_lon180]
                data2_nl = deepcopy(data2)      # new data array
                data2_nl[:n_neg_lon] = data2[i_lon180-1:]
                data2_nl[n_neg_lon:] = data2[:i_lon180-1]
                plt.pcolormesh(grid_lon_nl, grid_lat, data2_nl, cmap=cmap, vmin=vmin, 
                               vmax=vmax, edgecolor=edgecolor, lw=lw, **kwargs)
                lon_min, lon_max = grid_lon_nl[[0,-1], 0]   # new longitude bounds
            else:
                plt.pcolormesh(grid_lon, grid_lat, data2, cmap=cmap, vmin=vmin, 
                               vmax=vmax, edgecolor=edgecolor, lw=lw, **kwargs)
            plt.gca().axis('equal')
            plt.gca().set_adjustable('box')
            plt.xlabel('Longitude [°]')
            plt.ylabel('Latitude [°]')
            plt.xlim(lon_min, lon_max)
            plt.ylim(lat_min, lat_max)
        if cbar:
            if cb_title == 'auto':
                cb_title = r'Reflectance @ $\lambda$' + ' = {0:.2f} µm'.format(lam)
            cb = plt.colorbar()
            cb.set_label(cb_title)
        if grid:
            ax = plt.figure(Nfig).get_axes()[0]
            lonlim = ax.get_xlim()
            latlim = ax.get_ylim()
            lon_sgn = np.sign(lonlim[1] - lonlim[0])
            lat_sgn = np.sign(latlim[1] - latlim[0])
            lon_grid = np.arange(np.round(lonlim[0]/10)*10, np.round(lonlim[1]/10)*10+lon_sgn, 
                        10 * lon_sgn)   # 10° grid in longitude
            lat_grid = np.arange(np.round(latlim[0]/10)*10, np.round(latlim[1]/10)*10+lat_sgn, 
                        10 * lat_sgn)   # 10° grid in latitude
            plt.grid()
            if polar:
                ax.set_rticks(lat_grid)
            else:
                ax.set_xticks(lon_grid)
                ax.set_yticks(lat_grid)
        plt.title(title)
        plt.tight_layout()
    # Output
    if out:
        mask2 = np.clip(mask, 0, 1)
        return data2, mask2, grid_lat, grid_lon, mask_obs

##----------------------------------------------------------------------------------------
## Sauvegarde résultats
def save_map_omega_list(omega_list, lat_min=-90, lat_max=90, lon_min=0, lon_max=360,
                        pas_lat=0.1, pas_lon=0.1, lam=1.085, data_list=None, data_desc='', 
                        mask_list=None, negative_values=False, sav_filename='auto', ext='',
                        base_folder='../data/OMEGA/sav_map_list_v2/', sub_folder=''):
    """Save the output of the omega_plots.show_omega_list_v2() function with the requested
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
    negative_values : bool, optional (default False)
        Set if the negative values are considered as relevant data or not.
    sav_filename : str, optional (default 'auto')
        The saving file name.
        | If 'auto' -> Automatically generated.
    ext : str, optional (default '')
        Extension to add at the end of the filename (useful in case of automatic generation).
    base_folder : str, optional (default '../data/OMEGA/sav_map_list_v2/')
        The base folder to save the data.
    sub_folder : str, optional (default '')
        The subfolder to save the data.
        Final path = "base_folder / sub_folder / sav_filename"
    """
    # Initialization filename
    if sav_filename == 'auto':
        sav_filename = ('res_show_omega_list_v2__lat{0:0>2d}-{1:0>2d}_pas{2:0}_'
                        + 'lon{3:0>3d}-{4:0>3d}_pas{5:0}__{6:s}.pkl').format(
                            lat_min, lat_max, pas_lat, lon_min, lon_max, pas_lon, ext)
    sav_filename2 = os.path.join(base_folder, sub_folder, sav_filename)
    if data_list is None:
        data_desc = 'Reflectance @ λ = {0:0} µm'.format(lam)
    elif data_desc == '':
        data_desc = 'unknown input data'
    # Compute the data sampling
    data, mask, grid_lat, grid_lon, mask_obs = show_omega_list_v2(omega_list,
                lam, lat_min, lat_max, lon_min, lon_max, pas_lat, pas_lon,
                data_list=data_list, mask_list=mask_list, negative_values=negative_values,
                plot=False, out=True)
    # Sav file
    input_params = {
        'omega_list' : od.get_names(omega_list),
        'lat_min' : lat_min,
        'lat_max' : lat_max,
        'lon_min' : lon_min,
        'lon_max' : lon_max,
        'pas_lat' : pas_lat,
        'pas_lon' : pas_lon,
        'data'    : data_desc,
        'filename': sav_filename,
        'datetime': datetime.datetime.now().strftime('%d/%m/%Y %H:%M')}
    save_file = {
        'data' : data,
        'mask' : mask,
        'grid_lat' : grid_lat,
        'grid_lon' : grid_lon,
        'mask_obs' : mask_obs,
        'infos' : input_params
                }
    uf.save_pickle(save_file, sav_filename2, True)
    
def load_map_omega_list(filename):
    """Load and return the result of omega_plots.show_omega_list_v2() previously saved
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
    """
    loaded_dict = uf.load_pickle(filename, True)
    data, mask, grid_lat, grid_lon, mask_obs, infos = loaded_dict.values()
    return data, mask, grid_lat, grid_lon, mask_obs, infos

def show_omega_list_v2_man(data, grid_lat, grid_lon, infos, cmap='Greys_r', vmin=None, vmax=None, 
                           title='auto', Nfig=None, polar=False, cbar=True, cb_title='auto',
                           grid=True, negatives_longitudes=False,
                           edgecolor='face', lw=0.1, **kwargs):
    """Display an composite map from a list OMEGA/MEx observations, previously sampled on 
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
    edgecolor : {'none', None, 'face', color', color sequence}, optional (default 'face')
        The color of the edges, see documentation of plt.pcolormesh for more details.
        > Added in version 2.2.8 to fix display due for new version of matplotlib.
    lw : float, optional (default 0.1)
        The line width of the edges (if diplayed).
    **kwargs:
        Optional arguments for the plt.pcolormesh() function.
    """
    lat_min, lat_max = infos['lat_min'], infos['lat_max']
    lon_min, lon_max = infos['lon_min'], infos['lon_max']
    if title == 'auto':
        title = 'Composite map from OMEGA/MEx observations' 
    fig = plt.figure(Nfig)
    Nfig = fig.number   # get the actual figure number if Nfig=None
    if polar:
        ax = plt.axes(polar=True)
        plt.pcolormesh(grid_lon*np.pi/180, grid_lat, data, cmap=cmap, 
                    vmin=vmin, vmax=vmax, edgecolor=edgecolor, lw=lw, **kwargs)
        ax.set_yticklabels([])  # remove the latitude values in the plot
        plt.xlim(0, 2*np.pi)
        if np.abs(lat_max) >= np.abs(lat_min):
            latlim = (lat_max, lat_min)
        else:
            latlim = (lat_min, lat_max)
        if latlim[0] > 0:   # Northern hemisphere
            ax.set_theta_offset(-np.pi/2)   # longitude origin at the bottom
        else:               # Southern hemisphere
            ax.set_theta_offset(np.pi/2)    # longitude origin at the top
            ax.set_theta_direction(-1)      # clockwise theta
        plt.ylim(latlim)
    else:
        if negatives_longitudes and (lon_max > 180):
            n_neg_lon = np.sum(grid_lon[:,0] > 180) # nb of negative longitudes (>180°)
            i_lon180 = np.where(grid_lon[:,0] > 180)[0][0] # 1st index of lon > 180°
            grid_lon_nl = deepcopy(grid_lon)        # new longitude grid [-180°, 180°]
            grid_lon_nl[:n_neg_lon] = grid_lon[i_lon180-1:-1] - 360
            grid_lon_nl[n_neg_lon:] = grid_lon[:i_lon180]
            data_nl = deepcopy(data)      # new data array
            data_nl[:n_neg_lon] = data[i_lon180-1:]
            data_nl[n_neg_lon:] = data[:i_lon180-1]
            plt.pcolormesh(grid_lon_nl, grid_lat, data_nl, cmap=cmap, vmin=vmin, 
                           vmax=vmax, edgecolor=edgecolor, lw=lw, **kwargs)
            lon_min, lon_max = grid_lon_nl[[0,-1], 0]   # new longitude bounds
        else:
            plt.pcolormesh(grid_lon, grid_lat, data, cmap=cmap, vmin=vmin, 
                           vmax=vmax, edgecolor=edgecolor, lw=lw, **kwargs)
        plt.gca().axis('equal')
        plt.gca().set_adjustable('box')
        plt.xlabel('Longitude [°]')
        plt.ylabel('Latitude [°]')
        plt.xlim(lon_min, lon_max)
        plt.ylim(lat_min, lat_max)
    if cbar:
        if cb_title == 'auto':
            cb_title = infos['data']
        cb = plt.colorbar()
        cb.set_label(cb_title)
    if grid:
        ax = plt.figure(Nfig).get_axes()[0]
        lonlim = ax.get_xlim()
        latlim = ax.get_ylim()
        lon_sgn = np.sign(lonlim[1] - lonlim[0])
        lat_sgn = np.sign(latlim[1] - latlim[0])
        lon_grid = np.arange(np.round(lonlim[0]/10)*10, np.round(lonlim[1]/10)*10+lon_sgn, 
                    10 * lon_sgn)   # 10° grid in longitude
        lat_grid = np.arange(np.round(latlim[0]/10)*10, np.round(latlim[1]/10)*10+lat_sgn, 
                    10 * lat_sgn)   # 10° grid in latitude
        plt.grid()
        if polar:
            ax.set_rticks(lat_grid)
        else:
            ax.set_xticks(lon_grid)
            ax.set_yticks(lat_grid)
    plt.title(title)
    plt.tight_layout()

##----------------------------------------------------------------------------------------
## Affichage spectres spécifiques extraits plots interactifs
def plot_psp(sp1_id, *args, sp2_id=(None, None), Nfig=None, sp_dict=picked_spectra, **kwargs):
    """Plot previously picked spectra from interactive plots.
    If two spectra id are given, the ration sp1/sp2 is showed.

    Parameters
    ==========
    sp1_id : tuple of int (nfig, sp_nb)
        nfig : The figure number of the selected spectra.
        sp_nb : The number of the spectra in this figure (starting at 1).
    *args : 
        Optional arguments for the plt.plot() function.
    sp2_id : tuple of int (nfig, sp_nb), optional (default (None, None))
        nfig : The figure number of the selected spectra.
        sp_nb : The number of the spectra in this figure (starting at 1).
    Nfig : int or str or None, optional (default None)
        The target figure ID.
    sp_dict : dict, optional (default picked_spectra)
        The dictionary containing the picked spectra from interactive figures.
        Default is the current one.
    **kwargs:
        Optional arguments for the plt.plot() function.
    """
    nfig1, n_sp1 = sp1_id
    nfig2, n_sp2 = sp2_id
    if (n_sp2 is None) or (nfig2 is None):
        lam = sp_dict[nfig1][0]
        sp = sp_dict[nfig1][n_sp1]
        ylabel = 'Reflectance'
    else:
        lam1, lam2 = sp_dict[nfig1][0], sp_dict[nfig2][0]
        sp_1, sp_2 = sp_dict[nfig1][n_sp1], sp_dict[nfig2][n_sp2]
        lam = od.shared_lam([lam1, lam2])
        mask_lam1 = uf.where_closer_array(lam, lam1)
        mask_lam2 = uf.where_closer_array(lam, lam2)
        sp = sp_1[mask_lam1] / sp_2[mask_lam2]
        ylabel = 'Ratioed reflectance'
    plt.figure(Nfig)
    plt.plot(lam, sp, *args, **kwargs)
    plt.xlabel('λ [µm]')
    plt.ylabel(ylabel)
    plt.tight_layout()

# ##----------------------------------------------------------------------------------------
# ## Affichage cartes composites - sauvegardées high-level
# def proj_grid_v3(Lat, Lon, data, lat_min=-90, lat_max=90, lon_min=0, lon_max=360,
                 # pas_lat=0.1, pas_lon=0.1):
    # """Sample the data from the input OMEGA/MEx observation on a given lat/lon grid.

    # Parameters
    # ==========
    # Lat : 2D array
        # The array of latitudes for each pixel of data.
    # Lon : 2D array
        # The array of longitudes for each pixel of data.
    # data : 2D array
        # The initial array of values from an OMEGA/MEx observation.
        # e.g.: Refelectance at selected wvl, spectra, or derived data such as IBD map.
    # lat_min : float, optional (default -90)
        # The minimal latitude of the grid.
    # lat_max : float, optional (default 90)
        # The maximum latitude of the grid.
    # lon_min : float, optional (default 0)
        # The minimal longitude of the grid.
    # lon_max : float, optional (default 360)
        # The maximal longitude of the grid.
    # pas_lat : float, optional (default 0.1)
        # The latitude intervals of the grid.
    # pas_lon : float, optional (default 0.1)
        # The longitude intervals of the grid.

    # Returns
    # =======
    # grid_data : 2D array (dim : Nlon x Nlat)
        # The data values, sampled on the new lat/lon grid.
    # mask : 2D array
        # The array indicating where the new grid has been filled by the OMEGA data.
    # grid lat : 2D array
        # The new latitude grid.
    # grid lon : 2D array
        # The new longitude grid.
    # """
    # # Initialisation
    # lat_array = np.arange(lat_min, lat_max+pas_lat, pas_lat)
    # lon_array = np.arange(lon_min, lon_max+pas_lon, pas_lon)
    # Nlon, Nlat = len(lon_array)-1, len(lat_array)-1
    # grid_lat, grid_lon = np.meshgrid(lat_array, lon_array)
    # grid_data = np.zeros((Nlon, Nlat))
    # mask = np.zeros((Nlon, Nlat))
    # lat2 = np.floor(np.clip(Lat, lat_min, lat_max-0.1*pas_lat) /pas_lat) * pas_lat
    # lon2 = np.floor(np.clip(Lon, lon_min, lon_max-0.1*pas_lon) /pas_lon) * pas_lon
    # nx, ny = lat2.shape
    # # Sampling on the new grid
    # for j in range(ny):
        # for i in range(nx):
            # longi, lati = lon2[i,j], lat2[i,j]
            # i_lon = int(longi/pas_lon - lon_min/pas_lon)
            # j_lat = int(lati/pas_lat - lat_min/pas_lat)
            # data_tmp = data[i,j]
            # if (not np.isnan(data_tmp)) & (data_tmp > 0):   # Filtrage régions sans données
                # grid_data[i_lon,j_lat] += data_tmp
                # mask[i_lon,j_lat] += 1
    # grid_data[grid_data==0] = np.nan
    # grid_data2 = grid_data / mask       # Normalisation
    # mask2 = np.clip(mask, 0, 1)
    # return grid_data2, mask2, grid_lat, grid_lon

# def check_list_data_mask(data_list, mask_list, disp=True):
    # """Check the compatibility between data_list and mask_list.
    # Raise ValueError if uncompatibility.

    # Parameters
    # ==========
    # data_list : 3D array
        # List of high-level map associated to a sample of OMEGA/MEx observations.
    # mask_list : 3D array
        # List of masks to remove the corrupted pixels of each OMEGA/MEx observation.
    # disp : bool, optional (default True)
        # Enable the display of the result of the test.
    # """
    # if len(mask_list) != len(data_list):
        # raise ValueError("data_list and mask_list must have the same size")
    # else:
        # for i in range(len(data_list)):
            # if mask_list[i].shape != data_list[i].shape:
                # raise ValueError("The shapes of items {0} of data_list and mask_list does not match.".format(i))
    # if disp:
        # print("\033[01;32mCompatibility between data_list and mask_list OK\033[0m")

# def show_omega_list_v3(data_list, geom_list, lat_min=-90, lat_max=90, lon_min=0, lon_max=360,
                       # pas_lat=0.1, pas_lon=0.1, cmap='Greys_r', vmin=None, vmax=None, 
                       # title='auto', Nfig=None, polar=False, cbar=True, cb_title='auto',
                       # mask_list=None, plot=True, grid=True, out=False, 
                       # negatives_longitudes=False, **kwargs):
    # """Display an composite map from a list OMEGA/MEx observations, sampled on a new lat/lon grid.
    
    # Specificities of v3 : use previously saved high-level data from OMEGA/MEx observations.

    # Parameters
    # ==========
    # data_list : 3D array
        # List of 2D array maps corresponding to severals OMEGA/MEx observations.
    # geom_list : 1D array of dict
        # List of dictionaries containing geometric informations corresponding to the 
        # OMEGA/MEx observations that are in `data_list` in the **same order**.
    # lat_min : float, optional (default -90)
        # The minimal latitude of the grid.
    # lat_max : float, optional (default 90)
        # The maximum latitude of the grid.
    # lon_min : float, optional (default 0)
        # The minimal longitude of the grid.
    # lon_max : float, optional (default 360)
        # The maximal longitude of the grid.
    # pas_lat : float, optional (default 0.1)
        # The latitude intervals of the grid.
    # pas_lon : float, optional (default 0.1)
        # The longitude intervals of the grid.
    # cmap : str, optional (default 'Greys_r')
        # The matplotlib colormap.
    # vmin : float or None, optional (default None)
        # The lower bound of the coloscale.
    # vmax : float or None, optional (default None)
        # The upper bound of the colorscale.
    # title : str, optional (default 'auto')
        # The title of the figure.
    # Nfig : int or str or None, optional (default None)
        # The target figure ID.
    # polar : bool, optional (default False)
        # If True -> Use a polar projection for the plot.
    # cbar : bool, optional (default True)
        # If True -> Diplay the colorbar.
    # cb_title : str, optional (default 'auto')
        # The title of the colorbar.
    # mask_list : 3D array
        # 1D array of the same dimension of `omega_list` containing the masks to remove the
        # corrupted pixels of each observaiton, in the **same order** that the observations of 
        # `data_list`.
        # Each mask is a 2D array, filled with 1 for good pixels and NaN for bad ones.
    # plot : bool, optional (default True)
        # If True -> Diplay the final figure.
    # grid : bool, optional (default True)
        # Enable the display of the lat/lon grid.
    # out : bool, optional (default False)
        # If True -> Return output.
    # negatives_longitudes : bool, optional (default False)
        # Argument for non-polar plots.
        # | True -> longitudes between 0° and 360°.
        # | False -> longitudus between -180° and 180°.
    # **kwargs:
        # Optional arguments for the plt.pcolormesh() function.

    # Returns (if out=True)
    # =======
    # data : 2D array (dim : Nlon x Nlat)
        # The omega reflectance at lam, sampled on the new lat/lon grid.
    # mask : 2D array
        # The array indicating where the new grid has been filled by the OMEGA data.
    # grid lat : 2D array
        # The new latitude grid.
    # grid lon : 2D array
        # The new longitude grid.
    # mask_obs : 2D array of str
        # The array indicating which observations have been used to fill each grid position.
    # """
    # # Init grids
    # lat_array = np.arange(lat_min, lat_max+pas_lat, pas_lat)
    # lon_array = np.arange(lon_min, lon_max+pas_lon, pas_lon)
    # Nlon, Nlat = len(lon_array)-1, len(lat_array)-1
    # grid_lat, grid_lon = np.meshgrid(lat_array, lon_array)
    # data, mask = np.zeros((Nlon, Nlat)), np.zeros((Nlon, Nlat))
    # mask_obs = np.ndarray((Nlon, Nlat), dtype=object)
    # mask_obs.fill('')
    # # Test compatibility
    # # > add test avec geom
    # if not (mask_list is None):
        # check_list_data_mask(data_list, mask_list, disp=True)
    # # Loop on observations -> sampling and projection on the same grid
    # for i in tqdm(range(len(data_list))):
        # if mask_list is None:
            # data_tmp = data_list[i]     # Data without mask
        # else:
            # data_tmp = data_list[i] * mask_list[i]  # Data with mask
        # data0, mask0 = proj_grid_v3(geom_list[i]['lat'], geom_list[i]['lon'], data_tmp, 
                                    # lat_min, lat_max, lon_min, lon_max, pas_lat, pas_lon)[:2]
        # data += np.nan_to_num(data0)    # Conversion NaN -> 0 pour somme des images
        # mask += mask0
        # mask_obs[mask0 == 1] += (geom_list[i]['obsname'] + ',')
    # data[mask == 0] = np.nan
    # data2 = data/mask   # Normalisation
    # # Affichage figure
    # if plot:
        # if title == 'auto':
            # title = 'Composite map from OMEGA/MEx observations' 
        # fig = plt.figure(Nfig)
        # Nfig = fig.number   # get the actual figure number if Nfig=None
        # if polar:
            # ax = plt.axes(polar=True)
            # plt.pcolormesh(grid_lon*np.pi/180, grid_lat, data2, cmap=cmap, 
                        # vmin=vmin, vmax=vmax, **kwargs)
            # ax.set_yticklabels([])  # remove the latitude values in the plot
            # plt.xlim(0, 2*np.pi)
            # if np.abs(lat_max) >= np.abs(lat_min):
                # latlim = (lat_max, lat_min)
            # else:
                # latlim = (lat_min, lat_max)
            # if latlim[0] > 0:   # Northern hemisphere
                # ax.set_theta_offset(-np.pi/2)   # longitude origin at the bottom
            # else:               # Southern hemisphere
                # ax.set_theta_offset(np.pi/2)    # longitude origin at the top
                # ax.set_theta_direction(-1)      # clockwise theta
            # plt.ylim(latlim)
        # else:
            # if negatives_longitudes and (lon_max > 180):
                # n_neg_lon = np.sum(grid_lon[:,0] > 180) # nb of negative longitudes (>180°)
                # i_lon180 = np.where(grid_lon[:,0] > 180)[0][0] # 1st index of lon > 180°
                # grid_lon_nl = deepcopy(grid_lon)        # new longitude grid [-180°, 180°]
                # grid_lon_nl[:n_neg_lon] = grid_lon[i_lon180-1:-1] - 360
                # grid_lon_nl[n_neg_lon:] = grid_lon[:i_lon180]
                # data2_nl = deepcopy(data2)      # new data array
                # data2_nl[:n_neg_lon] = data2[i_lon180-1:]
                # data2_nl[n_neg_lon:] = data2[:i_lon180-1]
                # plt.pcolormesh(grid_lon_nl, grid_lat, data2_nl, cmap=cmap, vmin=vmin, 
                               # vmax=vmax, **kwargs)
                # lon_min, lon_max = grid_lon_nl[[0,-1], 0]   # new longitude bounds
            # else:
                # plt.pcolormesh(grid_lon, grid_lat, data2, cmap=cmap, vmin=vmin, 
                               # vmax=vmax, **kwargs)
            # plt.gca().axis('equal')
            # plt.gca().set_adjustable('box')
            # plt.xlabel('Longitude [°]')
            # plt.ylabel('Latitude [°]')
            # plt.xlim(lon_min, lon_max)
            # plt.ylim(lat_min, lat_max)
        # if cbar:
            # if cb_title == 'auto':
                # cb_title = r'Reflectance @ $\lambda$' + ' = {0:.2f} µm'.format(lam)
            # cb = plt.colorbar()
            # cb.set_label(cb_title)
        # if grid:
            # ax = plt.figure(Nfig).get_axes()[0]
            # lonlim = ax.get_xlim()
            # latlim = ax.get_ylim()
            # lon_sgn = np.sign(lonlim[1] - lonlim[0])
            # lat_sgn = np.sign(latlim[1] - latlim[0])
            # lon_grid = np.arange(np.round(lonlim[0]/10)*10, np.round(lonlim[1]/10)*10+lon_sgn, 
                        # 10 * lon_sgn)   # 10° grid in longitude
            # lat_grid = np.arange(np.round(latlim[0]/10)*10, np.round(latlim[1]/10)*10+lat_sgn, 
                        # 10 * lat_sgn)   # 10° grid in latitude
            # plt.grid()
            # if polar:
                # ax.set_rticks(lat_grid)
            # else:
                # ax.set_xticks(lon_grid)
                # ax.set_yticks(lat_grid)
        # plt.title(title)
        # plt.tight_layout()
    # # Output
    # if out:
        # mask2 = np.clip(mask, 0, 1)
        # return data2, mask2, grid_lat, grid_lon, mask_obs

# def save_map_omega_list_v3(data_list, geom_list, lat_min=-90, lat_max=90, lon_min=0, lon_max=360,
                           # pas_lat=0.1, pas_lon=0.1, data_desc='', mask_list=None, sav_filename='auto', ext='',
                           # base_folder='../data/OMEGA/sav_map_list_v2/', sub_folder=''):
    # """Save the output of the omega_plots.show_omega_list_v2() function with the requested
    # parameters as a dictionary.

    # Parameters
    # ==========
    # data_list : 3D array
        # List of 2D array maps corresponding to severals OMEGA/MEx observations.
    # geom_list : 1D array of dict
        # List of dictionaries containing geometric informations corresponding to the 
        # OMEGA/MEx observations that are in `data_list` in the **same order**.
    # lat_min : float, optional (default -90)
        # The minimal latitude of the grid.
    # lat_max : float, optional (default 90)
        # The maximum latitude of the grid.
    # lon_min : float, optional (default 0)
        # The minimal longitude of the grid.
    # lon_max : float, optional (default 360)
        # The maximal longitude of the grid.
    # pas_lat : float, optional (default 0.1)
        # The latitude intervals of the grid.
    # pas_lon : float, optional (default 0.1)
        # The longitude intervals of the grid.
    # data_desc : str, optional (default '')
        # Description of the data contained in data_list (if used).
    # mask_list : 3D array
        # 1D array of the same dimension of `omega_list` containing the masks to remove the
        # corrupted pixels of each observaiton, in the **same order** than the observations of 
        # `omega_list`.
        # Each mask is a 2D array, filled with 1 for good pixels and NaN for bad ones.
    # sav_filename : str, optional (default 'auto')
        # The saving file name.
        # | If 'auto' -> Automatically generated.
    # ext : str, optional (default '')
        # Extension to add at the end of the filename (useful in case of automatic generation).
    # base_folder : str, optional (default '../data/OMEGA/sav_map_list_v2/')
        # The base folder to save the data.
    # sub_folder : str, optional (default '')
        # The subfolder to save the data.
        # Final path = "base_folder / sub_folder / sav_filename"
    # """
    # # Initialization filename
    # if sav_filename == 'auto':
        # sav_filename = ('res_show_omega_list_v3__lat{0:0>2d}-{1:0>2d}_pas{2:0}_'
                        # + 'lon{3:0>3d}-{4:0>3d}_pas{5:0}__{6:s}.pkl').format(
                            # lat_min, lat_max, pas_lat, lon_min, lon_max, pas_lon, ext)
    # sav_filename2 = os.path.join(base_folder, sub_folder, sav_filename)
    # if data_desc == '':
        # data_desc = 'unknown data'
    # # Compute the data sampling
    # data, mask, grid_lat, grid_lon, mask_obs = show_omega_list_v3(data_list,
                # geom_list, lat_min, lat_max, lon_min, lon_max, pas_lat, pas_lon,
                # mask_list=mask_list, plot=False, out=True)
    # # Sav file
    # input_params = {
        # 'omega_list' : np.array([geom_obs['obsname'] for geom_obs in geom_list]),
        # 'lat_min' : lat_min,
        # 'lat_max' : lat_max,
        # 'lon_min' : lon_min,
        # 'lon_max' : lon_max,
        # 'pas_lat' : pas_lat,
        # 'pas_lon' : pas_lon,
        # 'data'    : data_desc,
        # 'filename': sav_filename,
        # 'datetime': datetime.datetime.now().strftime('%d/%m/%Y %H:%M')}
    # save_file = {
        # 'data' : data,
        # 'mask' : mask,
        # 'grid_lat' : grid_lat,
        # 'grid_lon' : grid_lon,
        # 'mask_obs' : mask_obs,
        # 'infos' : input_params
                # }
    # uf.save_pickle(save_file, sav_filename2, True)
    
##----------------------------------------------------------------------------------------
## End of code
##----------------------------------------------------------------------------------------
