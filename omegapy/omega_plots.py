#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## omega_plots.py
## Created by Aurélien STCHERBININE
## Last modified by Aurélien STCHERBININE : 06/03/2020

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
# Local
from . import useful_functions as uf
from . import omega_data as od
from .omega_data import OMEGAdata

if __name__ == '__main__':
    # Activation of the interactive mode
    plt.ion()
    # De-activation of default matplotlib keybindings for keyboard arrows
    if 'left' in mpl.rcParams['keymap.back']:
        mpl.rcParams['keymap.back'].remove('left')
    if 'right' in mpl.rcParams['keymap.forward']:
        mpl.rcParams['keymap.forward'].remove('right')

# Name of the current file
py_file = 'omega_plots.py'

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
    cmap : str (optional, default 'Greys_r')
        The matplotlib colormap.
    vmin : float or None (optional, default None)
        The lower bound of the coloscale.
    vmax : float or None (optional, default None)
        The upper bound of the colorscale.
    cb_title : str (optional, default '')
        The title of the colorbar.
    Nfig : int or str or None (optional, default None)
        The target figure ID.
    """
    fig = plt.figure(Nfig)
    plt.imshow(cube[:,:,i_lam], cmap=cmap, vmin=vmin, vmax=vmax,
               aspect='equal', origin='lower', interpolation=None)
    cb = plt.colorbar()
    cb.set_label(cb_title)
    plt.tight_layout()

def show_omega(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None,
               title='auto', xlim=(None, None), ylim=(None, None), Nfig=None):
    """Display an OMEGA/MEx observation in a rectangular pixel grid.

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
    """
    if ((lam_unit == 'm') or isinstance(lam, float)) and (lam < 10):
        i_lam = uf.where_closer(lam, omega.lam)
    else:
        i_lam = deepcopy(lam)
    lam = omega.lam[i_lam]
    if refl:
        cube = omega.cube_rf
        cb_title = 'Reflectance'
    else:
        cube = omega.cube_i
        cb_title = r'Radiance [W.m$^{-2}$.sr$^{-1}$.µm$^{-1}$]'
    if title == 'auto':
        title = ('OMEGA/MEx observation {0}\n'.format(omega.name) + 
                r'$\lambda$' + ' = {0:.2f} µm'.format(lam))
    show_cube(cube, i_lam, cmap, vmin, vmax, cb_title, Nfig)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.title(title)
    plt.tight_layout()

def show_omega_v2(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None,
                  alpha=None, title='auto', lonlim=(None, None), latlim=(None, None), Nfig=None,
                  polar=False, cbar=True):
    """Display an OMEGA/MEx observation with respect of the lat/lon coordinates of the pixels,
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
    """
    if ((lam_unit == 'm') or isinstance(lam, float)) and (lam < 10):
        i_lam = uf.where_closer(lam, omega.lam)
    else:
        i_lam = deepcopy(lam)
    lam = omega.lam[i_lam]
    if refl:
        cube = omega.cube_rf
        cb_title = 'Reflectance'
    else:
        cube = omega.cube_i
        cb_title = r'Radiance [W.m$^{-2}$.sr$^{-1}$.µm$^{-1}$]'
    if title == 'auto':
        title = ('OMEGA/MEx observation {0}\n'.format(omega.name) + 
                r'$\lambda$' + ' = {0:.2f} µm'.format(lam))
    fig = plt.figure(Nfig)
    if polar:
        ax = plt.axes(polar=True)
        plt.pcolormesh(omega.lon_grid*np.pi/180, omega.lat_grid, cube[:,:,i_lam], cmap=cmap, 
                       alpha=alpha, vmin=vmin, vmax=vmax)
        ax.set_yticklabels([])  # remove the latitude values in the plot
        ax.set_theta_offset(-np.pi/2)   # longitude origin at the bottom
        if latlim[0] is None:
            if np.max(omega.lat) > 0:
                latlim = (90, np.min(omega.lat_grid)-1)
            else:
                latlim = (-90, np.max(omega.lat_grid)+1)
    else:
        plt.pcolormesh(omega.lon_grid, omega.lat_grid, cube[:,:,i_lam], cmap=cmap, alpha=alpha,
                       vmin=vmin, vmax=vmax)
        plt.gca().axis('equal')
        plt.xlabel('Longitude [°]')
        plt.ylabel('Latitude [°]')
    if cbar:
        cb = plt.colorbar()
        cb.set_label(cb_title)
    plt.xlim(lonlim)
    plt.ylim(latlim)
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

def show_omega_interactif_v2(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', 
                          vmin=None, vmax=None, title='auto', autoyscale=True,
                          alpha=None, lonlim=(None, None), latlim=(None, None),
                          polar=False):
    """Affichage interactif d'un cube de données.
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
    """
    # Initialisation
    if refl:
        yaxis = 'Reflectance'
        cube = omega.cube_rf
    else:
        yaxis = r'Radiance [W.m$^{-2}$.sr$^{-1}$.µm$^{-1}$]'
        cube = omega.cube_i
    ny, nx, nlam = cube.shape
    if polar:
        lon, lat = omega.lon*np.pi/180, omega.lat
    else:
        lon, lat = omega.lon, omega.lat
    bij = np.zeros((ny, nx), dtype=int)
    for j in range(ny):
        for i in range(nx):
            bij[j,i] = 10000*j + i
    # fig1, ax1 = plt.subplots(1,1)
    fig1 = plt.figure()
    nfig = fig1.number
    # ax1.scatter(lon, lat, c=bij, marker='s', s=1, picker=True, alpha=0)
    show_omega_v2(omega, lam, refl, lam_unit, cmap, vmin, vmax, alpha, title, 
                    lonlim, latlim, nfig, polar)
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
        lati = omega.lat[ycoord, xcoord]
        longi = omega.lon[ycoord, xcoord]
        line = plt.plot(omega.lam, cube[ycoord, xcoord], 
                label='lat = {0:.2f}° | lon = {1:.2f}° | pixel coord = ({2:d}, {3:d})'.format(
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
            ymin, ymax = vmin, vmax
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
        # longi, lati = artist.get_offsets()[ind]
        # if not ctrl:        # Si Ctrl enfoncée, pas le plot précédent est conservé
            # plt.clf()
        # print(longi, xcoord, lati, ycoord)
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
## Affichage IBD
def show_ibd_v2(omega, ibd, cmap='viridis', vmin=None, vmax=None, alpha=None, title='auto', 
                cb_title = 'IBD', lonlim=(None, None), latlim=(None, None), Nfig=None, 
                polar=False, cbar=True):
    """Affichage ibd avec pcolormesh.
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
    """
    if title == 'auto':
        title = ('OMEGA/MEx observation {0}'.format(omega.name))
    fig = plt.figure(Nfig)
    if polar:
        ax = plt.axes(polar=True)
        plt.pcolormesh(omega.lon_grid*np.pi/180, omega.lat_grid, ibd, cmap=cmap, 
                       alpha=alpha, vmin=vmin, vmax=vmax)
        ax.set_yticklabels([])  # remove the latitude values in the plot
        ax.set_theta_offset(-np.pi/2)   # longitude origin at the bottom
    else:
        plt.pcolormesh(omega.lon_grid, omega.lat_grid, ibd, cmap=cmap, alpha=alpha,
                       vmin=vmin, vmax=vmax)
        plt.gca().axis('equal')
        plt.xlabel('Longitude [°]')
        plt.ylabel('Latitude [°]')
    if cbar:
        cb = plt.colorbar()
        cb.set_label(cb_title)
    plt.xlim(lonlim)
    plt.ylim(latlim)
    plt.title(title)
    plt.tight_layout()

##----------------------------------------------------------------------------------------
## Projection grille
def proj_grid(omega, data, lat_min=-90, lat_max=90, lon_min=0, lon_max=360,
              pas_lat=0.1, pas_lon=0.1):
    """Sample the data from the input OMEGA/MEx observation on a given lat/lon grid.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA/MEx observation
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

    Returns
    =======
    grid_data : 2D array (dim : Nlon x Nlat)
        The omega reflectance at lam, sampled on the new lat/lon grid.
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
            if (not np.isnan(data_tmp)) & (data_tmp > 0):   # Filtrage régions sans données
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
    disp : bool (optional, default True)
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

def show_omega_list_v2(omega_list, lam, lat_min=-90, lat_max=90, lon_min=0, lon_max=360,
                       pas_lat=0.1, pas_lon=0.1, cmap='Greys_r', vmin=None, vmax=None, 
                       title='auto', Nfig=None, polar=False, cbar=True, cb_title='auto',
                       data_list=None, plot=True, out=False, **kwargs):
    """Display an composite map from a list OMEGA/MEx observations, sampled on a new lat/lon grid.

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
    """
    # Sampling on same grid
    lat_array = np.arange(lat_min, lat_max+pas_lat, pas_lat)
    lon_array = np.arange(lon_min, lon_max+pas_lon, pas_lon)
    Nlon, Nlat = len(lon_array)-1, len(lat_array)-1
    grid_lat, grid_lon = np.meshgrid(lat_array, lon_array)
    data, mask = np.zeros((Nlon, Nlat)), np.zeros((Nlon, Nlat))
    if data_list is None:
        for omega in tqdm(omega_list):
            i_lam = uf.where_closer(lam, omega.lam)
            data0, mask0 = proj_grid(omega, omega.cube_rf[:,:,i_lam], lat_min, lat_max,
                                     lon_min, lon_max, pas_lat, pas_lon)[:2]
            data += np.nan_to_num(data0)    # Conversion NaN -> 0 pour somme des images
            mask += mask0
    else:
        check_list_data_omega(omega_list, data_list, disp=True)
        for i, omega in enumerate(tqdm(omega_list)):
            data0, mask0 = proj_grid(omega, data_list[i], lat_min, lat_max,
                              lon_min, lon_max, pas_lat, pas_lon)[:2]
            data += np.nan_to_num(data0)    # Conversion NaN -> 0 pour somme des images
            mask += mask0
    data[mask == 0] = np.nan
    data2 = data/mask   # Normalisation
    # Affichage figure
    if plot:
        if title == 'auto':
            title = 'Composite map from OMEGA/MEx observations' 
        fig = plt.figure(Nfig)
        if polar:
            ax = plt.axes(polar=True)
            plt.pcolormesh(grid_lon*np.pi/180, grid_lat, data2, cmap=cmap, 
                        vmin=vmin, vmax=vmax, **kwargs)
            ax.set_yticklabels([])  # remove the latitude values in the plot
            ax.set_theta_offset(-np.pi/2)   # longitude origin at the bottom
            plt.xlim(0, 2*np.pi)
            if (lat_max == 90) & (lat_min >= 0):
                latlim = (lat_max, lat_min)
            else:
                latlim = (lat_min, lat_max)
            plt.ylim(latlim)
        else:
            plt.pcolormesh(grid_lon, grid_lat, data2, cmap=cmap, vmin=vmin, 
                        vmax=vmax, **kwargs)
            plt.gca().axis('equal')
            plt.xlabel('Longitude [°]')
            plt.ylabel('Latitude [°]')
            plt.xlim(lon_min, lon_max)
            plt.ylim(lat_min, lat_max)
        if cbar:
            if cb_title == 'auto':
                cb_title = r'Reflectance at $\lambda$' + ' = {0:.2f} µm'.format(lam)
            cb = plt.colorbar()
            cb.set_label(cb_title)
        plt.title(title)
        plt.tight_layout()
    # Output
    if out:
        mask2 = np.clip(mask, 0, 1)
        return data2, mask2, grid_lat, grid_lon

##----------------------------------------------------------------------------------------
## End of code
##----------------------------------------------------------------------------------------
