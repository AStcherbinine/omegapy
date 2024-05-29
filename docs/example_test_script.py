#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## example_test_script.py
## Created by Aurélien STCHERBININE
## Last modified by Aurélien STCHERBININE : 29/05/2024

##----------------------------------------------------------------------------------------
"""Example test for OMEGA-Py: import, apply corrections, compute band depth and display
figures for OMEGA observation 0979_2.
"""
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
## Package importation
import omegapy.omega_data as od
import omegapy.omega_plots as op
import omegapy.useful_functions as uf

import numpy as np
import matplotlib.pyplot as plt
plt.ion()   # Activation of Matplotlib interactive mode

##----------------------------------------------------------------------------------------
## Parameters to adjust
# Number of simultaneous processes uses to perform the thermal correction
npool = 10

# Colormap for reflectance maps
cmap_refl = 'Greys_r'

# Colormap for BD maps
cmap_bd = 'Blues'

########
## Alternative with cmocean (if installed)
# import cmocean.cm as cmo
# cmap_bd = cmo.ice
########

## If needed, set the paths with the functions below
## (see https://astcherbinine.github.io/omegapy/configuration/#windows-or-if-you-have-troubles-using-the-environment-variables)
# od.set_omega_bin_path("/path/to/binary/files/folder/")
# od.set_omega_py_path("/path/to/omegapy-made/files/folder/")
##

##----------------------------------------------------------------------------------------
## Data importation & correction
# Load the data cube
omega0 = od.OMEGAdata('0979_2')

# Apply thermal and atmospheric corrections
omega = od.corr_therm_atm(
    omega0, 
    npool = npool,    # Adjust npool according to your system
    )

##----------------------------------------------------------------------------------------
## OMEGA mask to hide bad pixels / calibration lines
mask = od.omega_mask(
    omega,
    hide_128 = True,
    emer_lim = 10,
    inci_lim = 70,
    tempc_lim = -194,
    limsat_c = 500
    )

##----------------------------------------------------------------------------------------
## Compute band depth at 1.5μm
bd_15 = od.BD_omega(
    omega, 
    [1.5, 1.51], 
    1.30, 
    1.71, 
    norm = True
    )

# Mask on band depth values for overplotting
mask_bd15 = np.ones(bd_15.shape)    # Initialisation with array of 1
mask_bd15[bd_15 < 0.1] = np.nan     # NaN for the pixels with no ice that we want to hide

##----------------------------------------------------------------------------------------
## Display figures
# Surface reflectance
op.show_omega_v2(
    omega, 
    lam = 1.085, 
    polar = True, 
    vmin = 0, 
    vmax = 0.6,
    Nfig = "reflectance",
    cmap = cmap_refl,
    )

# 1.5μm Band Depth
op.show_data_v2(
    omega, 
    data = bd_15, 
    cb_title = "1.5μm BD", 
    polar = True, 
    cmap = cmap_bd, 
    vmin = 0, 
    vmax = 0.75,
    Nfig = "BD15",
    )

# Overplotting 1.5μm BD over surface reflectance
op.show_omega_v2(
    omega,
    lam = 1.085,
    polar = True,
    vmin = 0,
    vmax = 0.5,
    cbar = False,
    Nfig = "overplot",
    cmap = cmap_refl,
    )

op.show_data_v2(
    omega,
    data = bd_15,
    mask = mask_bd15,
    polar = True,
    cmap = cmap_bd,
    cb_title = "1.5 μm BD",
    vmin = 0,
    vmax = 0.75,
    Nfig = "overplot",
    title = 'Overplot map for ORB0979_2',
    )

plt.figure('reflectance').savefig('omegapy_ORB0979_2_reflectance.png')
plt.figure('BD15').savefig('omegapy_ORB0979_2_BD15.png')
plt.figure('overplot').savefig('omegapy_ORB0979_2_overplot_reflectance_BD15.png')

# Interactive display of the observation (@ λ = 1.085 µm)
op.show_omega_interactif_v2(
    omega, 
    lam = 1.085, 
    cmap = cmap_refl, 
    vmin = 0, 
    vmax = 0.5, 
    polar = True,
    mask = mask,
    title = 'Interactive map',
    autoyscale = False,
    ylim_sp = (0, 0.6),
    )

plt.show()

##----------------------------------------------------------------------------------------
# Search for the index of λ = 1.085 µm in the wavelength array
i_lam = uf.where_closer(1.085, omega.lam)

##----------------------------------------------------------------------------------------
## Print completed
print("\n\033[32;1m[OMEGA-Py test example completed]\033[0m")

##----------------------------------------------------------------------------------------
## End of code
##----------------------------------------------------------------------------------------
