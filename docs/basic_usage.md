# Basic usage

~~~python
# Package importation
import omegapy.omega_data as od
import omegapy.omega_plots as op
import omegapy.useful_functions as uf

# OMEGA file importation (ORB0979_3)
omega = od.OMEGAdata('0979_3')

# Atmospheric correction
omega_corr_atm = od.corr_atm(omega)

# Simultaneous Atmospheric & Thermal corrections (for the use of the L-channel)
# > Use the `npool` argument to control the number of simultaneous processes 
# > used to compute the thermal correction 
# > (e.g., npool=15 is usually a nice choice if your system can handle it)
# > Note: multiprocessing is currently not available for Windows
omega_corr_therm_atm = od.corr_therm_atm(omega, npool=1)

# Thermal correction only
omega_corr_therm = od.corr_therm(omega, npool=1)

# OMEGA mask to hide bad pixels / calibration lines
mask = od.omega_mask(
    omega_corr_therm_atm, 
    hide_128=True, 
    emer_lim=10, 
    inci_lim=70, 
    tempc_lim=-194, 
    limsat_c=500
    )

# Interactive display of the observation (@ λ = 1.085 µm)
op.show_omega_interactif_v2(
    omega_corr_therm_atm, 
    lam=1.085, 
    cmap='Greys_r', 
    vmin=0, 
    vmax=0.5, 
    polar=True,
    mask=mask
    )

# Search for the index of λ = 1.085 µm in the wavelength array
i_lam = uf.where_closer(1.085, omega.lam)
~~~

