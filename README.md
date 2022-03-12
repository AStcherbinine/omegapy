![version](https://img.shields.io/badge/version-2.2.8-blue)
![pythonversion](https://img.shields.io/badge/Python-3.7+-blue)

<p align="center">
<img width="250" height="250" src="https://github.com/AStcherbinine/omegapy/blob/master/docs/logo_omegapy_small2.png">
</p>

# OMEGA-Py : Python tools for OMEGA data

Importation and display of OMEGA/MEx observations in Python 3, based on the IDL *SOFT10* routines developped in the IAS planetary team.

> **Disclaimer:** This module is not the official software distributed by the OMEGA team.

## Installation & Update
### Method 1: from PyPI (recommended)
**Installation:** `pip3 install omegapy`

**Update:** `pip3 install omegapy --upgrade` 


### Method 2: from the GitHub repository (development version)
**Installation:** Clone the repository and install with pip:

~~~bash
git clone https://github.com/AStcherbinine/omegapy.git
cd omegapy
pip3 install .
~~~

**Update:** Go to the previously cloned repository, pull the last updates, and install them with pip:
~~~bash
cd omegapy
git pull
pip3 install .
~~~

## Configuration
You have to configure the default path of the folders containing the OMEGA binary and omegapy-made files
using the environment variables `OMEGA_BIN_PATH` (for the binary .QUB and .NAV files)
and `OMEGA_PY_PATH` (for the omegapy-made files).

To do so, add the following lines to your `~/.bashrc` :
~~~bash
export OMEGA_BIN_PATH="/path/to/binary/files/folder/"
export OMEGA_PY_PATH="/path/to/omegapy-made/files/folder/"
~~~
*Adapt the path to suit your own architecture.*

## Basic usage
~~~python
# package importation
import omegapy.omega_data as od
import omegapy.omega_plots as op
import omegapy.useful_functions as uf
# OMEGA file importation (ORB0964_2)
omega = od.OMEGAdata('0964_2')
# Atmospheric correction
omega_corr_atm = od.corr_atm(omega_corr)
# Simultaneous Atmospheric & Thermal corrections (for the use of the L-channel)
# > Use the `npool` argument to control the number of simultaneous processes used to compute the thermal correction 
# > (e.g., npool=15 is usually a nice choice if your system can handle it)
omega_corr_therm_atm = od.corr_therm_atm(omega_corr, npool=1)
# Thermal correction only
omega_corr_therm = od.corr_therm(omega, npool=1)
# Interactive display of the observation (@ λ = 1.085 µm)
op.show_omega_interactif_v2(omega_corr_therm_atm, lam=1.085, cmap='Greys_r', vmin=0, vmax=0.5, polar=True)
# Search for the index of λ = 1.085 µm in the wavelength array
i_lam = uf.where_closer(1.085, omega.lam)
~~~

See [`docs/*.md`](https://github.com/AStcherbinine/omegapy/blob/master/docs/) or the interactive IPython help for more details.

## Credits

© Aurélien Stcherbinine (2020–2022)

Institut d'Astrophysique Spatiale (IAS), Université Paris-Saclay, CNRS, Orsay, France

LATMOS/IPSL, UVSQ Université Paris-Saclay, Sorbonne Université, CNRS, Guyancourt, France


## License
This package is released under a MIT open source license. See [`LICENSE`](https://github.com/AStcherbinine/omegapy/blob/master/LICENSE) for more details.
