![version](https://img.shields.io/badge/version-1.2-blue)
![pythonversion](https://img.shields.io/badge/Python-3.7-blue)
![idlversion](https://img.shields.io/badge/IDL-8.3-blue)

# OMEGA Py : Python OMEGA analysis tools

Importation and display of OMEGA/MEx observations in Python 3, based on IDL routines developped in the IAS planetary team.

## Installation
Clone the repository (or download the sources), and install with pip:

~~~bash
git clone https://git.ias.u-psud.fr/astcherb1/omegapy
cd omegapy
pip install .
~~~

## Update
Go to the previously cloned repository, pull the last updates, and install with pip:
~~~bash
cd omegapy
git pull
pip install .
~~~

## Usage
~~~python
# package importation
import omegapy.omega_data as od
import omegapy.omega_plots as op
# OMEGA file importation (ORB0964_2)
omega = od.OMEGAdata('0964_2')
# Thermal correction
omega_corr_therm = od.corr_therm(omega)
# Atmospheric correction
omega_corr_therm_atm = od.corr_atm(omega_corr_therm)
# Interactive display of the observation (@ λ = 1.085 µm)
op.show_omega_interactif_v2(omega_corr_therm_atm, lam=1.085, cmap='Greys_r', vmin=0, vmax=0.5, polar=True)
~~~

See [`docs/*.md`](docs/) or the interactive IPython help for more details.

## Credits
Developped at the Institut d'Astrophysique Spatiale (IAS), Université Paris-Saclay, Orsay, France  

© Aurélien Stcherbinine (2020)
