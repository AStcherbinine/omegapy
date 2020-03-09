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
In [1]: import omegapy.omega_data as od

In [2]: import omegapy.omega_plots as op
~~~

See `doc.md` or the interactive IPython help for more details.

--------------------
Developped at the Institut d'Astrophysique Spatiale (IAS), Université Paris-Saclay, Orsay, France  

© Aurélien Stcherbinine (2020)
