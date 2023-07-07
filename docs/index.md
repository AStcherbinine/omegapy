# OMEGA-Py: Python tools for OMEGA data
![version](https://img.shields.io/badge/version-2.2.7-blue)
![pythonversion](https://img.shields.io/badge/Python-3.7+-blue)

OMEGA-Py is a Python module dedicated to the use and analysis of data from the OMEGA instrument onboard the ESA Mars-Express orbiter.

> **Disclaimer:** This module is not the official software distributed by the OMEGA team.

#### Features:
 - Importation of raw PSA-format data.
 - Data correction from instrumental effects.
 - Thermal and atmospheric corrections.
 - Visualization of the data with interactive tools.

## Installation
### Method 1: from PyPI (recommended)
~~~bash
pip3 install omegapy
~~~

### Method 2: from the GitHub repository (development version)
~~~bash
git clone https://github.com/AStcherbinine/omegapy.git
cd omegapy
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


## Reading of binaries

## Data correction
### Atmospheric correction

### Thermal correction

### Aerosols correction
**TODO**

## Data visualization
![Interactive display example](img/exemple_affichage_interactif_ORB0979_3_alb226.png)

## Data handling

## Futures improvement
 - Compatibility with files downloaded from the PDS (lowercase letters) and not only PSA (uppercase letters)
 - Automatic download of files from the PSA FTP
 - Use of a custom atmospheric spectrum for the atmospheric correction
 - Optimization of the customization of display functions
 - Add slider to change the displayed wavelength for the reflectance for interactive plots
 - Simplify the use of the thermal and atmospheric correction functions -> only one with multiple arguments
 - Add filtering options in the `find_cube` function

------------
## Credits

© Aurélien Stcherbinine (2020–2022)

Institut d'Astrophysique Spatiale (IAS), Université Paris-Saclay, CNRS, Orsay, France

LATMOS/IPSL, UVSQ Université Paris-Saclay, Sorbonne Université, CNRS, Guyancourt, France


## License
This package is released under a MIT open source license. See [`LICENSE`](https://github.com/AStcherbinine/omegapy/blob/master/LICENSE) for more details.
