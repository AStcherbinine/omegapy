<!--![version](https://img.shields.io/badge/version-2.4-blue)-->
<!--![pythonversion](https://img.shields.io/badge/Python-3.7+-blue)-->
<!--[![DOI](https://zenodo.org/badge/349763849.svg)](https://zenodo.org/badge/latestdoi/349763849)-->

<p align="center">
<img width="250" height="250" src="logo_omegapy_small2.png">
</p>

# OMEGA-Py : Python tools for OMEGA data

Importation and display of OMEGA/MEx observations in Python 3, based on the IDL *SOFT10* routines developped in the IAS planetary team.

!!! warning "Disclaimer"
    This module is not the official software distributed by the OMEGA team.

<!-- > **Disclaimer:** This module is not the official software distributed by the OMEGA team.-->

-------------
## Main features
 - Importation of raw PSA-format data.
 - Data correction from instrumental effects.
 - Thermal and atmospheric corrections.
 - Visualization of the data with interactive tools.

## Futures improvement
 - Compatibility with files downloaded from the PDS (lowercase letters) and not only PSA (uppercase letters)
 - Automatic download of files from the PSA FTP
 - Use of a custom atmospheric spectrum for the atmospheric correction
 - Optimization of the customization of display functions
 - Add slider to change the displayed wavelength for the reflectance for interactive plots
 - Simplify the use of the thermal and atmospheric correction functions -> only one with multiple arguments
 - Add filtering options in the `find_cube` function

