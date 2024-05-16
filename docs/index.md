<!--![version](https://img.shields.io/badge/version-2.4-blue)-->
<!--![pythonversion](https://img.shields.io/badge/Python-3.7+-blue)-->
<!--[![DOI](https://zenodo.org/badge/349763849.svg)](https://zenodo.org/badge/latestdoi/349763849)-->

<p align="center">
<img width="250" height="250" src="logo_omegapy_small2.png">
</p>

# OMEGA-Py : Python tools for OMEGA data

<!--Importation and display of OMEGA/MEx observations in Python 3, based on the IDL *SOFT10* routines developped in the IAS planetary team.-->

<!--!!! warning "Disclaimer"-->
<!--    This module is not the official software distributed by the OMEGA team.-->

<!-- > **Disclaimer:** This module is not the official software distributed by the OMEGA team.-->

*OMEGA-Py* is a Python 3 module dedicated to the scientific use of data
provided by the *Observatoire pour la Minéralogie, l'Eau, les Glaces et
l'Activité* (OMEGA) instrument onboard the ESA Mars Express (MEx) orbiter[^1].
It has been developed as an alternative to the historical *SOFT 10* IDL routines
of the official OMEGA software provided by the instrument team[^2].

The module notably includes a re-implementation of the most recent release of
the IDL OMEGA software, but also contains several additional data reduction
functions such as build-in atmospheric and thermal corrections (using
previously published methods) and graphics tools including interactive
visualization of the data or generation of composite OMEGA maps.

[^1]: J.-P. Bibring, A. Soufflot, M. Berthé, et al. (2004). 
OMEGA : Observatoire pour la Minéralogie, l'Eau, les Glaces et l'Activité.
*ESA Publication Division, 1240*, 37

[^2]: [ftp://psa.esac.esa.int/pub/mirror/MARS-EXPRESS/OMEGA/MEX-M-OMEGA-2-EDR-FLIGHT-EXT7-V1.0/SOFTWARE/](ftp://psa.esac.esa.int/pub/mirror/MARS-EXPRESS/OMEGA/MEX-M-OMEGA-2-EDR-FLIGHT-EXT7-V1.0/SOFTWARE/)

## Why this module?
The accessibility of data returned by space missions is a crucial point to ensure the
development of open science.
While the OMEGA dataset is public with yearly releases as of 2024, the legacy pipeline uses
a proprietary software and several crucial data reduction algorithms are not public,
thus severely hindering its use.
Since the beginning of the science phase of the OMEGA instrument in 2004, the instrument
team has provided 10 releases of the IDL software (`SOFT01` to `SOFT10`) to read the
level 1B binary files that can be downloaded from the ESA [PSA](https://archives.esac.esa.int/psa/#!Table%20View/OMEGA=instrument) 
and generate level 2A data with reflectance spectra.

However, the presence of an IDL solution only may raise some concerns:

 * The cost of an IDL license, as it is a proprietary language, makes it not accessible to everyone.
 * As the community (and especially the youngest generation) is moving to use mostly Python
   instead of IDL, the requirement to use the IDL language to access OMEGA data can limit its
   accessibility.

In addition, over the past years, the OMEGA dataset had a reputation in the community
for being challenging and requiring a lot of investment to use.
With *OMEGA-Py* we aim to tackle this reputation by providing a free all-in-one toolbox
to load, correct, analyze, and visualize the OMEGA data, and thus make the unique OMEGA dataset
rich of 20 years of observations easily accessible to the community and especially to
the younger generation of scientists and students.


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
 - Simplify the use of the thermal and atmospheric correction functions --> only one with multiple arguments
 - Add more filtering options in the `find_cube` function

