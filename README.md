![version](https://img.shields.io/badge/version-3.0.1-blue)
![pythonversion](https://img.shields.io/badge/Python-3.7+-blue)
[![DOI](https://zenodo.org/badge/349763849.svg)](https://zenodo.org/doi/10.5281/zenodo.7818828)


<p align="center">
<img width="250" height="250" src="https://github.com/AStcherbinine/omegapy/blob/master/docs/logo_omegapy_small2.png">
</p>

# OMEGA-Py : Python tools for OMEGA data

Importation and display of OMEGA/MEx observations in Python 3, based on the IDL *SOFT10* routines developped in the IAS planetary team.

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

## Documentation
Full documentation of the OMEGA-Py module is available at [astcherbinine.github.io/omegapy](https://astcherbinine.github.io/omegapy/)

Planetary Data Workshop 2023: [abstract](https://github.com/AStcherbinine/omegapy/blob/master/docs/Stcherbinine_PDW2023_7007_omegapy.pdf) & [slides](https://github.com/AStcherbinine/omegapy/blob/master/docs/PDW_Flagstaff_Stcherbinine_omegapy_upload.pdf)

## Community guidelines
To contribute to OMEGA-Py, report an issue or seek support, please refer to the community guidelines
page of the documentation available [here](https://astcherbinine.github.io/omegapy/community/).

## Citing OMEGA-Py
If you are using OMEGA-Py in your research, please cite it according to the guidelines available [here](https://astcherbinine.github.io/omegapy/credits/).

## Credits

© Aurélien Stcherbinine (2020–2024)

Institut d'Astrophysique Spatiale (IAS), Université Paris-Saclay, CNRS, Orsay, France

LATMOS/IPSL, UVSQ Université Paris-Saclay, Sorbonne Université, CNRS, Guyancourt, France

Department of Astronomy and Planetary Science, Northern Arizona University, Flagstaff, AZ, USA

Institut de Recherche en Astrophysique et Planétologie (IRAP), CNES, Université Toulouse III,
CNRS, Toulouse, France


## License
This package is released under a MIT open source license. See [`LICENSE`](https://github.com/AStcherbinine/omegapy/blob/master/LICENSE) for more details.
