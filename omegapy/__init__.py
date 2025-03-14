#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''OMEGA-Py :

Submodules
==========
omega_data :
    Importation and correction (thermal + atmospheric) of OMEGA/MEx observations.
omega_plots :
    Display of OMEGA/MEx observations.
useful_functions :
    Useful python functions, some of them are used in omega_data and omega_plots.
'''

name = 'omegapy'

from . import omega_data, omega_plots, useful_functions

import json
import os

package_path = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(package_path, 'package.json'), 'r') as f:
    data_json = json.load(f)

__version__ = data_json['version']
