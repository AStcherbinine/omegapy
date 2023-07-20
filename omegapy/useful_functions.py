#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## useful_functions.py
## Created by Aurélien STCHERBININE
## Last modified by Aurélien STCHERBININE : 19/07/2023

##-----------------------------------------------------------------------------------
"""Useful generics functions.
"""
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
## Packages
# Global
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from scipy.optimize import curve_fit
import pickle
import datetime
import scipy.constants as const
import glob

##-----------------------------------------------------------------------------------
## Fitting
def f_lin(x, a, b):
    """Linear function: returns f(x) = a*x + b
    
    Parameters
    ----------
    x : float or ndarray
    a : float
        The line slope.
    b : float
        The origin ordinate.
    
    Returns
    -------
    f(x) = a*x + b : float or ndarray
    """
    return a*x + b

def reg_lin(X, Y, **kwargs):
    """Return the result of the linear regression ( f(x) = a*x + b )
    on the input values.
    
    Parameters
    ----------
    X : ndarray
        The X-values.
    Y : ndarray
        The Y-values.
    **kwargs
        Optional keyword arguments to pass to the `scipy.optimize.curve_fit` function.

    Returns
    -------
    a : float
        Slope of the fitted line.
    b : float
        Origin ordinate of the fitted line.
    """
    a, b = curve_fit(f_lin, X, Y, **kwargs)[0]
    return a, b

def planck(lam, T):
    """Return the Black body radiance, associated to the input wavelength and
    temperature. According to the Planck's law.

    Parameters
    ----------
    lam : float or array-like
        The wavelength (in m).
    T : float
        The temperature (in K).
    
    Returns
    -------
    B_lam : float or array-like
        The spectral radiance (in W.m-2.sr-1.m-1).
    """
    h = const.h
    c = const.c
    kB = const.k
    B_lam = (2*h*c*c) / (lam**5) / (np.exp(h*c / (lam*kB*T)) - 1)
    return B_lam

def fit_black_body(lam, sp, T_bounds=(0, 1e6)):
    """Return the temperature associated to the fitted black body thermical
    spectrum.

    Parameters
    ----------
    lam : array-like
        The wavelength array (in m).
    sp : array-like
        The spectral radiance (in W.m-2.sr-1.m-1) to be fitted.
    T_bounds : 2-tuple, default (0, 1e6)
        The bounds for the temperature fitting.

    Returns
    -------
    T : float
        The temperature of the fitted Planck's law radiance (in K).
    """
    T = curve_fit(planck, lam, sp, bounds=T_bounds)[0][0]
    return T

def degre2(x, a, b, c):
    """2nd degree polynomial.

    Parameters
    ----------
    x : array-like or float
    a : float
    b : float
    c : float

    Returns
    -------
    y : float
        `y = a*x**2 + b*x + c`
    """
    return a*x*x + b*x + c

def degre3(x, a, b, c, d):
    """3rd degree polynomial.

    Parameters
    ----------
    x : array-like or float
    a : float
    b : float
    c : float
    d : float

    Returns
    -------
    y : float
        `y = a*x**3 + b*x**2 + c*x + d`
    """
    return a*x*x*x + b*x*x + c*x + d

##-----------------------------------------------------------------------------------
## Filters
def median_filter(sp, n):
    """Apply a median filter on the values of the spectrum, by replacing each value
    by the median of the values in a 2n+1 wide window, centered on the considered value.

    Parameters
    ----------
    sp : ndarray
        Array of transmittance values.
    n : int
        The len of the window the moving median is 2n+1.

    Returns
    -------
    sp_med : ndarray
        Filtered transmittance array.
    """
    if n==0:
        return sp
    elif n < 0:
        raise ValueError('n must be >= 0')
    sp_med = deepcopy(sp)
    for i in range(n):
        if np.isnan(sp[i]):
            sp_med[i] = np.nan
        else:
            sp_med[i] = np.nanmedian(sp_med[:2*i+1])
    for i in range(n, len(sp)-n):
        if np.isnan(sp[i]):
            sp_med[i] = np.nan
        else:
            sp_med[i] = np.nanmedian(sp_med[i-n:i+n+1])
    for i in range(len(sp)-n, len(sp)):
        if np.isnan(sp[i]):
            sp_med[i] = np.nan
        else:
            sp_med[i] = np.nanmedian(sp_med[-2*(len(sp)-i):])
    return sp_med

def moving_average(sp, n):
    """Apply a moving average filter on the values of the spectrum, by replacing each value
    by the average of the values in a 2n+1 wide window, centered on the considered value.

    Parameters
    ----------
    sp : ndarray
        Array of the transmittance values.
    n : int
        The len of the window of the moving average is 2n+1.

    Returns
    -------
    sp_med : ndarray
        Filtered transmittance array.
    """
    if n==0:
        return sp
    elif n < 0:
        raise ValueError('n must be >= 0')
    sp_moy = deepcopy(sp)
    for i in range(n):
        if np.isnan(sp[i]):
            sp_moy[i] = np.nan
        else:
            sp_moy[i] = np.nanmean(sp_moy[:2*i+1])
    for i in range(n, len(sp)-n):
        if np.isnan(sp[i]):
            sp_moy[i] = np.nan
        else:
            sp_moy[i] = np.nanmean(sp_moy[i-n:i+n+1])
    for i in range(len(sp)-n, len(sp)):
        if np.isnan(sp[i]):
            sp_moy[i] = np.nan
        else:
            sp_moy[i] = np.nanmean(sp_moy[-2*(len(sp)-i):])
    return sp_moy

##-----------------------------------------------------------------------------------
## Search
def where_closer(value, array):
    """Return the index of the closest value to `value` in `array`.
    
    Parameters
    ----------
    value : float
        Searched value.
    array : ndarray
        The array.

    Returns
    -------
    i : int
        The index of the closer value to `value` in `array`.
    """
    array2 = np.abs(array - value)
    i_closer = np.where(array2 == np.nanmin(array2))[0][0]
    return i_closer

def where_closer_array(values, array):
    """Return the list of the indexes of the closest values to `values` in `array`.

    Parameters
    ----------
    values : ndarray
        Array of searched values.
    array : ndarray
        The array.

    Returns
    -------
    I : ndarray
        Array of the index of the closest values in `array`.
    """
    i_closer = []
    for val in values:
        i_closer.append(where_closer(val, array))
    return np.array(i_closer)

##-----------------------------------------------------------------------------------
## Filename search
def myglob(basename, exclude=[]):
    """Return the absolute path according to the input `basename`.
    If multiple files corresponds to `basename`, the user will be asked
    to choose one.

    --------------------------------------------
    | `int` --> Select the corresponding filename.</br>
    | `q`/`quit`/`exit` --> Return `None`.</br>
    | `a`/`all` --> Return the list of all filenames.</br>
    --------------------------------------------

    Parameters
    ----------
    basename : str
        The basename of the target file.
    exclude : array-like of str, default []
        List of sub-strings to exclude from the results.

    Returns
    -------
    fname : str
        The absolute path of the selected file.
    """
    fnames = glob.glob(basename)
    if not isinstance(exclude, (list, np.ndarray)):
        raise ValueError("exclude parameter must be a list or numpy.ndarray")
    if len(exclude) > 0:
        fnames2 = []
        for name in fnames:
            test = True
            for excl in exclude:
                if excl in name:
                    test = False
                    continue
            if test:
                fnames2.append(name)
        fnames = fnames2
    fnames.sort()
    if fnames == []:
        # raise ValueError("No such file found.")
        print("\033[1;33mNo such file found.\033[0m")
        return None
    elif len(fnames) == 1:
        return fnames[0]
    else:
        dico = {}
        print('\033[1m{0} files found :\033[0m'.format(len(fnames)))
        for i, fname in enumerate(fnames):
            dico[str(i+1)] = fname
            print('{0:>2d} : \033[3m{1}\033[0m'.format(i+1, fname))
        print('\n\033[1mEnter the corresponding number to select one filename :\033[0m')
        while True:
            try:
                # n = input('Selection : ')
                n = input('>>> ')
                if n in dico.keys():
                    return dico[n]
                elif n=='q' or n=='quit' or n=='exit':
                    return None
                elif n=='a' or n=='all':
                    return fnames
                else:
                    print('Error, please enter an integer between 1 and ' 
                        + '{0}'.format(len(fnames)))
            except KeyboardInterrupt:
                return None

##-----------------------------------------------------------------------------------
## Sorting
def sort_dict(dico):
    """Sort a dictionary by its keys values.

    Parameters
    ----------
    dico : dict
        The input unsorted dictionary.

    Returns
    -------
    dico_sorted : dict
        The sorted dictionary.
    """
    # Conversion en np.arrays
    values = np.array(list(dico.values()))
    keys = np.array(list(dico.keys()))
    # Tri par valeurs de clé croissantes
    i_ord = np.argsort(keys)
    keys2 = deepcopy(keys[i_ord])
    values2 = deepcopy(values[i_ord])
    # Sauvegarde dans un nouveau dictionnaire 
    dico_sorted = {}
    for i in range(len(keys2)):
        dico_sorted[keys2[i]] = values2[i]
    return dico_sorted

##-----------------------------------------------------------------------------------
## Saving / Loading
def save_pickle(obj, target_path, disp=True):
    """Save an object at the selected path using the pickle module.

    Parameters
    ----------
    obj : Object
        The object to save.
    target_path : str
        The saving path name.
    disp : bool, default True
        Control the display.</br>
            | `True` --> Print the saving filename.</br>
            | `False` --> Nothing printed.
    """
    with open(target_path, 'wb') as output:
        pickle.dump(obj, output)
    if disp:
        print('\033[01;34mSaved as \033[0;03m' + target_path + '\033[0m')

def load_pickle(filename, disp=True):
    """Load and return a previously saved object with pickle.

    Parameters
    ----------
    filename : str
        The file path.
    disp : bool, default True
        Control the display.</br>
            | `True` --> Print the loading filename.</br>
            | `False` --> Nothing printed.

    Returns
    -------
    obj : Object
        The loaded object.
    """
    filename2 = myglob(filename)
    with open(filename2, 'rb') as input_file:
        obj = pickle.load(input_file)
        if disp:
            print('\033[03m' + filename2 + '\033[0;01;34m loaded\033[0m')
        return obj

##-----------------------------------------------------------------------------------
## Test existence avant sauvegarde
def test_security_overwrite(path):
    """Test if a file already exists, and if yes ask the user if he wants to
    ovewrite it or not.

    Parameters
    ----------
    path : str
        The target file path.

    Returns
    -------
    overwrite : bool
        | `True` --> No existent file, or overwriting allowed.</br>
        | `False` --> Existent file, no overwriting.
    """
    erase = 'n'
    if glob.glob(path) != []:
        try:
            erase = input('Do you really want to erase and replace \033[3m' + path +
                        '\033[0m ? (y/N) ')
        except KeyboardInterrupt:
            erase = 'n'
        if erase != 'y' :
            print("\033[1mFile preserved\033[0m")
            return False
        else:
            return True
    else:
        return True

##-----------------------------------------------------------------------------------
## End of code
##-----------------------------------------------------------------------------------
