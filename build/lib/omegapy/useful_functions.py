#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## useful_functions.py
## Created by Aurélien STCHERBININE
## Last modified by Aurélien STCHERBININE : 29/05/2024

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
def myglob(basename, exclude=[], recursive=False):
    """Return the absolute path according to the input `basename`.
    If multiple files corresponds to `basename`, the user will be asked
    to choose one.

    --------------------------------------------
    * `int` --> Select the corresponding filename.</br>
    * `q`/`quit`/`exit` --> Return `None`.</br>
    * `a`/`all` --> Return the list of all filenames.</br>

    --------------------------------------------

    Parameters
    ----------
    basename : str
        The basename of the target file.
    exclude : array-like of str, default []
        List of sub-strings to exclude from the results.
    recursive : bool, default False
        If recursive is True, the pattern `**` will match any files and
        zero or more directories and subdirectories.

    Returns
    -------
    fname : str
        The absolute path of the selected file.
    """
    fnames = glob.glob(basename, recursive=recursive)
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
## Fonction similaires IDL
def idl_spline(X, Y, T, sigma = 1.0):
    """Performs a cubic spline interpolation.

    Parameters
    ----------
    X : ndarray
        The abcissa vector. Values MUST be monotonically increasing.
    Y : ndarray
        The vector of ordinate values corresponding to X.
    T : ndarray
        The vector of abcissae values for which the ordinate is
        desired. The values of T MUST be monotonically increasing.
    Sigma : float, default 1.0
        The amount of "tension" that is applied to the curve. The
		default value is 1.0. If sigma is close to 0, (e.g., .01),
		then effectively there is a cubic spline fit. If sigma
		is large, (e.g., greater than 10), then the fit will be like
		a polynomial interpolation.
    
    Returns
    -------
    spl : ndarray
	    Vector of interpolated ordinates.</br>
	    Result(i) = value of the function at T(i).
    """
    n = min(len(X), len(Y))
    if n <= 2:
        print('X and Y must be arrays of 3 or more elements.')
    if sigma != 1.0:
        sigma = min(sigma, 0.001)
    yp = np.zeros(2*n)
    delx1 = X[1]-X[0]
    dx1 = (Y[1]-Y[0])/delx1
    nm1 = n-1
    nmp = n+1
    delx2 = X[2]-X[1]
    delx12 = X[2]-X[0]
    c1 = -(delx12+delx1)/(delx12*delx1)
    c2 = delx12/(delx1*delx2)
    c3 = -delx1/(delx12*delx2)
    slpp1 = c1*Y[0]+c2*Y[1]+c3*Y[2]
    deln = X[nm1]-X[nm1-1]
    delnm1 = X[nm1-1]-X[nm1-2]
    delnn = X[nm1]-X[nm1-2]
    c1 = (delnn+deln)/(delnn*deln)
    c2 = -delnn/(deln*delnm1)
    c3 = deln/(delnn*delnm1)
    slppn = c3*Y[nm1-2]+c2*Y[nm1-1]+c1*Y[nm1]
    sigmap = sigma*nm1/(X[nm1]-X[0])
    dels = sigmap*delx1
    exps = np.exp(dels)
    sinhs = 0.5*(exps-1/exps)
    sinhin = 1/(delx1*sinhs)
    diag1 = sinhin*(dels*0.5*(exps+1/exps)-sinhs)
    diagin = 1/diag1
    yp[0] = diagin*(dx1-slpp1)
    spdiag = sinhin*(sinhs-dels)
    yp[n] = diagin*spdiag
    delx2 = X[1:]-X[:-1]
    dx2 = (Y[1:]-Y[:-1])/delx2
    dels = sigmap*delx2
    exps = np.exp(dels)
    sinhs = 0.5*(exps-1/exps)
    sinhin = 1/(delx2*sinhs)
    diag2 = sinhin*(dels*(0.5*(exps+1/exps))-sinhs)
    diag2 = np.concatenate([np.array([0]), diag2[:-1]+diag2[1:]])
    dx2nm1 = dx2[nm1-1]
    dx2 = np.concatenate([np.array([0]), dx2[1:]-dx2[:-1]])
    spdiag = sinhin*(sinhs-dels)
    for i in range(1, nm1):
        diagin = 1/(diag2[i]-spdiag[i-1]*yp[i+n-1])
        yp[i] = diagin*(dx2[i]-spdiag[i-1]*yp[i-1])
        yp[i+n] = diagin*spdiag[i]
    diagin = 1/(diag1-spdiag[nm1-1]*yp[n+nm1-1])
    yp[nm1] = diagin*(slppn-dx2nm1-spdiag[nm1-1]*yp[nm1-1])
    for i in range(n-2, -1, -1):
        yp[i] = yp[i]-yp[i+n]*yp[i+1]
    m = len(T)
    subs = np.repeat(nm1, m)
    s = X[nm1]-X[0]
    sigmap = sigma*nm1/s
    j = 0
    for i in range(1, nm1+1):
        while T[j] < X[i]:
            subs[j] = i
            j += 1
            if j == m:
                break
        if j == m:
            break
    subs1 = subs-1
    del1 = T-X[subs1]
    del2 = X[subs]-T
    dels = X[subs]-X[subs1]
    exps1 = np.exp(sigmap*del1)
    sinhd1 = 0.5*(exps1-1/exps1)
    exps = np.exp(sigmap*del2)
    sinhd2 = 0.5*(exps-1/exps)
    exps = exps1*exps
    sinhs = 0.5*(exps-1/exps)
    spl = (yp[subs]*sinhd1+yp[subs1]*sinhd2)/sinhs+((Y[subs]-yp[subs])*del1+(Y[subs1]-yp[subs1])*del2)/dels
    if m == 1:
        return spl[0]
    else:
        return spl

##-----------------------------------------------------------------------------------
## End of code
##-----------------------------------------------------------------------------------
