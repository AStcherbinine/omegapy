#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## omega_data.py
## Created by Aurélien STCHERBININE
## Last modified by Aurélien STCHERBININE : 20/01/2023

##----------------------------------------------------------------------------------------
"""Importation and correction of OMEGA/MEx observations from binaries files.
"""
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
## Packages
# Global
import numpy as np
from copy import deepcopy
from tqdm import tqdm
import scipy.constants as const
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.io import readsav
from scipy import interpolate
from scipy import ndimage
import datetime
import pickle
import os
import glob
import pandas as pd
import multiprocessing as mp
import itertools
import ctypes
# Local
from . import useful_functions as uf

# Name of the current file
_py_file = 'omega_data.py'
_Version = 2.2

# Path of the package files
package_path = os.path.abspath(os.path.dirname(__file__))
# Path of the directory containing the OMEGA binary files (.QUB and .NAV)
_omega_bin_path = os.getenv('OMEGA_BIN_PATH', default='/data2/opt/geomeg/data/product/')
# Path of the directory containing the OMEGA python-made files
_omega_py_path = os.getenv('OMEGA_PY_PATH', default='/data/mex-omegj/data1/omega_python/omegapy/')
# Warnings for non-defined path
if os.getenv('OMEGA_BIN_PATH') is None:
    print("\033[33mWarning: $OMEGA_BIN_PATH not defined, set to '/data2/opt/geomeg/data/product/' by default.\033[0m")
if os.getenv('OMEGA_PY_PATH') is None:
    print("\033[33mWarning: $OMEGA_PY_PATH not defined, set to '/data/mex-omegj/data1/omega_python/omegapy/' by default.\033[0m")

##-----------------------------------------------------------------------------------
## Fonctions internes privées pour importation données OMEGAdata
def _read_header(filename):
    """Lecture du header d'un fichier .NAV ou .QUB  d'une observation OMEGA/MEx
    (format PDS).

    Parameters
    ==========
    filename : str
        The path of the target file.

    Returns
    =======
    header_dict : dict
        The data contained in the header.
        type : {'keyword' : 'value', ...}
    """
    # Initialisation
    nav_file = open(filename, 'rb')
    header_dict = {}
    ktemp = ''
    # lecture du header ligne par ligne
    for line in nav_file:
        lined = line.decode('utf8').replace('\r\n', '') # décodage en UTF-8
        i_eq = lined.find('=')
        if i_eq != -1:                              # Cas "keyword = value"
            keyword = lined[:i_eq].replace(' ', '') # keyword : retrait espaces
            value = lined[i_eq+2:].lstrip()         # valeur : retrait espaces après le '='
            header_dict[keyword] = value
            ktemp = keyword
        elif (lined == '') or (lined[:2] == '/*'):  # ligne vide ou commentaire
            continue
        elif lined == 'END':                        # fin header
            break
        else:                                       # description sur plusieurs lignes
            header_dict[ktemp] += lined.lstrip()
    # Fermeture du fichier
    nav_file.close()
    # Uniformisation des clés
    keys = header_dict.keys()
    if not 'SPACECRAFT_POINTING_MODE' in keys:
        if 'SC_POINTING_MODE' in keys:
            header_dict['SPACECRAFT_POINTING_MODE'] = deepcopy(header_dict['SC_POINTING_MODE'])
            del header_dict['SC_POINTING_MODE']
        elif 'SPACECRACRAFT_POINTING_MODE' in keys:
            header_dict['SPACECRAFT_POINTING_MODE'] = deepcopy(header_dict['SPACECRACRAFT_POINTING_MODE'])
            del header_dict['SPACECRACRAFT_POINTING_MODE']
        else:
            header_dict['SPACECRAFT_POINTING_MODE'] = '"N/A"'
    if header_dict['SPACECRAFT_POINTING_MODE'] == '"UNK"':
        header_dict['SPACECRAFT_POINTING_MODE'] = '"UNKNOWN"'
    # Sortie
    return header_dict

def _compute_local_time(lon, subsol_lon):
    """Compute the local time at each point of an OMEGA observation.

    Parameters
    ==========
    lon : 2D array
        The longitude of each pixel (deg).
    subsol_lon : float
        The longitude of the sub-solar point at observation time.

    Returns
    =======
    loct : 2D array of floats
        The array of the local time for each pixel of the observation.
    """
    # Initialisation
    loct = -np.ones(lon.shape)
    # Note : sub_solar_longitude = longitude du midi au moment de l'orbite
    # Calcul écart au point sub-solaire
    delta_lon = subsol_lon - lon    # en longitude
    delta_loct = np.abs(delta_lon * 12/180)     # conversion en heure
    mask0 = (delta_loct > 12)
    delta_loct[mask0] = 24 - delta_loct[mask0]  # écart à midi < 12h (en absolu)
    # Calcul heure locale en chaque point
    mask = (delta_lon > 180) | ((delta_lon < 0) & (delta_lon > -180))
    loct[mask] = 12 + delta_loct[mask]                  # Aprem -> on rajoute des heures
    loct[mask==False] = 12 - delta_loct[mask==False]    # Matin -> on retire des heures
    # Output
    return loct

def _utc_to_my(dt):
    """Convert a UTC datetime to the corresponding Martian Year (MY).
    
    Martian Years are numbered according to the calendar proposed by R. Todd Clancy 
    (Clancy et al., Journal of Geophys. Res 105, p 9553, 2000): Martian Year 1 begins 
    (at a time such that Ls=0) on April 11th, 1955.

    Parameters
    ==========
    dt : datetime.datetime
        The UTC datetime object.

    Returns
    =======
    my : int
        The corresponding Martian Year.
    """
    datetime_my1 = datetime.datetime(1955, 4, 11)   # Start MY 1
    my_sol_duration = 668.6     # Nb of Martian sols during a MY
    sol_sec_duration = 88775.245    # Duration of a sol in seconds
    my = int( (dt - datetime_my1).total_seconds() // (my_sol_duration * sol_sec_duration)) + 1
    return my

def _read_cube(filename_qub, disp=True):
    """Python implementation of the `readcube.pro` routine from the SOFT09 OMEGA pipeline.
    Extract and return the binary data from an OMEGA/MEx .QUB file.

    Parameters
    ==========
    filename_qub : str
        The path to the .QUB file.
    disp : bool, optional (default True)
        Enable or disable the display of informations during the file reading.
        | True -> Enable display.

    Returns
    =======
    idat : 3D ndarray

    sdat0 : 2D ndarray

    sdat1 : 3D ndarray

    info : 1D ndarray

    """
    # Initialisation
    info = np.zeros(6)
    # Lecture header
    hd_qub = _read_header(filename_qub)
    test = int.from_bytes(str.encode(hd_qub['COMMAND_DESC'][34]), byteorder='big') - 48
    if test > 9:
        test -= 7
    flagc = np.zeros(3, dtype=np.int64)
    if test // 8:
        flagc[2] = 1
    if (test // 4) and 1:
        flagc[1] = 1
    if (test // 2) and 1:
        flagc[0] = 1
    lrec = np.int64(hd_qub['RECORD_BYTES'])
    dqual = np.float64(hd_qub['DATA_QUALITY_ID'])
    info[5] = dqual + 0.001
    nrec = np.int64(hd_qub['LABEL_RECORDS'])
    # ^IMAGE ? (nrec + nax) ?
    tps_exp = np.array(hd_qub['EXPOSURE_DURATION'][1:14].split(','), dtype=np.float64)
    info[:3] = tps_exp * flagc
    info[3] = np.int64(hd_qub['DOWNTRACK_SUMMING'])
    info[4] = np.float64(hd_qub['INST_CMPRS_RATE'])
    npixel, npara, nscan = np.array(hd_qub['CORE_ITEMS'][1:-1].split(','), dtype=np.int64)
    cbyte = np.int64(hd_qub['CORE_ITEM_BYTES'])
    sax0, sax1, sax2 = np.array(hd_qub['SUFFIX_ITEMS'][1:6].split(','), dtype=np.int64)
    sbyte = np.int64(hd_qub['SUFFIX_BYTES'])
    # print infos
    if disp:
        print('core:   {0:>8d}{1:>8d}{2:>8d}  cbyte:{3:>8d}'.format(npixel, npara, nscan, cbyte))
        print('suffix: {0:>8d}{1:>8d}{2:>8d}  sbyte:{3:>8d}'.format(sax0, sax1, sax2, sbyte))
    #----------------------------
    # Lecture données cube
    # nscan
    #  |-- npara  +  npixel x sax1 (int 32bits)
    #       |-- npixel (int 16bits)  +  1 (int 32 bits)
    class F_line(ctypes.Structure):
        _fields_ = [('C_line', ctypes.c_int16 * npixel), ('S_line', ctypes.c_int32)]

    class SS_line(ctypes.Structure):
        _fields_ = [('SS_line', ctypes.c_int32 * npixel)]

    class F_frame(ctypes.Structure):
        _fields_ = [('F_line', F_line * npara), ('SS_frame', SS_line * sax1)]

    class F_Qube(ctypes.Structure):
        _fields_ = [('F_frame', F_frame * nscan)]

    fqube = F_Qube()
    skip = lrec * nrec  # taille header (en bytes)
    with open(filename_qub, 'rb') as cube_file:
        data_hd = cube_file.read(skip)  # bytes header
        data_qub = cube_file.readinto(fqube)     # bytes data

    idat = np.ndarray((nscan, npara, npixel), np.int32)
    sdat0 = np.ndarray((nscan, npara), np.int32)
    sdat1 = np.ndarray((nscan, sax1, npixel), np.int32)

    fqube2 = np.ctypeslib.as_array(fqube.F_frame)
    idat[:,:,:] = fqube2['F_line']['C_line']
    sdat0[:,:] = fqube2['F_line']['S_line']
    if sax1 > 0:
        sdat1[:,:,:] = fqube2['SS_frame']['SS_line']
    # On remet dans le même ordre que readomega
    # idat : (npixel, npara, nscan)
    # sdat0: (sax0, npara, nscan) = (npara, nscan)
    # sdat1: (npixel, sax1, nscan)
    idat = np.swapaxes(idat, 0, 2)
    sdat0 = np.transpose(sdat0)
    sdat1 = np.swapaxes(sdat1, 0, 2)
    return idat, sdat0, sdat1, info

class CubeError(Exception):
    """Exception raised if encounter an issue when reading the OMEGA/MEx cube
    binaries data.
    """
    def __init__(self, message):
        self.message = message

def _readomega(cube_id, disp=True, corrV=True, corrL=True, data_path='_omega_bin_path'):
    """Python implementation of the `readomega.pro` routine from the SOFT09 OMEGA pipeline.
    
    Parameters
    ==========
    cube_id : str
        The observation ID (format ORBXXXX_X).
    disp : bool, optional (default True)
        Enable or disable the display of informations during the file reading.
        | True -> Enable display.
    corrV : bool, optional (default True)
        If True, compute the correction on the visible channel (Vis).
    corrL : bool, optional (default True)
        If True, compute the correction on the long-IR channel (L).
    data_path : str, optional (default _omega_bin_path)
        The path of the directory containing the data (.QUB) and 
        navigation (.NAV) files.

    Returns
    =======
    out_data : dict
        {'ldat', 'jdat', 'wvl', 'ic', 'specmars'}
    out_geom : dict
        {'latitude', 'longitude', 'emergence', 'incidence', 'altitude', 'ut_time',
         'temperature, 'saturation_C', 'saturation_vis', 'lat_grid', 'lon_grid'}
    """
    # Default path
    if data_path == "_omega_bin_path":
        data_path = _omega_bin_path
    # Filename
    nomgeo0 = cube_id + '.NAV'
    nomfic0 = cube_id + '.QUB'
    nomfic = os.path.join(data_path, nomfic0)
    nomgeo = os.path.join(data_path, nomgeo0)
    if disp:
        print('\n\033[1mComputing OMEGA observation {0:s}\033[0m'.format(cube_id))
    # Orbit number in base 10
    orbnum = int.from_bytes(str.encode(nomfic0[3]), byteorder='big') - 48
    if orbnum > 9:
        orbnum -= 7
    orbnum = 1000 * orbnum + int(nomfic0[4:7])
    # Read cube data
    idat, sdat0, sdat1, info = _read_cube(nomfic, disp)
    # Stop if only one line in the cube
    # if (len(idat.shape) != 3) or (len(idat) == 1):
    if idat.shape[2] == 1:
        raise CubeError('Only one line in cube {0:s}'.format(nomfic0))
    # Extract data from info array
    exposure = info[0:3]
    summation = info[3]
    bits_per_data = info[4]
    data_quality = int(info[5])

    #--------------------------
    # Test presence geometry cube
    if not os.path.exists(nomgeo):
        raise CubeError('No corresponding NAV cube {0:s}'.format(nomgeo0))

    #--------------------------
    # Data from the geometry file
    #--------------------------
    # Data from the .NAV header
    #--------------------------
    hd_nav = _read_header(nomgeo)
    npixel, npara, nscan = np.array(hd_nav['CORE_ITEMS'][1:-1].split(','), dtype=np.int64)
    lrec = np.int64(hd_nav['RECORD_BYTES'])
    nrec = np.int64(hd_nav['LABEL_RECORDS'])
    try:
        nau = np.int64(np.float64(hd_nav['SOLAR_DISTANCE']) / 14960)
    except KeyError:
        nau = 15000
        print('\033[1;33mNo Solar distance -> 1.5 AU\033[0m')
    #--------------------------
    # > Lien dossier à adapter
    specmars = np.loadtxt(os.path.join(package_path, 'OMEGA_dataref', 'specsol_0403.dat')) # Spectre solaire Terre
    dmars = nau*1e-4
    specmars /= dmars**2    # Spectre solaire au niveau de Mars
    #--------------------------
    # Lecture géometrie
    class F_line_nav(ctypes.Structure):
        _fields_ = [('C_line', ctypes.c_int32 * npixel)]

    class F_frame_nav(ctypes.Structure):
        _fields_ = [('F_line', F_line_nav * npara)]

    class F_Qube_nav(ctypes.Structure):
        _fields_ = [('F_frame', F_frame_nav * nscan)]

    fqube_nav = F_Qube_nav()
    skip = lrec * nrec  # taille header (en bytes)
    geocube = np.ndarray((nscan, npara, npixel), np.int32)
    with open(nomgeo, 'rb') as nav_cube:
        data_hd = nav_cube.read(skip)  # bytes header
        data_qub = nav_cube.readinto(fqube_nav)     # bytes data

    fqube_nav2 = np.ctypeslib.as_array(fqube_nav.F_frame)
    geocube[:,:,:] = fqube_nav2['F_line']['C_line']
    # On remet dans le même sens que readomega
    geocube = np.swapaxes(geocube, 0, 2)

    trans = np.pi / 180 * 1e-4
    # extraction données
    ecl = np.cos(geocube[:, 2, :] * trans)
    longitude = geocube[:, 6, :] * 1e-4
    latitude = geocube[:, 7, :] * 1e-4
    altitude = geocube[:, 12, :] * 1e-3
    emergence = geocube[:, 3, :] * 1e-4
    incidence = geocube[:, 2, :] * 1e-4
    ut_time = geocube[:, 1, :]
    lon_grid = np.swapaxes(geocube[:, 13:17, :], 1, 2) * 1e-4
    lat_grid = np.swapaxes(geocube[:, 17:21, :], 1, 2) * 1e-4
    
    #--------------------------
    # Preliminary pipeline tool
    #--------------------------
    npix, _, nbal = idat.shape

    fond2 = np.loadtxt(os.path.join(package_path, 'OMEGA_dataref', 'fond2.dat'), dtype=int)
    if exposure[1] > 4:
        fond2 *= 2
    fondcur = np.transpose(fond2[128:] * np.ones((nbal, 128)))
    
    # Init jdat
    jdat = deepcopy(idat).astype(np.float64)
    idat_leq1 = np.where(idat <= 1)
    jdat[idat_leq1] = 1e-5
    # Nombre pixels à 0 IR
    pix0IR = len(np.where(idat[:, 0:256, :] <= 0)[0])
    if disp:
        print('        0 or less IR: {:>8d}'.format(pix0IR))

    hkmin = np.min(sdat1[14, 1, :])
    if hkmin < 6:
        hkmin = 6
    indj = np.where(sdat1[14, 1, :] >= hkmin)[0]
    i6 = indj[0]
    indj = np.where(sdat1[14, 1, :] >= hkmin+1)[0]
    balHK = indj[0] - i6
    indi = np.where((sdat1[14, 1, :] >= 6) & (sdat0[10, :] >= 100))[0]
    balsm = balHK*8

    if (balHK > 0) and (nbal > indi[0]+balsm):
        ndeb = indi[0]
        nf = len(indi) - 1
        nfin = indi[nf]
        if sdat0[26, nbal-1] < 20:
            sdat0[:, nbal-1] = sdat0[:, nbal-2]
        for k in range(256):
            b = sdat0[k, indi].astype(np.float64)
            a = np.concatenate([
                2*b[0] - b[:balsm][::-1], b, 2*b[nf] - b[nf-balsm+1:] ])
            c = ndimage.uniform_filter(a, balsm, mode='nearest')[balsm:balsm+nf+1]
            tck = interpolate.splrep(indi, c, k=3)  # Cubic spline interpolation
            d = (interpolate.splev(ndeb + np.arange(nbal-ndeb), tck) 
                    - sdat0[k, ndeb:nbal])
            jdat[:, k, ndeb:nbal] += d

    jdat[:, 128:256, :] -= fondcur
    # Valeur min jdat = 1e-5 (pas 0)
    jdat_eq0 = np.where(jdat < 1e-5)
    jdat[jdat_eq0] = 1e-5
    # Normalisation sommation
    if summation != 1:
        jdat /= summation
    # Correction linéarité voie C
    linearC = np.loadtxt(os.path.join(package_path, 'OMEGA_dataref', 'linearC.dat'))
    jdat[:, 0:128, :] = linearC[np.clip((jdat[:, 0:128, :] + 0.5).astype(int), None, 4095)]
    # Correspondance longueur d'onde
    wvl = np.loadtxt(os.path.join(package_path, 'OMEGA_dataref', 'lambda_0403.dat'))
    # Importation fichiers
    bound = np.loadtxt(os.path.join(package_path, 'OMEGA_dataref', 'boundcur.dat'), unpack=True)
    if exposure[1] < 4:
        mtf = np.loadtxt(os.path.join(package_path, 'OMEGA_dataref', 'mtf120315_25.dat'))
        rap = np.loadtxt(os.path.join(package_path, 'OMEGA_dataref', 'rapcur_25.dat'), unpack=True)
    else:
        mtf = np.loadtxt(os.path.join(package_path, 'OMEGA_dataref', 'mtf120315_50.dat'))
        rap = np.loadtxt(os.path.join(package_path, 'OMEGA_dataref', 'rapcur_50.dat'), unpack=True)
    ib0 = np.where(orbnum > bound[0, :])[0]
    mtf[ib0] *= rap[0, ib0]
    ib1 = np.where(orbnum > bound[1, :])[0]
    mtf[ib1] *= rap[1, ib1]
    ib2 = np.where(orbnum > bound[2, :])[0]
    mtf[ib2] *= rap[2, ib2]

    #--------------------------
    # Correction voie L
    #--------------------------
    if corrL:
        orbmax = int(os.path.getsize(os.path.join(package_path, 'OMEGA_dataref', 'rapmtflcur.bin')) / 512 - 1)
        orbcur = deepcopy(orbnum)
        if orbcur > orbmax:
            orbcur = orbmax
            print("\033[1;33mWarning: L channel corrected as for orbit {:d}\033[0m".format(orbcur))
        offs = 512 * orbcur

        class Rapmtflcur(ctypes.Structure):
            _fields_ = [('line', ctypes.c_float * 128)]

        raplb = Rapmtflcur()
        with open(os.path.join(package_path, 'OMEGA_dataref', 'rapmtflcur.bin'), 'rb') as rapmtf_file:
            _ = rapmtf_file.read(offs)
            nbytes = rapmtf_file.readinto(raplb)
        rapl = np.ctypeslib.as_array(raplb.line)
        mtf[128:256] *= rapl

    if exposure[1] > 5:
        mtf[0:256] *= (exposure[1] / 5)
    jdat_tmp = np.swapaxes(jdat, 1, 2)  # Réordonnement axes pour calcul
    jdat_tmp[:, :, :256] /= mtf[:256]
    jdat = np.swapaxes(jdat_tmp, 1, 2)  # On remet dans le bon ordre

    ic = np.where(mtf < 10000)[0]

    #--------------------------
    # Correction voie Vis
    #--------------------------
    if corrV:
        #--------------------------
        # Bit error correction
        #--------------------------
        vis = idat[:, 256:, :].astype(np.float64)
        pixels, _, lines = idat.shape
        exptime = info[2] / 1000

        if pixels == 128:
            level = info[3] * 4095 * 2
        else:
            level = 4095 * info[3]
        # VIS -> pixels négatifs
        i_vis_neg = np.where(vis <= 0)
        counter_neg = len(i_vis_neg[0])
        vis[i_vis_neg] = 0.001
        # VIS -> pixels supérieurs à valeur max
        i_vis_pos = np.where(vis > level)
        counter_pos = len(i_vis_pos[0])
        vis[i_vis_pos] = level
        # VIS -> pixels saturés
        i_vis_sat = np.where((vis <= level) & (vis > (0.8*level)))
        counter_sat = len(i_vis_sat[0])
        vis[i_vis_sat] = 0.8 * level
        # VIS -> spikes
        plan3 = deepcopy(vis[:, 3, :])
        i_spike3 = np.where(plan3 > (0.5 * (vis[:, 2, :] + vis[:, 3, :]) + 50))
        counter_spike3 = len(i_spike3[0])
        plan3[i_spike3] = 0.5 * (vis[:, 2, :] + vis[:, 3, :])[i_spike3]
        vis[:, 3, :] = plan3
        counter_spikes = 0
        # Despiking si plus de 10 lignes
        if lines >= 10:
            # Version python avec filtre sur le cube en une fois
            # Note : on retire les 3 points extrêmes qui ne sont pas pris en compte avec IDL
            median_lines = ndimage.median_filter(vis, size=(1,1,7), mode='nearest')
            i_spikes = np.where(np.abs(median_lines[:,:,3:-3] - vis[:,:,3:-3]) > 200)
            counter_spikes = len(i_spikes[0])
            vis[:,:,3:-3][i_spikes] = median_lines[:,:,3:-3][i_spikes]
        
        if disp:
            print(' negative pixels VIS: {:>8d}'.format(counter_neg))
            print('anomalous pixels VIS: {:>8d}'.format(counter_pos))
            print('saturated pixels VIS: {:>8d}'.format(counter_sat))
            print('          spikes VIS: {:>8d}'.format(counter_spikes))

        #--------------------------
        # Bias correction
        #--------------------------
        percen = 3.4e-6
        fper = 0
        pend = 5
        if pixels == 128:
            fper = 1
        elif pixels == 64:
            fper = 2
        elif pixels == 32:
            fper = 4
        elif pixels == 16:
            fper = 8

        lv = wvl[256:]
        xa, xb = lv[0], lv[95]
        ya, yb = (percen*fper), (percen*fper + (percen*fper)/100*pend)
        a = (yb - ya) / (xb - xa)
        b = ya - a*xa
        Strl_nc = lv * a + b
        # Dimensions Strl_nc : (96) -> (96, 596) pour calcul en tableaux
        Strl_nc = np.repeat(np.array([Strl_nc]), lines, axis=0).T

        strL2 = np.sum(vis, axis=0)
        strL1 = np.sum(vis, axis=(0,1))

        vis = (vis - 50*strL2*percen*fper) - 0.75*strL1*Strl_nc

        #--------------------------
        # Smear correction
        #--------------------------
        bands = 96
        # Correction for the CCD cleaning time
        ft = 0.87
        time_int = info[2] - ft

        # Check for saturation in the 90:95 channels
        vis_tmp = np.swapaxes(vis, 0, 1)    # Inversions axes pour calcul (lam, npix, nscan)
        means = np.mean(vis_tmp[0:85], axis=0)
        i_sat = np.where(vis_tmp[85:] > means)
        vis_tmp[i_sat] = 0
        vis = np.swapaxes(vis_tmp, 0, 1)        # On remet dans l'ordre

        # Smearing term
        smear = np.ndarray(vis.shape, dtype=np.float64)
        for i in range(bands):
            smear[:, i, :] = np.sum(vis[:, i:, :], axis=1) * ft / time_int / bands

        vis -= smear

        #--------------------------
        # Flat cube
        #--------------------------
        # Lecture fichier binaire flat
        class FlatVisLine(ctypes.Structure):
            _fields_ = [('Line', ctypes.c_uint32 * 128)]

        class FlatVisFrame(ctypes.Structure):
            _fields_ = [('Frame', FlatVisLine * 96)]

        Sliceb = FlatVisFrame()
        with open(os.path.join(package_path, 'OMEGA_dataref', 'flatVIS050701.bin'), 'rb') as flat_file:
            nbytes = flat_file.readinto(Sliceb)
        Slice = np.ctypeslib.as_array(Sliceb.Frame)
        Slice = Slice['Line'].T
        # Indices début / fin
        istart = int((128 - pixels) / 2)
        iend = int(istart + pixels)
        # Inversion axes pour calcul
        Slice_tmp = Slice.T
        vis_tmp = np.swapaxes(vis, 0, 2)
        # Correction
        vis_tmp = vis_tmp * 2**15 / Slice_tmp[:, istart:iend]
        # On remet dans l'ordre
        vis = np.swapaxes(vis_tmp, 0, 2)

        #--------------------------
        # Radiometric calibration
        #--------------------------
        f = 1
        if pixels == 128:
            f = 2   # Internal summation
        f2 = 1
        if exptime == 0.1:
            f2 = 2
        elif exptime == 0.2:
            f2 = 4
        #--------------------------
        def ordcorr(lam, sp):
            """
            Parameters
            ==========
            lam : ndarray
                VNIR wavelength in microns.
            sp : ndarray
                The spectrum to correct.

            Returns
            =======
            I_corr : ndarray
                The corrected spectrum.
            """
            #_______________________________________________________________________
            # Assignment of the variable k[3, kn], from file Ceoff_II_96_ias.txt
            # k[0, :] = spectral channels with the second order contamination
            # k[1, :] = wavelengths correspondent to each spectral channel
            # k[2, :] = values of the coefficients for each channel
            nk = 17 # Total number of channels with the contamination
            k = np.ndarray((3, nk))
            k[0] = np.arange(nk) + 79
            k[1] = np.array([952.100, 959.600, 967.000, 974.300, 981.700, 988.900, 
                    996.300, 1003.70, 1010.90, 1018.20, 1025.50, 1032.90, 1040.30, 
                    1047.94, 1055.34, 1062.52, 1069.87]) / 1000
            k[2] = np.array([0.00478825, 0.00764561, 0.0117481, 0.0186357, 0.0291402, 
                    0.0430787, 0.0561666,  0.0694787,  0.0831372,  0.0993437, 0.113755, 
                    0.124017, 0.123925, 0.118724, 0.110151, 0.106578, 0.106578])
            #_______________________________________________________________________
            # Formula : I_corr = I_raw-kII*I_1_raw
            # I_corr = corrected final Intensity
            # I_raw = raw spectrum measured by VNIR (initial intensity)
            # KII = coefficients for the second order
            # I_1_raw = raw intensity responsable of the second order contribution
            I_raw = np.zeros(96)
            I_1_raw = np.zeros(96)
            kII = np.zeros(96)
            I_corr = np.zeros(96)
            #_______________________________________________________________________
            # Assignment of I_raw
            # lam = deepcopy(wvl[256:352])
            I_raw = deepcopy(sp)
            #_______________________________________________________________________
            # Assignment of I_1_raw
            ch = k[0]   # Spectral channels concerning the II ord
            l2 = k[1]   # Lambda of the II order
            l1 = l2/2   # Lambda of the first order
            nl = len(l1)
            ch_start = 16   # First spectral channel of the first order wavelength
            ch_end = 23     # Last spectral channel of the first order wavelength
            # Cubic spline interpolation
            spl = interpolate.InterpolatedUnivariateSpline(lam[ch_start:ch_end+1], I_raw[ch_start:ch_end+1], k=3)
            I_1 = spl(l1)
            I_1_raw[int(k[0, 0]):int(k[0, nl-1])+1] = I_1
            #_______________________________________________________________________
            # Assignment of kII
            kII[int(k[0, 0]):int(k[0, nl-1])+1] = k[2]
            #_______________________________________________________________________
            # Computation of I_corr
            I_corr = I_raw - kII * I_1_raw
            #_______________________________________________________________________
            # Output
            return I_corr
        #--------------------------
        for i in range(lines):
            for m in range(pixels):
                I_corr = ordcorr(wvl[256:352], vis[m, : , i])
                jdat[m, 256:352, i] = I_corr / (f2 * mtf[256:352] * f * summation)

    #--------------------------
    # Reflectance factor calculation
    #--------------------------
    ldat = deepcopy(jdat)
    for n in range(352):
        ldat[:, n, :] = ldat[:, n, :] / (specmars[n] * ecl[:, :])
    albedo = deepcopy(ldat[:, 11, :])

    #--------------------------
    # Saturation
    #--------------------------
    # Voie C
    saturation_c = sdat0[39, :] - (idat[:, 39, :] / summation)
    # Voie visible
    saturation_vis = idat[:, 299, :] / summation
    # Temperature voie C
    temperature = sdat1[0, 2, :] * 0.001
    # Temperature voie L
    # temperature = sdat1[1, 2, :] * 0.001

    #--------------------------
    # Output data
    #--------------------------
    out_data = {
        'ldat'  : np.swapaxes(ldat, 0, 2),
        'jdat'  : np.swapaxes(jdat, 0, 2),
        'wvl'   : wvl,
        'ic'    : ic,
        'specmars'  : specmars
    }
    out_geom = {
        'latitude'  : latitude.T,
        'longitude' : longitude.T,
        'emergence' : emergence.T,
        'incidence' : incidence.T,
        'altitude'  : altitude.T,
        'ut_time'   : ut_time.T,
        'temperature'   : temperature,
        'saturation_c'  : saturation_c.T,
        'saturation_vis': saturation_vis.T,
        'lat_grid'  : np.moveaxis(lat_grid, [0,1,2], [1,0,2]),
        'lon_grid'  : np.moveaxis(lon_grid, [0,1,2], [1,0,2])
    }
    return out_data, out_geom

##-----------------------------------------------------------------------------------
## Class OMEGAdata - Importation and 1st correction of OMEGA data cubes
class OMEGAdata:
    """Importation of OMEGA/MEx observation.

    Parameters
    ==========
    obs : str
        The name of the OMEGA observation.
    empty : bool, optional (default False)
        If True, return an empty OMEGAdata object.
    data_path : str, optional (default _omega_bin_path)
        The path of the directory containing the data (.QUB) and 
        navigation (.NAV) files.
    corrV : bool, optional (default True)
        If True, compute the correction on the visible channel (Vis).
    corrL : bool, optional (default True)
        If True, compute the correction on the long-IR channel (L).
    disp : bool, optional (default True)
        Enable or disable the display of informations during the file reading.
        | True -> Enable display.

    Attributes
    ==========
    name : str
        The observation ID.
    lam : 1D array
        The wavelength array (in µm).
    cube_i : 3D array
        The hyperspectral data cube in physical units (W.m-2.sr-1.µm-1).
        dim : [X, Y, wvl]
    cube_rf : 3D array
        The I/F ratio hyperspectral data cube.
        dim : [X, Y, wvl]
    ls : float
        The Solar longitude of the observation (deg).
    lat : 2D array
        The latitude of each pixel (deg).
    lon : 2D array
        The longitude of each pixel (deg).
    alt : 2D array
        The elevation of the pixel footprint center point (km).
    loct : 2D array of floats
        The array of the local time for each pixel of the observation.
    my : int
        The Martian Year number at the time of the observation.
    emer : 2D array
        The angle of emergent line (from the surface) (deg).
    inci : 2D array
        The incidence angle at the surface (deg).
    specmars : 1D array
        The Solar radiation spectrum on Mars (W.m-2.sr-1.µm-1).
    utc : datetime.datetime
        The average UTC time of the observation.
    ic : dict
        The index of the used spectral pixels for each channel.
    lat_grid : 2D array
        The latitude grid of the observation (from the edge of the pixels).
    lon_grid : 2D array
        The longitude grid of the observation (from the edge of the pixels).
    surf_temp : 2D array
        The retrieved surface temperature of each pixel (from the thermal correction).
    sensor_temp_c : 1D array
        The temperature of the sensor for each line of the (for the C-channel).
    saturation_c : 2D array
        Information about the saturation of the C-channel.
    saturation_vis : 2D array
        Information about the saturation of the Vis-channel.
    summation : int
        The downtrack summing.
    bits_per_data : float
        The compression rate in bits per data.
    data_quality : int
        Information about the data quality, from 0 to 5 depending on missing lines and
        compression errors. (See SOFT09_readme.txt for more details.)
    lrec : int
        The number of bytes in each physical record in the data product file.
    nrec : int
        The number of physical records that make up the PDS product label.
    sol_dist_au : float
        The distance between the center of the observation and the Sun (a.u.).
    npixel : int
        The number of pixels of the length of scan (can be 16, 32, 64, or 128 pixels).
    nscan : int
        Number of scanned pixel lines.
    npara : int
        Number of parameters describing the geometry of the observation.
    point_mode : str
        The pointing mode of the instrument.
    target : str
        The name of the target of the observation ('MARS', 'PHOBOS' or 'DEIMOS') .
    mode_channel : int
        Information about the presence of each channel.
    orient : array
        The vector orientation of the spacecraft.
    subsol_lat : float
        Latitude of the sub-solar point at observation time (deg).
    subsol_lon : float
        Longitude of the sub-solar point at observation time (deg).
    min_lat : float
        Southernmost latitude of the observation (deg).
    max_lat : float
        Northernmost latitude of the observation (deg).
    min_lon : float
        Easternmost longitude of the observation (deg).
    max_lon : float
        Westernmost longitude of the observation (deg).
    slant : float
        Distance from the spacecraft to the center of the observation along the line of sight (km).
    quality : int
        The quality level of the cube.
        (0: corrupted | 1: good | 128 : corrupted mode 128)
    therm_corr : bool
        | True -> Thermal correction applied.
        | False -> No thermal correction.
    therm_corr_infos : dict
        Information about the thermal correction (date, method).
    atm_corr : bool
        | True -> Atmospheric correction applied.
        | False -> No atmospheric correction.
    atm_corr_infos : dict
        Information about the atmospheric correction (date, method).
    version : int
        The major release version of the omegapy.omega_data.py file used.
    add_infos : str
        Additional informations about the observation.
        Show in the OMEGAdata representation.
    """

    def __init__(self, obs='', empty=False, data_path="_omega_bin_path", corrV=True, corrL=True, disp=True):
        # Default paths
        if data_path == "_omega_bin_path":
            data_path = _omega_bin_path
        # Infos
        self.version = int(_Version)
        self.therm_corr = False
        self.atm_corr = False
        self.therm_corr_infos = {'datetime': None, 'method': None}
        self.atm_corr_infos = {'datetime': None, 'method': None}
        self.quality = 1
        self.add_infos = ''

        if not empty:
            obs_name = uf.myglob(os.path.join(data_path, '*' + obs + '*.QUB'))
            if obs_name is None:
                print("\033[1;33mAborted\033[0m")
                empty = True

        if empty:
            # Data
            self.name = ''
            self.lam = np.array([])
            self.cube_i = np.array([[[]]])
            self.cube_rf = np.array([[[]]])
            self.ls = np.nan
            self.lat = np.array([[]])
            self.lon = np.array([[]])
            self.alt = np.array([[]])
            self.loct = np.array([[]])
            self.emer = np.array([[]])
            self.inci = np.array([[]])
            self.specmars = np.array([])
            self.utc = datetime.datetime.now()
            self.my = np.nan
            self.orbit = None
            self.surf_temp = np.array([[]])
            self.ic = {'V' : np.arange(265, 333),
                       'C' : np.arange(8, 123),
                       'L' : np.arange(137, 256)}
            self.lon_grid = np.array([[]])
            self.lat_grid = np.array([[]])
            self.sensor_temp_c = np.array([])
            self.saturation_c = np.array([[]])
            self.saturation_vis = np.array([[]])
            self.surf_temp = np.array([[]])
            self.summation = None
            self.bits_per_data = None
            self.data_quality = None
            self.mode_channel = None
            # Nav
            self.lrec = None
            self.nrec = None
            self.sol_dist_au = None
            self.npixel = None
            self.npara = None
            self.nscan = None
            self.point_mode = None
            self.orient = np.array([])
            self.subsol_lon = None
            self.subsol_lat = None
            self.min_lat = None
            self.max_lat = None
            self.min_lon = None
            self.max_lon = None
            self.slant = None
            self.target = None

        else:
            # obs_name = uf.myglob(os.path.join(data_path, '*' + obs + '*.QUB'))
            # if obs_name is None:
                # print("\033[1;33mAborted\033[0m")
                # return None
            nomfic0 = os.path.split(obs_name)[1][:-4]    # Récupération nom + décodage UTF-8
            data_dict, geom_dict = _readomega(nomfic0, disp=disp, corrV=corrV, corrL=corrL, 
                                              data_path=data_path)

            if disp:
                print("\n\033[01;34mComputing data extraction and correction...\033[0m", end=' ')
            # Extract values
            ldat = data_dict['ldat']
            jdat = data_dict['jdat']
            wvl = data_dict['wvl']
            ic = data_dict['ic']
            specmars = data_dict['specmars']
            lat = geom_dict['latitude']
            lon = geom_dict['longitude']
            alt = geom_dict['altitude']
            emer = geom_dict['emergence']
            inci = geom_dict['incidence']
            utc = geom_dict['ut_time']
            temperature = geom_dict['temperature']
            saturation_c = geom_dict['saturation_c']
            saturation_vis = geom_dict['saturation_vis']
            lon_px = geom_dict['lon_grid']
            lat_px = geom_dict['lat_grid']
            # Correction of OMEGA data (same as clean_spec.pro)
            ic_C= ic[(ic >= 8) & (ic <= 122)]        # IR short (voie C)
            ic_L = ic[(ic >= 137) & (ic <= 255)]     # IR long (voie L)
            ic_V = ic[(ic >= 265) & (ic <= 332)]     # visible
            ic2 = np.concatenate([ic_V, ic_C, ic_L])
            lam = wvl[ic2]
            specmars = specmars[ic2]
            lam_mask_all = deepcopy(lam)
            lam_mask_all[:] = True
            lam_mask = {'all' : lam_mask_all}
            orbit_nb = int(nomfic0[3:-2])
            # Cube in physical units (W.m-2.sr-1.µm-1)
            cube_i = jdat[:, ic2, :]
            # Cube of reflectance factor I/F
            cube_rf = ldat[:, ic2, :]
            # Cube as [X, Y, lam]
            cube_i2 = np.swapaxes(cube_i, 1, 2)
            cube_rf2 = np.swapaxes(cube_rf, 1, 2)
            # Observation UTC date & time
            Y, M, D, h, m, s = np.median(utc[:,:6], axis=0).astype(np.int64)
            utc_dt = datetime.datetime(Y, M, D, h, m, s)
            # Longitude pixels grid
            ny, nx = lon.shape
            lon_grid = np.zeros((ny+1, nx+1))
            lon_grid[1:,1:] = lon_px[:,:,0]
            lon_grid[1:,0] = lon_px[:,0,1]
            lon_grid[0,1:] = lon_px[0,:,3]
            lon_grid[0,0] = lon_px[0,0,2]
            # Longitude pixels grid
            lat_grid = np.zeros((ny+1, nx+1))
            lat_grid[1:,1:] = lat_px[:,:,0]
            lat_grid[1:,0] = lat_px[:,0,1]
            lat_grid[0,1:] = lat_px[0,:,3]
            lat_grid[0,0] = lat_px[0,0,2]
            # Storage as class arguments
            self.lam = lam.astype(np.float64)
            self.cube_i = cube_i2.astype(np.float64)
            self.cube_rf = cube_rf2.astype(np.float64)
            self.lat = lat.astype(np.float64)
            self.lon = lon.astype(np.float64)
            self.alt = alt.astype(np.float64)
            self.emer = emer.astype(np.float64)
            self.inci = inci.astype(np.float64)
            self.specmars = specmars.astype(np.float64)
            self.name = nomfic0
            self.orbit = orbit_nb
            self.utc = utc_dt
            self.ic = {'V' : ic_V.astype(int), 
                       'C' : ic_C.astype(int), 
                       'L' : ic_L.astype(int)}
            self.lam_ma = lam_mask
            self.lat_grid = lat_grid
            self.lon_grid = lon_grid
            self.sensor_temp_c = temperature.astype(np.float64)
            self.saturation_c = saturation_c.astype(np.float64)
            self.saturation_vis = saturation_vis.astype(np.float64)
            #--------------------------
            # Data from the .QUB header
            #--------------------------
            hd_qub = _read_header(obs_name[:-4] + '.QUB')
            self.summation = np.int64(hd_qub['DOWNTRACK_SUMMING'])
            self.bits_per_data = np.float64(hd_qub['INST_CMPRS_RATE'])
            self.data_quality = np.int64(hd_qub['DATA_QUALITY_ID'])
            mode_channel_tmp = hd_qub['COMMAND_DESC'][34:36]
            if mode_channel_tmp == 'EF':
                self.mode_channel = 1
            elif mode_channel_tmp == '80':
                self.mode_channel = 2
            elif mode_channel_tmp == 'C7':
                self.mode_channel = 3
            else:
                self.mode_channel = mode_channel_tmp
            #--------------------------
            # Data from the .NAV header
            #--------------------------
            hd_nav = _read_header(obs_name[:-4] + '.NAV')
            npixel, npara, nscan = np.array(hd_nav['CORE_ITEMS'][1:-1].split(','), dtype=np.int64)
            self.lrec = np.int64(hd_nav['RECORD_BYTES'])
            self.nrec = np.int64(hd_nav['LABEL_RECORDS'])
            self.sol_dist_au = np.float64(hd_nav['SOLAR_DISTANCE']) / 14960e4
            npixel, npara, nscan = np.array(hd_nav['CORE_ITEMS'][1:-1].split(','), dtype=np.int64)
            self.npixel = npixel
            self.npara = npara
            self.nscan = nscan
            self.point_mode = hd_nav['SPACECRAFT_POINTING_MODE'][1:-1]
            self.orient = np.array(hd_nav['SPACECRAFT_ORIENTATION'][1:-1].split(','), dtype=np.int64)
            self.ls = np.float64(hd_nav['SOLAR_LONGITUDE'])
            self.subsol_lon = np.float64(hd_nav['SUB_SOLAR_LONGITUDE'])
            self.subsol_lat = np.float64(hd_nav['SUB_SOLAR_LATITUDE'])
            self.min_lat = np.float64(hd_nav['MINIMUM_LATITUDE'])
            self.max_lat = np.float64(hd_nav['MAXIMUM_LATITUDE'])
            self.min_lon = np.float64(hd_nav['WESTERNMOST_LONGITUDE'])
            self.max_lon = np.float64(hd_nav['EASTERNMOST_LONGITUDE'])
            self.slant = np.float64(hd_nav['SLANT_DISTANCE'])
            self.target = hd_nav['TARGET_NAME']
            #--------------------------
            temp_init = np.zeros(self.lat.shape)
            temp_init[:] = np.nan
            self.surf_temp = temp_init
            #--------------------------
            # Local time & MY
            self.loct = _compute_local_time(self.lon, self.subsol_lon)
            self.my = _utc_to_my(self.utc)
            #--------------------------
            # Cube quality
            OBC = readsav(os.path.join(package_path, 'OMEGA_dataref', 'OBC_OMEGA_OCT2017.sav'))
            good_orbits_OBC = np.array(OBC['good_orbits'][0], dtype=int)
            corrupted_orbits_csv = pd.read_csv(os.path.join(package_path, 'OMEGA_dataref', 'corrupted_obs.csv'), comment='#',
                                                skipinitialspace=True)
            corrupted_orbits = np.array(corrupted_orbits_csv['corrupted_obs'], dtype=str)
            corrupted_orbits_comments = np.array(corrupted_orbits_csv['comment'], dtype=str)
            if (npixel==128) & (orbit_nb >= 513):
                self.quality = 128
                self.add_infos = 'Corrupted 128 pixels cube'
            if orbit_nb not in good_orbits_OBC:
                self.quality = 0
                self.add_infos = 'Corrupted orbit'
            if nomfic0 in corrupted_orbits:
                self.quality = 0
                i_obs = int(np.where(corrupted_orbits==nomfic0)[0])
                self.add_infos = corrupted_orbits_comments[i_obs]
            #--------------------------
            # End of data extraction & correction
            if disp:
                print("\033[01;32m[done]\033[0m")


    def __copy__(self):
        new_omega = OMEGAdata(empty=True)
        # Data
        new_omega.name = self.name
        new_omega.lam = self.lam
        new_omega.cube_i = self.cube_i
        new_omega.cube_rf = self.cube_rf
        new_omega.ls = self.ls
        new_omega.lat = self.lat
        new_omega.lon = self.lon
        new_omega.alt = self.alt
        new_omega.loct = self.loct
        new_omega.my = self.my
        new_omega.emer = self.emer
        new_omega.inci = self.inci
        new_omega.specmars = self.specmars
        new_omega.utc = self.utc
        new_omega.orbit = self.orbit
        new_omega.ic = self.ic
        new_omega.lam_ma = self.lam_ma
        new_omega.lon_grid = self.lon_grid
        new_omega.lat_grid = self.lat_grid
        new_omega.sensor_temp_c = self.sensor_temp_c
        new_omega.saturation_c = self.saturation_c
        new_omega.saturation_vis = self.saturation_vis
        new_omega.surf_temp = self.surf_temp
        new_omega.summation = self.summation
        new_omega.bits_per_data = self.bits_per_data
        new_omega.data_quality = self.data_quality
        new_omega.mode_channel = self.mode_channel
        # Nav
        new_omega.lrec = self.lrec
        new_omega.nrec = self.nrec
        new_omega.sol_dist_au = self.sol_dist_au
        new_omega.npixel = self.npixel
        new_omega.npara = self.npara
        new_omega.nscan = self.nscan
        new_omega.point_mode = self.point_mode
        new_omega.orient = self.orient
        new_omega.subsol_lon = self.subsol_lon
        new_omega.subsol_lat = self.subsol_lat
        new_omega.min_lat = self.min_lat
        new_omega.max_lat = self.max_lat
        new_omega.min_lon = self.min_lon
        new_omega.max_lon = self.max_lon
        new_omega.slant = self.slant
        new_omega.target = self.target
        # Infos
        new_omega.quality = self.quality
        new_omega.therm_corr = self.therm_corr
        new_omega.atm_corr = self.atm_corr
        new_omega.therm_corr_infos = self.therm_corr_infos
        new_omega.atm_corr_infos = self.atm_corr_infos
        new_omega.version = self.version
        new_omega.add_infos = self.add_infos
        return new_omega

    def __deepcopy__(self, memo):
        new_omega = OMEGAdata(empty=True)
        memo[id(self)] = new_omega
        # Data
        new_omega.name = deepcopy(self.name, memo)
        new_omega.lam = deepcopy(self.lam, memo)
        new_omega.cube_i = deepcopy(self.cube_i, memo)
        new_omega.cube_rf = deepcopy(self.cube_rf, memo)
        new_omega.ls = deepcopy(self.ls, memo)
        new_omega.lat = deepcopy(self.lat, memo)
        new_omega.lon = deepcopy(self.lon, memo)
        new_omega.alt = deepcopy(self.alt, memo)
        new_omega.loct = deepcopy(self.loct, memo)
        new_omega.my = deepcopy(self.my, memo)
        new_omega.emer = deepcopy(self.emer, memo)
        new_omega.inci = deepcopy(self.inci, memo)
        new_omega.specmars = deepcopy(self.specmars, memo)
        new_omega.utc = deepcopy(self.utc, memo)
        new_omega.orbit = deepcopy(self.orbit, memo)
        new_omega.ic = deepcopy(self.ic, memo)
        new_omega.lam_ma = deepcopy(self.lam_ma, memo)
        new_omega.lon_grid = deepcopy(self.lon_grid, memo)
        new_omega.lat_grid = deepcopy(self.lat_grid, memo)
        new_omega.sensor_temp_c = deepcopy(self.sensor_temp_c, memo)
        new_omega.saturation_c = deepcopy(self.saturation_c, memo)
        new_omega.saturation_vis = deepcopy(self.saturation_vis, memo)
        new_omega.surf_temp = deepcopy(self.surf_temp, memo)
        new_omega.summation = deepcopy(self.summation, memo)
        new_omega.bits_per_data = deepcopy(self.bits_per_data, memo)
        new_omega.data_quality = deepcopy(self.data_quality, memo)
        new_omega.mode_channel = deepcopy(self.mode_channel, memo)
        # Nav
        new_omega.lrec = deepcopy(self.lrec, memo)
        new_omega.nrec = deepcopy(self.nrec, memo)
        new_omega.sol_dist_au = deepcopy(self.sol_dist_au, memo)
        new_omega.npixel = deepcopy(self.npixel, memo)
        new_omega.npara = deepcopy(self.npara, memo)
        new_omega.nscan = deepcopy(self.nscan, memo)
        new_omega.point_mode = deepcopy(self.point_mode, memo)
        new_omega.orient = deepcopy(self.orient, memo)
        new_omega.subsol_lon = deepcopy(self.subsol_lon, memo)
        new_omega.subsol_lat = deepcopy(self.subsol_lat, memo)
        new_omega.min_lat = deepcopy(self.min_lat, memo)
        new_omega.max_lat = deepcopy(self.max_lat, memo)
        new_omega.min_lon = deepcopy(self.min_lon, memo)
        new_omega.max_lon = deepcopy(self.max_lon, memo)
        new_omega.slant = deepcopy(self.slant, memo)
        new_omega.target = deepcopy(self.target, memo)
        # Infos
        new_omega.quality = deepcopy(self.quality, memo)
        new_omega.therm_corr = deepcopy(self.therm_corr, memo)
        new_omega.atm_corr = deepcopy(self.atm_corr, memo)
        new_omega.therm_corr_infos = deepcopy(self.therm_corr_infos, memo)
        new_omega.atm_corr_infos = deepcopy(self.atm_corr_infos, memo)
        new_omega.version = deepcopy(self.version, memo)
        new_omega.add_infos = deepcopy(self.add_infos, memo)
        return new_omega

    def __eq__(self, other):
        if isinstance(other, OMEGAdata):
            return ((self.name == other.name) and
                    (self.ls == other.ls) and
                    (self.my == other.my) and
                    (self.therm_corr == other.therm_corr) and 
                    (self.atm_corr == other.atm_corr) and
                    (self.cube_rf == other.cube_rf).all())
        else:
            return False

    def __repr__(self):
        description = """
        OMEGA/MEx observation {0} – (v{1})

        Ls = {2:.1f}° – MY {3:.0f}
        
        Cube quality : {4}
        Thermal correction : {5}
        Atmospheric correction : {6}

        \033[3m{7}\033[0m""".format(self.name, self.version, self.ls, self.my, self.data_quality, 
                                    self.therm_corr, self.atm_corr, self.add_infos)
        return description
    
    def get_header_qub(self, data_path='_omega_bin_path'):
        """Return the data from the header of the .QUB file, as a dictionary.

        See the OMEGA ECAID for informations about the header entries.
        
        Parameters
        ==========
        data_path : str, optional (default _omega_bin_path)
            The path of the directory containing the data (.QUB) files.

        Returns
        =======
        hd_qub : dict
            Dictionary containing the data from the ORBXXXX_X.QUB file.
        """
        # Default path
        if data_path == "_omega_bin_path":
            data_path = _omega_bin_path
        qub_path = os.path.join(data_path, self.name+'.QUB')
        hd_qub = _read_header(qub_path)
        return hd_qub

    def get_header_nav(self, data_path='_omega_bin_path'):
        """Return the data from the header of the .NAV file, as a dictionary.

        See the OMEGA ECAID for informations about the header entries.
        
        Parameters
        ==========
        data_path : str, optional (default _omega_bin_path)
            The path of the directory containing the navigation (.NAV) files.

        Returns
        =======
        hd_nav : dict
            Dictionary containing the data from the ORBXXXX_X.NAV file.
        """
        # Default path
        if data_path == "_omega_bin_path":
            data_path = _omega_bin_path
        nav_path = os.path.join(data_path, self.name+'.NAV')
        hd_nav = _read_header(nav_path)
        return hd_nav

##-----------------------------------------------------------------------------------
## Recherche observation
def find_cube(lon0, lat0, cmin=0, cmax=10000, out=False, data_path='_omega_bin_path'):
    """Display the available OMEGA/MEx cubes with observations of the target
    latitude and longitude, Python translation of the IDL procedure `findcub.pro`.

    Parameters
    ==========
    lon0 : float
        The target longitude (in degrees).
    lat0 : float
        The target latitude (in degrees).
    cmin : float, optional (default 0)
        The minimum orbit number.
    cmax : float, optional (default 10000)
        The maximum orbit number.
    out : bool, optional (default False)
        If True -> return output
    data_path : str, optional (default _omega_bin_path)
        The path of the directory containing the data (.QUB) and 
        navigation (.NAV) files.

    Returns (If out == True)
    =======
    cub_list : array-like
        List of matching observations.
        Format : (orbit, x, y, dmin, altMEx, inci, emer, phas, loct, Ls, MY)
    """
    # Default path
    if data_path == "_omega_bin_path":
        data_path = _omega_bin_path
    #-----------------------------
    # Internal function testin
    def testin(x0, y0, x1, y1):
        """Internal function for find_cube: test if the point of coordinates 
        (x0, y0) is include in the (x1, y1) grid.

        Parameters
        ==========
        x0 : float
            X-coordinate of the point.
        y0 : float
            Y-coordinate of the point.
        x1 : 1D array
            X-coordinates of the observation.
        y1 : 1D array
            Y-coordinates of the observations.

        Returns
        =======
        res : bool
            | True if the point (x0, y0) is in the observation grid.
            | False if it's not.
        """
        nb = len(x1)
        x2 = np.concatenate([x1, [x1[0]]])
        y2 = np.concatenate([y1, [y1[0]]])
        dx = x2 - x0
        dy = y2 - y0
        atot = 0
        for n in range(nb):
            ps = dx[n] * dx[n+1] + dy[n] * dy[n+1]
            pv = dx[n] * dy[n+1] - dx[n+1] * dy[n]
            atot = atot + np.arctan2(pv, ps)
        if np.abs(atot) > 3:
            return True
        else:
            return False
    #-----------------------------
    # Initialization
    res_folder = os.path.join(package_path, 'res_findcube')
    trans = np.pi / 180
    x0 = np.cos(lon0*trans) * np.cos(lat0*trans)
    y0 = np.sin(lon0*trans) * np.cos(lat0*trans)
    z0 = np.sin(lat0*trans)

    x1, y1 = np.zeros((2, 10))
    nomc = []

    nomout = os.path.join(res_folder, 'orbites_lg{lon:.0f}lt{lat:.0f}.dat'.format(lon=lon0, lat=lat0))
    path_cubliste = os.path.join(res_folder, 'cubelist')
    path_cubindex = os.path.join(package_path, 'OMEGA_dataref', 'cubindex.ref')
    with open(path_cubliste, 'w') as f:
        f.write('# long: {lon:7.3f}  lat: {lat:7.3f}\n\n'.format(lon=lon0, lat=lat0))
        f.write('#{0:^9s} {1:^6s}{2:^6s}{3:^8s}{4:^9s}{5:^7s}{6:^8s}{7:^8s}{8:^8s}{9:^8s}{10:^4s}\n'.format(
                'orbit', 'x', 'y', 'dmin', 'altMEx', 'inci', 'emer', 'phas', 'loct', 'Ls', 'MY'))
    with open(nomout, 'w') as f:
        f.write('')
    nhits = 0
    geocube = 0
    cubindex = open(path_cubindex, 'rb')
    # Search for matching observations
    # South pole
    if lat0 <= -60:
        for ncube in range(10000):
            nomcube = cubindex.readline().decode('utf8').replace('\n', '')
            norb = int(nomcube[3:7])
            if norb == 0:
                break   # End of file
            lon1 = np.fromstring(cubindex.readline().decode('utf8').replace('\n', ''), sep=' ')
            lat1 = np.fromstring(cubindex.readline().decode('utf8').replace('\n', ''), sep=' ')
            if (norb < cmin) or (norb > cmax):
                continue
            if np.min(lat1) > -60:
                continue
            x1 = np.cos(lon1*trans) * np.cos(lat1*trans)
            y1 = np.sin(lon1*trans) * np.cos(lat1*trans)
            if testin(x0, y0, x1, y1):
                nomc.append(nomcube)
                nhits += 1
    # North pole
    elif lat0 >= 60:
        for ncube in range(10000):
            nomcube = cubindex.readline().decode('utf8').replace('\n', '')
            norb = int(nomcube[3:7])
            if norb == 0:
                break   # End of file
            lon1 = np.fromstring(cubindex.readline().decode('utf8').replace('\n', ''), sep=' ')
            lat1 = np.fromstring(cubindex.readline().decode('utf8').replace('\n', ''), sep=' ')
            if (norb < cmin) or (norb > cmax):
                continue
            if np.max(lat1) < 60:
                continue
            x1 = np.cos(lon1*trans) * np.cos(lat1*trans)
            y1 = np.sin(lon1*trans) * np.cos(lat1*trans)
            if testin(x0, y0, x1, y1):
                nomc.append(nomcube)
                nhits += 1
    # Intermediate region
    else:
        for ncube in range(10000):
            nomcube = cubindex.readline().decode('utf8').replace('\n', '')
            norb = int(nomcube[3:7])
            if norb == 0:
                break   # End of file
            lon1 = np.fromstring(cubindex.readline().decode('utf8').replace('\n', ''), sep=' ')
            lat1 = np.fromstring(cubindex.readline().decode('utf8').replace('\n', ''), sep=' ')
            if (norb < cmin) or (norb > cmax):
                continue
            x1 = np.cos(lon1*trans) * np.cos(lat1*trans)
            y1 = np.sin(lon1*trans) * np.cos(lat1*trans)
            z1 = np.sin(lat1*trans)
            ps = x0*x1 + y0*y1 + z0*z1
            if np.max(ps) < 0.85:
                continue
            lon2 = lon1 - lon0
            i1 = np.where(lon2 < -180)[0]
            if len(i1) > 0:
                lon2[i1] += 360
            i2 = np.where(lon2 > 180)[0]
            if len(i2) > 0:
                lon2[i2] -= 360
            if testin(0, lat0, lon2, lat1):
                nomc.append(nomcube)
                nhits += 1
    cubindex.close()
    # Find position & infos for each observation
    print('{0:^10s} {1:^6s}{2:^6s}{3:^8s}{4:^9s}{5:^7s}{6:^8s}{7:^8s}{8:^8s}{9:^8s}{10:^4s}'.format(
            'orbit', 'x', 'y', 'dmin', 'altMEx', 'inci', 'emer', 'phas', 'loct', 'Ls', 'MY'))
    for n in range(nhits):
        testfile = os.path.join(data_path, nomc[n]+'.NAV')
        if os.path.exists(testfile) == False:
            print('{0:8s}{1:s}'.format(nomc[n], '\033[3m   No corresponding .NAV file\033[0m'))
            continue
        #--------------------------
        # Data from the geometry .NAV file
        #--------------------------
        hd_nav = _read_header(testfile)
        npixel, npara, nscan = np.array(hd_nav['CORE_ITEMS'][1:-1].split(','), dtype=np.int64)
        lrec = np.int64(hd_nav['RECORD_BYTES'])
        nrec = np.int64(hd_nav['LABEL_RECORDS'])
        solong = np.float64(hd_nav['SOLAR_LONGITUDE'])
        sslong = np.float64(hd_nav['SUB_SOLAR_LONGITUDE'])
        #--------------------------
        # Lecture géometrie
        class F_line_nav(ctypes.Structure):
            _fields_ = [('C_line', ctypes.c_int32 * npixel)]

        class F_frame_nav(ctypes.Structure):
            _fields_ = [('F_line', F_line_nav * npara)]

        class F_Qube_nav(ctypes.Structure):
            _fields_ = [('F_frame', F_frame_nav * nscan)]

        fqube_nav = F_Qube_nav()
        skip = lrec * nrec  # taille header (en bytes)
        geocube = np.ndarray((nscan, npara, npixel), np.int32)
        with open(testfile, 'rb') as nav_cube:
            data_hd = nav_cube.read(skip)  # bytes header
            data_qub = nav_cube.readinto(fqube_nav)     # bytes data

        fqube_nav2 = np.ctypeslib.as_array(fqube_nav.F_frame)
        geocube[:,:,:] = fqube_nav2['F_line']['C_line']
        # On remet dans le même sens que readomega
        geocube = np.swapaxes(geocube, 0, 2)
        #--------------------------
        longa = geocube[:, 6, :] * 1e-4
        lata = geocube[:, 7, :] * 1e-4
        xa = np.cos(longa*trans) * np.cos(lata*trans)
        ya = np.sin(longa*trans) * np.cos(lata*trans)
        za = np.sin(lata*trans)
        xr = y0 * za - z0 * ya
        yr = z0 * xa - x0 * za
        zr = x0 * ya - y0 * xa
        dist = np.sqrt(xr*xr + yr*yr + zr*zr) * 3393
        distmin = np.min(dist)
        [i0], [j0] = np.where(dist == distmin)
        inci = geocube[i0, 2, j0] * 1e-4
        emer = geocube[i0, 3, j0] * 1e-4
        phas = geocube[i0, 10, j0] * 1e-4
        slant = geocube[i0, 11, j0] * 1e-3
        alt = geocube[i0, 12, j0] * 1e-3
        possible_geom_corruption = False
        try:
            Y, M, D, h, m, s = geocube[:6, 1, j0]
            utc_dt = datetime.datetime(Y, M, D, h, m, s)
        except:     # Possible corruption of some geometry lines
            Y, M, D, h, m, s = np.median(geocube[:6, 1, :], axis=1).astype(np.int64)
            utc_dt = datetime.datetime(Y, M, D, h, m, s)
            possible_geom_corruption = True
        my = _utc_to_my(utc_dt)
        loct = _compute_local_time(longa, sslong)[i0, j0]
        obs_output = '{0:9s}{1:6d}{2:6d}{3:8.2f}{4:9.1f}{5:8.2f}{6:8.2f}{7:8.2f}{8:8.2f}{9:8.2f}{10:4d}'.format(
                        nomc[n], i0, j0, distmin, slant, inci, emer, phas, loct, solong, my)
        if possible_geom_corruption:
            obs_output = '\033[3m' + obs_output + '\033[0m'
        print(obs_output)

        with open(path_cubliste, 'a') as f_cublist:
            f_cublist.write('{0:9s}{1:6d}{2:6d}{3:8.2f}{4:9.1f}{5:8.2f}{6:8.2f}{7:8.2f}{8:8.2f}{9:8.2f}{10:4d}\n'.format(
                    nomc[n], i0, j0, distmin, slant, inci, emer, phas, loct, solong, my))
        with open(nomout, 'a') as f_listeobs:
            f_listeobs.write('{0:s}\n'.format(nomc[n]))
    # Output (if out==True)
    if out:
        cub_list = np.genfromtxt(path_cubliste, skip_header=3,
                                 dtype=None, encoding='utf8')
        return cub_list

##-----------------------------------------------------------------------------------
## Sauvegarde / Importation
def save_omega(omega, savname='auto', folder='', base_folder='_omega_py_path',
               pref ='', suff='', disp=True):
    """Save an OMEGA object at the selected path using the pickle module.

    Final_path = base_folder + folder + savname

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA/MEx observation object.
    savname : str, optional (default 'auto')
        The saving filename.
        | If 'auto' -> savname = 'pref_omega.name_ext.pkl'
    folder : str, optional (default '')
        The subfolder to save the data.
    base_folder : str, optional (default _omega_py_path)
        The base folder path.
    pref : str, optional (default '')
        The prefix of the savname.
    suff : str, optional (default '')
        The suffix of the savname.
    disp : bool
        Control the display.
            | True -> Print the saving filename.
            | False -> Nothing printed.
    """
    # Default path
    if base_folder == "_omega_py_path":
        base_folder = _omega_py_path
    # Initialisation nom fichier auto
    if savname == 'auto':
        if (len(suff)>0) and (suff[0] != '_'):
            suff = '_' + suff
        if (len(pref)>0) and (pref[-1] != '_'):
            pref = pref + '_'
        savname = '{pref}{name}{suff}.pkl'.format(name=omega.name, pref=pref, suff=suff)
    # Chemin sav fichier
    target_path = os.path.join(base_folder, folder, savname)
    # Sauvegarde pickle
    with open(target_path, 'wb') as output:
        pickle.dump(omega, output)
    if disp:
        print('\033[01;34mSaved as \033[0;03m' + target_path + '\033[0m')

def load_omega(filename, disp=True):
    """Load and return a previously saved OMEGAdata object (with save_omega()).

    Parameters
    ==========
    filename : str
        The file path.
    disp : bool
        Control the display.
            | True -> Print the loading filename.
            | False -> Nothing printed.

    Returns
    =======
    omega : OMEGAdata 
        The loaded object of OMEGA/MEx observation.
    """
    filename2 = uf.myglob(filename)
    with open(filename2, 'rb') as input_file:
        omega = pickle.load(input_file)
        if disp:
            print('\033[03m' + filename2 + '\033[0;01;34m loaded\033[0m')
        return omega

def load_omega_list(basename, disp=True):
    """Load a list of saved OMEGAdata objects, using load_omega().

    Parameters
    ==========
    basename : str
        The file path basename.
    disp : bool
        Control the display.
            | True -> Print the loading filename.
            | False -> Nothing printed.

    Returns
    =======
    omega_list : ndarray of OMEGAdata objects
        The array of loaded objects of OMEGA/MEx observation.
    """
    path_list = glob.glob(basename)
    omega_list = []
    for i in range(len(path_list)):
        omega_list.append(load_omega(path_list[i], disp))
    return np.array(omega_list)

def autosave_omega(omega, folder='auto', base_folder='_omega_py_path', security=True, disp=True):
    """Save an OMEGA object at the selected path using the pickle module, with automatic
    configuration of the target name.

    Final_path = base_folder + folder + name{_corr_therm_atm}.pkl

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA/MEx observation object.
    folder : str, optional (default 'auto')
        The subfolder to save the data.
        | If 'auto' -> folder = 'vX', where X is the major release version of the used code.
    base_folder : str, optional (default _omega_py_path)
        The base folder path.
    security : bool, optional (default True)
        Enable / disable checking before overwriting a file.
        | True -> Check if the target file already exists before overwriting on it.
                  And if is the case, you will be asked for a confirmation.
        | False -> Didn't care about the already existing files.
    disp : bool
        Control the display.
            | True -> Print the saving filename.
            | False -> Nothing printed.
    """
    # Default path
    if base_folder == "_omega_py_path":
        base_folder = _omega_py_path
    # Initialisation nom fichier auto
    if omega.therm_corr and omega.atm_corr:
        suff = '_corr_therm_atm'
    elif omega.therm_corr:
        suff = '_corr_therm'
    elif omega.atm_corr:
        suff = '_corr_atm'
    else:
        suff = ''
    savname = '{name}{suff}.pkl'.format(name=omega.name, suff=suff)
    if folder == 'auto':
        folder = 'v' + str(int(omega.version))
    # Chemin sav fichier
    target_path = os.path.join(base_folder, folder, savname)
    # Testing existent file
    if security:
        write = uf.test_security_overwrite(target_path)
    else:
        write = True
    # Sauvegarde pickle
    if write:
        with open(target_path, 'wb') as output:
            pickle.dump(omega, output)
        if disp:
            print('\033[01;34mSaved as \033[0;03m' + target_path + '\033[0m')

def autoload_omega(obs_name, folder='auto', version=_Version, base_folder='_omega_py_path',
                   therm_corr=None, atm_corr=None, disp=True, bin_folder='_omega_bin_path'):
    """Load and return a previously saved OMEGAdata object using pickle (with autosave_omega()).

    Parameters
    ==========
    obs_name : str
        The observation ID.
    folder : str, optional (default 'auto')
        The subfolder where the data is.
        | If 'auto' -> folder = 'vX', where X is the major release version of the used code.
    version : float, optional (default _Version)
        The version of the target file (if folder is 'auto').
    base_folder : str, optional (default _omega_py_path)
        The base folder path.
    therm_corr : bool or None, optional (default None)
        | True -> Only results with thermal correction.
        | False -> Only results without thermal correction.
        | None -> Both with and without thermal correction.
    atm_corr : bool or None, optional (default None)
        | True -> Only results with atmospheric correction.
        | False -> Only results without atmospheric correction.
        | None -> Both with and without atmospheric correction.
    disp : bool
        Control the display.
            | True -> Print the loading filename.
            | False -> Nothing printed.
    bin_folder : str, optional (default _omega_bin_path)
        The path of the directory containing the data (.QUB) and 
        navigation (.NAV) files.

    Returns
    =======
    omega : OMEGAdata 
        The loaded object of OMEGA/MEx observation.
    """
    # Default paths
    if base_folder == "_omega_py_path":
        base_folder = _omega_py_path
    if bin_folder == "_omega_bin_path":
        bin_folder = _omega_bin_path
    # Initialisation
    ext = ''
    excl = []
    if therm_corr:
        ext += '_therm'
    elif therm_corr == False:
        excl.append('therm')
    if atm_corr:
        ext += '_atm'
    elif atm_corr == False:
        excl.append('atm')
    filename = '*{name}*{corr_ext}*.pkl'.format(name=obs_name, corr_ext=ext)
    if folder == 'auto':
        Mversion = int(version)
        folder = 'v' + str(Mversion)
    filename2 = uf.myglob(os.path.join(base_folder, folder, filename), exclude=excl)
    if filename2 is None:
        if (therm_corr in [None, False]) and (atm_corr in [None, False]):
            obs_name_bin = glob.glob(os.path.join(bin_folder, '*' + obs_name + '*.QUB'))
            if len(obs_name_bin) == 0 :
                return None
            else:
                print('\033[1mMatching binary files:\033[0m')
                return OMEGAdata(obs_name, data_path=bin_folder)
        else:
            return None
    else:
        with open(filename2, 'rb') as input_file:
            omega = pickle.load(input_file)
            if disp:
                print('\033[03m' + filename2 + '\033[0;01;34m loaded\033[0m')
            return omega

##-----------------------------------------------------------------------------------
## Correction thermique 1 - Détermitation réflectance à 5µm à partir de celle à 2.4µm
## puis CN -- Méthode historique (Calvin & Erard 1997)
_omega_tmp = None
def _corr_therm_sp(args):
    """Remove the thermal component in an OMEGA spectrum, using the historical method
    based on the reflectance determination using reference spectra from Calvin & Erard (1997).
    
    Important note : this function is dedicated to internal use in the corr_therm() function,
        as it use non-local objects related to the multiprocessing.

    Parameters
    ==========
    args : tuple of 3 elements
      | x : int
      |     The x-coordinate of the pixel.
      | y : int
      |     The y-coordinate of the pixel.
      | disp : bool, optional (default True)
      |     If True display the fitted temperature/reflectance in the console.

    Returns
    =======
    sp_rf_corr : 1D array
        The reflectance spectrum, corrected from the thermal component.
    T_fit : float
        The retrieved surface temperature (in K).
    x : int
        The x-coordinate of the pixel.
    y : int
        The y-coordinate of the pixel.
    """
    # Extraction arguments
    x, y, disp = args
    # Extraction données
    global _omega_tmp
    omega = _omega_tmp
    lam = omega.lam
    sp_sol = omega.specmars
    sp_rf = omega.cube_rf[y, x]
    sp_i = omega.cube_i[y, x]
    ecl = np.cos(omega.inci[y, x] * np.pi/180)
    # spectels #97-#112 des spectres de ref <-> 2.3-2.5µm
    fref = os.path.join(package_path, 'OMEGA_dataref', 'refclair_sombr_omega_CL.dat') # from Erard and Calvin (1997)
    sp_clair, sp_sombre = np.loadtxt(fref, unpack=True)
    i_lam1, i_lam2 = uf.where_closer_array([2.3, 2.5], lam)
    alb_clair = np.average(sp_clair[97:113])
    alb_sombre = np.average(sp_sombre[97:113])
    alb_C = np.average(sp_rf[i_lam1:i_lam2+1])
    # Simulation d'un spectre : thèse D. Jouglet, p159
    # CL des spectres de référence `clair` et `sombre` (Calvin & Erard 1997)
    coeff = (alb_C - alb_sombre) / (alb_clair - alb_sombre)
    sp_simu = coeff * sp_clair + (1-coeff) * sp_sombre
    # Sélection des spectels 5.03-5.09µm (4 derniers voie L)
    sp_simu_5m = sp_simu[252:256]
    i_lam3, i_lam4 = uf.where_closer_array([5.03, 5.09], lam)
    i_lam4 += 1
    def simu_sp_5microns(i_lams, T):
        i1, i2 = i_lams.astype(int)
        lam2 = lam[i1:i2] * 1e-6    # Conversion en m
        sp_sol2 = sp_sol[i1:i2]
        Blam = uf.planck(lam2, T) * 1e-6    # Loi de Planck en W.m-2.sr-1.µm-1
        sp_simu2 = sp_simu_5m * sp_sol2 * ecl + (1-sp_simu_5m) * Blam
        return sp_simu2
    try:
        # Fit de la température
        T_fit = curve_fit(simu_sp_5microns, (i_lam3, i_lam4), sp_i[i_lam3:i_lam4], 
                            bounds=(0,400))[0][0]
        # Réflectance
        refl = np.average(sp_simu[252:256])
        if disp:
            print('Temperature = {0:.3f} K   |   Reflectance = {1:.5f}'.format(T_fit, refl))
        # Correction thermique spectre
        Blam = uf.planck(lam*1e-6, T_fit) * 1e-6   # En W.m-2.sr-1.µm-1
        sp_rf_corr = (sp_i - Blam) / (sp_sol*ecl - Blam)
    except ValueError:
        # Si fit impossible (infs ou NaN dans le spectre) -> NaN partout
        sp_rf_corr = deepcopy(sp_rf)
        sp_rf_corr.fill(np.nan)
        T_fit = np.nan
    # Output
    return sp_rf_corr, T_fit, x, y

def corr_therm(omega, npool=1):
    """Remove the thermal component in the OMEGA hyperspectral cube.
    
    Parallelization is implemented using the multiprocessing module. The number of
    process to run is controlled by the npool argument.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA observation data.
    npool : int, optional (default 1)
        Number of parallelized worker process to use.

    Returns
    =======
    omega_corr : OMEGAdata
        The input OMEGA observation, where the reflectance is corrected from
        the thermal component.
    """
    # Test correction
    if omega.therm_corr:
        print("\033[1;33mThermal correction already applied\033[0m")
        return deepcopy(omega)
    # Initialisation
    global _omega_tmp
    _omega_tmp = deepcopy(omega)
    omega_corr = deepcopy(omega)
    ny, nx, nlam = omega.cube_i.shape
    rf_corr = np.zeros((ny,nx,nlam), dtype=np.float64)
    surf_temp = np.zeros((ny,nx), dtype=np.float64)
    # Itérateur
    it = [(x, y, False) for x, y in itertools.product(range(nx), range(ny))]
    # Correction thermique
    # chunksize = len(it) // npool    # Approx size of each process
    chunksize = 1
    # pool = mp.Pool(npool)
    with mp.Pool(npool) as pool:
        for res in tqdm(pool.imap_unordered(_corr_therm_sp, it, chunksize), total=len(it), desc='Thermal correction'):
            sp_rf_corr, T_fit, x, y = res
            rf_corr[y,x] = sp_rf_corr
            surf_temp[y,x] = T_fit
        pool.close()
    _omega_tmp = None
    omega_corr.cube_rf = rf_corr
    omega_corr.surf_temp = surf_temp
    # Update infos
    omega_corr.therm_corr = True
    omega_corr.therm_corr_infos['datetime'] = datetime.datetime.now()
    omega_corr.therm_corr_infos['method'] = '(M1) Calvin & Erard'
    # Sortie
    # tfin = time.time()
    # print('Duration : {0:.0f} min {1:.2f} sec'.format((tfin-tini)//60, (tfin-tini)%60))
    return omega_corr
    
##-----------------------------------------------------------------------------------
## Correction thermique 2 - Fit CN et réflectance en même temps
def corr_therm2_sp(omega, x, y, disp=True):
    """Remove the thermal component in an OMEGA spectrum, with simultaneous retriving
    of reflectance and temperature.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA observation data.
    x : int
        The x-coordinate of the pixel.
    y : int
        The y-coordinate of the pixel.
    disp : bool, optional (default True)
        If True display the fitted temperature/reflectance in the console.

    Returns
    =======
    lam : 1D array
        The wavelength array (in µm).
    sp_rf_corr : 1D array
        The reflectance spectrum, corrected from the thermal component.
    T_fit : float
        The retrieved surface temperature (in K).
    """
    # Test correction
    if omega.therm_corr:
        print("\033[1;33mThermal correction already applied\033[0m")
        return omega.lam, omega.sp_rf[y, x]
    # Extraction données
    lam = omega.lam
    sp_rf = omega.cube_rf[y, x]
    sp_i = omega.cube_i[y, x]
    sp_sol = omega.specmars
    ecl = np.cos(omega.inci[y, x] * np.pi/180)
    # Sélection des spectels 5.03-5.09µm (4 derniers voie L)
    i_lam3, i_lam4 = uf.where_closer_array([5.03, 5.09], lam)
    i_lam4 += 1
    def simu_sp_5microns2(i_lams, T, refl):
        i1, i2 = i_lams.astype(int)
        lam2 = lam[i1:i2] * 1e-6    # Conversion en m
        sp_sol2 = sp_sol[i1:i2]
        Blam = uf.planck(lam2, T) * 1e-6    # Loi de Planck en W.m-2.sr-1.µm-1
        sp_simu2 = refl * sp_sol2 * ecl + (1-refl) * Blam
        return sp_simu2
    try:
        # Fit de la température et réflectance
        T_fit, refl = curve_fit(simu_sp_5microns2, (i_lam3, i_lam4), sp_i[i_lam3:i_lam4], #p0=(280, 0.5),
                                bounds=([0, 0], [400, 0.5]))[0]
        if disp:
            print('Temperature = {0:.3f} K   |   Reflectance = {1:.5f}'.format(T_fit, refl))
        # Correction thermique spectre
        Blam = uf.planck(lam*1e-6, T_fit) * 1e-6   # En W.m-2.sr-1.µm-1
        sp_rf_corr = (sp_i - Blam) / (sp_sol*ecl - Blam)
    except ValueError:
        # Si fit impossible (infs ou NaN dans le spectre) -> NaN partout
        sp_rf_corr = deepcopy(sp_rf)
        sp_rf_corr.fill(np.nan)
        T_fit = np.nan
    return lam, sp_rf_corr, T_fit

def corr_therm2(omega):
    """Remove the thermal component in the OMEGA hyperspectral cube, 
    with simultaneous retriving of reflectance and temperature.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA observation data.

    Returns
    =======
    omega_corr : OMEGAdata
        The input OMEGA observation, where the reflectance is corrected from
        the thermal component.
    """
    # Test correction
    if omega.therm_corr:
        print("\033[1;33mThermal correction already applied\033[0m")
        return deepcopy(omega)
    # Initialisation
    omega2 = deepcopy(omega)
    omega_corr = deepcopy(omega)
    ny, nx, nlam = omega2.cube_i.shape
    # Correction spectres
    for x in tqdm(range(nx)):
        for y in tqdm(range(ny)):
            sp_rf_corr, surf_temp = corr_therm2_sp(omega2, x, y, disp=False)[1:]
            omega_corr.cube_rf[y,x] = sp_rf_corr
            omega_corr.surf_temp[y,x] = surf_temp
    # Sortie
    omega_corr.therm_corr = True
    omega_corr.therm_corr_infos['datetime'] = datetime.datetime.now()
    omega_corr.therm_corr_infos['method'] = '(M2) Simultaneous refl & temp'
    return omega_corr

##-----------------------------------------------------------------------------------
## Correction atmosphérique 1 - 1.93µm & 2.01µm au même niveau de réflectance
def corr_atm_sp(lam, sp_rf, tr_atm):
    """Remove the atmospheric component in an OMEGA spectrum.

    Parameters
    ==========
    lam : 1D array
        The wavelength array.
    sp_rf : 1D array
        The reflectance spectrum.
    tr_atm : 1D array
        Atmospheric transmission spectrum.

    Returns
    =======
    sp_rf_corr : 1D array
        The reflectance spectrum, corrected from the atmospheric component.
    """
    # TODO > retirer /0

    # Détermination exposant
    i_lam1, i_lam2 = uf.where_closer_array([1.93, 2.01], lam)
    expo = np.log(sp_rf[i_lam1] / sp_rf[i_lam2]) / np.log(tr_atm[i_lam1] / tr_atm[i_lam2])
    print(expo)
    # Correction
    sp_rf_corr = sp_rf * tr_atm**(-expo)
    # Sortie
    return sp_rf_corr
    
def corr_atm(omega):
    """Remove the atmospheric component in the OMEGA hyperspectral cube.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA observation data.

    Returns
    =======
    omega_corr : OMEGAdata
        The input OMEGA observation, where the reflectance is corrected from
        the atmospheric component.
    """
    # Test correction
    if omega.atm_corr:
        print("\033[1;33mAtmospheric correction already applied\033[0m")
        return deepcopy(omega)
    # Initialisation
    omega2 = deepcopy(omega)
    omega_corr = deepcopy(omega)
    ny, nx, nlam = omega2.cube_rf.shape
    lam = omega2.lam
    cube_rf = omega2.cube_rf
    ic_CL = np.concatenate([omega.ic['C'], omega.ic['L']])
    nV = len(omega.ic['V'])
    # Chargement données atmosphère
    atmorap = np.loadtxt(os.path.join(package_path, 'OMEGA_dataref', 'omega_atmorap_CL.dat'))
    tr_atm = np.ones(nlam)
    tr_atm[nV:] = atmorap[ic_CL]    # données atm uniquement voies C & L
    # Détermination exposant
    i_lam1, i_lam2 = uf.where_closer_array([1.93, 2.01], lam)
    expo = np.log(cube_rf[:,:,i_lam1] / cube_rf[:,:,i_lam2]) / np.log(tr_atm[i_lam1] / tr_atm[i_lam2])
    # Correction spectres
    for x in tqdm(range(nx), desc='Atmospheric correction'):
        for y in range(ny):
            sp_rf_corr = cube_rf[y,x] * tr_atm**(-expo[y,x])
            omega_corr.cube_rf[y,x] = sp_rf_corr
    # Sortie
    omega_corr.atm_corr = True
    omega_corr.atm_corr_infos['datetime'] = datetime.datetime.now()
    omega_corr.atm_corr_infos['method'] = 'M1 : same reflectance level at 1.93μm and 2.01μm'
    return omega_corr
    
##-----------------------------------------------------------------------------------
## Correction atmosphérique 2 - Spectre plat entre 1.97µm et 2.00µm
def f_min(expo, sp_rf, sp_atm):
    sp_rf_corr = sp_rf * sp_atm**(-expo)
    ecart = 0
    for i in range(len(sp_rf)):
        ecart += (sp_rf_corr[i] - np.average(sp_rf_corr))**2
    return ecart

def corr_atm2_sp(lam, sp_rf, tr_atm):
    """Remove the atmospheric component in an OMEGA spectrum.

    Parameters
    ==========
    lam : 1D array
        The wavelength array.
    sp_rf : 1D array
        The reflectance spectrum.
    tr_atm : 1D array
        Atmospheric transmission spectrum.

    Returns
    =======
    sp_rf_corr : 1D array
        The reflectance spectrum, corrected from the atmospheric component.
    """
    # Détermination exposant
    i_lams = uf.where_closer_array([1.97, 1.98, 2.00], lam)
    sp_rf2 = sp_rf[i_lams]
    sp_atm2 = tr_atm[i_lams]
    expo0 = 1
    res_opt = minimize(f_min, expo0, args=(sp_rf2, sp_atm2))
    # if res_opt.success:
    expo = res_opt.x[0]
    print(expo)
    # Correction
    sp_rf_corr = sp_rf * tr_atm**(-expo)
    # Sortie
    return sp_rf_corr
    
def corr_atm2(omega):
    """Remove the atmospheric component in the OMEGA hyperspectral cube.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA observation data.

    Returns
    =======
    omega_corr : OMEGAdata
        The input OMEGA observation, where the reflectance is corrected from
        the atmospheric component.
    """
    # Test correction
    if omega.atm_corr:
        print("\033[1;33mAtmospheric correction already applied\033[0m")
        return deepcopy(omega)
    # Initialisation
    omega2 = deepcopy(omega)
    omega_corr = deepcopy(omega)
    ny, nx, nlam = omega2.cube_rf.shape
    lam = omega2.lam
    cube_rf = omega2.cube_rf
    ic_CL = np.concatenate([omega.ic['C'], omega.ic['L']])
    nV = len(omega.ic['V'])
    # Chargement données atmosphère
    atmorap = np.loadtxt(os.path.join(package_path, 'OMEGA_dataref', 'omega_atmorap_CL.dat'))
    tr_atm = np.ones(nlam)
    tr_atm[nV:] = atmorap[ic_CL]    # données atm uniquement voies C & L
    # Détermination exposant
    i_lams = uf.where_closer_array([1.97, 1.98, 2.0], lam)
    cube_rf2 = cube_rf[:,:,i_lams]
    sp_atm2 = tr_atm[i_lams]
    expo0 = 1
    # Correction spectres
    for x in tqdm(range(nx)):
        for y in range(ny):
            expo = minimize(f_min, expo0, args=(cube_rf2[y,x], sp_atm2)).x[0]
            omega_corr.cube_rf[y,x] = cube_rf[y,x] * tr_atm**(-expo)
    # Sortie
    omega_corr.atm_corr = True
    omega_corr.atm_corr_infos['datetime'] = datetime.datetime.now()
    omega_corr.atm_corr_infos['method'] = 'M2 : flattest spectra between 1.97µm and 2.00µm'
    return omega_corr
    
##-----------------------------------------------------------------------------------
## **Voie L** -- Correction thermique 1 et atmosphérique 1 simultanées
def _corr_therm_atm_sp(args):
    """Remove the thermal and atmospheric component in an OMEGA spectrum, using the historical method
    based on the reflectance determination using reference spectra from Calvin & Erard (1997).
    
    Important note : this function is dedicated to internal use in the corr_therm() function,
        as it use non-local objects related to the multiprocessing.

    Parameters
    ==========
    args : tuple of 3 elements
      | x : int
      |     The x-coordinate of the pixel.
      | y : int
      |     The y-coordinate of the pixel.
      | disp : bool, optional (default True)
      |     If True display the fitted temperature/reflectance in the console.

    Returns
    =======
    sp_rf_corr : 1D array
        The reflectance spectrum, corrected from the thermal and atmospheric component.
    T_fit : float
        The retrieved surface temperature (in K).
    x : int
        The x-coordinate of the pixel.
    y : int
        The y-coordinate of the pixel.
    """
    # Extraction arguments
    x, y, disp = args
    # Extraction données
    global _omega_tmp
    omega = _omega_tmp
    lam = omega.lam
    sp_sol = omega.specmars
    sp_rf = omega.cube_rf[y, x]
    sp_i = omega.cube_i[y, x]
    ecl = np.cos(omega.inci[y, x] * np.pi/180)
    # spectels #97-#112 des spectres de ref <-> 2.3-2.5µm
    fref = os.path.join(package_path, 'OMEGA_dataref', 'refclair_sombr_omega_CL.dat') # from Erard and Calvin (1997)
    sp_clair, sp_sombre = np.loadtxt(fref, unpack=True)
    i_lam1, i_lam2 = uf.where_closer_array([2.3, 2.5], lam)
    alb_clair = np.average(sp_clair[97:113])
    alb_sombre = np.average(sp_sombre[97:113])
    alb_C = np.average(sp_rf[i_lam1:i_lam2+1])
    # Simulation d'un spectre : thèse D. Jouglet, p159
    # CL des spectres de référence `clair` et `sombre` (Calvin & Erard 1997)
    coeff = (alb_C - alb_sombre) / (alb_clair - alb_sombre)
    sp_simu = coeff * sp_clair + (1-coeff) * sp_sombre
    # Sélection des spectels 5.03-5.09µm (4 derniers voie L)
    sp_simu_5m = sp_simu[252:256]
    i_lam3, i_lam4 = uf.where_closer_array([5.03, 5.09], lam)
    i_lam4 += 1
    #--- Atmosphère
    nlam = len(lam)
    ic_CL = np.concatenate([omega.ic['C'], omega.ic['L']])
    nV = len(omega.ic['V'])
    # Chargement données atmosphère
    atmorap = np.loadtxt(os.path.join(package_path, 'OMEGA_dataref', 'omega_atmorap_CL.dat'))
    tr_atm = np.ones(nlam)
    tr_atm[nV:] = atmorap[ic_CL]    # données atm uniquement voies C & L
    # Détermination exposant
    i_lam1, i_lam2 = uf.where_closer_array([1.93, 2.01], lam)
    expo = np.log(sp_rf[i_lam1] / sp_rf[i_lam2]) / np.log(tr_atm[i_lam1] / tr_atm[i_lam2])
    #---
    def simu_sp_5microns(i_lams, T):
        i1, i2 = i_lams.astype(int)
        lam2 = lam[i1:i2] * 1e-6    # Conversion en m
        sp_sol2 = sp_sol[i1:i2]
        Blam = uf.planck(lam2, T) * 1e-6    # Loi de Planck en W.m-2.sr-1.µm-1
        sp_simu2 = sp_simu_5m * sp_sol2 * ecl + (1-sp_simu_5m) * Blam
        return sp_simu2
    try:
        # Fit de la température
        T_fit = curve_fit(simu_sp_5microns, (i_lam3, i_lam4), sp_i[i_lam3:i_lam4], 
                            bounds=(0,400))[0][0]
        # Réflectance
        refl = np.average(sp_simu[252:256])
        if disp:
            print('Temperature = {0:.3f} K   |   Reflectance = {1:.5f}'.format(T_fit, refl))
        # Correction thermique spectre
        Blam = uf.planck(lam*1e-6, T_fit) * 1e-6   # En W.m-2.sr-1.µm-1
        # sp_rf_corr = (sp_i - Blam) / (sp_sol*ecl - Blam)
        sp_rf_corr = (sp_i*tr_atm**(-expo) - Blam) / (sp_sol*ecl - Blam)
    except ValueError:
        # Si fit impossible (infs ou NaN dans le spectre) -> NaN partout
        sp_rf_corr = deepcopy(sp_rf)
        sp_rf_corr.fill(np.nan)
        T_fit = np.nan
    # Output
    return sp_rf_corr, T_fit, x, y

def corr_therm_atm(omega, npool=1):
    """Remove the thermal and atmospheric component in the OMEGA hyperspectral cube.
    
    Parallelization is implemented using the multiprocessing module. The number of
    process to run is controlled by the npool argument.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA observation data.
    npool : int, optional (default 1)
        Number of parallelized worker process to use.

    Returns
    =======
    omega_corr : OMEGAdata
        The input OMEGA observation, where the reflectance is corrected from
        the thermal and atmospheric component.
    """
    # Test correction
    if omega.therm_corr and omega.atm_corr:
        print("\033[1;33mThermal & atmosphecir corrections already applied\033[0m")
        return deepcopy(omega)
    # Initialisation
    global _omega_tmp
    _omega_tmp = deepcopy(omega)
    omega_corr = deepcopy(omega)
    ny, nx, nlam = omega.cube_i.shape
    rf_corr = np.zeros((ny,nx,nlam), dtype=np.float64)
    surf_temp = np.zeros((ny,nx), dtype=np.float64)
    # Itérateur
    it = [(x, y, False) for x, y in itertools.product(range(nx), range(ny))]
    # Correction thermique
    # chunksize = len(it) // npool    # Approx size of each process
    chunksize = 1
    # pool = mp.Pool(npool)
    with mp.Pool(npool) as pool:
        for res in tqdm(pool.imap_unordered(_corr_therm_atm_sp, it, chunksize), total=len(it), desc='Thermal correction'):
            sp_rf_corr, T_fit, x, y = res
            rf_corr[y,x] = sp_rf_corr
            surf_temp[y,x] = T_fit
        pool.close()
    _omega_tmp = None
    omega_corr.cube_rf = rf_corr
    omega_corr.surf_temp = surf_temp
    # Update infos
    omega_corr.therm_corr = True
    omega_corr.therm_corr_infos['datetime'] = datetime.datetime.now()
    omega_corr.therm_corr_infos['method'] = '(M1) Calvin & Erard - Simultaneous atm (L channel)'
    omega_corr.atm_corr = True
    omega_corr.atm_corr_infos['datetime'] = datetime.datetime.now()
    omega_corr.atm_corr_infos['method'] = 'M1 : same reflectance level at 1.93μm and 2.01μm - Simultaneous therm (L channel)'
    # Sortie
    # tfin = time.time()
    # print('Duration : {0:.0f} min {1:.2f} sec'.format((tfin-tini)//60, (tfin-tini)%60))
    return omega_corr
    
##-----------------------------------------------------------------------------------
## Correction mode 128
def corr_mode_128(omega):
    """Correction corrupted pixels mode 128.
    """
    omega_corr = deepcopy(omega)
    ic128 = deepcopy(omega.ic)
    npixel = omega.npixel
    nscan = omega.nscan
    if npixel != 128:
        print('\033[1mNot a 128 pixel cube\033[0m')
    elif (npixel==128) & (omega.orbit >= 513):
        print('\033[33mCorrupted 128 pixel cube\033[0m')
        omega128_interp = readsav(os.path.join(package_path, 'OMEGA_dataref', 'omega128_interpol.sav'))
        if str.encode(omega.name[3:]) in omega128_interp['cublist']:
            i_omega = np.where(omega128_interp['cublist'] == str.encode(omega.name[3:]))[0][0]
            cubtype = omega128_interp['cubstatus'][i_omega]
            if cubtype <= 0:
                print('\033[01;33;41mNo correction available (good, corrupted or unusual cube)\033[0m')
            else:
                if cubtype == 1:
                    print('Parity 1 : spectel 28 corrupted in even lines')
                    firsteven, firstodd = 28, 12
                elif cubtype == 2:
                    print('Parity 2 : spectel 28 corrupted in odd lines')
                    firsteven, firstodd = 12, 28
                even = 2 * (np.arange((nscan-2)//2, dtype=int) + 1)     # even lines
                odd  = 2 * np.arange((nscan-1)//2, dtype=int) + 1       # odd lines
                cube_i = omega.cube_i
                cube_rf = omega.cube_rf
                cube_i_corr = deepcopy(cube_i)
                cube_rf_corr = deepcopy(cube_rf)
                for w in range(11): # loop on spectral position
                    cube_rf_corr[even, 80:95, firsteven+32*w:firsteven+3+32*w] = 0.5 * (
                        cube_rf[even+1, 80:95, firsteven+32*w:firsteven+3+32*w] + 
                        cube_rf[even-1, 80:95, firsteven+32*w:firsteven+3+32*w] )
                    cube_rf_corr[odd, 80:95, firstodd+32*w:firstodd+3+32*w] = 0.5 * (
                        cube_rf[odd+1, 80:95, firstodd+32*w:firstodd+3+32*w] + 
                        cube_rf[odd-1, 80:95, firstodd+32*w:firstodd+3+32*w] )
                    cube_rf_corr[0, 80:95, firsteven+32*w:firsteven+3+32*w] = (
                        cube_rf[1, 80:95, firsteven+32*w:firsteven+3+32*w] )
                    if (nscan/2)*2 == nscan:
                        cube_rf_corr[nscan-1, 80:95, firstodd+32*w:firstodd+3+32*w] = (
                            cube_rf[nscan-2, 80:95, firstodd+32*w:firstodd+3+32*w])
                    else:
                        cube_rf_corr[nscan-1, 80:95, firsteven+32*w:firsteven+3+32*w] = (
                            cube_rf[nscan-2, 80:95, firsteven+32*w:firsteven+3+32*w] )
                omega_corr.cube_i = cube_i_corr
                omega_corr.cube_rf = cube_rf_corr
        else:
            print('\033[01;33;41mCube not in list\033[0m')
        # lam_mask['mode128'] = 
        # omega_corr.add_infos += '\nWarning: Corrupted 128 pixel cube'
    return omega_corr
    
##-----------------------------------------------------------------------------------
## Correction & sauvegarde
def corr_save_omega(obsname, folder='auto', base_folder='_omega_py_path', security=True,
                    overwrite=True, compress=True, npool=1):
    """Correction and saving of OMEGA/MEx observations.
    
    Parallelization is implemented using the multiprocessing module. The number of
    process to run is controlled by the npool argument.

    Parameters
    ==========
    obsname : str
        The name of the OMEGA observation.
    folder : str, optional (default 'auto')
        The subfolder to save the data.
        | If 'auto' -> folder = 'vX', where X is the major release version of the used code.
    base_folder : str, optional (default _omega_py_path)
        The base folder path.
    security : bool, optional (default True)
        Enable / disable checking before overwriting a file.
        | True -> Check if the target file already exists before overwriting on it.
                  And if is the case, you will be asked for a confirmation.
        | False -> Didn't care about the already existing files.
    overwrite : bool, optional (default True)
        If security is False, default choice for overwriting on existent file.
    compress : bool, optional (default True)
        If True, the radiance cube after correction is removed (i.e. set to None)
        in order to reduce the size of the saved file.
    npool : int, optional (default 1)
        Number of parallelized worker process to use.
    """
    if folder == 'auto':
        folder = 'v' + str(int(_Version))
    omega = OMEGAdata(obsname)
    name = omega.name
    # Default path
    if base_folder == "_omega_py_path":
        base_folder = _omega_py_path
    # path synthax
    basename = os.path.join(base_folder, folder, name, '{0}.pkl')
    # Testing existent file
    if os.path.exists(basename.format('_corr_therm_atm')):
        exists = True
    else:
        exists = False
    if security:
        overwrite = uf.test_security_overwrite(basename.format('*'))
    if (not exists) or (exists and overwrite):
        save_omega(omega, folder=folder, base_folder=base_folder)
        print('\n\033[01mThermal correction\033[0m')
        omega_corr = corr_therm(omega, npool)
        if compress:
            omega_corr.cube_i = None
        save_omega(omega_corr, folder=folder, base_folder=base_folder, suff='corr_therm')
        print('\n\033[01mAtmospheric correction\033[0m')
        omega_corr_atm = corr_atm(omega_corr)
        save_omega(omega_corr_atm, folder=folder, base_folder=base_folder, suff='corr_therm_atm')
    else:
        print('\n\033[01;34mExistent files preserved for {0} - v{1}\033[0m\n'.format(name, _Version))

def corr_save_omega_list(liste_obs, folder='auto', base_folder='_omega_py_path',
                         security=True, overwrite=True, compress=True, npool=1):
    """Correction and saving of a list of OMEGA/MEx observations.
    
    Parallelization is implemented using the multiprocessing module. The number of
    process to run is controlled by the npool argument.

    Parameters
    ==========
    liste_obs : list of str
        The list of the name of the OMEGA observations.
    folder : str, optional (default 'auto')
        The subfolder to save the data.
        | If 'auto' -> folder = 'vX', where X is the major release version of the used code.
    base_folder : str, optional (default _omega_py_path)
        The base folder path.
    security : bool, optional (default True)
        Enable / disable checking before overwriting a file.
        | True -> Check if the target file already exists before overwriting on it.
                  And if is the case, you will be asked for a confirmation.
        | False -> Do not care about the already existing files.
    overwrite : bool, optional (default True)
        If security is False, default choice for overwriting on existent file.
    compress : bool, optional (default True)
        If True, the radiance cube after correction is removed (i.e. set to None)
        in order to reduce the size of the saved file.
    npool : int, optional (default 1)
        Number of parallelized worker process to use.
    """
    # Default path
    if base_folder == "_omega_py_path":
        base_folder = _omega_py_path
    N = len(liste_obs)
    if folder == 'auto':
        folder = 'v' + str(int(_Version))
    for i, obsname in enumerate(liste_obs):
        print('\n\033[01mComputing observation {0} / {1} : {2}\033[0m\n'.format(i+1, N, obsname))
        corr_save_omega(obsname, folder, base_folder, security, overwrite, compress, npool)
    print("\n\033[01;32m Done\033[0m\n")

##-----------------------------------------------------------------------------------
## Correction & sauvegarde - v2 (L)
def corr_save_omega2(obsname, folder='auto', base_folder='_omega_py_path', security=True,
                    overwrite=True, compress=True, npool=1):
    """Correction and saving of OMEGA/MEx observations.
    
    Parallelization is implemented using the multiprocessing module. The number of
    process to run is controlled by the npool argument.

    Parameters
    ==========
    obsname : str
        The name of the OMEGA observation.
    folder : str, optional (default 'auto')
        The subfolder to save the data.
        | If 'auto' -> folder = 'vX', where X is the major release version of the used code.
    base_folder : str, optional (default _omega_py_path)
        The base folder path.
    security : bool, optional (default True)
        Enable / disable checking before overwriting a file.
        | True -> Check if the target file already exists before overwriting on it.
                  And if is the case, you will be asked for a confirmation.
        | False -> Didn't care about the already existing files.
    overwrite : bool, optional (default True)
        If security is False, default choice for overwriting on existent file.
    compress : bool, optional (default True)
        If True, the radiance cube after correction is removed (i.e. set to None)
        in order to reduce the size of the saved file.
    npool : int, optional (default 1)
        Number of parallelized worker process to use.
    """
    if folder == 'auto':
        folder = 'v' + str(int(_Version))
    omega = OMEGAdata(obsname)
    name = omega.name
    # Default path
    if base_folder == "_omega_py_path":
        base_folder = _omega_py_path
    # path synthax
    basename = os.path.join(base_folder, folder, name, '{0}.pkl')
    # Testing existent file
    if os.path.exists(basename.format('_corr_therm_atm')):
        exists = True
    else:
        exists = False
    if security:
        overwrite = uf.test_security_overwrite(basename.format('*'))
    if (not exists) or (exists and overwrite):
        # save_omega(omega, folder=folder, base_folder=base_folder)
        print('\n\033[01mThermal & atmospheric corrections\033[0m')
        omega_corr = corr_therm_atm(omega, npool)
        if compress:
            omega_corr.cube_i = None
        save_omega(omega_corr, folder=folder, base_folder=base_folder, suff='corr_therm_atm')
    else:
        print('\n\033[01;34mExistent files preserved for {0} - v{1}\033[0m\n'.format(name, _Version))

def corr_save_omega2_list(liste_obs, folder='auto', base_folder='_omega_py_path',
                         security=True, overwrite=True, compress=True, npool=1):
    """Correction and saving of a list of OMEGA/MEx observations.
    
    Parallelization is implemented using the multiprocessing module. The number of
    process to run is controlled by the npool argument.

    Parameters
    ==========
    liste_obs : list of str
        The list of the name of the OMEGA observations.
    folder : str, optional (default 'auto')
        The subfolder to save the data.
        | If 'auto' -> folder = 'vX', where X is the major release version of the used code.
    base_folder : str, optional (default _omega_py_path)
        The base folder path.
    security : bool, optional (default True)
        Enable / disable checking before overwriting a file.
        | True -> Check if the target file already exists before overwriting on it.
                  And if is the case, you will be asked for a confirmation.
        | False -> Do not care about the already existing files.
    overwrite : bool, optional (default True)
        If security is False, default choice for overwriting on existent file.
    compress : bool, optional (default True)
        If True, the radiance cube after correction is removed (i.e. set to None)
        in order to reduce the size of the saved file.
    npool : int, optional (default 1)
        Number of parallelized worker process to use.
    """
    # Default path
    if base_folder == "_omega_py_path":
        base_folder = _omega_py_path
    N = len(liste_obs)
    if folder == 'auto':
        folder = 'v' + str(int(_Version))
    for i, obsname in enumerate(liste_obs):
        print('\n\033[01mComputing observation {0} / {1} : {2}\033[0m\n'.format(i+1, N, obsname))
        corr_save_omega(obsname, folder, base_folder, security, overwrite, compress, npool)
    print("\n\033[01;32m Done\033[0m\n")

##-----------------------------------------------------------------------------------
## Liste observations depuis fichier CSV
def import_list_obs_csv(filename):
    """Import a list of observations ID from a csv file generated by JMars.

    Parameters
    ==========
    filename : str
        The target path of the csv file.

    Returns
    =======
    liste_obs : array of str
        The list of observations ID from the csv file.
    """
    df = pd.read_csv(filename)
    columns_id = df.columns
    liste_obs = np.array(df['product_id_trunc'])
    return liste_obs

##-----------------------------------------------------------------------------------
## Setters
def set_omega_bin_path(new_path):
    """Set the global private _omega_bin_path variable to new_path.

    Parameters
    ==========
    new_path : str
        The new path of the OMEGA binary files (.QUB and .NAV).
    """
    if not isinstance(new_path, str):
        raise ValueError('new_path must be a str')
    global _omega_bin_path
    _omega_bin_path = new_path

def set_omega_py_path(new_path):
    """Set the global private _omega_py_path variable to new_path.

    Parameters
    ==========
    new_path : str
        The new path of the OMEGA python-made files.
    """
    if not isinstance(new_path, str):
        raise ValueError('new_path must be a str')
    global _omega_py_path
    _omega_py_path = new_path

##-----------------------------------------------------------------------------------
## Getters
def get_omega_bin_path():
    """Return the vavue of the global private _omega_bin_path variable.

    Returns
    =======
    omega_bin_path : str
        The path of the OMEGA binary files (.QUB and .NAV).
    """
    return deepcopy(_omega_bin_path)

def get_omega_py_path():
    """Return the vavue of the global private _omega_py_path variable.

    Returns
    =======
    omega_py_path : str
        The new path of the OMEGA python-made files.
    """
    return deepcopy(_omega_py_path)

def get_names(omega_list):
    """Return the array of the observation ID of each OMEGA/MEx observation in omega_list.

    Parameters
    ==========
    omega_list : array of OMEGAdata
        The input array of OMEGA observations.
    
    Returns
    =======
    names : ndarray
        The array of the omega_list observations ID.
    """
    names = []
    for omega in omega_list:
        names.append(omega.name)
    return names

def get_ls(omega_list):
    """Return the array of the Solar longitude of each OMEGA/MEx observation in omega_list.

    Parameters
    ==========
    omega_list : array of OMEGAdata
        The input array of OMEGA observations.
    
    Returns
    =======
    ls : ndarray
        The array of the omega_list Ls.
    """
    ls = []
    for omega in omega_list:
        ls.append(omega.ls)
    return ls

##-----------------------------------------------------------------------------------
## Update cube quality
def update_cube_quality(obs_name='ORB*.pkl', folder='auto', version=_Version, 
                        base_folder='_omega_py_path'):
    """Update the quality attribute of previously saved OMEGAdata objects.

    Parameters
    ==========
    obs_name : str, optional (default 'ORB*.pkl')
        The files basename.
    folder : str, optional (default 'auto')
        The subfolder where the data is.
        | If 'auto' -> folder = 'vX', where X is the major release version of the used code.
    version : float, optional (default _Version)
        The version of the target file (if folder is 'auto').
        Default is the current code version.
    base_folder : str, optional (default _omega_py_path)
        The base folder path.
    """
    # Default path
    if base_folder == "_omega_py_path":
        base_folder = _omega_py_path
    # Initialisation
    if obs_name[-4] != '.pkl':
        obs_name += '.pkl'
    if folder == 'auto':
        folder = 'v' + str(int(version))
    basename = uf.myglob(os.path.join(base_folder, folder, obs_name))
    # Load list corrupted obs
    OBC = readsav(os.path.join(package_path, 'OMEGA_dataref', 'OBC_OMEGA_OCT2017.sav'))
    good_orbits_OBC = np.array(OBC['good_orbits'][0], dtype=int)
    corrupted_orbits_csv = pd.read_csv(os.path.join(package_path, 'OMEGA_dataref', 'corrupted_obs.csv'), 
                                       comment='#', skipinitialspace=True)
    corrupted_orbits = np.array(corrupted_orbits_csv['corrupted_obs'], dtype=str)
    corrupted_orbits_comments = np.array(corrupted_orbits_csv['comment'], dtype=str)
    # Loop on obs in the selected folder
    fnames = glob.glob(basename)
    if fnames == []:
        print("\033[1;33mNo such file found.\033[0m")
    else:
        for fname in tqdm(fnames):
            omega = load_omega(fname, disp=False)
            omega.quality = 1
            if (omega.npixel==128) & (omega.orbit >= 513):
                omega.quality = 128
                omega.add_infos = 'Corrupted 128 pixels cube'
            if omega.orbit not in good_orbits_OBC:
                omega.quality = 0
                omega.add_infos = 'Corrupted orbit'
            if omega.name in corrupted_orbits:
                omega.quality = 0
                i_obs = int(np.where(corrupted_orbits==omega.name)[0])
                omega.add_infos = corrupted_orbits_comments[i_obs]
            save_omega(omega, fname, '', '', '', '', False)
        print('\033[1m{0} files updated\033[0m'.format(len(fnames)))

##-----------------------------------------------------------------------------------
## Importation liste OMEGAdata avec filtrage automatisé
def load_omega_list2(liste_obs, therm_corr=True, atm_corr=True, **kwargs):
    """Load a list of saved OMEGAdata objects, using load_omega().

    Parameters
    ==========
    liste_obs : array of str
        List of OMEGA/MEx observations ID.
    therm_corr : bool or None, optional (default None)
        | True -> Only results with thermal correction.
        | False -> Only results without thermal correction.
        | None -> Both with and without thermal correction.
    atm_corr : bool or None, optional (default None)
        | True -> Only results with atmospheric correction.
        | False -> Only results without atmospheric correction.
        | None -> Both with and without atmospheric correction.
    **kwargs:
        Optional arguments for autoload_omega().

    Returns
    =======
    omega_list : list of OMEGAdata objects
        The list of loaded objects of OMEGA/MEx observation.
    """
    omega_list = []
    Nabs = 0
    OBC = readsav(os.path.join(package_path, 'OMEGA_dataref', 'OBC_OMEGA_OCT2017.sav'))
    good_orbits_OBC = np.array(OBC['good_orbits'][0], dtype=int)
    for i, obsname in enumerate(tqdm(liste_obs)):
        omega = autoload_omega(obsname, therm_corr=therm_corr, atm_corr=atm_corr, disp=False,
                               **kwargs)
        if omega is None:
            Nabs += 1
            continue
        if not omega.orbit in good_orbits_OBC:
            continue
        if omega.quality == 0:
            continue
        if omega.target != 'MARS':
            continue
        if omega.mode_channel != 1:
            continue
        if omega.data_quality == 0:
            continue
        if omega.point_mode == 'N/A':
            continue
        if (int(omega.name[-1]) == 0) and (omega.npixel == 64) and (omega.bits_per_data == 1):
            continue
        # if omega.npixel == 16:
            # continue
        omega_list.append(omega)
    Ntot = len(liste_obs)
    Nacc = len(omega_list)
    Nrej = Ntot - Nacc - Nabs
    print('\n\033[1m{0} observations in list_obs\n'.format(Ntot) +
          '{0} loaded, {1} rejected, {2} not found\033[0m\n'.format(Nacc, Nrej, Nabs))
    return omega_list

def test_cube(obs):
    """Test the quality of an OMEGA/MEx observation from the header informations
    witout open it.

    Parameters
    ==========
    obs : str
        The name of the OMEGA observation.

    Returns
    =======
    test_quality : bool
        | True -> Accepted observation.
        | False -> Rejected observation.
    """
    # Recherhe nom de fichier
    data_path = _omega_bin_path
    obs_name = uf.myglob(os.path.join(data_path, '*' + obs + '*.QUB'))
    if obs_name is None:
        print("\033[1;33mAborted\033[0m")
        return False
    nomfic0 = os.path.split(obs_name)[1][:-4]    # Récupération nom + décodage UTF-8
    numCube = int(nomfic0[-1])
    # Lecture header fichier .QUB
    hd_qub = _read_header(obs_name[:-4] + '.QUB')
    summation = np.int64(hd_qub['DOWNTRACK_SUMMING'])
    bits_per_data = np.float64(hd_qub['INST_CMPRS_RATE'])
    data_quality = np.int64(hd_qub['DATA_QUALITY_ID'])
    mode_channel_tmp = hd_qub['COMMAND_DESC'][34:36]
    if mode_channel_tmp == 'EF':
        mode_channel = 1
    elif mode_channel_tmp == '80':
        mode_channel = 2
    elif mode_channel_tmp == 'C7':
        mode_channel = 3
    else:
        mode_channel = mode_channel_tmp
    # Lecture header fichier .NAV
    if glob.glob(obs_name[:-4] + '.NAV') == []:
        return False    # Pas de fichier .NAV
    hd_nav = _read_header(obs_name[:-4] + '.NAV')
    npixel, npara, nscan = np.array(hd_nav['CORE_ITEMS'][1:-1].split(','), dtype=np.int64)
    point_mode = hd_nav['SPACECRAFT_POINTING_MODE'][1:-1]
    target = hd_nav['TARGET_NAME']
    # Test si cube OK
    if target != 'MARS':
        return False
    elif mode_channel != 1:
        return False
    elif data_quality == 0:
        return False
    elif point_mode == 'N/A':
        return False
    elif (numCube == 0) and (npixel == 64) and (bits_per_data == 1):
        return False
    else:
        return True

##-----------------------------------------------------------------------------------
## List available good observations & generate csv file
def compute_list_good_observations(savfilename='liste_good_obs.csv', 
                                   folder='../data/OMEGA/liste_obs', security=True):
    """Scan the available OMEGA/MEx data cubes and list the observations considered as 
    good quality.
    The results are saved in the specified csv file.

    Parameters
    ==========
    savfilename : str
        The name of the csv file to save the data.
    folder : str
        The name of the folder where the saved file will be located.
        Final saved file path = folder + savfilename
    security : bool, optional (default True)
        Enable / disable checking before overwriting a file.
        | True -> Check if the target file already exists before overwriting on it.
                  And if is the case, you will be asked for a confirmation.
        | False -> Didn't care about the already existing files.
    """
    # Test existence fichier de sauvegarde
    sav_file_path = os.path.join(folder, savfilename)
    if security:
        test_overwrite = uf.test_security_overwrite(sav_file_path)
        if not test_overwrite:
            return None
    # Liste observations disponibles
    bin_obs_list = glob.glob(os.path.join(_omega_bin_path, 'ORB*.QUB'))
    bin_obs_list.sort()
    # Initialisation
    gobs = open(sav_file_path, 'w', encoding='utf-8')
    gobs.write('obsname, Ls [°], lat_min [°], lat_max [°], lon_min [°], lon_max [°], '
               + 'UTC date/time, Npixel, Nscan\n')
    Nacc = 0
    # Test qualité de chaque observation
    for obs_name in tqdm(bin_obs_list):
        nomfic0 = os.path.split(obs_name)[1][:-4]    # Récupération nom + décodage UTF-8
        numCube = nomfic0[-1]
        # Lecture header fichier .QUB
        hd_qub = _read_header(obs_name[:-4] + '.QUB')
        summation = np.int64(hd_qub['DOWNTRACK_SUMMING'])
        bits_per_data = np.float64(hd_qub['INST_CMPRS_RATE'])
        data_quality = np.int64(hd_qub['DATA_QUALITY_ID'])
        mode_channel_tmp = hd_qub['COMMAND_DESC'][34:36]
        if mode_channel_tmp == 'EF':
            mode_channel = 1
        elif mode_channel_tmp == '80':
            mode_channel = 2
        elif mode_channel_tmp == 'C7':
            mode_channel = 3
        else:
            mode_channel = mode_channel_tmp
        # Lecture header fichier .NAV
        if glob.glob(obs_name[:-4] + '.NAV') == []:
            continue
        hd_nav = _read_header(obs_name[:-4] + '.NAV')
        npixel, npara, nscan = np.array(hd_nav['CORE_ITEMS'][1:-1].split(','), dtype=np.int64)
        point_mode = hd_nav['SPACECRAFT_POINTING_MODE'][1:-1]
        target = hd_nav['TARGET_NAME']
        test = True
        # Test si cube OK
        if target != 'MARS':
            continue
        elif mode_channel != 1:
            continue
        elif data_quality == 0:
            continue
        elif point_mode == 'N/A':
            continue
        elif (numCube == '0') and (npixel == 64) and (bits_per_data == 1):
            continue
        else:
            # Si OK -> sauvegarde des infos dans le fichier
            gobs.write(('{obsname:s}, {ls:s}, {lat_min:s}, {lat_max:s}, {lon_min:s},'
                    +'{lon_max:s}, {utc:s}, {npixel:d}, {nscan:d}\n').format(obsname = nomfic0, 
                                    ls = hd_nav['SOLAR_LONGITUDE'],
                                    lat_min = hd_nav['MINIMUM_LATITUDE'],
                                    lat_max = hd_nav['MAXIMUM_LATITUDE'],
                                    lon_min = hd_nav['WESTERNMOST_LONGITUDE'],
                                    lon_max = hd_nav['EASTERNMOST_LONGITUDE'],
                                    utc = hd_nav['START_TIME'][:16],
                                    npixel = npixel,
                                    nscan = nscan))
            Nacc += 1
    gobs.close()
    # Résultats
    Ntot = len(bin_obs_list)
    Nrej = Ntot - Nacc
    print('\n\033[1m{0} observations found\n'.format(Ntot) +
            '{0} accepted, {1} rejected\033[0m'.format(Nacc, Nrej))
    print('\n\033[01;34mResults saved in \033[0;03m' + sav_file_path + '\033[0m')

##-----------------------------------------------------------------------------------
## UTC to MY
def utc_to_my(dt):
    """Convert a UTC datetime to the corresponding Martian Year (MY).
    
    Martian Years are numbered according to the calendar proposed by R. Todd Clancy 
    (Clancy et al., Journal of Geophys. Res 105, p 9553, 2000): Martian Year 1 begins 
    (at a time such that Ls=0) on April 11th, 1955.

    Parameters
    ==========
    dt : datetime.datetime
        The UTC datetime object.

    Returns
    =======
    my : int
        The corresponding Martian Year.
    """
    datetime_my1 = datetime.datetime(1955, 4, 11)   # Start MY 1
    my_sol_duration = 668.6     # Nb of Martian sols during a MY
    sol_sec_duration = 88775.245    # Duration of a sol in seconds
    my = int( (dt - datetime_my1).total_seconds() // (my_sol_duration * sol_sec_duration)) + 1
    return my

##----------------------------------------------------------------------------------------
## Shared wavelength array
def shared_lam(lam_list):
    """Return a list of wavelength shared by all the input wavelength arrays.

    Parameters
    ==========
    lam_list : list of 1D np.array
        The list of wavelength array.

    Returns
    =======
    lam2 : 1D np.array
        The wavelength array that contains only wavelength shared by all the arrays of
        lam_list.
    """
    lam0 = deepcopy(lam_list[0])
    lam2 = []
    for lami in lam0:
        test = True
        for lam_array in lam_list:
            if not (lami in lam_array):
                test = False
                break
        if test:
            lam2.append(lami)
    lam2 = np.array(lam2)
    return lam2

def shared_lam_omegalist(omega_list):
    """Return a list of wavelength shared by all the wavelength arrays of the input
    OMEGA/MEx observations.

    Parameters
    ==========
    omega_list : list of OMEGAdata
        The list of OMEGA/MEx observations.

    Returns
    =======
    lam2 : 1D np.array
        The wavelength array that contains only wavelength shared by all the arrays of
        lam_list.
    """
    lam0 = deepcopy(omega_list[0].lam)
    lam2 = []
    for lami in lam0:
        test = True
        for omega in omega_list:
            if not (lami in omega.lam):
                test = False
                break
        if test:
            lam2.append(lami)
    lam2 = np.array(lam2)
    return lam2

##-----------------------------------------------------------------------------------
## End of code
##-----------------------------------------------------------------------------------
