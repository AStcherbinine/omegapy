#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## omega_data.py
## Created by Aurélien STCHERBININE
## Last modified by Aurélien STCHERBININE : 09/03/2020

##----------------------------------------------------------------------------------------
"""Importation of OMEGA observations in the OMEGAdata class.
Using IDL routines containing in ./omega_routines/.
"""
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
## Packages
# Global
import numpy as np
from copy import deepcopy
from tqdm import tqdm
import pidly
import scipy.constants as const
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.io import readsav
from datetime import datetime
import pickle
import os
import glob
import pandas as pd
# Local
from . import useful_functions as uf

# Name of the current file
py_file = 'omega_data.py'
Version = 1.2

##----------------------------------------------------------------------------------------
## Class OMEGAdata - Importation of OMEGA data cubes
class OMEGAdata:
    """Importation of OMEGA/MEx observation, using the readomega_vpy.pro IDL routine.

    Parameters
    ==========
    obs : str
        The name of the OMEGA observation.
    empty : bool, optional (default False)
        If True, return an empty OMEGAdata object.
    data_path : str, optional (default '/data2/opt/geomeg/data/product/')
        The path of the directory containing the data (.QUB) and 
        navigation (.NAV) files.

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
    emer : 2D array
        The angle of emergent line (from the surface) (deg).
    inci : 2D array
        The incidence angle at the surface (deg).
    specmars : 1D array
        The Solar radiation spectrum on Mars (W.m-2.sr-1.µm-1).
    utc : datetime.datetime
        The average UCT time of the observation.
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
    version : float
        The version of the omega_data.py file used.
    attributes : list
        The list of the attributes of the object.
    **TO BE COMPLETED**
    """

    def __init__(self, obs='', empty=False, data_path='/data2/opt/geomeg/data/product/'):
        self.attributes = ['name', 'lam', 'cube_i', 'cube_rf', 'ls', 'lat', 'lon', 'alt',
                           'emer', 'inci', 'specmars', 'utc', 'quality', 'therm_corr',
                           'atm_corr', 'therm_corr_infos', 'atm_corr_infos', 
                           'version', 'attributes']
        # Infos
        self.version = Version
        self.therm_corr = False
        self.atm_corr = False
        self.therm_corr_infos = {'datetime': None, 'method': None}
        self.atm_corr_infos = {'datetime': None, 'method': None}
        self.quality = 1
        self.add_infos = ''

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
            self.emer = np.array([[]])
            self.inci = np.array([[]])
            self.specmars = np.array([])
            self.utc = datetime.now()
            self.orbit = None
            self.surf_temp = np.array([[]])
            self.ic = {'V' : np.arange(265, 333),
                       'C' : np.arange(8, 123),
                       'L' : np.arange(137, 256)}
            self.lon_grid = np.array([[]])
            self.lat_grid = np.array([[]])
            self.surf_temp = np.array([[]])
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

        else:
            idl = pidly.IDL()
            idl("cd, 'omega_routines'")
            obs_name = uf.myglob(data_path + '*' + obs + '*.QUB')
            if obs_name is None:
                print("\033[1;33mAborted\033[0m")
                return None
            nomfic0 = obs_name[obs_name.rfind('/')+1:-4]    # Récupération nom + décodage UTF-8
            idl("readomega_vpy, '{0}', ldat, jdat, wvl, ic, specmars, ".format(nomfic0)
                + "latitude, longitude, emergence, incidence, altitude, solarlongi, ut_time, geocube")

            print("\n\033[01;34mComputing data extraction and correction...\033[0m", end=' ')
            # Extract values from the idl session
            data_dict = idl.ev_list(['ldat', 'jdat', 'wvl', 'ic', 'specmars'], use_cache=True)
            ldat = data_dict['ldat']
            jdat = data_dict['jdat']
            wvl = data_dict['wvl']
            ic = data_dict['ic']
            specmars = data_dict['specmars']
            infos_dict = idl.ev_list(['latitude', 'longitude', 'emergence', 'incidence',
                                    'altitude', 'solarlongi', 'ut_time', 'geocube'], use_cache=True)
            lat = infos_dict['latitude']
            lon = infos_dict['longitude']
            alt = infos_dict['altitude']
            emer = infos_dict['emergence']
            inci = infos_dict['incidence']
            ls = infos_dict['solarlongi']
            utc = infos_dict['ut_time']
            geocube = infos_dict['geocube']
            # Close the idl session
            idl.close()
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
            cube_i2 = np.moveaxis(cube_i, [0,1,2], [0,2,1])
            cube_rf2 = np.moveaxis(cube_rf, [0,1,2], [0,2,1])
            # Observation UTC date & time
            Y, M, D, h, m, s = np.average(utc[:,:6], axis=0).astype(np.int64)
            utc_dt = datetime(Y, M, D, h, m, s)
            # Longitude pixels grid
            ny, nx = lon.shape
            lon_px = np.moveaxis(geocube[:,13:17,:], [0,1,2], [0,2,1]) * 1e-4
            lon_grid = np.zeros((ny+1, nx+1))
            lon_grid[1:,1:] = lon_px[:,:,0]
            lon_grid[1:,0] = lon_px[:,0,1]
            lon_grid[0,1:] = lon_px[0,:,3]
            lon_grid[0,0] = lon_px[0,0,2]
            # Longitude pixels grid
            lat_px = np.moveaxis(geocube[:,17:21,:], [0,1,2], [0,2,1]) * 1e-4
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
            self.ls = ls.astype(np.float64)
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
            #--------------------------
            # Data from the .NAV file
            #--------------------------
            nav = read_header(obs_name[:-4] + '.NAV')
            npixel, npara, nscan = np.array(nav['CORE_ITEMS'][1:-1].split(','), dtype=np.int64)
            self.lrec = np.int64(nav['RECORD_BYTES'])
            self.nrec = np.int64(nav['LABEL_RECORDS'])
            self.sol_dist_au = np.float64(nav['SOLAR_DISTANCE']) / 14960e4
            npixel, npara, nscan = np.array(nav['CORE_ITEMS'][1:-1].split(','), dtype=np.int64)
            self.npixel = npixel
            self.npara = npara
            self.nscan = nscan
            self.point_mode = nav['SPACECRAFT_POINTING_MODE'][1:-1]
            self.orient = np.array(nav['SPACECRAFT_ORIENTATION'][1:-1].split(','), dtype=np.int64)
            # self.ls = np.float64(nav['SOLAR_LONGITUDE'])
            self.subsol_lon = np.float64(nav['SUB_SOLAR_LONGITUDE'])
            self.subsol_lat = np.float64(nav['SUB_SOLAR_LATITUDE'])
            self.min_lat = np.float64(nav['MINIMUM_LATITUDE'])
            self.max_lat = np.float64(nav['MAXIMUM_LATITUDE'])
            self.min_lon = np.float64(nav['WESTERNMOST_LONGITUDE'])
            self.max_lon = np.float64(nav['EASTERNMOST_LONGITUDE'])
            self.slant = np.float64(nav['SLANT_DISTANCE'])
            #--------------------------
            temp_init = np.zeros(self.lat.shape)
            temp_init[:] = np.nan
            self.surf_temp = temp_init
            #--------------------------
            # Cube quality
            OBC = readsav('../data/OMEGA_dataref/OBC_OMEGA_OCT2017.sav')
            good_orbits_OBC = np.array(OBC['good_orbits'][0], dtype=int)
            corrupted_orbits_csv = pd.read_csv('../data/OMEGA_dataref/corrupted_obs.csv', comment='#',
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
        new_omega.emer = self.emer
        new_omega.inci = self.inci
        new_omega.specmars = self.specmars
        new_omega.utc = self.utc
        new_omega.orbit = self.orbit
        new_omega.ic = self.ic
        new_omega.lam_ma = self.lam_ma
        new_omega.lon_grid = self.lon_grid
        new_omega.lat_grid = self.lat_grid
        new_omega.surf_temp = self.surf_temp
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
        new_omega.emer = deepcopy(self.emer, memo)
        new_omega.inci = deepcopy(self.inci, memo)
        new_omega.specmars = deepcopy(self.specmars, memo)
        new_omega.utc = deepcopy(self.utc, memo)
        new_omega.orbit = deepcopy(self.orbit, memo)
        new_omega.ic = deepcopy(self.ic, memo)
        new_omega.lam_ma = deepcopy(self.lam_ma, memo)
        new_omega.lon_grid = deepcopy(self.lon_grid, memo)
        new_omega.lat_grid = deepcopy(self.lat_grid, memo)
        new_omega.surf_temp = deepcopy(self.surf_temp, memo)
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
                    (self.therm_corr == other.therm_corr) and 
                    (self.atm_corr == other.atm_corr) and
                    (self.cube_i == other.cube_i).all())
        else:
            return False

    def __repr__(self):
        description = """
        OMEGA/MEx observation {0} - (v{1})
        
        Cube quality : {2}
        Thermal correction : {3}
        Atmospheric correction : {4}

        \033[3m{5}\033[0m""".format(self.name, self.version, self.quality, self.therm_corr, 
                                    self.atm_corr, self.add_infos)
        # Ajout Ls ?
        # Ajout UTC ? (utc.strftime('%d-%m-%Y %H:%M'))
        return description

##----------------------------------------------------------------------------------------
## Lecture fichier navigation .NAV
def read_header(filename):
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
    if header_dict['SPACECRAFT_POINTING_MODE'] == '"UNK"':
        header_dict['SPACECRAFT_POINTING_MODE'] = '"UNKNOWN"'
    # Sortie
    return header_dict

##----------------------------------------------------------------------------------------
## Recherche observation
def find_cube(lat, lon, cmin=0, cmax=10000, out=False):
    """Display the available OMEGA/MEx cubes with observations of the target
    latitude and longitude, using the IDL procedure `findcub.pro`.

    Parameters
    ==========
    lat : float
        The target latitude (in degrees).
    lon : float
        The target longitude (in degrees).
    cmin : float, optional (default 0)
        The minimum orbit number.
    cmax : float, optional (default 10000)
        The maximum orbit number.
    out : bool, optional (default False)
        If True -> return output

    Returns (If out == True)
    =======
    cub_list : array-like
        List of matching observations.
        Format : (orbit, x, y, dmin, altMEx, inci, emer, phas, Ls)
    """
    idl = pidly.IDL()
    idl("cd, 'omega_routines'")
    idl.pro('findcub', lon, lat, cmin, cmax)
    idl.close()
    if out:
        cub_list = np.genfromtxt('omega_routines/cubliste', skip_header=2, skip_footer=1,
                                 dtype=None, encoding='utf8')
        return cub_list

##----------------------------------------------------------------------------------------
## Sauvegarde / Importation
def save_omega(omega, savname='auto', folder='', base_folder='../data/OMEGA/',
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
    base_folder : str, optional (default '../data/OMEGA/')
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
    # Vérification syntaxe chemin
    if (len(base_folder) > 0) and (base_folder[-1] != '/'):
        base_folder += '/'
    if (len(folder) > 0) and (folder[-1] != '/'):
        folder += '/'
    # Initialisation nom fichier auto
    if savname == 'auto':
        if (len(suff)>0) and (suff[0] != '_'):
            suff = '_' + suff
        if (len(pref)>0) and (pref[-1] != '_'):
            pref = pref + '_'
        savname = '{pref}{name}{suff}.pkl'.format(name=omega.name, pref=pref, suff=suff)
    # Chemin sav fichier
    target_path = base_folder + folder + savname
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

##----------------------------------------------------------------------------------------
## Correction thermique 1 - Détermitation réflectance à 5µm à partir de celle à 2.4µm
## puis CN -- Méthode historique (Calvin & Erard 1997)
def corr_therm_sp(omega, x, y, disp=True):
    """Remove the thermal component in an OMEGA spectrum, using the historical method
    based on the reflectance determination using reference spectra from Calvin & Erard (1997).

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
    # spectels #97-#112 des spectres de ref <-> 2.3-2.5µm
    fref = '../data/OMEGA_dataref/refclair_sombr_omega_CL.dat' # from Erard and Calvin (1997)
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
        Blam = uf.planck(lam2, T) * 1e-6    # Loi de Planck en *µm-1
        sp_simu2 = sp_simu_5m * sp_sol2 * ecl + (1-sp_simu_5m) * Blam
        return sp_simu2
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
    return lam, sp_rf_corr, T_fit

def corr_therm(omega):
    """Remove the thermal component in the OMEGA hyperspectral cube.

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
            sp_rf_corr, surf_temp = corr_therm_sp(omega2, x, y, disp=False)[1:]
            omega_corr.cube_rf[y,x] = sp_rf_corr
            omega_corr.surf_temp[y,x] = surf_temp
    # Sortie
    omega_corr.therm_corr = True
    omega_corr.therm_corr_infos['datetime'] = datetime.now()
    omega_corr.therm_corr_infos['method'] = '(M1) Calvin & Erard'
    return omega_corr
    
##----------------------------------------------------------------------------------------
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
    # Fit de la température et réflectance
    T_fit, refl = curve_fit(simu_sp_5microns2, (i_lam3, i_lam4), sp_i[i_lam3:i_lam4], #p0=(280, 0.5),
                            bounds=([0, 0], [400, 0.5]))[0]
    if disp:
        print('Temperature = {0:.3f} K   |   Reflectance = {1:.5f}'.format(T_fit, refl))
    # Correction thermique spectre
    Blam = uf.planck(lam*1e-6, T_fit) * 1e-6   # En W.m-2.sr-1.µm-1
    sp_rf_corr = (sp_i - Blam) / (sp_sol*ecl - Blam)
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
    omega_corr.therm_corr_infos['datetime'] = datetime.now()
    omega_corr.therm_corr_infos['method'] = '(M2) Simultaneous refl & temp'
    return omega_corr

##----------------------------------------------------------------------------------------
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
    atmorap = np.loadtxt('../data/OMEGA_dataref/omega_atmorap_CL.dat')
    tr_atm = np.ones(nlam)
    tr_atm[nV:] = atmorap[ic_CL]    # données atm uniquement voies C & L
    # Détermination exposant
    i_lam1, i_lam2 = uf.where_closer_array([1.93, 2.01], lam)
    expo = np.log(cube_rf[:,:,i_lam1] / cube_rf[:,:,i_lam2]) / np.log(tr_atm[i_lam1] / tr_atm[i_lam2])
    # Correction spectres
    for x in tqdm(range(nx)):
        for y in tqdm(range(ny)):
            sp_rf_corr = cube_rf[y,x] * tr_atm**(-expo[y,x])
            omega_corr.cube_rf[y,x] = sp_rf_corr
    # Sortie
    omega_corr.atm_corr = True
    omega_corr.atm_corr_infos['datetime'] = datetime.now()
    omega_corr.atm_corr_infos['method'] = 'M1 : same reflectance level at 1.93μm and 2.01μm'
    return omega_corr
    
##----------------------------------------------------------------------------------------
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
    atmorap = np.loadtxt('../data/OMEGA_dataref/omega_atmorap_CL.dat')
    tr_atm = np.ones(nlam)
    tr_atm[nV:] = atmorap[ic_CL]    # données atm uniquement voies C & L
    # Détermination exposant
    i_lams = uf.where_closer_array([1.97, 1.98, 2.0], lam)
    cube_rf2 = cube_rf[:,:,i_lams]
    sp_atm2 = tr_atm[i_lams]
    expo0 = 1
    # Correction spectres
    for x in tqdm(range(nx)):
        for y in tqdm(range(ny)):
            expo = minimize(f_min, expo0, args=(cube_rf2[y,x], sp_atm2)).x[0]
            omega_corr.cube_rf[y,x] = cube_rf[y,x] * tr_atm**(-expo)
    # Sortie
    omega_corr.atm_corr = True
    omega_corr.atm_corr_infos['datetime'] = datetime.now()
    omega_corr.atm_corr_infos['method'] = 'M2 : flattest specra between 1.97µm and 2.00µm'
    return omega_corr
    
##----------------------------------------------------------------------------------------
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
        omega128_interp = readsav('../data/OMEGA_dataref/omega128_interpol.sav')
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
    
##----------------------------------------------------------------------------------------
## Correction cube OMEGA
def corr_omega(omega):
    """Remove the thermal and atmospheric component in the OMEGA hyperspectral cube.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA observation data.

    Returns
    =======
    omega_corr : OMEGAdata
        The input OMEGA observation, where the reflectance is corrected from
        the thermal and atmospheric component.
    """
    # Initialisation
    omega2 = deepcopy(omega)
    omega_corr = deepcopy(omega)
    ny, nx, nlam = omega2.cube_i.shape
    # Correction spectres
    for x in tqdm(range(nx)):
        for y in tqdm(range(ny)):
            # Correction thermique
            lam, sp_rf_corr = corr_therm_sp(omega2, x, y, disp=False)
            # Correction atmosphérique
            # sp_rf_corr2 = corr_atm_sp()
            omega_corr.cube_rf[y,x] = sp_rf_corr
    # Sortie
    omega_corr.therm_corr = True
    omega_corr.atm_corr = True
    return omega_corr

##----------------------------------------------------------------------------------------
## Sauvegarde
def test_security_overwrite(path):
    """Test if a file already exists, and if yes ask the user if he wants to
    ovewrite it or not.

    Parameters
    ==========
    path : str
        The target file path.

    Returns
    =======
    overwrite : bool
        | True -> No existent file, or overwriting allowed.
        | False -> Existent file, no overwriting.
    """
    erase = 'n'
    if glob.glob(path) != []:
        try:
            erase = input('Do you really want to erase and replace "' + path +
                        '" ? (y/N) ')
        except KeyboardInterrupt:
            erase = 'n'
        if erase != 'y' :
            print("file preserved")
            return False
        else:
            return True
    else:
        return True

def corr_save_omega(obsname, folder='auto', base_folder='../data/OMEGA/', security=True,
                    overwrite=True, compress=True):
    """Correction and saving of OMEGA/MEx observations.

    Parameters
    ==========
    obsname : str
        The name of the OMEGA observation.
    folder : str, optional (default 'auto')
        The subfolder to save the data.
        | If 'auto' -> folder = 'vX.X', where X.X is the Version of the current code.
    base_folder : str, optional (default '../data/OMEGA/')
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
    """
    if folder == 'auto':
        folder = 'v' + str(Version)
    omega = OMEGAdata(obsname)
    name = omega.name
    # path synthax
    if (base_folder != '') and (base_folder[-1] != '/'):
        base_folder += '/'
    if (folder != '') and (folder[-1] != '/'):
        folder += '/'
    basename = base_folder + folder + name + '{0}.pkl'
    # Testing existent file
    if os.path.exists(basename.format('_corr_therm_atm')):
        exists = True
    else:
        exists = False
    if security:
        overwrite = test_security_overwrite(basename.format('*'))
    if (not exists) or (exists and overwrite):
        save_omega(omega, folder=folder, base_folder=base_folder)
        print('\n\033[01mThermal correction\033[0m')
        omega_corr = corr_therm(omega)
        if compress:
            omega_corr.cube_i = None
        save_omega(omega_corr, folder=folder, base_folder=base_folder, suff='corr_therm')
        print('\n\033[01mAtmospheric correction\033[0m')
        omega_corr_atm = corr_atm(omega_corr)
        save_omega(omega_corr_atm, folder=folder, base_folder=base_folder, suff='corr_therm_atm')
    else:
        print('\n\033[01;34mExistent files preserved for {0} - v{1}\033[0m\n'.format(name, Version))

def corr_save_omega_list(liste_obs, folder='auto', base_folder='../data/OMEGA/',
                         security=True, overwrite=True, compress=True):
    """Correction and saving of a list of OMEGA/MEx observations.

    Parameters
    ==========
    liste_obs : list of str
        The list of the name of the OMEGA observations.
    folder : str, optional (default 'auto')
        The subfolder to save the data.
        | If 'auto' -> folder = 'vX.X', where X.X is the Version of the current code.
    base_folder : str, optional (default '../data/OMEGA/')
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
    """
    N = len(liste_obs)
    if folder == 'auto':
        folder = 'v' + str(Version)
    for i, obsname in enumerate(liste_obs):
        print('\n\033[01mComputing observation {0} / {1} : {2}\033[0m\n'.format(i+1, N, obsname))
        corr_save_omega(obsname, folder, base_folder, security, overwrite, compress)
    print("\n\033[01;32m Done\033[0m\n")

##----------------------------------------------------------------------------------------
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

##----------------------------------------------------------------------------------------
## End of code
##----------------------------------------------------------------------------------------
