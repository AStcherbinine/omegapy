# OMEGA Py documentation - v1.2

## `omegapy.omega_data`

Importation of OMEGA observations in the OMEGAdata class.
Using IDL routines containing in *./omega_routines/*.

`class omegapy.omega_data.OMEGAdata(obs='', empty=False, data_path='/data2/opt/geomeg/data/product/')`

`omegapy.omega_data.find_cube(lat, lon, cmin=0, cmax=10000, out=False)`

`omegapy.omega_data.save_omega(omega, savname='auto', folder='', base_folder='../data/OMEGA/', pref ='', suff='', disp=True)`

`omegapy.omega_data.load_omega(filename, disp=True)`

`omegapy.omega_data.load_omega_list(basename, disp=True)`

`omegapy.omega_data.import_list_obs_csv(filename)`

`omegapy.omega_data.corr_therm(omega)`

`omegapy.omega_data.corr_therm2(omega)`

`omegapy.omega_data.corr_atm(omega)`

`omegapy.omega_data.corr_atm2(omega)`

`omegapy.omega_data.corr_save_omega(obsname, folder='auto', base_folder='../data/OMEGA/', security=True, overwrite=True, compress=True)`

`omegapy.omega_data.corr_save_omega_list(liste_obs, folder='auto', base_folder='../data/OMEGA/', security=True, overwrite=True, compress=True)`


### OMEGAdata class
~~~python
class omegapy.omega_data.OMEGAdata(obs='', empty=False, data_path='/data2/opt/geomeg/data/product/'):
    Importation of OMEGA/MEx observation, using the readomega_vpy.pro IDL routine.

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
    
    Methods
    =======
    __init__(self, obs='', empty=False, data_path='/data2/opt/geomeg/data/product/')
    
    __copy__(self)

    __deepcopy(self, memo)

    __eq__(self, other)

    __repr__(self)
~~~

### Observation search
~~~python
omegapy.omega_data.find_cube(lat, lon, cmin=0, cmax=10000, out=False):
    Display the available OMEGA/MEx cubes with observations of the target
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
~~~

### OMEGAdata files handling
~~~python
omegapy.omega_data.save_omega(omega, savname='auto', folder='', base_folder='../data/OMEGA/',
               pref ='', suff='', disp=True):
    Save an OMEGA object at the selected path using the pickle module.

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
~~~

~~~python
omegapy.omega_data.load_omega(filename, disp=True):
    Load and return a previously saved OMEGAdata object (with save_omega()).

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
~~~

~~~python
omegapy.omega_data.load_omega_list(basename, disp=True):
    Load a list of saved OMEGAdata objects, using load_omega().

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
~~~

~~~python
omegapy.omega_data.import_list_obs_csv(filename):
    Import a list of observations ID from a csv file generated by JMars.

    Parameters
    ==========
    filename : str
        The target path of the csv file.

    Returns
    =======
    liste_obs : array of str
        The list of observations ID from the csv file.
~~~

### Thermal correction
Correction thermique 1 - Détermitation réflectance à 5µm à partir de celle à 2.4µm puis CN -- Méthode historique (Calvin & Erard 1997)

~~~python
omegapy.omega_data.corr_therm(omega):
    Remove the thermal component in the OMEGA hyperspectral cube.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA observation data.

    Returns
    =======
    omega_corr : OMEGAdata
        The input OMEGA observation, where the reflectance is corrected from
        the thermal component.
~~~

Correction thermique 2 - Fit CN et réflectance en même temps

~~~python
omegapy.omega_data.corr_therm2(omega):
    Remove the thermal component in the OMEGA hyperspectral cube, 
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
~~~

### Atmospheric correction
Correction atmosphérique 1 - 1.93 µm & 2.01 µm au même niveau de réflectance

~~~python
omegapy.omega_data.corr_atm(omega):
    Remove the atmospheric component in the OMEGA hyperspectral cube.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA observation data.

    Returns
    =======
    omega_corr : OMEGAdata
        The input OMEGA observation, where the reflectance is corrected from
        the atmospheric component.
~~~

Correction atmosphérique 2 - Spectre plat entre 1.97 µm et 2.00 µm

~~~python
omegapy.omega_data.corr_atm2(omega):
    Remove the atmospheric component in the OMEGA hyperspectral cube.

    Parameters
    ==========
    omega : OMEGAdata
        The OMEGA observation data.

    Returns
    =======
    omega_corr : OMEGAdata
        The input OMEGA observation, where the reflectance is corrected from
        the atmospheric component.
~~~

### Correction & Saving
Import an OMEGA/MEx observation, apply a thermal and atmospheric correction (M1) and save the OMEGAdata object at each step.

~~~
omegapy.omega_data.corr_save_omega(obsname, folder='auto', base_folder='../data/OMEGA/', security=True,
                    overwrite=True, compress=True):
    Correction and saving of OMEGA/MEx observations.

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
~~~

~~~
omegapy.omega_data.corr_save_omega_list(liste_obs, folder='auto', base_folder='../data/OMEGA/',
                         security=True, overwrite=True, compress=True):
    Correction and saving of a list of OMEGA/MEx observations.

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
~~~
