# OMEGA-Py documentation - v2.0

## `omegapy.omega_data`

Importation and correction of OMEGA/MEx observations from binaries files.

`class OMEGAdata(obs='', empty=False, data_path=_omega_bin_path, corrV=True, corrL=True, disp=True)`

`find_cube(lon0, lat0, cmin=0, cmax=10000, out=False)`

`autosave_omega(omega, folder='auto', base_folder=_omega_py_path, security=True, disp=True)`

`autoload_omega(obs_name, folder='auto', version=_Version, base_folder=_omega_py_path, therm_corr=None, atm_corr=None, disp=True)`

`save_omega(omega, savname='auto', folder='', base_folder=_omega_py_path, pref ='', suff='', disp=True)`

`load_omega(filename, disp=True)`

`load_omega_list(basename, disp=True)`

`load_omega_list2(liste_obs, therm_corr=True, atm_corr=True, **kwargs)`

`import_list_obs_csv(filename)`

`corr_therm(omega, npool=1)`

`corr_therm2(omega)`

`corr_atm(omega)`

`corr_atm2(omega)`

`corr_save_omega(obsname, folder='auto', base_folder=_omega_py_path, security=True, overwrite=True, compress=True, npool=1)`

`corr_save_omega_list(liste_obs, folder='auto', base_folder=_omega_py_path, security=True, overwrite=True, compress=True, npool=1)`

`set_omega_bin_path(new_path)`

`set_omega_py_path(new_path)`

`get_omega_bin_path()`

`get_omega_py_path()`

`get_names(omega_list)`

`get_ls(omega_list)`

`update_cube_quality(obs_name='ORB*.pkl', folder='auto', version=_Version, base_folder=_omega_py_path)`

`test_cube(obs)`

`compute_list_good_observations(savfilename='liste_good_obs.csv', folder='../data/OMEGA/liste_obs', security=True)`

`utc_to_my(dt)`

`shared_lam(lam_list)`

`shared_lam_omegalist(omega_list)`


### OMEGAdata class
~~~python
class omegapy.omega_data.OMEGAdata(obs='', empty=False, data_path=_omega_bin_path, corrV=True, corrL=True, disp=True):
    Importation of OMEGA/MEx observation.

    Parameters
    ==========
    obs : str
        The name of the OMEGA observation.
    empty : bool, optional (default False)
        If True, return an empty OMEGAdata object.
    data_path : str, optional (default _omega_py_path)
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
        The average UCT time of the observation.
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
        Latitude of the sub-solar point at observation time.
    subsol_lon : float
        Longitude of the sub-solar point at observation time.
    min_lat : float
        Southernmost latitude of the observation (deg).
    max_lat : float
        Northernmost latitude of the observation (deg).
    min_lon : float
        Easternmost longitude of the observation (deg).
    max_lon : float
        Westernmost longitude of the observation (dego.
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
    
    Methods
    =======
    __init__(self, obs='', empty=False, data_path=_omega_bin_path)
    
    __copy__(self)

    __deepcopy(self, memo)

    __eq__(self, other)

    __repr__(self)
~~~

### Observation search
~~~python
omegapy.omega_data.find_cube(lon0, lat0, cmin=0, cmax=10000, out=False):
    Display the available OMEGA/MEx cubes with observations of the target
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

    Returns (If out == True)
    =======
    cub_list : array-like
        List of matching observations.
        Format : (orbit, x, y, dmin, altMEx, inci, emer, phas, loct, Ls, MY)
~~~

### OMEGAdata files handling
~~~python
omegapy.omega_data.autosave_omega(omega, folder='auto', base_folder=_omega_py_path, security=True, disp=True):
    Save an OMEGA object at the selected path using the pickle module, with automatic
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
~~~

~~~python
omegapy.omega_data.autoload_omega(obs_name, folder='auto', version=_Version, base_folder=_omega_py_path,
                   therm_corr=None, atm_corr=None, disp=True):
    Load and return a previously saved OMEGAdata object using pickle (with autosave_omega()).

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

    Returns
    =======
    omega : OMEGAdata 
        The loaded object of OMEGA/MEx observation.
~~~

~~~python
omegapy.omega_data.save_omega(omega, savname='auto', folder='', base_folder=_omega_py_path,
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
omegapy.omega_data.load_omega_list2(liste_obs, therm_corr=True, atm_corr=True, **kwargs)
    Load a list of saved OMEGAdata objects, using load_omega().

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
omegapy.omega_data.corr_therm(omega, npool=1):
    Remove the thermal component in the OMEGA hyperspectral cube.
    
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

~~~python
omegapy.omega_data.corr_save_omega(obsname, folder='auto', base_folder=_omega_py_path, security=True,
                    overwrite=True, compress=True, npool=1):
    Correction and saving of OMEGA/MEx observations.
    
    Parallelization is implemented using the multiprocessing module. The number of
    process to run is controlled by the npool argument.

    Parameters
    ==========
    obsname : str
        The name of the OMEGA observation.
    folder : str, optional (default 'auto')
        The subfolder to save the data.
        | If 'auto' -> folder = 'vX.X', where X.X is the given value of code version.
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
~~~

~~~
omegapy.omega_data.corr_save_omega_list(liste_obs, folder='auto', base_folder=_omega_py_path,
                         security=True, overwrite=True, compress=True, npool=1):
    Correction and saving of a list of OMEGA/MEx observations.
    
    Parallelization is implemented using the multiprocessing module. The number of
    process to run is controlled by the npool argument.

    Parameters
    ==========
    liste_obs : list of str
        The list of the name of the OMEGA observations.
    folder : str, optional (default 'auto')
        The subfolder to save the data.
        | If 'auto' -> folder = 'vX.X', where X.X is the Version of the current code.
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
~~~

~~~python
omegapy.omega_data.set_omega_bin_path(new_path):
    Set the global private _omega_bin_path variable to new_path.

    Parameters
    ==========
    new_path : str
        The new path of the OMEGA binary files (.QUB and .NAV).
~~~

~~~python
omegapy.omega_data.set_omega_py_path(new_path):
    Set the global private _omega_py_path variable to new_path.

    Parameters
    ==========
    new_path : str
        The new path of the OMEGA python-made files.
~~~

~~~python
omegapy.omega_data.get_omega_bin_path():
    Return the vavue of the global private _omega_bin_path variable.

    Returns
    =======
    omega_bin_path : str
        The path of the OMEGA binary files (.QUB and .NAV).
~~~

~~~python
omegapy.omega_data.get_omega_py_path():
    Return the vavue of the global private _omega_py_path variable.

    Returns
    =======
    omega_py_path : str
        The new path of the OMEGA python-made files.
~~~

~~~python
omegapy.omega_data.get_names(omega_list):
    Return the array of the observation ID of each OMEGA/MEx observation in omega_list.

    Parameters
    ==========
    omega_list : array of OMEGAdata
        The input array of OMEGA observations.
    
    Returns
    =======
    names : ndarray
        The array of the omega_list observations ID.
~~~

~~~python
omegapy.omega_data.get_ls(omega_list):
    Return the array of the Solar longitude of each OMEGA/MEx observation in omega_list.

    Parameters
    ==========
    omega_list : array of OMEGAdata
        The input array of OMEGA observations.
    
    Returns
    =======
    ls : ndarray
        The array of the omega_list Ls.
~~~


### Update quality
~~~Python
omegapy.omega_data.update_cube_quality(obs_name='ORB*.pkl', folder='auto', version=_Version, 
                        base_folder=_omega_py_path):
    Update the quality attribute of previously saved OMEGAdata objects.

    Parameters
    ==========
    obs_name : str, optional (default 'ORB*.pkl')
        The files basename.
    folder : str, optional (default 'auto')
        The subfolder where the data is.
        | If 'auto' -> folder = 'vX.X', where X.X is the given value of code version.
    version : float, optional (default _Version)
        The version of the target file (if folder is 'auto').
        Default is the current code version.
    base_folder : str, optional (default _omega_py_path)
        The base folder path.
~~~


### Testing observations quality
~~~python
omegapy.omega_data.test_cube(obs):
    Test the quality of an OMEGA/MEx observation from the header informations
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
~~~

~~~python
omegapy.omega_data.compute_list_good_observations(savfilename='liste_good_obs.csv', 
                                   folder='../data/OMEGA/liste_obs', security=True):
    Scan the available OMEGA/MEx data cubes and list the observations considered as 
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
~~~

### Conversion UTC to MY
~~~python
omegapy.omega_data.utc_to_my(dt):
    Convert a UTC datetime to the corresponding Martian Year (MY).
    
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
~~~

### Shared wavelength array between different observations
~~~python
omegapy.omega_data.shared_lam(lam_list):
    Return a list of wavelength shared by all the input wavelength arrays.

    Parameters
    ==========
    lam_list : list of 1D np.array
        The list of wavelength array.

    Returns
    =======
    lam2 : 1D np.array
        The wavelength array that contains only wavelength shared by all the arrays of
        lam_list.
~~~

~~~python
omegapy.omega_data.shared_lam_omegalist(omega_list):
    Return a list of wavelength shared by all the wavelength arrays of the input
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
~~~
