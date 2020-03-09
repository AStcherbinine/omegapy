# OMEGAPy documentation - v1.2

## omegapy.omega_data

Importation of OMEGA observations in the OMEGAdata class.
Using IDL routines containing in ./omega_routines/.

### OMEGAdata class
~~~python
omegapy.omega_data.OMEGAdata(obs='', empty=False, data_path='/data2/opt/geomeg/data/product/'):
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
    OMEGAdata.__init__(self, obs='', empty=False, data_path='/data2/opt/geomeg/data/product/')
    
    OMEGAdata.__copy__(self)

    OMEGAdata.__deepcopy(self, memo)

    OMEGAdata.__eq__(self, other)

    OMEGAdata.__repr__(self)
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
