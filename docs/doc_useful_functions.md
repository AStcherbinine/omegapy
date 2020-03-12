# OMEGA-Py documentation - v1.2

## `omegapy.useful_functions`

Some useful generic functions.

`where_closer(value, array)`

`where_closer_array(values, array)`

`myglob(basename)`

`sort_dict(dico)`

`save_pickle(obj, target_path, disp=True)`

`load_pickle(filename, disp=True)`

`reg_lin(X, Y)`

`planck(lam, T)`

`fit_black_body(lam, sp, T_bounds=(0, 1e6))`

`moyenne_glissante(sp, n)`

`filtre_median(sp, n)`


### Index search in an array

~~~python
omegapy.useful_functions.where_closer(value, array):
    Renvoie l'indice de la valeur la plus proche de celle recherchée dans array.
    
    Parameters
    ==========
    values : float
        Searched value.
    array : ndarray
        The array.

    Returns
    =======
    i : int
        The index of the closer value to value in array.
~~~

~~~python
omegapy.useful_functions.where_closer_array(values, array):
    Renvoie la liste des indices des valeurs les plus proches de celles recherchées
    dans array.

    Parameters
    ==========
    values : ndarray
        Array of searched values.
    array : ndarray
        The array.

    Returns
    =======
    I : ndarray
        Array of the index of the closer values in the array.
~~~

### File search
~~~python
omegapy.useful_functions.myglob(basename):
    Return the absolute path according to the input basename.
    If mulitple files corresponds to the basename, the user will be asked
    to choose one.

    --------------------------------------------
    | int -> Select the corresponding filename.
    | q/quit/exit -> Return None.
    | a/all -> Return the list of all filenames.
    --------------------------------------------

    Parameters
    ==========
    basename : str
        The basename of the target file.

    Returns
    =======
    fname : str
        The absolute path of the selected file.
~~~

### Sorting
~~~python
omegapy.useful_functions.sort_dict(dico):
    Sort a dictionary by its keys values.

    Parameters
    ==========
    dico : dict
        The input unsorted dictionary.

    Returns
    =======
    dico_sorted : dict
        The sordet dictionary.
~~~

### Saving / loading objects
~~~python
omegapy.useful_functions.save_pickle(obj, target_path, disp=True):
    Save an object at the selected path using the pickle module.

    Parameters
    ==========
    obj : Object
        The object to save.
    target_path : str
        The saving path name.
    disp : bool
        Control the display.
            | True -> Print the saving filename.
            | False -> Nothing printed.
~~~

~~~python
omegapy.useful_functions.load_pickle(filename, disp=True):
    Load and return a previously saved object with pickle.

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
    obj : Object
        The loaded object.
~~~

### Curve fitting
~~~python
omegapy.useful_functions.reg_lin(X, Y):
    Renvoie le résultat de la régression linéaire ( f(x) = a*x + b ) sur les valeurs 
    en entrée.
    
    Parameters
    ==========
    X : ndarray
        The X-values.
    Y : ndarray
        The Y-values.

    Returns
    =======
    a : float
        Slope of the fitted line.
    b : float
        Origin ordinate of the fitted line.
~~~

~~~python
omegapy.useful_functions.planck(lam, T):
    Return the Black body radiance, associated to the input wavelength and
    temperature. According to the Planck's law.

    Parameters
    ==========
    lam : float or array-like
        The wavelength (in m).
    T : float
        The temperature (in K).
    
    Returns
    =======
    B_lam : float or array-like
        The spectral radiance (in W.m-2.sr-1.m-1).
~~~

~~~python
omegapy.useful_functions.fit_black_body(lam, sp, T_bounds=(0, 1e6)):
    Return the temperature associated to the fitted black body thermical
    spectrum.

    Parameters
    ==========
    lam : array-like
        The wavelength array (in m).
    sp : array-like
        The spectral radiance (in W.m-2.sr-1.m-1) to be fitted.
    bounds : 2-tuple, optional (default (0, 1e6))
        The bounds for the temperature fitting.

    Returns
    =======
    T : float
        The temperature of the fitted Planck's law radiance (in K).
~~~

### Smoothing
~~~python
omegapy.useful_functions.moyenne_glissante(sp, n):
    Applique un filtre de moyenne glissante sur les valeurs du spectre en remplaçant 
    chaque valeur par la moyenne des valeurs dans un intervalle de dimension 2n+1 centré 
    sur la valeur en question.

    Parameters
    ==========
    sp : ndarray
        Array of the transmittance values.
    n : int
        The len of the window of the moving average is 2n+1.

    Returns
    =======
    sp_med : ndarray
        Filtered transmittance array.
~~~

~~~python
omegapy.useful_functions.filtre_median(sp, n):
    Applique un filtre médian sur les valeurs du spectre en remplaçant chaque valeur 
    par la médiane des valeurs dans un intervalle de dimension 2n+1 centré sur la valeur
    en question.

    Parameters
    ==========
    sp : ndarray
        Array of transmittance values.
    n : int
        The len of the window the moving median is 2n+1.

    Returns
    =======
    sp_med : ndarray
        Filtered transmittance array.
~~~
