## Check if the module has been installed

After installing the module with `pip` (as described [here](../installation)), you can check that the module has been actually
installed in your environment by opening a Python terminal and trying to import the module:

~~~python
import omegapy
import omegapy.omega_data as od
import omegapy.omega_plots as op
import omegapy.useful_functions as uf
~~~

If these commands do not raise any error, the *OMEGA-Py* module is installed in your Python environment.

!!! failure
    If the module is not available, here are a few things that you can try:

     * Restart your terminal or application to reload the environment.
     * Check that you are in the expected Python environment.
     * Re-install the module:
       ~~~bash
        pip3 install omegapy
       ~~~

    If it is still not working, does it work if you try to install an other module?
    Search for how to install a Python module with your local setup.


## Check the module version

You can also check that you have the latest version of the module installed:

![omegapy-version](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fraw.githubusercontent.com%2FAStcherbinine%2Fomegapy%2Fmaster%2Fomegapy%2Fpackage.json&query=%24.version&prefix=v&logoColor=white&label=omegapy&labelColor=grey&color=blue&link=https%3A%2F%2Fpypi.org%2Fproject%2Fomegapy%2F)

**From a Python console**
~~~python
omegapy.__version__
~~~

**From a terminal with pip**
~~~bash
pip freeze | grep omegapy
~~~

!!! failure
    If you don't have the last version, update it with
    ~~~bash
    pip3 install omegapy --upgrade
    ~~~


## Check the paths configuration

Please refer to [this page](../configuration) for how to configure the default paths for the OMEGA binary and the 
omegapy-made files.

You can check the values stored in them by displaying the output of
[`od.get_omega_bin_path()`](../reference/omega_data/#omega_data.get_omega_bin_path)
(OMEGA *.QUB* & *.NAV* binary files) and
[`od.get_omega_py_path()`](../reference/omega_data/#omega_data.get_omega_py_path)
(omegapy-made files).

!!! failure
    If these functions do not return the expected path, please refer to
    [this paragraph](../configuration/#windows-or-if-you-have-troubles-using-the-environment-variables)
    and use the 
    [`od.set_omega_bin_path()`](../reference/omega_data/#omega_data.set_omega_bin_path)
    and
    [`od.set_omega_py_path()`](../reference/omega_data/#omega_data.set_omega_py_path)
    functions.


## Test to process and display an OMEGA observation

The final test would be to actually download and process an OMEGA observation.

For instance, download the data from cube 2 of orbit 0979:

 * [ORB0979_2.QUB](https://archives.esac.esa.int/psa/ftp/MARS-EXPRESS/OMEGA/MEX-M-OMEGA-2-EDR-FLIGHT-V1.0/DATA/ORB09/ORB0979_2.QUB)
 * [ORB0979_2.NAV](https://archives.esac.esa.int/psa/ftp/MARS-EXPRESS/OMEGA/MEX-M-OMEGA-2-EDR-FLIGHT-V1.0/DATA/GEM09/ORB0979_2.NAV)

Then download the [`example_test_script.py`](example_test_script.py) Python file and run it
to load, process, and display this OMEGA observation.

??? abstract "Source code for [`example_test_script.py`](example_test_script.py)"
    ~~~python title="example_test_script.py"
    --8<-- "docs/example_test_script.py"
    ~~~

This example script will:

 * load the binary files;
 * apply the atmospheric and thermal corrections to the OMEGA observation;
 * compute the 1.5μm band depth criteria (see the [band depth example](../example1_band_depth));
 * display and save the reflectance, 1.5μm BD, and overplotted maps 
   (see the [overplotting example](../example2_overplot));
 * display the interactive visualization tool (if running the code in an IPython shell).

You can then check that you have reproduced the figures shown in the examples.

!!! tip "Notes"
    If running the code directly from a bash terminal (with `python3 example_test_script.py`),
    the data will be processed, and the figures saved in your current directory.
    But you will not have access to the interactive visualization of the cube.

    To use the interactive aspect of *OMEGA-Py*, open an IPython terminal to run the code
    with `#!python run example_test_script.py` (or execute it from your favorite IDE with
    a IPython console).


!!! failure
    * If encountering an error in the data loading/processing, check your installation with the above steps.
    * If encountering an error with displaying the interactive map, check that you have an appropriate `matplotlib`
      backend installed and activated (e.g., `PyQt5`) and
      refer to [this paragraph](../data_visualization/#interactive-visualization).
