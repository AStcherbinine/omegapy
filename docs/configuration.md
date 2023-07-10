# Configuration
You have to configure the default path of the folders containing the OMEGA binary and omegapy-made files
using the environment variables `OMEGA_BIN_PATH` (for the binary .QUB and .NAV files)
and `OMEGA_PY_PATH` (for the omegapy-made files).

## Linux
To do so, add the following lines to your `~/.bashrc` (or `~/.zshrc`, depending on which shell you are using):
~~~bash
export OMEGA_BIN_PATH="/path/to/binary/files/folder/"
export OMEGA_PY_PATH="/path/to/omegapy-made/files/folder/"
~~~
*Adapt the path to suit your own architecture.*

**Tip:** You can check that these variables are properly set up by typing the following command in a new terminal: `echo $OMEGA_BIN_PATH` and `echo $OMEGA_PY_PATH`.
It should print the path you previously set.

## MacOS
Similar to the Linux procedure, except that the `~/.bashrc` file may not be loaded automatically by default.
In that case, use instead `~/.bash_profile`.

**Note for more recent MacOS versions:** The default shell has been changed from bash to zsh in the more recent versions of MacOS. Thus, if you are using a zsh shell, edit the `~/.zshrc` file instead of `~/.bash_profile` or `~/.bashrc`.

## Windows (or if you have troubles using the environment variables)
If you are using Windows, you cannot easily set these environment variables.
Lucky you, there is a solution!

Note that it also apply if you are using another OS but had troubles setting the environment variables as described above (i.e., you are seeing these warnings when loading omegapy: `Warning: $OMEGA_BIN_PATH not defined` and/or `Warning: $OMEGA_PY_PATH not defined`).

In that case, you can set these path directly with Python using the `omega_data.set_omega_bin_path()` and `omega_data.set_omega_py_path()` functions.
Assuming you have already load `omegapy.omega_data` as `od`, simply execute:
~~~python
od.set_omega_bin_path("/path/to/binary/files/folder/")
od.set_omega_py_path("/path/to/omegapy-made/files/folder/")
~~~
*Adapt the path to suit your own architecture.*

You will have to run these commands everytime you start a new Python console, so I suggest to put these lines at the beginning of your script, just after the omegapy import.

