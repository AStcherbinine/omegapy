# OMEGA-Py documentation - v1.2

## [`omegapy.omega_data`](doc_omega_data.md)

Importation of OMEGA observations in the OMEGAdata class.
Using IDL routines containing in *./omega_routines/*.

`class OMEGAdata(obs='', empty=False, data_path=omega_bin_path)`

`find_cube(lat, lon, cmin=0, cmax=10000, out=False)`

`save_omega(omega, savname='auto', folder='', base_folder=omega_py_path, pref ='', suff='', disp=True)`

`load_omega(filename, disp=True)`

`load_omega_list(basename, disp=True)`

`import_list_obs_csv(filename)`

`corr_therm(omega)`

`corr_therm2(omega)`

`corr_atm(omega)`

`corr_atm2(omega)`

`corr_save_omega(obsname, folder='auto', base_folder=omega_py_path, security=True, overwrite=True, compress=True)`

`corr_save_omega_list(liste_obs, folder='auto', base_folder=omega_py_path, security=True, overwrite=True, compress=True)`

`set_omega_bin_path(new_path)`

`set_omega_py_path(new_path)`

`get_omega_bin_path()`

`get_omega_py_path()`

## [`omegapy.omega_plots`](doc_omega_plots.md)

Display of OMEGAdata cubes.

`show_omega(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, title='auto', xlim=(None, None), ylim=(None, None), Nfig=None)`

`show_omega_v2(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, alpha=None, title='auto', lonlim=(None, None), latlim=(None, None), Nfig=None, polar=False, cbar=True)`

`show_omega_interactif(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, title='auto', autoyscale=True, xlim=(None, None), ylim=(None, None))`

`show_omega_interactif_v2(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, title='auto', autoyscale=True, alpha=None, lonlim=(None, None), latlim=(None, None), polar=False)`

`show_ibd_v2(omega, ibd, cmap='viridis', vmin=None, vmax=None, alpha=None, title='auto', cb_title = 'IBD', lonlim=(None, None), latlim=(None, None), Nfig=None, polar=False, cbar=True)`

