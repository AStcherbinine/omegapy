# OMEGA-Py documentation - v1.2

## [`omegapy.omega_data`](doc_omega_data.md)

Importation of OMEGA observations in the OMEGAdata class.
Using IDL routines containing in *./omega_routines/*.

`class OMEGAdata(obs='', empty=False, data_path=omega_bin_path)`

`find_cube(lat, lon, cmin=0, cmax=10000, out=False)`

`autosave_omega(omega, folder='auto', base_folder=omega_py_path, security=True, disp=True)`

`autoload_omega(obs_name, folder='auto', version=Version, base_folder=omega_py_path, disp=True)`

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

`update_cube_quality(obs_name='ORB*.pkl', folder='auto', version=Version, base_folder=omega_py_path)`


## [`omegapy.omega_plots`](doc_omega_plots.md)

Display of OMEGAdata cubes.

`show_omega(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, title='auto', xlim=(None, None), ylim=(None, None), Nfig=None)`

`show_omega_v2(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, alpha=None, title='auto', lonlim=(None, None), latlim=(None, None), Nfig=None, polar=False, cbar=True)`

`show_omega_interactif(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, title='auto', autoyscale=True, xlim=(None, None), ylim=(None, None))`

`show_omega_interactif_v2(omega, lam=1.085, refl=True, lam_unit='m', data=None, cmap='Greys_r', cb_title='IBD', title='auto', vmin=None, vmax=None, autoyscale=True, ylim_sp=(None, None), alpha=None, lonlim=(None, None), latlim=(None, None), polar=False, cbar=True, grid=True)`

`show_ibd_v2(omega, ibd, cmap='viridis', vmin=None, vmax=None, alpha=None, title='auto', cb_title = 'IBD', lonlim=(None, None), latlim=(None, None), Nfig=None, polar=False, cbar=True)`

`show_omega_list_v2(omega_list, lam=1.085, lat_min=-90, lat_max=90, lon_min=0, lon_max=360, pas_lat=0.1, pas_lon=0.1, cmap='Greys_r', vmin=None, vmax=None, title='auto', Nfig=None, polar=False, cbar=True, cb_title='auto', data_list=None, plot=True, grid=True, out=False, **kwargs)`


## [`omegapy.useful_functions`](doc_useful_functions.md)

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
