# OMEGA-Py documentation - v2.0

## [`omegapy.omega_data`](doc_omega_data.md)

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


## [`omegapy.omega_plots`](doc_omega_plots.md)

Display of OMEGAdata cubes.

`show_omega(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, title='auto', xlim=(None, None), ylim=(None, None), Nfig=None)`

`show_omega_v2(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, alpha=None, title='auto', lonlim=(None, None), latlim=(None, None), Nfig=None, polar=False, cbar=True, grid=True, mask=None, negatives_longitudes='auto')`

`show_omega_interactif(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, title='auto', autoyscale=True, xlim=(None, None), ylim=(None, None))`

`show_omega_interactif_v2(omega, lam=1.085, refl=True, lam_unit='m', data=None, cmap='Greys_r', cb_title='data', title='auto', vmin=None, vmax=None, autoyscale=True, ylim_sp=(None, None), alpha=None, lonlim=(None, None), latlim=(None, None), polar=False, cbar=True, grid=True, mask=None, lam_mask=None, negatives_longitudes='auto')`

`show_data_v2(omega, data, cmap='viridis', vmin=None, vmax=None, alpha=None, title='auto', cb_title = 'data', lonlim=(None, None), latlim=(None, None), Nfig=None, polar=False, cbar=True, grid=True, mask=None, negatives_longitudes='auto')`

`show_omega_list_v2(omega_list, lam=1.085, lat_min=-90, lat_max=90, lon_min=0, lon_max=360, pas_lat=0.1, pas_lon=0.1, cmap='Greys_r', vmin=None, vmax=None, title='auto', Nfig=None, polar=False, cbar=True, cb_title='auto', data_list=None, mask_list=None, negative_values=False, plot=True, grid=True, out=False, negatives_longitudes=False, **kwargs)`

`save_map_omega_list(omega_list, lat_min=-90, lat_max=90, lon_min=0, lon_max=360, pas_lat=0.1, pas_lon=0.1, lam=1.085, data_list=None, data_desc='', mask_list=None, negative_values=False, sav_filename='auto', ext='', base_folder='../data/OMEGA/sav_map_list_v2/', sub_folder='')`

`load_map_omega_list(filename)`

`show_omega_list_v2_man(data, grid_lat, grid_lon, infos, cmap='Greys_r', vmin=None, vmax=None, title='auto', Nfig=None, polar=False, cbar=True, cb_title='auto', grid=True, negatives_longitudes=False, **kwargs)`

`plot_psp(sp1_id, *args, sp2_id=(None, None), Nfig=None, sp_dict=picked_spectra, **kwargs)`


## [`omegapy.useful_functions`](doc_useful_functions.md)

Some useful generic functions.

`where_closer(value, array)`

`where_closer_array(values, array)`

`myglob(basename, exclude=[])`

`sort_dict(dico)`

`save_pickle(obj, target_path, disp=True)`

`load_pickle(filename, disp=True)`

`test_security_overwrite(path)`

`reg_lin(X, Y)`

`planck(lam, T)`

`fit_black_body(lam, sp, T_bounds=(0, 1e6))`

`moyenne_glissante(sp, n)`

`filtre_median(sp, n)`
