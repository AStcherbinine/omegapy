# OMEGA-Py documentation - v1.2

## [`omegapy.omega_data`](doc_omega_data.md)

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


## [`omegapy.omega_plots`](doc_omega_plots.md)

Display of OMEGAdata cubes.

`omegapy.omega_plots.show_omega(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, title='auto', xlim=(None, None), ylim=(None, None), Nfig=None)`

`omegapy.omega_plots.show_omega_v2(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, alpha=None, title='auto', lonlim=(None, None), latlim=(None, None), Nfig=None, polar=False, cbar=True)`

`omegapy.omega_plots.show_omega_interactif(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, title='auto', autoyscale=True, xlim=(None, None), ylim=(None, None))`

`omegapy.omega_plots.show_omega_interactif_v2(omega, lam, refl=True, lam_unit='m', cmap='Greys_r', vmin=None, vmax=None, title='auto', autoyscale=True, alpha=None, lonlim=(None, None), latlim=(None, None), polar=False)`

`omegapy.omega_plots.show_ibd_v2(omega, ibd, cmap='viridis', vmin=None, vmax=None, alpha=None, title='auto', cb_title = 'IBD', lonlim=(None, None), latlim=(None, None), Nfig=None, polar=False, cbar=True)`

