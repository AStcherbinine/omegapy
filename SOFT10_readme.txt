OMEGA software (IDL version)
release SOFT08
april 2012
********************************************************
changes from SOFT01 to SOFT02
- new version of readomega which includes flat-fielding and
  second order corrections for the VIS channel
- new photometric function (mtf050110_25.dat and mtf050110_50.dat)
  with some adjustments for the VIS channel

Instrumental problems:
- additional information on bad spectels
- additional information on a possible digital saturation
  of the dark in the L channel
- complementary remarks on the co-registration of channels, with
  a warning on the validity range of the L channel
- numerical perturbation for the 128 pixel modes after orbit 511
********************************************************
changes from SOFT02 to SOFT03
- more detailed explanation on omega_path
- improved pipeline for VIS channel deconvolution. The impact of 
  the 2 msec read time of the CCD is now taken into account. 
  the offset correction, bias correction, flat and transfer function
  in the visible are also modified (new photometric function: mtf050701)
  new routines: smear_corr.pro, offset_corr.pro

! This results in a much improved overlap with SWIR-C at 1 µm.
  it should however be noted that the visible channel has an effective
  spatial resolution which is close to 1 pixel along the line, 
  but close to 4 pixels (5 mrad) along the swath. This prevents
  meaningful comparison of VIS and IR spectra for spectral features
  of very small patches (a few pixels), as they are much weaker
  in the lower resolution VIS channel.
  
- updated information on instrumental effects
********************************************************
changes from SOFT03 to SOFT04

- a new version of readomega is provided, which handles the evolution
  of IR spectels due to radiation damage. Unreliable spectels are
  given very low values. An evolving photometric function is used
  for spectels which are still stable, but with a reduced efficiency 
  
  (see corresponding section "time evolution of bad and hot IR spectels")
  
  a sligthly revised version of the photometric function (for spectel 69) 
  is now implemented so as to be compatible with the new version of
  readomega
  
- the new readomega also removes an estimate of the thermal background 
  from optics
  
- a routine "readpart.pro" is provided so as to read part of a large
  cube when memory size becomes a problem
  the syntax and result are similar to readomega, but the user is
  requested to provide a number of scans to skip and a number of
  scans to read in the cube. Entering "0" as the number of scans
  to read will result in reading all remaining scans
********************************************************
changes from SOFT04 to SOFT05
- recovery procedure for the lack of dark data for the last scan 
  in summed modes (the dark from the previous scan is considered)
  This resulted in unusable spectra for the last line of summed cubes 
- update of the source files for following spectel evolution 
  ("bound" and "rap") 
  
********************************************************
changes from SOFT05 to SOFT06
- IMPORTANT !!! this release includes a first cut at the issue
  of the L channel evolution with time. Specific photometric functions
  for the L channel are derived for each specific orbit number depending
  on the evolution of the internal calibration signal. All orbits should
  now be usable for issues such as temperature measurements which require
  a reliable photometric response for this channel. Some caution is still
  required for transition regions between states of the L channel
  (this can be seen when the photometric function changes rapidly from 
  one observation to the next). 
  the present release provides a first estimate based on internal 
  calibration values up to orbit 6410 (file rapmtflcur.bin). 
  If readomega or readpart is applied to later orbits 
  (which can only be the case for coI's), the following warning:  
  "warning: L channel corrected as for orbit 6410" 
  will be issued, so as to indicate that the L channel correction 
  at that stage is not yet based on actual analysis of the internal
  calibration results. 
   
- update of the source files for following spectel evolution 
  ("bound" and "rap") have been replaced by a currently applicable version
  (boundcur.dat, rapcur_25.dat, rapcur_50.dat). Those in SOFT06 correspond
  to spectel evolution until 31/12/2008
*******************************************************
changes from SOFT06 to SOFT07

- implementation of the DATA_QUALITY keyword
  a data_quality variable is provided as an output of readomega (readpart). 
  It incidates possible problems within the data cube. 
  A data cube is transmitted as sets of 64 spectra ("slices"). There are 
  500 to 1500 slices in a typical data cube. Each slice is transmitted
  as 6 TM packets. A typical cube corresponds to 3000 to 9000 TM packets
  The values of this index are the following:
     5 : perfect
     4 : one data gap   (gap of 40 TM packets or more)
     3 : missing data   (up to 9 TM packets are missing or corrupted)
     2 : acceptable     (< 3 gaps, < 10 isolated missing packets)
     1 : poor           (3 to 5 gaps or 10 to 100 isolated missing packets)
     0 : bad            (> 5 gaps or > 100 missing/corrupted packets)

- update of the source files for following spectel evolution 
  the "currently applicable version" (boundcur.dat, rapcur_25.dat, 
  rapcur_50.dat). Those in SOFT07, including "rapmtflcur.bin",
  relevant to L channel correction, correspond to spectel evolution 
  and flight calibration levels until orbit 7790 (29/01/2010)
**********************************************************************
changes from SOFT07 to SOFT08
- IMPORTANT !!!
the major change is linked to the loss of the C cooler after orbit 8485. 
For all subsequent observations, the C channel is not switched on, and
many observations use only the VIS channel. The exposure time has been 
set to 0 for channels which are not switched on (all corresponding data
are also set to 0)
- implementation of a new visible mtf (mtf120315_25.dat and 
  mtf120315_50.dat) provided by F. Altieri
- extension of the L channel correction to orbit 9123
**********************************************************************
changes from SOFT08 to SOFT09
- adjustment of the determination of the "ic" array (valid spectels)
  so as to eliminate all C channel spectels after orbit 8500
- extension of the L channel correction to orbit 11252
**********************************************************************
changes from SOFT09 to SOFT10
-extension of the L channel correction to orbit 13131
**********************************************************************

 
unzipping SOFT10.ZIP will create a new subdirectory, SOFT10, so that
further releases (SOFT11, ...) will not erase previous versions
(this can be useful for comparing the results of the different updates). 

So as to use the new version, the user should copy the full content of 
SOFT10 to the OMEGA working directory, from where IDL should be started.

IMPORTANT ! The first thing a user should do after unzipping a new
software release and copying the contents to the working directory 
is to modify the file:
         omega_path
in the working directory (from which IDL is run)
replacing the first line with the path to the directory
where the user puts data files (.QUB), and the second line
with the path to the directory where the user puts the geometry files (.NAV)
the two directories can be the same
a "/" must be appended for UNIX, LINUX...
a "\" must be appended for windows

An OMEGA observation corresponds to data obtained with the same observation
parameters (scan length, integration time for the VIS and IR channels, 
downtrack summing, compression). There are in general several OMEGA 
observations for the same orbit, as the scan length is changed depending
on the spacecraft altitude so as to best match the sampling rate and
the drift velocity on the surface.

The OMEGA data files and geometry files are named from the corresponding
OMEGA observation: ORBNNNN_S where NNNN is the orbit number (format I4.4)
and S is the rank of the OMEGA observation on that orbit (starting with 0).
The extensions are .QUB for the data, .NAV for the geometry
example: ORB0018_1.QUB and ORB0018_1.NAV correspond to the second 
observation on orbit 18 (scan length: 32 pixels). 

The main procedure for reading level 1b OMEGA data is readomega.pro

After entering IDL and typing ".run readomega", one has simply to enter 
the name of the requested OMEGA observation in answer to the query

example:

IDL (enter)
.run readomega (enter)
OMEGA observation:
obsname (enter)      /* obsname example: ORB0072_2 */

**** NOTE **** THE NAME MUST BE ENTERED WITHOUT THE EXTENSION

readomega first checks that obsname.QUB exists in the directory
indicated by the path for data files in "omega_path". 

If the cube is not found, it prints:
 "file OMEGA_PATH_DATA/obsname.QUB not found"
and goes back to the query. Otherwise, readomega reads the
pds cube, then it look for the corresponding geometry cube 
obsname.NAV in the directory indicated by the path for geometry files 
in the "omega_path" file. If there is no such file, it prints: 
"no corresponding NAV cube". If the NAV file exists, readomega reads 
the geometry cube. 

the resulting variables are:

idat(npix,352,nscan): integer 16 bits, raw OMEGA data
      npix: 16, 32, 64 or 128 pixels
      nscans: number of scans in the cube
      the 352 wavelengths are in the following order:
      IR C (128 values), then IR L (128 values), then Visible (96 values)

wvl(352) OMEGA wavelengths (file 'lambda_0304.dat'). these wavelengths
      are valid at the center of the IR and VIS field. The wavelength
      range changes only slightly across the field of view, in particular
      in the IR.

jdat(npix, 352, nscan): floating point, OMEGA data physical units 
      (W/m2/steradian/µm)
      There is a second order contribution to the last 35 spectels of 
      the VIS channel, which is corrected. Nevertherless, these spectels
      are not as reliable as spectels 256 to 315. The L channel (128:255)
      shows variations in terms of internal calibration signal, which
      do not correlate simply with DN levels from Mars or Phobos

note: the last 4 scans (16 pixels), 2 scans (32 pixels)
      or 1 scan (64, 128 pixels) of idat and jdat
      have only IR data (spectels 0 to 255)

There is calibration data at the beginning of each cube for the visible
channel, at the beginning of the ORBNNNN_0.QUB cubes for the IR channels
The number of calibration scans depend on the pixel length (16 to 128)
and the summation (1, 2 or 4 for 128 pixel modes)

      for every cube, the first scan (128 pixels x 4), 3 scans (128 x 2) 
      7 scans (128 x 1), 14 scans (64 pixels), 
      28 scans (32 pixels), 56 scans (16 pixels) of the visible channel
      (spectels 256 and above) correspond to an internal calibration

      for cubes with names NNNN_0 (first cube in a sequence), the first
      6 scans (128 x 4), 12 scans (128 x 2), 24 scans (128 x 1)
      48 scans (64 pixels), 96 scans (32 pixels) or 192 scans (16 pixels)
      of the IR channel (spectels 0 to 255) correspond to an internal 
      calibration (closed shutter, lamp on at 6 different levels, 
      in order 0,4,3,2,1,0

specmars(352): floating point, solar spectrum at the distance of Mars
      at the time of observations. jdat(i,*,n)/specmars(*) provides 
      the I/F values for the spectral range dominated by reflected 
      sunlight. It is meaningless for the spectral range dominated
      by thermal emission (usually above 3.5 µm)

sdat0(352,nscan): 
       0:256 contains the corresponding IR dark
       raw data: sdat0(0:255,n)-idat(i,0:255,n)
       this can be useful when there are problems
       with the IR (digital saturation of the dark, detector saturation) 

       detector saturation in the IR (only observed in the C channel
       with 5 msec integration time) can be checked by plotting:
       sdat0(0:127,n)-idat(i,0:127,n). If there is a plateau at values
       close to 300 DN around spectel 41 (the highest DN level), these
       values are saturated.

       Perturbations of the dark are most prominent with 16 pixels, 
       5 msec integration time at high DN levels from channels 20 to 50.
       A simple way to check for such effects is to plot sdat0(41,*)  
       (dark level for spectel 41, the highest DN level
       in the signal). If it is not nearly constant, and there is a dip 
       when the signal is high, the dark is perturbed. A non perturbed
       dark needs then to be recovered from a region of the cube with
       lower signal values (or another cube from the same orbit). 
       The proper level of the data can then be recovered by:
 
       idat(i,0:127,n)=idat -sdat0 (0:127,n) + sdat0 (0:127, unperturbed)


       With 2.5 msec modes, improper adjustment of the polarization 
       parameter can set part of the L channel dark current (orbits 12xx)
       (spectels 128 to 255) beyond the 4095 quantization limit.
       the bottom of the 3 µm band is perturbed
       This can be checked by plotting sdat0(128:255,n)
       so as to verify that it does not reach 4095 except for 
       channel 159 (dead).
       This can be recovered by searching for a non saturated dark
       in other cubes, then scaling it (0 level = 4290) and shifting it
       so as to match the part of the dark below 4095
      

sdat1(npix,7,nscan):
       housekeeping data as defined in the EAICD

exposure (float(3)): integration times for the IR (C channel), 
       the IR (L channel) and the VIS channel. The integration time
       is the same for IR C and IR L in all foreseen OMEGA modes.
       It can be either 2.5 msec or 5 msec. the visible integration
       time is in general set as 50 msec.      
 *** NEW ***
       the exposure time is set to 0. for channels which are not
       switched on. A VIS only observation will therefore 
       corresponds to 0., 0. and e.g. 50 msec.
    
       there are two photometric functions (352 values per files):
          mtfYYMM_25.dat (2.5 msec, release of year YY, month MM)
          mtfYYMM_50.dat ( 50 msec, release of year YY, month MM) 

       The readomega procedure selects the proper exposure time 
       for the IR and for the visible. In the present release, it uses the 
       mtf050110_25.dat nad mtf050110_50.dat files. It also takes into account
       the summation by two internal to the VIS channel for 128 pixel modes

IMPORTANT REMARK: new releases of the photometric function are likely to
       be necessary as the knowledge of the intrument improves.
       with this software release, relative variations of 5% or more
       between spectra are considered reliable. Smaller variations can
       be related to linearity effects, and should be considered with
       caution, in particular if they correlated with I/F (linearity). 
       The confidence on the absolute radiometric calibration is 
       at the 15% level. Three sets of photometric functions (040319, 
       050110 and 050701) are provided, which indicate the importance
       of updates. The 050701 set is used when computing jdat

bits_per_pixel : number of bits per pixel from data compression

summation : number of downtrack summations (only for modes with 128 pixels)
       this parameter can be 1.0, 2.0 or 4.0, in which case both idat
       and sdat0 have values can be a factor of 2 or 4 higher than
       in the nominal situation (summation 1.0). This factor is taken
       care of when computing jdat
      
The corresponding auxiliary pds cube ORBNNNN_X.NAV is read by readomega
if available (otherwise the user gets a "no corresponding geometry file"
message

products:

geocube: geometry cube with as many samples and lines as idat,
         and 51 bands. The definitions of each band are:
   0:   scan value (in DN) for each pixel
   1:   information on the time (hr:min:sec:msec)
   2:   incidence on the ellipsoid  (unit: 0.0001°)
   3:   emergence from the ellipsoid (unit: 0.0001°)
   4:   incidence wrt the local gravity field (atmospheric science)
   5:   emergence wrt the local gravity field
   6,7 : longitude and latitude of the center of the IR C pixel (0.0001°)
   8,9,10: incidence, emergence, phase for the IR C pixel (unit: 0.0001°)
   11 : distance in m
   12 : altitude (MOLA - ellipsoid) in m
        if > 65.535 km: altitude above the surface + 65536 m
        (for limb observations)
   13-16: longitudes of the corners of the IFOV (unit: 0.0001°)
   17-20: latitudes of the corners of the IFOV (unit: 0.0001°)
   21-35: same as 6-20 for the IR L pixel
   36-50: same as 6-20 for the Visible pixel

More details are given in the EAICD

dmars : sun - mars distance, in AU
specmars (352 values) : solar spectrum at OMEGA wavelengths at
        the heliocentric distance of Mars during observations
longi : longitude in ° of all OMEGA pixels (C channel)
lati : latitude in ° of all OMEGA pixels (C channel)
alt : MOLA altitude, in km, for all pixels in the cube (C channel)

IMPORTANT NOTE
   The geocube data set is very useful to obtain a first estimate
   of the region imaged by each OMEGA cube. Given the remaining uncertainties
   on position and attitude restitution, the observed variations in
   alignment of the C and L channel (from 1 to 3 pixels), and the 
   registration problems between spectels, the positionning error is 
   estimated as 3 mrad or better for most cubes. A specific effort using
   other data sets (in particular comparing the CO2 band strength with
   MOLA altimetry) is required if one requires a better positioning.

        another very useful check at medium to high incidences is 
        provided by the comparison between the cosine of the local 
        incidence (geocube(*,8,*) and the actual signal in the 
        continuum (e.g. spectel 11). This check of shadowing is 
        also useful so as to search for clouds (in particular in 
        Hellas) which destroy the otherwise good correlation
        between the expected illumination (assuming no rayleigh
        scattering) and actual signal 

        UPDATE ON INSTRUMENTAL PROBLEMS

- time evolution of hot and bad IR spectels
     since ground calibration, IR spectel 158 has been dead "cold" (0 value)
     and spectel 78 has been dead "hot" (always high values).
     since orbit 171, spectel 34 is also dead "hot".
     spectels 78, 171, and 34 from orbit 171 on are replaced by averages 
     of their two neighbours so as not to mess up spectral plots

     Cosmic ray degradation slowly increases the number of spectels
     which are very hot (500 DN or more above the dark), but not dead. 
     - 69, 88, 224 since the beginning
     - 188 since orbit 1147
     - 155 since orbit 1990
     the following spectels are moderately hot (100 to 500 DN above 
     the nominal dark)
     at orbit 2000:
     - 55, 66, 79, 85, 121, 127, 200, 222
     for any given cube, dead and hot pixels can be observed by
     plotting the dark: sdat0(0:255,n). They appear as spikes
     (upward for 158, downward for the others)

     very hot pixels are not reliable

     hot pixels have to be considered with caution, as they can 
     result in spikes in the decalibrated spectra. It is possible
     to recover most of the information for these spectral values
     with appropriate filtering tools or by using spectral ratios
     
     readomega now automatically handles hot and bad spectels as a
     function of orbit number. The first orbit after which a problem
     is observed is provided by the "bound060110.dat" file
     which generates the "bound" array (256 values, 1 per spectel)  
     orbit 0 means that the problem was observed during cruise,
     orbit 9999 means that the spectel has not (yet !) evolved.
     the photometric function (mtf) is multiplied by a factor
     1.e30 for spectels which have become unreliable. Otherwise,
     it is multiplied by a factor which is in general closer to 1
     for the 2.5 msec integration time than for 5 msec.
     An integer array (ic) is generated by readomega. It indicates
     which indices correspond to usable spectels on that date. 
     "plot,wvl(ic),jdat(pix,ic,line)" will therefore get rid of
     the spurious values due to unreliable spectels.

- time evolution of the L channel
     The L channel is not fully stable in the cross scan direction.
     The interval between the C and L channels varies between
     1.2 mrad and 3.1 mrad. There are mainly two states, which are
     correlated with the level of the calibration source in the L channel.
     when it is high, the separation is minimal. When it is low, 
     the separation is large. 

     There are (up to now) three series of orbit with a nominally high
     level of the calibration source:
      - orbits 98 to 511
      - orbits 939 to 1221
      - from orbit 2000 to 2600, the level slowly increases
        to reach a value only slightly lower than nominal

     A careful comparison of signals on bright terrains reveals a 
     likely impact on the L photometric function. 

     the new version of readomega provided by SOFT06 makes possible
     comparisons between spectra obtained at different levels of 
     the calibration source, as the photometric efficiency is modified
     accordingly.

- numerical perturbation problem for modes with 128 pixels
     Since orbit 511, pixels 80 to 95 are perturbed every 32 spectels
     across the 3 spectral ranges. For even lines, these spectels are:
     [12,13,14,15], [44,45,46,47] .....
     For odd lines, these spectels are:
     [28,29,30,31], [60,61,62,63] ....

     This can be (to some extent) corrected by replacing these 
     four values by the mean of the corresponding spectels on the 
     previous and next line, which have a different parity

         
