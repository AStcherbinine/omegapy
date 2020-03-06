
; 20/01/04

function ordcorr, l, spt

; ll: VNIR wavelength in micron
; spt: spettrum to be corrected


;___________________________________________________________________________________
; assignment of the variable k(3,kn), from file Coeff_II_96_ias.txt
; k (0,*) = spectral channels with the second order contamination
; k (1,*) = wavelengths correspondent to each spectral channel
; k (2,*) = values of the coefficients for each channel
nk = 17; total number of the channel with the contamination
k = fltarr(3,nk)
k(0,*) = findgen(nk)+79
k(1,*) = [952.100, 959.600, 967.000, 974.300, 981.700, 988.900, 996.300 ,$
		  1003.70, 1010.90, 1018.20, 1025.50, 1032.90, 1040.30, 1047.94, $
		  1055.34, 1062.52, 1069.87]/1000.

k(2,*) = [0.00478825, 0.00764561, 0.0117481, 0.0186357, 0.0291402, 0.0430787, $
          0.0561666,  0.0694787,  0.0831372,  0.0993437, 0.113755, 0.124017, 0.123925, $
           0.118724, 0.110151, 0.106578, 0.106578]
;__________________________________________________________________________________________
; MAIN
; Formula : I_corr = I_raw-kII*I_1_raw
; I_corr = corrected final Intensity
; I_raw = raw spectrum measured by VNIR (initial intensity)
; KII = coefficients for the second order
; I_1_raw = raw intensity responsable of the second order contribution

I_raw = fltarr(96)
I_1_raw = fltarr(96)
kII = fltarr(96)
I_corr= fltarr(96)

;_____________________________________
; Assignment of I_raw
I_raw = spt

;_________________________________
; Assignment of I_1_raw
ch = k(0,*); spectral channel concerning the II ord
l2 = k(1,*); lambda of the II ord
l1 = l2/2. ; lambda od first order
nl = n_elements(l1)
ch_start = 16; first spectral channel of the first order wavelength
ch_end   = 23; last spectral channel of the first order wavelength

I_1= spline(l(ch_start:ch_end), I_raw(ch_start:ch_end), l1)

I_1_raw(k(0,0):k(0,nl-1)) = I_1

;_____________________________________
; Assignment of kII
kII(k(0,0):k(0,nl-1)) = k(2,*)

;_____________________________________
; Computation of I_corr

I_corr = I_raw-kII*I_1_raw

return, I_corr

end
