function smear_corr_050701, spectrum, exposure_time,  smear

smear = fltarr(96)

spectrum1 = spectrum

; Correction for the CCD cleaning time
ft = 0.87
int = exposure_time - ft

;
bands = 96


; Check for saturation in the 90:95 channels
;
means = mean(spectrum1(0:84))
for i = 85, bands - 1 do begin
    if spectrum1(i) GT means then spectrum1(i) = 0
endfor

;;
; Smearing term
for i = 0, bands - 1 do begin
    smear(i) = total(spectrum1(i:bands-1)) * ft / int / bands
endfor


spectrum1 = spectrum1 - smear


return, spectrum1
end