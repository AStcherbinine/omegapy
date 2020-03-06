;SOFT09
;+ calculs fin + solarlongi
pro readomega_vpy, nomfic0, ldat, jdat, wvl, ic, specmars, latitude, longitude, emergence, incidence, altitude, solarlongi, ut_time, geocube

dirsoft = ''

openr,2,dirsoft+'omega_path'
datapath=''
geompath=''
readf,2,datapath
readf,2,geompath
close,2
;nomfic0=' '
;start:
;print,'OMEGA observation: '
;read,nomfic0

orbnum=(byte(nomfic0))(3)-48
if(orbnum gt 9) then orbnum=orbnum-7
orbnum=1000*orbnum+fix(strmid(nomfic0,4,3))
nomgeo=nomfic0+'.NAV'
nomfic0=nomfic0+'.QUB'
nomfic=datapath+nomfic0
openr,2,nomfic,ERROR=errflag
close,2
if (errflag ne 0) then begin
  print,'file ',nomfic,' not found'
  ;goto, start
  goto, fin
endif
readcube, nomfic, idat, sdat0, sdat1,info
if(size(idat))(0) ne 3 then begin
  print,'***** only one line in cube ',nomfic0,' ******'
  goto, fin
endif
exposure=0
exposure=info(0:2)
summation=info(3)
bits_per_data=info(4)
data_quality=fix(info(5))
nomgeo=geompath+nomgeo

openr,1,nomgeo,ERROR=errflag
if(errflag ne 0) then begin
  print,' no corresponding NAV cube'
  close,1
  dmars=1.52
  specmars=fltarr(352)
  openr,2,dirsoft+'specsol_0403.dat'
  readf,2,specmars
  close,2
  specmars=specmars/dmars/dmars     
  goto, nogeom
endif
trans= !pi/180.*1.e-4

geocube=0B
ilec=0B

read_geolbl,nomgeo, '', values, data

lrec   = data(0)
nrec   = data(1)
npixel = data(2)
npara  = data(3)
nscan  = data(4)
nau    = data(5)
solarlongi = data(8)
dmars=nau*1.e-4
close,2
openr,2,dirsoft+'specsol_0403.dat'
specmars=0B
specmars=fltarr(352)
readf,2,specmars
close,2
specmars=specmars/dmars/dmars

data = 0B
geocube=0B

geocube = lonarr(npixel, npara, nscan)
ilec = lonarr(npixel)

point_lun,1,lrec*nrec

for  k=0,nscan-1 do begin
   for j=0,npara-1 do begin
      readu, 1, ilec
      geocube(*,j,k) = ilec
   endfor
endfor
close, 1
if(geocube(1,1,0) gt 13) then geocube=swap_endian(temporary(geocube))
longi=reform(geocube(*,6,*)*1.e-4)
lati=reform(geocube(*,7,*)*1.e-4)
alt=reform(geocube(*,12,*)*1.e-3)
nogeom: 
close,/all
; preliminary pipeline tool

a=size(idat)
nbal=a(3)
npix=a(1)
jdat=0
jdat=float(idat)

fond2=intarr(256)
openr,2,dirsoft+'fond2.dat'
readf,2,fond2
close,2
if(exposure(1) gt 4.) then fond2=2*fond2
fondcur=fond2(128:255)#(0.*indgen(nbal)+1.)

jdat=0
jdat=float(idat)
i=where(idat le 1)
if(i(0) ne -1) then jdat(i)=1.e-5

pix0IR=0
i=where(idat(*,0:255,*) le 0)
if(i(0) ne -1) then pix0IR=(size(i))(3)
print,'       0 or less  IR: ',pix0IR

hkmin=min(sdat1(14,1,*))
if(hkmin lt 6) then hkmin=6
indj=where(sdat1(14,1,*) ge hkmin)
i6=indj(0)
indj=where(sdat1(14,1,*) ge hkmin+1)
balHK=indj(0)-i6
indi=where(sdat1(14,1,*) ge 6 and sdat0(10,*) gt 100)
balsm=balHK*8

if(balHK gt 0 and nbal gt indi(0)+balsm) then begin
  a=0.
  b=0.
  c=0.
  d=0.
  ndeb=indi(0)
  nf=(size(indi))(1)-1
  nfin=indi(nf)
  if(sdat0(26,nbal-1) lt 20) then sdat0(*,nbal-1)=sdat0(*,nbal-2)
  for k=0,255 do begin
    b=reform(float(sdat0(k,indi))) 
    a=[2*b(0)-rotate(b(0:balsm-1),2),b,2*b(nf)-rotate(b(nf-balsm+1:nf),2)]
    c=(smooth(a,balsm))(balsm:balsm+nf)
    d= -sdat0(k,ndeb:nbal-1)+spline(indi,c,ndeb+indgen(nbal-ndeb)) 
    for i=0,npix-1 do jdat(i,k,ndeb:nbal-1)= $
         jdat(i,k,ndeb:nbal-1)+d
  endfor
endif

for i=0,npix-1 do jdat(i,128:255,*)=jdat(i,128:255,*)-fondcur

i=where(jdat lt 1.e-5) 
if(i(0) ne -1) then jdat(i)=1.e-5

if(summation ne 1) then jdat=jdat/summation
linearC=0
linearC=fltarr(4096)
openr,2,dirsoft+'linearC.dat'
readf,2,linearC
close,2
jdat(*,0:127,*)=linearC(fix(jdat(*,0:127,*)+0.5))


wvl=0
wvl=fltarr(352)
openr,2,dirsoft+'lambda_0403.dat'
readf,2,wvl
close,2


mtf=fltarr(352)
rap=fltarr(3,256)
bound=intarr(3,256)
openr,2,dirsoft+'boundcur.dat'
readf,2,bound
close,2

if(exposure(1) lt 4.) then begin
  openr,2,dirsoft+'mtf120315_25.dat'
  readf,2,mtf
  close,2
  openr,2,dirsoft+'rapcur_25.dat'
  readf,2,rap
  close,2  
endif else begin
  openr,2,dirsoft+'mtf120315_50.dat'
  readf,2,mtf
  close,2
  openr,2,dirsoft+'rapcur_50.dat'
  readf,2,rap
  close,2  
endelse
ib=where(orbnum gt bound(0,*))
mtf(ib)=mtf(ib)*rap(0,ib)
ib=where(orbnum gt bound(1,*))
if(ib(0) ge 0) then mtf(ib)=mtf(ib)*rap(1,ib)
ib=where(orbnum gt bound(2,*))
if(ib(0) ge 0) then mtf(ib)=mtf(ib)*rap(2,ib)

; correction voie L
rapl=fltarr(128)
openr,2,'rapmtflcur.bin'
tmp=FSTAT(2)
orbmax=tmp.size/512-1
orbcur=orbnum
if(orbcur gt orbmax) then begin
  orbcur=orbmax
  print,orbcur,format='(" warning: L channel corrected as for orbit",I5)' 
endif
offs=long(512)*orbcur
point_LUN,2,offs
readu,2,rapl
close,2
mtf(128:255)=mtf(128:255)*rapl
if(exposure(1) gt 5.) then mtf(0:255)=mtf(0:255)/5.*exposure(1)

for n=0,255 do jdat(*,n,*)=jdat(*,n,*)/mtf(n)

ic=where(mtf lt 10000.)

;**************************************************************************************************************************
;bit error correction
;**************************************************************************************************************************

vis=float(idat(*,256:351,*))
siz = size(idat)
lines = siz(3)
pixels = siz(1)
exptime = info(2)/1000.
image_ratio = fltarr(pixels, lines)

level = 4095.*info(3)
if pixels EQ 128 then level = info(3)*4095*2.

i=where(vis le 0) 
counter_neg=0
if(i(0) ne -1) then begin
  counter_neg=(size(i))(3)
  vis(i)=0.001
endif

i=where(vis gt level)
counter_pos=0
if(i(0) ne -1) then begin 
  counter_pos=(size(i))(3)
  vis(i)=level
endif
i=where(vis le level and vis gt 0.8*level)
counter_sat=0
if(i(0) ne -1) then begin
  counter_sat=(size(i))(3)
  vis(i)=level*0.8
endif

plan3=reform(vis(*,3,*))
counter_spike3=0
i=where(plan3 gt 0.5*(vis(*,2,*)+vis(*,3,*))+50)
if(i(0) ne -1) then begin
  counter_spike3=(size(i))(3)
  plan3(i)=0.5*(reform(vis(*,2,*)+vis(*,3,*)))(i)
  vis(*,3,*)=plan3
endif

median_line = fltarr(lines)
column = fltarr(lines)

counter_spikes=long(0)

if(lines lt 10) then goto, nodespike

for k = 0, 95 do begin
    for m = 0, pixels - 1 do begin
       column(*) = vis(m, k, *)
       median_line = median(column, 7)
         for i = 0, lines-1 do begin
          if (abs(median_line(i) - vis(m, k, i)) GT 200) then begin
          vis(m, k, i) = median_line(i)
          counter_spikes = counter_spikes + 1
          endif
         endfor
    endfor
endfor
nodespike:

print, ' negative pixels VIS:', counter_neg
print, 'anomalous pixels VIS:', counter_pos
print, 'saturated pixels VIS:', counter_sat
print, '          spikes VIS:', counter_spikes

vis = float(fix(vis))

;******************************************************************************************************************************
; bias correction
;******************************************************************************************************************************



percen = 3.4e-6
fper = 0.
pend = 5.
if npixel eq 128 then fper = 1.
if npixel eq 64 then fper = 2.
if npixel eq 32 then fper = 4.
if npixel eq 16 then fper = 8.



strL1 =0.
strL2 = fltarr(96)
Strl_nc = 0.

lv = wvl(256:*)
lft1 = linfit([lv(0), lv(95)], [percen*fper, (percen*fper) +(percen*fper)/100.*pend ])
Strl_nc = (lv*lft1(1)+lft1(0))

for j = 0, lines-1 do begin 
      for k = 0, 96-1 do strL2(k) = total(vis(*, k, j))
      strL1 = total(vis(*, *, j))
      for i = 0, pixels -1 do begin
          vis(i, *, j) = (vis(i, *, j) - 50.*strL2*percen*fper)- 0.75*strL1*strl_nc;       
      endfor
endfor  

;******************************************************************************************************************************
; smear correction
;******************************************************************************************************************************

for i = 0, pixels-1 do begin
    for n = 0, lines-1 do begin
       vis(i, *, n) = smear_corr_050701( vis(i, *, n), info(2))

    endfor
endfor


;******************************************************************************************************************************
;flat cube
;******************************************************************************************************************************

slice = ulonarr(128, 96)

openr, 1, dirsoft+'flatVIS050701.bin'
readu, 1, slice
close, 1
if(slice(48,0) gt 40000) then slice=swap_endian(slice)


istart= (128 - pixels)/2.
iend  = istart + pixels - 1

for i = 0, lines-1 do begin

    vis(*, *, i) = vis(*, *, i) * 2.^15 / slice(istart:iend, *)

endfor

;***************************************************************************************************************************
;radiometric calibration
;***************************************************************************************************************************

f = 1.
if pixels eq 128 then  f = 2; internal summation
f2 = 1.
if exptime eq 0.1 then f2 = 2.
if exptime eq 0.2 then f2 = 4.

spectrum = fltarr(96)
 
for i = 0, lines - 1 do begin
    for m = 0, pixels - 1 do begin
       spectrum = vis(m, *, i)
       jdat(m, 256:351, i) = ordcorr(wvl(256:351), spectrum)/ ( f2 * mtf(256:351)) / f / summation
    endfor
    ;print, 'done %', 100*float(i)/(lines - 1)
endfor

;*************
;AJOUT MATHIEU
;*************

;useful quantities:
trans= !pi/180.*1.e-4
ecl=reform(cos(geocube(*,2,*)*trans))

longitude=reform(geocube(*,6,*)*1.e-4)
latitude=reform(geocube(*,7,*)*1.e-4)
altitude=reform(geocube(*,12,*)*1.e-3)
emergence=reform(geocube(*,3,*))*1.e-4
incidence=reform(geocube(*,2,*))*1.e-4
ut_time=reform(geocube(*,1,*))

;reflectance factor calculation
ldat=jdat
for n=0,351 do ldat(*,n,*)=ldat(*,n,*)/specmars(n)/ecl(*,*)
albedo=reform(ldat(*,11,*))
;atmospheric spectrum
atmorap=fltarr(256)  
openr,2,'atmorap.dat'
readf,2,atmorap
close,2


fin:
close,/all
end
