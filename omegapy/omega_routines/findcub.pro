;Version YL 25 mai 2011
;mathieu : ajout erreur cube .PRE non présent (omega_path -> dataom/
;mathieu : ajout lecture fichier findcube.pro
;mathieu : ajout affichage Ls

pro findcub,long0,lat0,cdeb,cfin

openr,2,'omega_path'
path=''
readf,2,path & readf,2,path
close,2

trans= !pi/180.
x0=cos(long0*trans)*cos(lat0*trans)
y0=sin(long0*trans)*cos(lat0*trans)
z0=sin(lat0*trans)
nomc=strarr(300)

long1=fltarr(10)
lat1=fltarr(10)
x1=long1
y1=lat1
nomcube=' '
close,/all


nomout='orbites_lg'+string(long0,format='(I0)')+'lt'+string(lat0,format='(I0)')+'.dat'
openw,4,nomout
openw,3,'cubliste'
printf,3,long0,lat0,format='(%" long: %7.3f  lat: %7.3f\n")'
openr,2,'cubindex.ref'
nhits=0
geocube=0

if(lat0 lt -60) then goto, sud
if(lat0 gt 60.) then goto, nord
for ncube=0,10000 do begin
  readf,2,nomcube
  norb=fix(strmid(nomcube,3,4))
  if(norb eq 0) then goto, fin
  readf,2,long1
  readf,2,lat1
  if(norb lt cdeb) then continue
  if(norb gt cfin) then goto,fin
  x1=cos(long1*trans)*cos(lat1*trans)
  y1=sin(long1*trans)*cos(lat1*trans)
  z1=sin(lat1*trans)
  ps=x0*x1+y0*y1+z0*z1
  if(max(ps) lt 0.85) then continue
  long2=long1-long0
  i=where(long2 lt -180.) 
  if(i(0) ne -1) then long2(i)=long2(i)+360.
  i=where(long2 gt 180.)
  if(i(0) ne -1) then long2(i)=long2(i)-360.
  if(testin(0,lat0,long2,lat1) eq 1) then begin
    nomc(nhits)=nomcube
    nhits=nhits+1
  endif
endfor
goto, fin
nord:
for ncube=0,10000 do begin
  readf,2,nomcube
  norb=fix(strmid(nomcube,3,4))
  if(norb eq 0) then goto, fin
  readf,2,long1
  readf,2,lat1
  if(norb lt cdeb) then continue
  if(norb gt cfin) then goto,fin
  if(max(lat1) lt 60.) then continue
  x1=cos(long1*trans)*cos(lat1*trans)
  y1=sin(long1*trans)*cos(lat1*trans)
  if(testin(x0,y0,x1,y1) eq 1) then begin
    nomc(nhits)=nomcube
    nhits=nhits+1
  endif
endfor
goto, fin
sud:
for ncube=0,10000 do begin
  readf,2,nomcube
  norb=fix(strmid(nomcube,3,4))
  if(norb eq 0) then goto, fin
  readf,2,long1
  readf,2,lat1
  if(norb lt cdeb) then continue
  if(norb gt cfin) then goto,fin
  if(min(lat1) gt -60.) then continue
  x1=cos(long1*trans)*cos(lat1*trans)
  y1=sin(long1*trans)*cos(lat1*trans)
  if(testin(x0,y0,x1,y1) eq 1) then begin
    nomc(nhits)=nomcube
    nhits=nhits+1
  endif
endfor
fin:
close,2
print,'   orbit    x    y    dmin   altMEx   inci   emer   phas    Ls'
for n=0,nhits-1 do begin
  ;testfile='dataom/'+nomc(n)+'.NAV'
  testfile=path+nomc(n)+'.NAV'
  a=FILE_SEARCH(testfile)
  ;if(a eq '') then testfile='dataom/'+nomc(n)+'.PRE'
  if(a eq '') then testfile=path+nomc(n)+'.PRE'
  ;a=FILE_SEARCH(testfile)  ;ajout mathieu
  ;if(a eq '') then begin
  ;    print,'donnees manquantes : ',testfile
  ;    goto,anomalie
  ;endif
  data=0B
  values=0B
  close,1
  read_geolbl,testfile,'',values,data
  lrec=fix(data(0))
  nrec=fix(data(1))
  npix=fix(data(2))
  npara=fix(data(3))
  nscan=fix(data(4))
  solong=float(data(8))
  geocube=lonarr(npix,npara,nscan)
  close,1
  openr,1,testfile
  point_lun,1,lrec*nrec
  readu,1,geocube
  close,1
  if(geocube(1,1,0) gt 13) then geocube=swap_endian(geocube)
  longa=geocube(*,6,*)*1.e-4
  lata=geocube(*,7,*)*1.e-4
  xa=cos(longa*trans)*cos(lata*trans)
  ya=sin(longa*trans)*cos(lata*trans)
  za=sin(lata*trans)  
  xr=y0*za-z0*ya
  yr=z0*xa-x0*za
  zr=x0*ya-y0*xa
  dist=sqrt(xr*xr+yr*yr+zr*zr)*3393.
  distmin=min(dist,ij)
  j0=long(ij/npix)
  i0=ij-j0*npix
  slant=geocube(i0,11,j0)*1.e-3
  inci=geocube(i0,2,j0)*1.e-4
  emer=geocube(i0,3,j0)*1.e-4
  phas=geocube(i0,10,j0)*1.e-4
print,nomc(n),i0,j0,distmin,slant,inci,emer,phas,solong,format='(A9,i5,i5,f8.2,f9.1,3f7.2,3f7.2)' 
printf,3,nomc(n),i0,j0,distmin,slant,inci,emer,phas,solong,format='(A9,i5,i5,f8.2,f9.1,3f7.2,3f7.2)'
printf,4,nomc(n)
anomalie:

endfor
;printf,3,'000000000'
close,3
close,4
end
