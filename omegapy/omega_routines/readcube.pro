pro readcube, nomfic, idat, sdat0, sdat1, info

close,3
close,4

idat=0
openr,3,nomfic,ERROR=err
line=' '
line1=' '
line2=' '
line3=' '
kax0=1
kax1=1
kax2=1
sax0=0
sax1=0
sax2=0
val=0
cbyte=2
axex=0
axey=1
axek=2
info=0
info=fltarr(6)
itype=1 ; SUN_INTEGER
orbnum=0

deb: 
  readf,3,line
  len=strlen(line)
  if strtrim(line,2) eq 'END' then goto, fin
  i=strpos(line,' = ')
  if i gt -1 then begin
    line1=strtrim(strmid(line,0,i),2)
    line2=strmid(line,i+3,len-i-2)
    if strmid(line1,8,4) eq 'DESC' then begin
      tst=byte(strmid(line2,34,1))-48
      if(tst gt 9) then tst=tst-7
      flagc=intarr(3)
      flagc(*)=0
      if(tst / 8) then flagc(2)=1
      if(tst /4 AND 1) then flagc(1)=1
      if(tst /2 AND 1) then flagc(0)=1
      goto,deb      
    endif
    if line1 eq 'RECORD_BYTES' then begin
      lrec=fix(line2)
      goto, deb
    endif
    if line1 eq 'DATA_QUALITY_ID' then begin
      dqual=fix(line2)
      info(5)=dqual*1.+0.001
      goto,deb
    endif
    if line1 eq 'LABEL_RECORDS' then begin 
      nrec=fix(line2)
      goto, deb
    endif
    if line1 eq '^IMAGE' then begin
      nrec=fix(line2)
      nax=2
      goto, deb
    endif
    if line1 eq 'EXPOSURE_DURATION' then begin
      line2=strtrim(line2,1)
      info(0)=float(strmid(line2,1,3))
      info(1)=float(strmid(line2,5,3))
      info(2)=float(strmid(line2,9,4))      
      info(0:2)=info(0:2)*flagc
      goto, deb
    endif
    if line1 eq 'DOWNTRACK_SUMMING' then begin
      info(3)=fix(line2)
      goto, deb
    endif
    if line1 eq 'INST_CMPRS_RATE' then begin
      info(4)=float(line2)
      goto, deb
    endif
    if line1 eq 'AXES' then nax=fix(line2)
    if line1 eq 'AXES_NAME' then begin
      if nax eq 2 then begin
	axex=0 & axek=1
	if line2 eq '(BAND,SAMPLE)' then begin
	  axex=1 & axek=0
	endif
      endif
      if nax eq 3 then begin
        axex=0 & axey=2 & axek=1
        if line2 eq '(BAND,SAMPLE,LINE)' then begin
          axex=1 & axey=2 & axek=0
        endif
        if line2 eq '(SAMPLE,LINE,BAND)' then begin
	  axex=0 & axey=1 & axek=2
        endif
      endif
      goto, deb
    endif
    if line1 eq 'LINES' then kax1=fix(line2)
    if line1 eq 'LINE_SAMPLES' then kax0=fix(line2)
    if line1 eq 'SAMPLE_TYPE' then begin
      line2=strtrim(line2,2)
      if line2 eq 'VAX_INTEGER' then itype=0
      goto, deb
    endif
    if line1 eq 'ORBIT_NUMBER' then begin
      orbnum=fix(line2)
      goto, deb
    endif
    if line1 eq 'CORE_ITEMS' then begin
      i=strpos(line2,',')
      kax0=fix(strmid(line2,1,i-1))
      line2=strmid(line2,i+1,strlen(line2)-i-1)
      i=strpos(line2,',')
      kax1=fix(strmid(line2,0,i))
      line2=strmid(line2,i+1,strlen(line2)-i-1)
      if nax gt 2 then begin
	i=strpos(line2,')')
        kax2=fix(strmid(line2,0,i))
      endif
      goto, deb
    endif
    if line1 eq 'CORE_ITEM_BYTES' then cbyte=fix(line2)	
    if line1 eq 'SUFFIX_ITEMS' then begin
      i=strpos(line2,',')
      sax0=fix(strmid(line2,1,i-1))
      line2=strmid(line2,i+1,strlen(line2)-i-1)
      i=strpos(line2,',')
      sax1=fix(strmid(line2,0,i))
      line2=strmid(line2,i+1,strlen(line2)-i-1)
      if nax gt 2 then begin
	i=strpos(line2,')')
        sax2=fix(strmid(line2,0,i))
      endif
    endif
    if line1 eq 'SUFFIX_BYTES' then sbyte=fix(line2)
  endif
goto,deb
fin:
print,' core:  ',kax0,kax1,kax2,' cbyte:',cbyte 
print,' suffix:',sax0,sax1,sax2,' sbyte:',sbyte

close,3
close,4
ilec=0
slec=0
idat=0
sdat0=0
sdat1=0

;ilec=intarr(kax0)
;idat=intarr(kax0,kax1,kax2)

;sdat0=lonarr(sax0,kax1,kax2)
;sdat1=lonarr(kax0,sax1,kax2)
;slec0=long(0)
;slec1=lonarr(kax0)

; New read, using assoc (S.E.)
        C_line = Make_array(kax0, Type = 2, /nozero) ; short int 
        F_line = {C_line: C_line, S_line: 0L}

        If sax1 NE 0 then begin
          SS_frame = reform(Make_array(kax0, sax1, Type = 3), kax0, sax1) ; long int
          F_frame = {F_line:replicate(F_line,kax1), SS_frame: SS_frame}
        endif else F_frame = {F_line:replicate(F_line,kax1)}
        SS_frame = 0B
        F_line = 0B
        F_Qube= replicate(F_frame,kax2)
        F_frame = 0B
       openr, 3, nomfic
       skip = lrec*nrec
       file = assoc(3, F_Qube, skip)
       element = file(0)
       close, 3
       file = 0B
       idat = element.F_line.C_line
       sdat0 = element.F_line.S_line
       If sax1 NE 0 then sdat1 = reform(element.SS_frame, kax0, sax1,kax2)
       element = 0B

if(sdat1(1,1,0) gt 255) then begin
  idat=swap_endian(temporary(idat))
  sdat0=swap_endian(temporary(sdat0))
  sdat1=swap_endian(temporary(sdat1))
endif

end
