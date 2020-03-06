;+
; NAME:
;     READ_GEOLBL
;
; PURPOSE:
;     Read the label of a geometrical data cube. In particular, you
;     can retrieve any parameters from it by specifying the keywords
;     you wish to get the values. Any call return the label size and
;     the cube dimension.
;
; CALLING SEQUENCE:
;     read_geolbl, filename, keywords, values, info
;
; INPUTS:
;    filename: geometrical data file name
;    keywords: string array containing a list of keyword

; KEYWORD PARAMETERS:
;    none
;
; OUTPUTS:
;    values: has the same dimension as the keywords argument.
;            It contains, if exist, the values of the requested
;            keywords.
;
;    info: integer array
;        info(0) = size of one label record
;        info(1) = number of label record
;        info(2) = number of pixel
;        info(3) = number of parameter
;        info(4) = number of scan
;
; EXAMPLE:
;     IDL>read_geolbl, 'TEST040519.PRE', ['PRODUCT_CREATION_TIME', 'AXIS_NAME'], values, info
;     IDL>print, values
;     2003-10-09T17:08:48.00 (SAMPLE,GEOMETRY,LINE)
;     IDL>info = long(info)
;     IDL>print, 'File size = ', ((info(0)*info(1) + info(2)*info(3)*info(4)*4 ) / 1024.), 'Ko'
;     File size =       16265.0Ko
;
; LIMITATIONS:
;
; MODIFICATION HISTORY:
;     October 2003, Nicolas Manaud.
;-


PRO read_geolbl, filename, keywords, values, info

line    = ' '
keyword = ' '
value   = ' '

npixel = 1
npara  = 1
nscan  = 1

openr, 100, filename

nb_item = n_elements(keywords)
values = strarr(nb_item)

max_line = 300
nau=15000

for n=0,max_line-1 do begin
  readf, 100, line
  len=strlen(line)
  if  strcompress(line, /REMOVE_ALL) eq 'END' then goto, quit
  i = strpos(line,' = ')
  if ( i gt -1 ) then begin
    keyword = strtrim(strmid(line,0,i),2)
    value   = strmid(line,i+3,len-i-2)

    for idx=0,nb_item-1 do if ( keywords(idx) eq keyword ) then values(idx) = value

    if keyword eq 'RECORD_BYTES'  then lrec=fix(value)
    if keyword eq 'LABEL_RECORDS' then nrec=fix(value)
    if keyword eq 'SOLAR_DISTANCE' then nau=fix(float(value)/14960.000)
    if keyword eq 'SUB_SOLAR_LONGITUDE' then slong=float(value)
    if keyword eq 'SUB_SOLAR_LATITUDE' then slat=float(value)
    if keyword eq 'SOLAR_LONGITUDE' then solarlong=float(value)
    
    if keyword eq 'CORE_ITEMS' then begin
      i      = strpos(value,',')
      npixel = fix(strmid(value,1,i-1))
      value  = strmid(value,i+1,strlen(value)-i-1)
      i      = strpos(value,',')
      npara  = fix(strmid(value,0,i))
      value  = strmid(value,i+1,strlen(value)-i-1)
      i      = strpos(value,')')
      nscan  = fix(strmid(value,0,i))
    endif
  endif
endfor
if(nau eq 15000) then print,'no solar distance -> 1.5 AU'


quit:

info = [ lrec, nrec, npixel, npara, nscan, nau,slong,slat,solarlong]

close,100

END
