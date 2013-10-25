pro dacget,filesdac,filesprm,data,narg=narg,flagsum=flagsum

if (n_elements(flagsum) ne 0) then flagsum0=flagsum else flagsum0=0

if (n_elements(filesdac) ne n_elements(filesprm)) then begin
  print,'Error: Number of files should be the same.'
  return
endif

if (product(file_test(filesdac))*product(file_test(filesprm)) eq 0) then begin
  print,'Error: Not all of specified files do exist'
  return
endif
mfx=n_elements(filesdac)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; check header
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

dacgetparam,filesprm[0],'mendian',mendian
dacgetparam,filesprm[0],'mtype',mtype
dacgetparam,filesprm[0],'mdimx',mdimx
mxar=intarr(mdimx)
for mdim=0,mdimx-1 do begin
  dacgetparam,filesprm[0],'m'+strtrim(mdim+1,2)+'x',mx0
  mxar[mdim]=mx0
endfor
marginar=intarr(mdimx)
for mdim=0,mdimx-1 do begin
  dacgetparam,filesprm[0],'margin'+strtrim(mdim+1,2),margin0
  marginar[mdim]=margin0
endfor

case mtype of
  4: mbyte=4
  5: mbyte=4
  6: mbyte=8
endcase

mmx=long(product(mxar,/integer))

mxar1=mxar
mxar2=mxar-2*marginar ; size of margin-cut array

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; get parallel information
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

mpear=intarr(mdimx,mfx)
mpexar=intarr(mdimx,mfx)
for mf=0,mfx-1 do begin
for mdim=0,mdimx-1 do begin
  dacgetparam,filesprm[mf],'mpe'+strtrim(mdim+1,2),mpe
  mpear[mdim,mf]=mpe
endfor
endfor

mpexar=intarr(mdimx,mfx)
mf=0
for mdim=0,mdimx-1 do begin
  dacgetparam,filesprm[mf],'mpe'+strtrim(mdim+1,2)+'x',mpex
  mpexar[mdim,mf]=mpex
endfor

mgxar2=mxar2*mpexar

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; only if "sequential output"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (mxar[mdimx-1] ne -1) then begin 
  if (n_elements(narg) ne 0) then $
    print,'Warning: file is for non-sequential output. narg is ignored.'
  goto, L100
endif

; count number of data sets
;
mmx=1L
if (mdimx ge 2) then mmx=long(product(mxar[0:mdimx-2],/integer))
datasize=mmx*mbyte

mf=0
openr,unit,filesdac[mf],/get_lun
fst=fstat(unit)
close,unit
free_lun,unit
filesize=fst.size
ndfx=filesize/datasize

; determine data set IDs to be read
;
if n_elements(narg) eq 0 then narg=-1
sz=size(narg)
allreadflag=0
narr=-1
if (sz[0] ge 1) then  begin
  whr=where(narg lt ndfx, mx)
  if (mx ge 1) then narr=narg[whr]
endif else begin
  if (narg ge 0) then begin
    if (narg lt ndfx) then narr=[narg]
  endif else begin
    allreadflag=1
    narr=indgen(ndfx)
  endelse
endelse
if (narr[0] eq -1) then begin
  print,'Error: No data of specified IDs'
  return
endif
ndx=n_elements(narr)

mxar1[mdimx-1]=ndx
mxar2[mdimx-1]=ndx
mgxar2[mdimx-1]=ndx

L100:

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; data reading
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

datag=make_array(mgxar2,type=mtype-1) ; for global data

for mf=0,mfx-1 do begin

  mxar20=mxar2

  ; read data
  ;
  if (mxar[mdimx-1] ne -1) then begin 
    ;;; non-sequential
    data=make_array(mmx,type=mtype-1)
    openr,unit,filesdac[mf],/get_lun
      readu,unit,data
    close,unit
    free_lun,unit
    data=reform(data,mxar)
  endif else begin 
    ;;;  sequential
    data0=make_array(mmx,type=mtype-1)
    data =make_array([mmx,ndx],type=mtype-1)
    openr,unit,filesdac[mf],/get_lun
    for nd=0,ndx-1 do begin
      n=narr[nd]
      pos=n*datasize
      point_lun,unit,pos
      readu,unit,data0
      data[*,nd]=data0
    endfor
    close,unit
    free_lun,unit
    data=reform(data,mxar1)
  endelse ; sequential or non-sequential

  ; cut margin
  ;
  execstring=''
  for mdim=0,mdimx-1 do begin
    mgn0=marginar[mdim]
    mgn1=mxar1[mdim]-marginar[mdim]-1
    addst=strtrim(mgn0,2)+':'+strtrim(mgn1,2)
    execstring=execstring+addst
    if (mdim lt mdimx-1) then execstring=execstring+','
  endfor
  execstring='data=data['+execstring+']'
  void=execute(execstring)
  if (mgn0 eq mgn1) then mxar20=[mxar20,1]
  data=reform(data,mxar20)

  ; endian change if necessary
  ;
  case mendian of
    1: begin  ; Data endian is Big
      data=swap_endian(data,/swap_if_little_endian)
      end
    0: begin  ; Data endian is Little
      data=swap_endian(data,/swap_if_big_endian)
      end
    else:
  endcase

  ; substitute to the grobal array
  ;
  tmpst=''
  for mdim=0,mdimx-1 do begin
    mgn0=mxar20[mdim]*mpear[mdim,mf]
    mgn1=mxar20[mdim]*(mpear[mdim,mf]+1)-1
    addst=strtrim(mgn0,2)+':'+strtrim(mgn1,2)
    tmpst=tmpst+addst
    if (mdim lt mdimx-1) then tmpst=tmpst+','
  endfor
  execstring='datag['+tmpst+']=data'
  if flagsum0 then execstring='datag['+tmpst+']=datag['+tmpst+']+data'
  void=execute(execstring)

endfor  ; mf

data=datag


;return
end
