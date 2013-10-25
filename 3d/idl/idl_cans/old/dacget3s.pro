pro dacget3s,files,data,narg=narg $
   ,ix=ix,jx=jx,kx=kx,nx=nx,narr=narr,ndx=ndx,margin=margin

if n_elements(narg) eq 0 then narg=-1
osb=!version.memory_bits/8 ; offset_byte

sz=size(narg)
allreadflag=0
if (sz[0] ge 1) then  begin
  narr=narg
endif else begin
  if (narg ge 0) then begin
    narr=[narg]
  endif else begin
    allreadflag=1
  endelse
endelse

mfx=n_elements(files)

if (mfx eq 1) then begin

file=files[0]

; check header
;
openr,unit,file,/get_lun

unity=0L
pos=long(osb)
point_lun,unit,pos
readu,unit,unity
if (unity gt 1000) then endian_flag=1 else endian_flag=0
; if (endian_flag = 1) then machine-arch and file-endian is different.

mver=0L
pos=pos+4L+osb*2
point_lun,unit,pos
readu,unit,mver
if endian_flag then byteorder,mver,/ntohl

;   for data type: integer=4, real=5, double=6
mtype=0L
pos=pos+4L+osb*2
point_lun,unit,pos
readu,unit,mtype
if endian_flag then byteorder,mtype,/ntohl

mndim=0L
pos=pos+4L+osb*2
point_lun,unit,pos
readu,unit,mndim
if endian_flag then byteorder,mndim,/ntohl

mdim=lonarr(mndim)
pos=pos+4L+osb*2
point_lun,unit,pos
readu,unit,mdim
if endian_flag then byteorder,mdim,/ntohl

ix=mdim[0] & jx=mdim[1] & kx=mdim[2]
if (mdim[mndim-1] eq -1) then begin
headersize=(4L+osb*2)*4+(mndim*4L+osb*2)
tmp=FSTAT(unit) & filesize=tmp.size
datasize=filesize-headersize
if (mtype eq 6) then dataunitsize=8L*ix*jx*kx else dataunitsize=4L*ix*jx*kx
dataperstep=dataunitsize+osb*2
ndx=datasize/dataperstep
endif else begin
  ndx=mdim[mndim-1]
endelse

if (allreadflag eq 1) then narr=indgen(ndx)

sz=size(narr)
ndmax=sz[1]

case mtype of
  4: begin
    data0=lonarr(ix,jx,kx)
    data=lonarr(ix,jx,kx,ndmax)
  end
  5: begin
    data0=fltarr(ix,jx,kx)
    data=fltarr(ix,jx,kx,ndmax)
  end
  6: begin
    data0=dblarr(ix,jx,kx)
    data=dblarr(ix,jx,kx,ndmax)
  end
endcase

for nd=0,ndmax-1 do begin
   nn=narr[nd]
   pos=headersize+nn*dataperstep      ; position in file

   if (pos lt filesize) then begin
     point_lun,unit,pos+osb*2
     readu,unit,data0
     data[*,*,*,nd]=data0
   endif else begin
     goto,label0
   endelse

endfor
label0:

ndx=nd
data=data[*,*,*,0:ndx-1]
nx=ndx

close,unit
free_lun,unit

if endian_flag then begin
case mtype of
  4: byteorder,data,/ntohl
  5: byteorder,data,/xdrtof
  6: byteorder,data,/xdrtod
endcase
endif

endif else begin

for mf=0,mfx-1 do begin

file=files[mf]

dacget3s,file,data0,narg=narg $
   ,ix=ix,jx=jx,kx=kx,nx=nx,narr=narr,ndx=ndx

if (mf eq 0) then begin
  igx=ix
  jgx=jx
  kgx=kx*mfx-2*margin*(mfx-1)
  sz=size(data0)
  mtype=sz[sz[0]+1]+1
  case mtype of
    4: begin
      data=lonarr(ix,jx,kgx,ndx)
    end
    5: begin
      data=fltarr(ix,jx,kgx,ndx)
    end
    6: begin
      data=dblarr(ix,jx,kgx,ndx)
    end
  endcase
endif

moff=mf*(kx-2*margin)+margin
data[*,*,moff:moff+kx-1-margin,*]=data0[*,*,margin:kx-1,*]

if (margin ge 1) then $
if (mf eq 0) then data[*,*,0:margin-1,*]=data0[*,*,0:margin-1,*]

endfor

ix=igx
jx=jgx
kx=kgx

endelse



return
end
