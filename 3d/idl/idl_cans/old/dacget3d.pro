pro dacget3d,files,data $
   ,ix=ix,jx=jx,kx=kx,margin=margin

mfx=n_elements(files)

if (mfx eq 1) then begin ; single file

osb=!version.memory_bits/8 ; offset_byte

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

case mtype of
  4: begin
    data=lonarr(ix,jx,kx)
  end
  5: begin
    data=fltarr(ix,jx,kx)
  end
  6: begin
    data=dblarr(ix,jx,kx)
  end
endcase

headersize=(4L+osb*2)*4+(mndim*4L+osb*2)
tmp=FSTAT(unit) & filesize=tmp.size

pos=headersize      ; position in file

   if (pos lt filesize) then begin
     point_lun,unit,pos+osb
     readu,unit,data
   endif else begin
     goto,label0
   endelse

label0:

close,unit
free_lun,unit

if endian_flag then begin
case mtype of
  4: byteorder,data,/ntohl
  5: byteorder,data,/xdrtof
  6: byteorder,data,/xdrtod
endcase
endif

endif else begin ; multi files

for mf=0,mfx-1 do begin

file=files[mf]

dacget3d,file,data0,ix=ix,jx=jx,kx=kx

if (mf eq 0) then begin
  igx=ix
  jgx=jx
  kgx=kx*mfx-2*margin*(mfx-1)
  sz=size(data0)
  mtype=sz[sz[0]+1]+1
  case mtype of
    4: begin
      data=lonarr(ix,jgx,kgx)
    end
    5: begin
      data=fltarr(ix,jgx,kgx)
    end
    6: begin
      data=dblarr(ix,jgx,kgx)
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
