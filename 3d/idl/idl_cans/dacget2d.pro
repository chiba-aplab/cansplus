pro dacget2d,files,data $
   ,ix=ix,jx=jx,margin=margin

mfx=n_elements(files)
if n_elements(msode0) eq 0 then msode0=4
msode=msode0

if (mfx eq 1) then begin ; single file

file=files[0]

; check header
;
openr,unit,file,/get_lun

L300:
unity=0L
pos=msode
point_lun,unit,pos
readu,unit,unity
if (unity eq 0) then begin
  msode=8
  goto,L300
endif
if (unity gt 1000) then endian_flag=1 else endian_flag=0
; if (endian_flag = 1) then machine-arch and file-endian is different.

mver=0L
pos=pos+4+msode*2
point_lun,unit,pos
readu,unit,mver
if endian_flag then byteorder,mver,/ntohl

;   for data type: integer=4, real=5, double=6
mtype=0L
pos=pos+4+msode*2
point_lun,unit,pos
readu,unit,mtype
if endian_flag then byteorder,mtype,/ntohl

mndim=0L
pos=pos+4+msode*2
point_lun,unit,pos
readu,unit,mndim
if endian_flag then byteorder,mndim,/ntohl

mdim=lonarr(mndim)
pos=pos+4+msode*2
point_lun,unit,pos
readu,unit,mdim
if endian_flag then byteorder,mdim,/ntohl

ix=mdim[0] & jx=mdim[1]

case mtype of
  4: begin
    data=lonarr(ix,jx)
  end
  5: begin
    data=fltarr(ix,jx)
  end
  6: begin
    data=dblarr(ix,jx)
  end
endcase

headersize=12L*4+mndim*4+8
tmp=FSTAT(unit) & filesize=tmp.size

pos=headersize      ; position in file

   if (pos lt filesize) then begin
     point_lun,unit,pos+msode
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

dacget2d,file,data0,ix=ix,jx=jx

if (mf eq 0) then begin
  igx=ix
  jgx=jx*mfx-2*margin*(mfx-1)
  sz=size(data0)
  mtype=sz[sz[0]+1]+1
  case mtype of
    4: begin
      data=lonarr(ix,jgx)
    end
    5: begin
      data=fltarr(ix,jgx)
    end
    6: begin
      data=dblarr(ix,jgx)
    end
  endcase
endif

moff=mf*(jx-2*margin)+margin
data[*,moff:moff+jx-1-margin,*]=data0[*,margin:jx-1,*]

if (margin ge 1) then $
if (mf eq 0) then data[*,0:margin-1,*]=data0[*,0:margin-1,*]

endfor

ix=igx
jx=jgx


endelse



return
end
