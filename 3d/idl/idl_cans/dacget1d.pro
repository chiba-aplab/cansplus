pro dacget1d,files,data $
   ,ix=ix,margin=margin

if n_elements(msode0) eq 0 then msode0=4
msode=msode0

mfx=n_elements(files)

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

ix=mdim[mndim-1]

case mtype of
  4: begin
    data=lonarr(ix)
  end
  5: begin
    data=fltarr(ix)
  end
  6: begin
    data=dblarr(ix)
  end
endcase

headersize=(4L+msode*2)*4+(4L*mndim+msode*2)
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

dacget1d,file,data0

if (mf eq 0) then begin
  sz=size(data0)
  ix=sz[1]
  igx=ix*mfx-2*margin*(mfx-1)
  mtype=sz[2]+1
  case mtype of
    4: begin
      data=lonarr(igx)
    end
    5: begin
      data=fltarr(igx)
    end
    6: begin
      data=dblarr(igx)
    end
  endcase
endif

moff=mf*(ix-2*margin)+margin
data[moff:moff+ix-1-margin]=data0[margin:ix-1]

if (margin ge 1) then $
  if (mf eq 0) then data[0:margin-1]=data0[0:margin-1]

endfor


endelse

ix=n_elements(data)


return
end
