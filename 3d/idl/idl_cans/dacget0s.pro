pro dacget0s,file,data,narg=narg $
   ,nx=nx,narr=narr,ndx=ndx

if n_elements(narg) eq 0 then narg=-1
if n_elements(msode0) eq 0 then msode0=4

msode=msode0

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

if (mdim[mndim-1] eq -1) then begin
headersize=(4L+msode*2)*4+(4L*mndim+msode*2)
tmp=FSTAT(unit) & filesize=tmp.size
datasize=filesize-headersize
if (mtype eq 6) then dataunitsize=8L else dataunitsize=4L
dataperstep=dataunitsize+msode*2
ndx=datasize/dataperstep
endif else begin
  ndx=mdim[mndim-1]
endelse

if (allreadflag eq 1) then narr=indgen(ndx)

sz=size(narr)
ndmax=sz[1]

case mtype of
  4: begin
    data0=0L
    data=lonarr(ndmax)
  end
  5: begin
    data0=0.e0
    data=fltarr(ndmax)
  end
  6: begin
    data0=0.d0
    data=dblarr(ndmax)
  end
endcase

for nd=0,ndmax-1 do begin
   nn=narr[nd]
   pos=headersize+nn*dataperstep      ; position in file

   if (pos lt filesize) then begin
     point_lun,unit,pos+msode
     readu,unit,data0
     data[nd]=data0
   endif else begin
     goto,label0
   endelse

endfor
label0:

nx=nd
data=data[0:nx-1]

close,unit
free_lun,unit

if endian_flag then begin
case mtype of
  4: byteorder,data,/ntohl
  5: byteorder,data,/xdrtof
  6: byteorder,data,/xdrtod
endcase
endif


return
end
