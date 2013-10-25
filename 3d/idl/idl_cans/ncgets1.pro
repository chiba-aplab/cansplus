pro ncgets1,files,narg,name,data $
   ,ix=ix,nx=nx,narr=narr,ndx=ndx

mfx=n_elements(files)

if (mfx eq 1) then begin

file=files[0]

sz=size(narg)
allreadflag=0
if (sz(0) ge 1) then  begin
  narr=narg
endif else begin
  if (narg ge 0) then begin
    narr=[narg]
  endif else begin
    allreadflag=1
  endelse
endelse

fid=ncdf_open(file,/nowrite)

did=ncdf_varid(fid,name)
var_info=ncdf_varinq(fid,did)
ncdf_diminq,fid,var_info.dim(0),namedim,ix
ncdf_diminq,fid,var_info.dim(1),namedim,nx

if (allreadflag eq 1) then begin
  ncdf_varget,fid,did,data
  narr=indgen(nx)
  ndx=nx
endif else begin
  whr=where(narr lt nx, ndx)
  if (ndx ge 1) then begin
    narr=narr(whr)
    data=fltarr(ix,ndx)
    for nd=0,ndx-1 do begin
      ncdf_varget,fid,did,data0,count=[ix,1],offset=[0,narr(nd)]
      data(*,nd)=data0
    endfor
  endif else begin
    print,'No data to be read'
    ndx=0
  endelse
endelse

ncdf_close,fid

endif else begin

for mf=0,mfx-1 do begin

file=files[mf]

ncgets1,file,narg,name,data0,ndx=nx

if (mf eq 0) then begin
  ncgetos,file,'margin',margin
  sz=size(data0)
  ix=sz[1]
  igx=ix*mfx-2*margin*(mfx-1)
  data=fltarr(igx,nx)
endif

moff=mf*(ix-2*margin)+margin
data[moff:moff+ix-1-margin,*]=data0[margin:ix-1,*]

if (mf eq 0) then data[0:margin-1,*]=data0[0:margin-1,*]

endfor

ix=igx

endelse

return
end
