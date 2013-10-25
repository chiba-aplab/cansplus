pro ncgetss,files,narg,name,t $
   ,nx=nx,narr=narr,ndx=ndx

if (n_elements(files) ge 2) then file=files[0] else file=files

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

nid=ncdf_dimid(fid,'nx')

ncdf_diminq,fid,nid,dum,nx

tid=ncdf_varid(fid,name)


ncdf_varget,fid,tid,t

if (allreadflag eq 1) then begin
  narr=indgen(nx)
  ndx=nx
endif else begin
  whr=where(narr lt nx, ndx)
  if (ndx ge 1) then begin
    narr=narr(whr)
  endif else begin
    print,'No data to be read'
    ndx=0
  endelse
endelse

ncdf_close,fid

if (ndx ge 1) then t=t(narr)

return
end
