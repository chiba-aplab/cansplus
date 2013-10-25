pro ncgets3,files,narg,name,data $
   ,nx=nx,narr=narr,ndx=ndx,ix=ix,jx=jx,kx=kx

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
ncdf_diminq,fid,var_info.dim(1),namedim,jx
ncdf_diminq,fid,var_info.dim(2),namedim,kx
ncdf_diminq,fid,var_info.dim(3),namedim,nx



if (allreadflag eq 1) then begin
  ncdf_varget,fid,did,data
  narr=indgen(nx)
  ndx=nx
endif else begin
  whr=where(narr lt nx, ndx)
  if (ndx ge 1) then begin
    narr=narr(whr)
    data=fltarr(ix,jx,kx,ndx)
    for nd=0,ndx-1 do begin
      ncdf_varget,fid,did,data0,count=[ix,jx,kx,1],offset=[0,0,0,narr(nd)]
      data(*,*,*,nd)=data0
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

ncgets3,file,narg,name,data0,ndx=nx,ix=ix,jx=jx,kx=kx

ncgetos,file,'ipe',ipe
ncgetos,file,'jpe',jpe
ncgetos,file,'kpe',kpe
ncgetos,file,'ipex',ipex
ncgetos,file,'jpex',jpex
ncgetos,file,'kpex',kpex

if (mf eq 0) then begin
  ncgetos,file,'margin',margin
  igx=ix*ipex-2*margin*(ipex-1)
  jgx=jx*jpex-2*margin*(jpex-1)
  kgx=kx*kpex-2*margin*(kpex-1)
  data=fltarr(igx,jgx,kgx,nx)
endif

ioff=ipe*(ix-2*margin)+margin
joff=jpe*(jx-2*margin)+margin
koff=kpe*(kx-2*margin)+margin
data[ioff:ioff+ix-1-margin,joff:joff+jx-1-margin,koff:koff+kx-1-margin,*] $
  =data0[margin:ix-1,margin:jx-1,margin:kx-1,*]

if (ipe eq 0) then $
   data[0:margin-1,joff:joff+jx-1-margin,koff:koff+kx-1-margin,*] $
   =data0[0:margin-1,margin:jx-1,margin:kx-1,*]
if (jpe eq 0) then $
   data[ioff:ioff+ix-1-margin,0:margin-1,koff:koff+kx-1-margin,*] $
   =data0[margin:ix-1,0:margin-1,margin:kx-1,*]
if (kpe eq 0) then $
   data[ioff:ioff+ix-1-margin,joff:joff+jx-1-margin,0:margin-1,*] $
   =data0[margin:ix-1,margin:jx-1,0:margin-1,*]

endfor

ix=igx
jx=jgx
kx=kgx

endelse

return
end
