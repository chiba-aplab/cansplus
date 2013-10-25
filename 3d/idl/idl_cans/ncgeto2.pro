pro ncgeto2,files,name,data $
   ,ix=ix,jx=jx,ipex=ipex,jpex=jpex

mfx=n_elements(files)

if (mfx eq 1) then begin ; single file

file=files[0]

fid=ncdf_open(file,/nowrite)
did=ncdf_varid(fid,name)
ncdf_varget,fid,did,data
ncdf_close,fid

endif else begin ; multi files

for mf=0,mfx-1 do begin

file=files[mf]

ncgeto2,file,name,data0

jpe=mf/ipex
ipe=mf mod ipex

if (mf eq 0) then begin
  ncgetos,file,'margin',margin
  sz=size(data0)
  ix=sz[1]
  igx=ix*ipex-2*margin*(ipex-1)
  jx=sz[2]
  jgx=jx*jpex-2*margin*(jpex-1)
  data=fltarr(igx,jgx)
endif

ioff=ipe*(ix-2*margin)+margin
joff=jpe*(jx-2*margin)+margin
data[ioff:ioff+ix-1-margin,joff:joff+jx-1-margin] $
   =data0[margin:ix-1,margin:jx-1]

if (ipe eq 0) then data[0:margin-1,joff:joff+jx-1-margin] $
  =data0[0:margin-1,margn:jx-1]
if (jpe eq 0) then data[ioff:ioff+ix-1-margin,0:margin-1] $
  =data0[margn:ix-1,0:margin-1]

endfor


endelse

sz=size(data)
ix=sz[1]
jx=sz[2]

return
end
