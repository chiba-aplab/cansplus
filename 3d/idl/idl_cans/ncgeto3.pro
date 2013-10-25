pro ncgeto3,files,name,data $
   ,ix=ix,jx=jx,kx=kx

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

ncgeto3,file,name,data0

if (mf eq 0) then begin
  ncgetos,file,'margin',margin
  sz=size(data0)
  ix=sz[1] & igx=ix
  jx=sz[2] & jgx=jx
  kx=sz[3]
  kgx=kx*mfx-2*margin*(mfx-1)
  data=fltarr(igx,jgx,kgx)
endif

moff=mf*(kx-2*margin)+margin
data[*,*,moff:moff+kx-1-margin]=data0[*,*,margin:kx-1]

if (mf eq 0) then data[*,*,0:margin-1]=data0[*,*,0:margin-1]

endfor


endelse

sz=size(data)
ix=sz[1]
jx=sz[2]
kx=sz[3]

return
end
