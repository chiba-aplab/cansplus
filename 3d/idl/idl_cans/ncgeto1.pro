pro ncgeto1,files,name,data $
   ,ix=ix

mfx=n_elements(files)

if (mfx eq 1) then begin ; single file

file=files[0]

fid=ncdf_open(file,/nowrite)
xid=ncdf_varid(fid,name)
ncdf_varget,fid,xid,data
ncdf_close,fid

endif else begin ; multi files

for mf=0,mfx-1 do begin

file=files[mf]

ncgeto1,file,name,data0

if (mf eq 0) then begin
  ncgetos,file,'margin',margin
  sz=size(data0)
  ix=sz[1]
  igx=ix*mfx-2*margin*(mfx-1)
  data=fltarr(igx)
endif

moff=mf*(ix-2*margin)+margin
data[moff:moff+ix-1-margin]=data0[margin:ix-1]

if (mf eq 0) then data[0:margin-1]=data0[0:margin-1]

endfor


endelse

ix=n_elements(data)


return
end
