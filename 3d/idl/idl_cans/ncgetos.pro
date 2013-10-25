pro ncgetos,files,name,value

if (n_elements(files) ge 2) then file=files[0] else file=files

fid=ncdf_open(file,/nowrite)
xid=ncdf_varid(fid,name)
ncdf_varget,fid,xid,value
ncdf_close,fid


return
end
