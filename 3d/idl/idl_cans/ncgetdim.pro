pro ncgetdim,file,name,ix

fid=ncdf_open(file,/nowrite)
iid=ncdf_dimid(fid,name)
ncdf_diminq,fid,iid,dum,ix


return
end
