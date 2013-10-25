pro ncputo1,idf,idi,name,x

ncdf_control,idf,/redef
idx=ncdf_vardef(idf,name,[idi])
ncdf_control,idf,/endef

xreal=float(x)
ncdf_varput,idf,idx,xreal


return
end
