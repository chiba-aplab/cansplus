pro ncdefss,idf,name,idda,idn

ncdf_control,idf,/redef
idda=ncdf_vardef(idf,name,idn)
ncdf_control,idf,/endef

return
end
