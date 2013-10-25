pro ncdefs3,idf,name,idda,idn,idi,idj,idk

ncdf_control,idf,/redef
idda=ncdf_vardef(idf,name,[idi,idj,idk,idn])
ncdf_control,idf,/endef

return
end
