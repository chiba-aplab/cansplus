pro ncopndm,idf,file,idn,idi,ix,idj,jx,idk,kx

idf=ncdf_create(file,/noclobber)

idn=ncdf_dimdef(idf,'nx',/unlimited)
idi=ncdf_dimdef(idf,'ix',ix)
idj=ncdf_dimdef(idf,'jx',jx)
idk=ncdf_dimdef(idf,'kx',kx)

ncdf_control,idf,/endef

return
end
