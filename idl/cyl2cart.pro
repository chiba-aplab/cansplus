pro cyl2cart, dat_car, x_car, y_car, z_car, dat_cyl, r_cyl, phi_cyl, z_cyl, phi_car

  info = size(dat_cyl)
  
  nr = info[1]
  nph = info[2]
  nz = info[3]
  tmp = dblarr(nr,nph+1,nz)
  tmp[*,0:nph-1,*] = dat_cyl
  tmp[*,nph,*] = dat_cyl[*,0,*]
  dat_cyl = tmp

  nx = nr*2
  ny = nr*2
  x_car = dblarr(nx)
  y_car = dblarr(ny)
  z_car = dblarr(nz)
  dat_car = dblarr(nx,ny,nz)

  rmax=max(r_cyl,min=rmin)
  zmax=max(z_cyl,min=zmin)
  x_car[nr:nx-1] = rmin+findgen(nr)*(rmax-rmin)/(nr-1.)
  x_car[0:nr-1]  = -reverse(x_car[nr:nx-1])
  y_car = x_car
  z_car = zmin+findgen(nz)*(zmax-zmin)/(nz-1.)

  r_tmp = dblarr(nr+1)
  r_tmp[0:nr-1] = r_cyl
  r_tmp[nr] = 2*r_cyl[nr-1]-r_cyl[nr-2]
  phi_tmp = dblarr(nph+1)
  phi_tmp[0:nph-1] = phi_cyl
  phi_tmp[nph] = 2*phi_cyl[nph-1]-phi_cyl[nph-2]
  z_tmp = dblarr(nz+1)
  z_tmp[0:nz-1] = z_cyl
  z_tmp[nz] = 2*z_cyl[nz-1]-z_cyl[nz-2]

  phi_car = dat_car
  for k=0,nz-1 do begin
    kz = value_locate(z_cyl,z_car[k])
    kz_int = kz+(z_car[k]-z_tmp[kz])/(z_tmp[kz+1]-z_tmp[kz])
  for j=0,ny-1 do begin
    sgn = signum(y_car[j])
  for i=0,nx-1 do begin
    r = sqrt(x_car[i]^2+y_car[j]^2)
    phi = acos(x_car[i]*sgn/r)-!pi*(min([0.0,sgn]))
    r = min([r,rmax])
    phi_car[i,j,k] = phi

    ir = value_locate(r_cyl,r)
    jph = value_locate(phi_cyl,phi)

    ir_int = ir+(r-r_tmp[ir])/(r_tmp[ir+1]-r_tmp[ir])
    jph_int = jph+(phi-phi_tmp[jph])/(phi_tmp[jph+1]-phi_tmp[jph])

    dat_car[i,j,k] = interpolate(dat_cyl,ir_int,jph_int,kz_int)
  endfor
  endfor
  endfor

end

