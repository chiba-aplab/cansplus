pro cyl2cart_2d, dat_car, x_car, z_car, dat_cyl, r_cyl, z_cyl

  info = size(dat_cyl)
  nr = info[1]
  nz = info[2]

  nx = nr*2
  x_car = fltarr(nx)
  z_car = fltarr(nz)
  dat_car = fltarr(nx,nz)

  rmax=max(r_cyl,min=rmin)
  zmax=max(z_cyl,min=zmin)
  x_car[nr:nx-1] = 0.0+(findgen(nr)+0.5)*(rmax-rmin)/(nr-1.)
  x_car[0:nr-1] = -reverse(x_car[nr:nx-1])
  z_car = zmin+findgen(nz)*(zmax-zmin)/(nz-1.)

  r_tmp = fltarr(nr+1)
  r_tmp[0:nr-1] = r_cyl
  r_tmp[nr] = 2*r_cyl[nr-1]-r_cyl[nr-2]
  z_tmp = fltarr(nz+1)
  z_tmp[0:nz-1] = z_cyl
  z_tmp[nz] = 2*z_cyl[nz-1]-z_cyl[nz-2]

  ir_int = fltarr(nx)
  for i=0,nx-1 do begin
    r = sqrt(x_car[i]^2)
    r = min([r,rmax])
    ir = value_locate(r_cyl,r)
    ir_int[i] = ir+(r-r_tmp[ir])/(r_tmp[ir+1]-r_tmp[ir])
  endfor
  kz_int = fltarr(nz)
  for k=0,nz-1 do begin
    kz = value_locate(z_cyl,z_car[k])
    kz_int[k] = kz+(z_car[k]-z_tmp[kz])/(z_tmp[kz+1]-z_tmp[kz])
  endfor

  dat_car = interpolate(dat_cyl,ir_int,kz_int,/grid)

end


pro cyl2cart_3d, dat_car, x_car, y_car, z_car, dat_cyl, r_cyl, phi_cyl, z_cyl, nzs, nze

  info = size(dat_cyl)
  
  nr = info[1]
  nph = info[2]
  nz = info[3]
  tmp = fltarr(nr,nph+1,nz)
  tmp[*,0:nph-1,*] = dat_cyl
  tmp[*,nph,*] = dat_cyl[*,0,*]

  nx = nr*2
  ny = nr*2
  x_car = fltarr(nx)
  y_car = fltarr(ny)
  z_car = fltarr(nz)
  dat_car = fltarr(nx,ny,nz)

  rmax=max(r_cyl,min=rmin)
  zmax=max(z_cyl,min=zmin)
  x_car[nr:nx-1] = 0.0+(findgen(nr)+0.5)*(rmax-rmin)/(nr-1.)
  x_car[0:nr-1]  = -reverse(x_car[nr:nx-1])
  y_car = x_car
  z_car = zmin+findgen(nz)*(zmax-zmin)/(nz-1.)

  r_tmp = fltarr(nr+1)
  r_tmp[0:nr-1] = r_cyl
  r_tmp[nr] = 2*r_cyl[nr-1]-r_cyl[nr-2]
  phi_tmp = fltarr(nph+1)
  phi_tmp[0:nph-1] = phi_cyl
  phi_tmp[nph] = 2*phi_cyl[nph-1]-phi_cyl[nph-2]
  z_tmp = fltarr(nz+1)
  z_tmp[0:nz-1] = z_cyl
  z_tmp[nz] = 2*z_cyl[nz-1]-z_cyl[nz-2]

  for k=nzs,nze do begin
    kz = value_locate(z_cyl,z_car[k])
    kz_int = kz+(z_car[k]-z_tmp[kz])/(z_tmp[kz+1]-z_tmp[kz])
  for j=0,ny-1 do begin
    sgn = signum(y_car[j])
  for i=0,nx-1 do begin
    r = sqrt(x_car[i]^2+y_car[j]^2)
    phi = acos(x_car[i]*sgn/r)-!pi*(min([0.0,sgn]))
    r = min([r,rmax])

    ir = value_locate(r_cyl,r)
    jph = value_locate(phi_cyl,phi)

    ir_int = ir+(r-r_tmp[ir])/(r_tmp[ir+1]-r_tmp[ir])
    jph_int = jph+(phi-phi_tmp[jph])/(phi_tmp[jph+1]-phi_tmp[jph])

    dat_car[i,j,k] = interpolate(tmp,ir_int,jph_int,kz_int)
  endfor
  endfor
  endfor

  dat_car = dat_car[*,*,nzs:nze]
  z_car = z_car[nzs:nze]

end


pro para_range, nst, ned, n1, n2, isize, irank

;;start and end of loop counter
  iwork1 = (n2-n1+1)/isize
  iwork2 = (n2-n1+1) mod (isize)
  nst = irank*iwork1+n1+min([irank,iwork2])
  ned = nst+iwork1-1
  if(iwork2 gt irank)then ned = ned+1
  
end


pro cyl2cart, para1, para2, para3, para4, para5, para6, para7, para8

case n_params() of
  6: cyl2cart_2d, para1, para2, para3, para4, para5, para6
  8: begin

     para1 = fltarr(n_elements(para6)*2,n_elements(para6)*2,n_elements(para8))
     para4 = fltarr(n_elements(para8))

     ntask = min([!cpu.tpool_nthreads-1,n_elements(para8)])
     if(float(!version.release) lt 8.6)then ntask=min([ntask,16])
     bridges = objarr(ntask)
     for l=0,ntask-1 do begin
       bridges[l] = obj_new('idl_idlbridge')
     endfor

     for l=0,ntask-1 do begin
       para_range,nzs,nze,0,n_elements(para8)-1,ntask,l
       bridges[l].execute,"@"+pref_get('IDL_STARTUP')
       bridges[l].setvar,"nzs",nzs
       bridges[l].setvar,"nze",nze
       bridges[l].setvar,"dat_cyl",para5
       bridges[l].setvar,"x_cyl",para6
       bridges[l].setvar,"y_cyl",para7
       bridges[l].setvar,"z_cyl",para8
       bridges[l].execute,"pushd, '"+getenv('PWD')+"'"
       bridges[l].execute,".compile cyl2cart"
       bridges[l].execute,"cyl2cart_3d,dat_car,x_car,y_car,z_car,dat_cyl,x_cyl,y_cyl,z_cyl,nzs,nze",/nowait
     endfor

     ostat = intarr(ntask)
     while(min(ostat) eq 0)do begin
       itask = where(ostat eq 0)
       foreach l, itask do begin
         if(bridges[l].status() eq 0) then begin
           nzs = bridges[l].getvar('nzs')
           nze = bridges[l].getvar('nze')
           para1[*,*,nzs:nze] = bridges[l].getvar('dat_car')    
           para2 = bridges[l].getvar('x_car')
           para3 = bridges[l].getvar('y_car')
           para4[nzs:nze] = bridges[l].getvar('z_car')
           ostat[l] = 1
         endif
       endforeach
     endwhile

     for l=0,ntask-1 do begin
       bridges[l].cleanup
     endfor

     end
  else: print,'Number of arguments is incorrect.'
endcase

end

