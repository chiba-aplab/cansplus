subroutine MC2_dxyz(ix,jx,kx,x,dx,y,dy,z,dz &
     ,ro,pr,vx,vy,vz,bx,by,bz,phi &
     ,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw &
     ,mdir)
  implicit none

  integer,intent(in) :: ix,jx,kx

  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: y,dy
  real(8),dimension(kx),intent(in) :: z,dz

  real(8),dimension(ix,jx,kx),intent(in) :: ro,pr,vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(in) :: bx,by,bz

  real(8),dimension(ix,jx,kx),intent(in) :: phi

  real(8),dimension(ix,jx,kx,2) :: row,prw,vxw,vyw,vzw
  real(8),dimension(ix,jx,kx,2) :: bxw,byw,bzw,phiw

  integer :: mdir


  call lr_state_MC2_dxyz(mdir,ix,jx,kx,dx,dy,dz,ro &
       ,row)
  call lr_state_MC2_dxyz(mdir,ix,jx,kx,dx,dy,dz,pr &
       ,prw)

  call lr_state_MC2_dxyz(mdir,ix,jx,kx,dx,dy,dz,vx &
       ,vxw)
  call lr_state_MC2_dxyz(mdir,ix,jx,kx,dx,dy,dz,vy &
       ,vyw)
  call lr_state_MC2_dxyz(mdir,ix,jx,kx,dx,dy,dz,vz &
       ,vzw)

  call lr_state_MC2_dxyz(mdir,ix,jx,kx,dx,dy,dz,bx &
       ,bxw)
  call lr_state_MC2_dxyz(mdir,ix,jx,kx,dx,dy,dz,by &
       ,byw)
  call lr_state_MC2_dxyz(mdir,ix,jx,kx,dx,dy,dz,bz &
       ,bzw)

  call lr_state_MC2_dxyz(mdir,ix,jx,kx,dx,dy,dz,phi &
       ,phiw)


  return
end subroutine MC2_dxyz
