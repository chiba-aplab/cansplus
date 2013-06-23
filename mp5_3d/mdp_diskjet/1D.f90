subroutine 1D(ix,jx,kx,x,dx,y,dy,z,dz &
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


  call lr_state_1D(mdir,ix,jx,kx,ro &
       ,row)
  call lr_state_1D(mdir,ix,jx,kx,pr &
       ,prw)

  call lr_state_1D(mdir,ix,jx,kx,vx &
       ,vxw)
  call lr_state_1D(mdir,ix,jx,kx,vy &
       ,vyw)
  call lr_state_1D(mdir,ix,jx,kx,vz &
       ,vzw)

  call lr_state_1D(mdir,ix,jx,kx,bx &
       ,bxw)
  call lr_state_1D(mdir,ix,jx,kx,by &
       ,byw)
  call lr_state_1D(mdir,ix,jx,kx,bz &
       ,bzw)

  call lr_state_1D(mdir,ix,jx,kx,phi &
       ,phiw)


  return
end subroutine 1D
