subroutine MP5toMC2(ix,jx,kx,x,dx,y,dy,z,dz &
     ,ro_input,ro,pr,vx,vy,vz,bx,by,bz,phi &
     ,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw &
     ,mdir,floor,ratio,xin)
  implicit none

  integer,intent(in) :: ix,jx,kx

  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: y,dy
  real(8),dimension(kx),intent(in) :: z,dz

  real(8),dimension(ix,jx,kx),intent(in) :: ro_input,ro,pr,vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(in) :: bx,by,bz

  real(8),dimension(ix,jx,kx),intent(in) :: phi

  real(8),dimension(ix,jx,kx,2) :: row,prw,vxw,vyw,vzw
  real(8),dimension(ix,jx,kx,2) :: bxw,byw,bzw,phiw

  integer :: mdir

  real(8),intent(in) :: floor,ratio,xin

  call switchMP5toMC2(mdir,ix,jx,kx,ro_input,ro &
       ,row,floor,ratio,xin,x)

  call switchMP5toMC2(mdir,ix,jx,kx,ro_input,pr &
       ,prw,floor,ratio,xin,x)

  call switchMP5toMC2(mdir,ix,jx,kx,ro_input,vx &
       ,vxw,floor,ratio,xin,x)

  call switchMP5toMC2(mdir,ix,jx,kx,ro_input,vy &
       ,vyw,floor,ratio,xin,x)

  call switchMP5toMC2(mdir,ix,jx,kx,ro_input,vz &
       ,vzw,floor,ratio,xin,x)

  call switchMP5toMC2(mdir,ix,jx,kx,ro_input,bx &
       ,bxw,floor,ratio,xin,x)

  call switchMP5toMC2(mdir,ix,jx,kx,ro_input,by &
       ,byw,floor,ratio,xin,x)

  call switchMP5toMC2(mdir,ix,jx,kx,ro_input,bz &
       ,bzw,floor,ratio,xin,x)

  call switchMP5toMC2(mdir,ix,jx,kx,ro_input,phi &
       ,phiw,floor,ratio,xin,x)

  return
end subroutine MP5toMC2
