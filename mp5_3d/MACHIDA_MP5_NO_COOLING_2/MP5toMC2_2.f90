subroutine MP5toMC2_2(ix,jx,kx,x,dx,y,dy,z,dz &
     ,ro,pr,vx,vy,vz,bx,by,bz,phi &
     ,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw &
     ,mdir,floor,limit)
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

  real(8),intent(in) :: floor,limit

  call switchMP5toMC2_2(mdir,ix,jx,kx,ro,ro &
       ,row,floor,limit,x,dx,y,dy,z,dz)

  call switchMP5toMC2_2(mdir,ix,jx,kx,ro,pr &
       ,prw,floor,limit,x,dx,y,dy,z,dz)

  call switchMP5toMC2_2(mdir,ix,jx,kx,ro,vx &
       ,vxw,floor,limit,x,dx,y,dy,z,dz)

  call switchMP5toMC2_2(mdir,ix,jx,kx,ro,vy &
       ,vyw,floor,limit,x,dx,y,dy,z,dz)

  call switchMP5toMC2_2(mdir,ix,jx,kx,ro,vz &
       ,vzw,floor,limit,x,dx,y,dy,z,dz)

  call switchMP5toMC2_2(mdir,ix,jx,kx,ro,bx &
       ,bxw,floor,limit,x,dx,y,dy,z,dz)

  call switchMP5toMC2_2(mdir,ix,jx,kx,ro,by &
       ,byw,floor,limit,x,dx,y,dy,z,dz)

  call switchMP5toMC2_2(mdir,ix,jx,kx,ro,bz &
       ,bzw,floor,limit,x,dx,y,dy,z,dz)

  call switchMP5toMC2_2(mdir,ix,jx,kx,ro,phi &
       ,phiw,floor,limit,x,dx,y,dy,z,dz)

  return
end subroutine MP5toMC2_2
