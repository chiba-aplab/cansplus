subroutine calPhysValueTest(margin,ix,jx,kx,x,dx,y,dy,z,dz,gm,ro,pr,vx,vy,vz,bx,by,bz &
     ,iCalNum,calPhysValues)
  use convert
  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  real(8),intent(in) :: gm

  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: y,dy
  real(8),dimension(kx),intent(in) :: z,dz

  real(8),dimension(ix,jx,kx),intent(in) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(in) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(in) :: bx,by,bz

  real(8),dimension(ix,jx,kx) :: rx,ry,rz,ee,kee
  real(8),dimension(ix,jx,kx) :: vol

  integer,intent(in) :: iCalNum
  real(8),dimension(iCalNum),intent(inout) :: calPhysValues
  real(8) :: roCal,rovxCal,eeCal,keeCal,volume

  integer :: i,j,k

  call convert_ptoc_m(ix,jx,kx,gm,ro,pr,vx,vy,vz,bx,by,bz &
     ,rx,ry,rz,ee)

  do k=1,kx
     do j=1,jx
        do i=1,ix
           kee(i,j,k) = 0.5d0*ro(i,j,k)*(vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2)
           vol(i,j,k) = 1.0d0
        enddo
     enddo
  enddo

  call calTest(margin,ix,jx,kx,x,dx,y,dy,z,dz,vol,volume)
  call calTest(margin,ix,jx,kx,x,dx,y,dy,z,dz,ro,roCal)
  call calTest(margin,ix,jx,kx,x,dx,y,dy,z,dz,rx,rovxCal)
  call calTest(margin,ix,jx,kx,x,dx,y,dy,z,dz,ee,eeCal)
  call calTest(margin,ix,jx,kx,x,dx,y,dy,z,dz,kee,keeCal)

  calPhysValues(1) = roCal
  calPhysValues(2) = rovxCal
  calPhysValues(3) = eeCal
  calPhysValues(4) = keeCal
  calPhysValues(5) = volume

  return
end subroutine calPhysValueTest
