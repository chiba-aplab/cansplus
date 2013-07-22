subroutine calPhysValueTestSur(mdir,mpid,margin,ix,jx,kx,x,dx,y,dy,z,dz,gm &
     ,ro,pr,vx,vy,vz,bx,by,bz &
     ,iCalNum,calPhysValues)
  use mpi_domain_xz
  use convert
  implicit none

  type(mpidomain),intent(in) :: mpid

  integer,intent(in) :: margin,ix,jx,kx
  integer,intent(in) :: mdir
  real(8),intent(in) :: gm

  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: y,dy
  real(8),dimension(kx),intent(in) :: z,dz

  real(8),dimension(ix,jx,kx),intent(in) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(in) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(in) :: bx,by,bz

  real(8),dimension(ix,jx,kx) :: rx,ry,rz,ee,kee
  real(8),dimension(ix,jx,kx) :: sur,fro,fee,frx

  integer,intent(in) :: iCalNum
  real(8),dimension(iCalNum),intent(inout) :: calPhysValues
  real(8) :: roCal,rovxCal,eeCal,keeCal,surface

  real(8) :: bxsq
  integer :: i,j,k

  call convert_ptoc_m(ix,jx,kx,gm,ro,pr,vx,vy,vz,bx,by,bz &
     ,rx,ry,rz,ee)

  do k=1,kx
     do j=1,jx
        do i=1,ix
           kee(i,j,k) = 0.5d0*ro(i,j,k)*(vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2)

           fro(i,j,k) = rx(i,j,k)
           bxsq = bx(i,j,k)**2

           frx(i,j,k) = rx(i,j,k)*vx(i,j,k) - bxsq
           fee(i,j,k) = -(vx(i,j,k)*(ee(i,j,k) + pr(i,j,k) - bxsq) &
                - bx(i,j,k)*(vy(i,j,k)*by(i,j,k) + vz(i,j,k)*bz(i,j,k)))

           sur(i,j,k) = 1.0d0

        enddo
     enddo
  enddo

  call calTestSurface(mdir,mpid,margin,ix,jx,kx,x,dx,y,dy,z,dz,sur,surface)
  call calTestSurface(mdir,mpid,margin,ix,jx,kx,x,dx,y,dy,z,dz,fro,roCal)
  call calTestSurface(mdir,mpid,margin,ix,jx,kx,x,dx,y,dy,z,dz,frx,rovxCal)
  call calTestSurface(mdir,mpid,margin,ix,jx,kx,x,dx,y,dy,z,dz,fee,eeCal)
  call calTestSurface(mdir,mpid,margin,ix,jx,kx,x,dx,y,dy,z,dz,kee,keeCal)

  calPhysValues(1) = roCal
  calPhysValues(2) = rovxCal
  calPhysValues(3) = eeCal
  calPhysValues(4) = keeCal
  calPhysValues(5) = surface

  return
end subroutine calPhysValueTestSur
