module bnd

  implicit none
  private

  public :: bnd__exec, bnd__absorb


contains


  subroutine bnd__exec(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta)

  use mpi_setup, only : mpid, mnull
  use boundary

  integer,intent(in) :: margin,ix,jx,kx
  real(8),intent(inout),dimension(ix,jx,kx) :: ro,pr
  real(8),intent(inout),dimension(ix,jx,kx) :: vx,vy,vz
  real(8),intent(inout),dimension(ix,jx,kx) :: bx,by,bz
  real(8),intent(inout),dimension(ix,jx,kx) :: phi,eta

!======================================================================
! inter-process communication by MPI
  call boundary__mpi(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi)

!======================================================================
! inner x-boundary
  call boundary__mpi_cyl(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi)

!======================================================================
! outer x-boudary
  if(mpid%r == mnull)then
     call bd_frex(1,margin,ro,ix,jx,kx)
     call bd_frex(1,margin,pr,ix,jx,kx)
     call bd_frex(1,margin,vx,ix,jx,kx)
     call bd_frex(1,margin,vy,ix,jx,kx)
     call bd_frex(1,margin,vz,ix,jx,kx)
     call bd_frex(1,margin,bx,ix,jx,kx)
     call bd_frex(1,margin,by,ix,jx,kx)
     call bd_frex(1,margin,bz,ix,jx,kx)
     call bd_frex(1,margin,phi,ix,jx,kx)
     call bd_frex(1,margin,eta,ix,jx,kx)
  end if

!----------------------------------------------------------------------|
! inner z-boundary
  if(mpid%d == mnull)then
     call bd_frez(0,margin,ro,ix,jx,kx)
     call bd_frez(0,margin,pr,ix,jx,kx)
     call bd_frez(0,margin,vx,ix,jx,kx)
     call bd_frez(0,margin,vy,ix,jx,kx)
     call bd_frez(0,margin,vz,ix,jx,kx)
     call bd_frez(0,margin,bx,ix,jx,kx)
     call bd_frez(0,margin,by,ix,jx,kx)
     call bd_frez(0,margin,bz,ix,jx,kx)
     call bd_frez(0,margin,phi,ix,jx,kx)
     call bd_frez(0,margin,eta,ix,jx,kx)
  end if

!----------------------------------------------------------------------|
! outer z-boundary
  if(mpid%t == mnull) then  
     call bd_frez(1,margin,ro,ix,jx,kx)
     call bd_frez(1,margin,pr,ix,jx,kx)
     call bd_frez(1,margin,vx,ix,jx,kx)
     call bd_frez(1,margin,vy,ix,jx,kx)
     call bd_frez(1,margin,vz,ix,jx,kx)
     call bd_frez(1,margin,bx,ix,jx,kx)
     call bd_frez(1,margin,by,ix,jx,kx)
     call bd_frez(1,margin,bz,ix,jx,kx)
     call bd_frez(1,margin,phi,ix,jx,kx)
     call bd_frez(1,margin,eta,ix,jx,kx)
  end if

  end subroutine bnd__exec


  subroutine bnd__absorb(ix,jx,kx,x,z,xin,ro,pr,vx,vy,vz,bx,by,bz,phi &
                        ,roi,pri,vxi,vyi,vzi,bxi,byi,bzi,phii)

  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: xin
  real(8),intent(in),dimension(ix) :: x
  real(8),intent(in),dimension(kx) :: z
  real(8),intent(in),dimension(ix,jx,kx) :: roi,pri
  real(8),intent(in),dimension(ix,jx,kx) :: vxi,vyi,vzi
  real(8),intent(in),dimension(ix,jx,kx) :: bxi,byi,bzi
  real(8),intent(in),dimension(ix,jx,kx) :: phii
  real(8),intent(inout),dimension(ix,jx,kx) :: ro,pr
  real(8),intent(inout),dimension(ix,jx,kx) :: vx,vy,vz
  real(8),intent(inout),dimension(ix,jx,kx) :: bx,by,bz
  real(8),intent(inout),dimension(ix,jx,kx) :: phi

  real(8) :: ai,ss,dx0
  integer :: i,j,k

  dx0 = 0.01d0
  do k=1,kx
    do j=1,jx
      do i=1,ix
         ss = sqrt(x(i)**2+z(k)**2)
         if(ss <= xin)then
            ai = 0.1d0*(1.0d0-tanh((ss-xin+5.0d0*dx0)/(2.0d0*dx0)))
            ro(i,j,k) = (1.0d0-ai)*ro(i,j,k)+ai*roi(i,j,k)
            pr(i,j,k) = (1.0d0-ai)*pr(i,j,k)+ai*pri(i,j,k)
            vx(i,j,k) = (1.0d0-ai)*vx(i,j,k)+ai*vxi(i,j,k)
            vy(i,j,k) = (1.0d0-ai)*vy(i,j,k)+ai*vyi(i,j,k)
            vz(i,j,k) = (1.0d0-ai)*vz(i,j,k)+ai*vzi(i,j,k)
            bx(i,j,k) = (1.0d0-ai)*bx(i,j,k)+ai*bxi(i,j,k)
            by(i,j,k) = (1.0d0-ai)*by(i,j,k)+ai*byi(i,j,k)
            bz(i,j,k) = (1.0d0-ai)*bz(i,j,k)+ai*bzi(i,j,k)
            phi(i,j,k) = (1.0d0-ai)*phi(i,j,k)+ai*phii(i,j,k)
         endif
      enddo
    enddo
  enddo

  end subroutine bnd__absorb


end module bnd
