subroutine bnd(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,x,z &
           ,xin,roi,pri,vxi,vyi,vzi,bxi,byi,bzi)

  use mpi_setup, only : mpid, mnull
  use boundary

  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  real(8),dimension(ix),intent(in) :: x
  real(8),dimension(kx),intent(in) :: z
  real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(inout) :: phi,eta
  real(8),dimension(ix,jx,kx),intent(in) :: roi,pri
  real(8),dimension(ix,jx,kx),intent(in) :: vxi,vyi,vzi
  real(8),dimension(ix,jx,kx),intent(in) :: bxi,byi,bzi
  real(8),intent(in) :: xin
  
!======================================================================
! inner x-boundary
  if(mpid%l == mnull)then
     call bd_synpx(0,margin,ro,ix,jx,kx)
     call bd_synpx(0,margin,pr,ix,jx,kx)
     call bd_synnx(0,margin,vx,ix,jx,kx)
     call bd_synnx(0,margin,vy,ix,jx,kx)
     call bd_synpx(0,margin,vz,ix,jx,kx)
     call bd_synnx(0,margin,bx,ix,jx,kx)
     call bd_synnx(0,margin,by,ix,jx,kx)
     call bd_synpx(0,margin,bz,ix,jx,kx)
     call bd_synpx(0,margin,phi,ix,jx,kx)
     call bd_synpx(0,margin,eta,ix,jx,kx)
  end if

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
  if(mpid%b == mnull)then
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

!----------------------------------------------------------------------|
! central boundary
  call absorb(ix,jx,kx,x,z,xin,ro,roi)
  call absorb(ix,jx,kx,x,z,xin,pr,pri)
  call absorb(ix,jx,kx,x,z,xin,vx,vxi)
  call absorb(ix,jx,kx,x,z,xin,vy,vyi)
  call absorb(ix,jx,kx,x,z,xin,vz,vzi)
  call absorb(ix,jx,kx,x,z,xin,bx,bxi)
  call absorb(ix,jx,kx,x,z,xin,by,byi)
  call absorb(ix,jx,kx,x,z,xin,bz,bzi)

end subroutine bnd


subroutine absorb(ix,jx,kx,x,z,xin,qq,qqi)

  implicit none

  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: xin
  real(8),intent(in),dimension(ix) :: x
  real(8),intent(in),dimension(kx) :: z
  real(8),intent(in),dimension(ix,jx,kx) :: qqi
  real(8),intent(out),dimension(ix,jx,kx) :: qq
  real(8) :: ai,ss,dx0
  integer :: i,j,k

  dx0 = 0.01d0
  do k=1,kx
    do j=1,jx
      do i=1,ix
         ss = sqrt(x(i)**2+z(k)**2)
         if(ss <= xin)then
         ai = 0.1d0*(1.0d0-tanh((ss-xin+5.0d0*dx0)/(2.0d0*dx0)))
         qq(i,j,k) = (1.0d0-ai)*qq(i,j,k)+ai*qqi(i,j,k)
         endif
      enddo
    enddo
  enddo
end subroutine absorb

