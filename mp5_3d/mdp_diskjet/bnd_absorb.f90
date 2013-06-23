!======================================================================|
subroutine bnd_absorb &
     & ( ix, jx, kx, x, z, xin &
     & , roi, pri, vxi, vyi, vzi, bxi, byi, bzi &
     & , ro,  pr,  vx,  vy,  vz,  bx,  by,  bz)
!======================================================================|
! NAME : bnd_absorb
!
! PURPOSE: Apply absorb boundary
!
! INPUTS:
!    ix,jx,kx     [integer] : grid number
!    x(ix), z(kx) [real(8)] : x-coord
!    xin          [real(8)] : x-ccord applied absorb boundary
!    roi (ix, jx, kx) [real(8)] : initial density
!    pri (ix, jx, kx) [real(8)] : initial pressure
!    vxi (ix, jx, kx) [real(8)] : initial x-velocity 
!    vyi (ix, jx, kx) [real(8)] : initial y-velocity
!    vzi (ix, jx, kx) [real(8)] : initial z-velocity
!    bxi (ix, jx, kx) [real(8)] : initial x-magnetic field
!    byi (ix, jx, kx) [real(8)] : initial y-magnetic field
!    bzi (ix, jx, kx) [real(8)] : initial z-magnetic field

! 
! INOUTPUTS:
!    ro (ix, jx, kx) [real(8)] : density
!    pr (ix, jx, kx) [real(8)] : pressure
!    vx (ix, jx, kx) [real(8)] : x-velocity 
!    vy (ix, jx, kx) [real(8)] : y-velocity
!    vz (ix, jx, kx) [real(8)] : z-velocity
!    bx (ix, jx, kx) [real(8)] : x-magnetic field
!    by (ix, jx, kx) [real(8)] : y-magnetic field
!    bz (ix, jx, kx) [real(8)] : z-magnetic field
!    
! OUTPUTS:
!
! INTERNAL VARIABLES
!
! UPDATE LOG:
! 2012/12/28 initial coding  (by H. ODA)
! 
!======================================================================!

  implicit none

  integer,intent(in) :: ix,jx,kx

  real(8),intent(in) :: xin

  real(8),dimension(ix),intent(in) :: x
  real(8),dimension(kx),intent(in) :: z

  real(8),dimension(ix,jx,kx),intent(in) :: roi,pri
  real(8),dimension(ix,jx,kx),intent(in) :: vxi,vyi,vzi
  real(8),dimension(ix,jx,kx),intent(in) :: bxi,byi,bzi

  real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz
!----------------------------------------------------------------------|

  call absorb(ix,jx,kx,x,z,xin,ro,roi)
  call absorb(ix,jx,kx,x,z,xin,pr,pri)

  call absorb(ix,jx,kx,x,z,xin,vx,vxi)
  call absorb(ix,jx,kx,x,z,xin,vy,vyi)
  call absorb(ix,jx,kx,x,z,xin,vz,vzi)

  call absorb(ix,jx,kx,x,z,xin,bx,bxi)
  call absorb(ix,jx,kx,x,z,xin,by,byi)
  call absorb(ix,jx,kx,x,z,xin,bz,bzi)

end subroutine bnd_absorb
