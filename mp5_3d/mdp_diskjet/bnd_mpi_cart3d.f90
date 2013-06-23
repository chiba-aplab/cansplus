!======================================================================|
subroutine bnd_mpi_cart3d &
     & ( ro, pr, vx, vy, vz, bx, by, bz, phi, eta &
     & , margin, ix, jx, kx ,rank)
!======================================================================!
! NAME : bnd_mpi_cart3d
!
! PURPOSE: Apply boundary condition
!
! INPUTS:
!    margin     [integer]             margin outside the boundary
!    ix, jx, kx [integer]             grid number
!    rank       [type(cans_mpi_rank)] structure for MPI rank
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
!    phi(ix, jx, kx) [real(8)] : 
!    eta(ix, jx, kx) [real(8)] : 
!    
! OUTPUTS
!
! INTERNAL VARIABLES
!
! NOTE: No boundary condition is required in the direction 
!       in which periods = .true. 
!
! UPDATE LOG:
! 2012/12/30 initial coding (by H. ODA)
!======================================================================!

!  use mpi_domain_xz
  use cans_type
  use cans_mpi

  implicit none

  include 'mpif.h'

  integer,intent(in) :: margin, ix, jx, kx
  !type(mpidomain) :: mpid
  type(cans_mpi_rank),intent(in) :: rank

  !real(8),dimension(ix) :: x

  real(8),dimension(ix,jx,kx),intent(inout) :: ro, pr
  real(8),dimension(ix,jx,kx),intent(inout) :: vx, vy, vz
  real(8),dimension(ix,jx,kx),intent(inout) :: bx, by, bz
  real(8),dimension(ix,jx,kx),intent(inout) :: phi, eta

  integer :: i,j,k
!----------------------------------------------------------------------|

!!$  call bd_pery(margin,ro,ix,jx,kx)
!!$  call bd_pery(margin,pr,ix,jx,kx)
!!$  
!!$  call bd_pery(margin,vx,ix,jx,kx)
!!$  call bd_pery(margin,vy,ix,jx,kx)
!!$  call bd_pery(margin,vz,ix,jx,kx)
!!$  
!!$  call bd_pery(margin,phi,ix,jx,kx)
!!$  call bd_pery(margin,eta,ix,jx,kx)
!!$  
!!$  call bd_pery(margin,bx,ix,jx,kx)
!!$  call bd_pery(margin,by,ix,jx,kx)
!!$  call bd_pery(margin,bz,ix,jx,kx)
!======================================================================!
  !-- set inner x-boundary --!
  !if(mpid%mpirank_2d(1) .eq.0)then
  if ( rank%dnn == mpi_proc_null ) then
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

  !-- set outer x-boudary --!
  !if(mpid%mpirank_2d(1) .eq. (mpid%mpisize_2d(1)-1))then
  if ( rank%unn == mpi_proc_null ) then
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

  !-- set  inner z-boundary --!
  !if(mpid%mpirank_2d(2) .eq. 0)then
  if ( rank%nnd == mpi_proc_null ) then
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

  !-- set outer z-boundary --!
  !if(mpid%mpirank_2d(2) .eq. (mpid%mpisize_2d(2)-1)) then  
  if ( rank%nnu == mpi_proc_null ) then
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

end subroutine bnd_mpi_cart3d
