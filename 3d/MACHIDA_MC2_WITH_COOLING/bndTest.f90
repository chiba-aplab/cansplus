subroutine bndTest(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,x)

  use mpi_domain_xz
  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  type(mpidomain) :: mpid

  real(8),dimension(ix) :: x

  real(8),dimension(ix,jx,kx) :: ro,pr
  real(8),dimension(ix,jx,kx) :: vx,vy,vz
  real(8),dimension(ix,jx,kx) :: bx,by,bz
  real(8),dimension(ix,jx,kx) :: phi,eta

  real(8) :: gm
  real(8) :: ro0,pr0
  real(8) :: vx0,vy0,vz0
  real(8) :: zbnd

  integer :: i,j,k

  gm = 5.0d0/3.0d0
  ro0 = 1.0d0
  pr0 = 1.0d0/gm
  vx0 = -1.0d0
  vy0 = 0.0d0
  vz0 = 0.0d0
  zbnd = 0.0d0

  call bd_pery(margin,ro,ix,jx,kx)
  call bd_pery(margin,pr,ix,jx,kx)
  
  call bd_pery(margin,vx,ix,jx,kx)
  call bd_pery(margin,vy,ix,jx,kx)
  call bd_pery(margin,vz,ix,jx,kx)
  
  call bd_pery(margin,phi,ix,jx,kx)
  call bd_pery(margin,eta,ix,jx,kx)
  
  call bd_pery(margin,bx,ix,jx,kx)
  call bd_pery(margin,by,ix,jx,kx)
  call bd_pery(margin,bz,ix,jx,kx)

!======================================================================
! inner x-boundary

  if(mpid%mpirank_2d(1) .eq.0)then

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
!!$     call bd_frex(0,margin,ro,ix,jx,kx)
!!$     call bd_frex(0,margin,pr,ix,jx,kx)
!!$     
!!$     call bd_frex(0,margin,vx,ix,jx,kx)
!!$     call bd_frex(0,margin,vy,ix,jx,kx)
!!$     call bd_frex(0,margin,vz,ix,jx,kx)
!!$     
!!$     call bd_frex(0,margin,bx,ix,jx,kx)
!!$     call bd_frex(0,margin,by,ix,jx,kx)
!!$     call bd_frex(0,margin,bz,ix,jx,kx)
!!$
!!$     call bd_frex(0,margin,phi,ix,jx,kx)
!!$     call bd_frex(0,margin,eta,ix,jx,kx)

  end if
!======================================================================
! outer x-boudary

  if(mpid%mpirank_2d(1) .eq. (mpid%mpisize_2d(1)-1))then

!!$     call bd_consx(1,margin,ro0,ro,ix,jx,kx)
!!$     call bd_consx(1,margin,pr0,pr,ix,jx,kx)
!!$
!!$     call bd_consx(1,margin,vx0,vx,ix,jx,kx)
!!$     call bd_consx(1,margin,vy0,vy,ix,jx,kx)
!!$     call bd_consx(1,margin,vz0,vz,ix,jx,kx)
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
! outer z-boundary
!


  if(mpid%mpirank_2d(2) .eq. (mpid%mpisize_2d(2)-1)) then  
     call bd_frez(1,margin,ro,ix,jx,kx)
     call bd_frez(1,margin,pr,ix,jx,kx)

!!$     call bd_consz(1,margin,zbnd,vx,ix,jx,kx)
!!$     call bd_consz(1,margin,zbnd,vy,ix,jx,kx)
!!$     call bd_consz(1,margin,zbnd,vz,ix,jx,kx)
     call bd_frez(1,margin,vx,ix,jx,kx)
     call bd_frez(1,margin,vy,ix,jx,kx)
     call bd_frez(1,margin,vz,ix,jx,kx)
     
     call bd_frez(1,margin,bx,ix,jx,kx)
     call bd_frez(1,margin,by,ix,jx,kx)
     call bd_frez(1,margin,bz,ix,jx,kx)

!!$  call bd_openbz_cyl(1,margin,bxs,bzs,bx,by,bz &
!!$       ,ix,jx,kx,x,dx,dz)

     call bd_frez(1,margin,phi,ix,jx,kx)
     call bd_frez(1,margin,eta,ix,jx,kx)
  end if

!----------------------------------------------------------------------|
! inner z-boundary
!

  if(mpid%mpirank_2d(2) .eq. 0)then

     call bd_frez(0,margin,ro,ix,jx,kx)
     call bd_frez(0,margin,pr,ix,jx,kx)

!!$     call bd_consz(0,margin,zbnd,vx,ix,jx,kx)
!!$     call bd_consz(0,margin,zbnd,vy,ix,jx,kx)
!!$     call bd_consz(0,margin,zbnd,vz,ix,jx,kx)
     call bd_frez(0,margin,vx,ix,jx,kx)
     call bd_frez(0,margin,vy,ix,jx,kx)
     call bd_frez(0,margin,vz,ix,jx,kx)
     
     call bd_frez(0,margin,bx,ix,jx,kx)
     call bd_frez(0,margin,by,ix,jx,kx)
     call bd_frez(0,margin,bz,ix,jx,kx)

!!$  call bd_openbz_cyl(0,margin,bxs,bzs,bx,by,bz &
!!$       ,ix,jx,kx,x,dx,dz)

     call bd_frez(0,margin,phi,ix,jx,kx)
     call bd_frez(0,margin,eta,ix,jx,kx)
  end if

  return
end subroutine bndTest
