subroutine bnd(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,x)
  use mpi_domain
  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  type(mpidomain) :: mpid

  real(8),dimension(ix) :: x

  real(8),dimension(ix,jx,kx) :: ro,pr
  real(8),dimension(ix,jx,kx) :: vx,vy,vz
  real(8),dimension(ix,jx,kx) :: bx,by,bz
  real(8),dimension(ix,jx,kx) :: phi,eta

  integer :: i,j,k

!======================================================================
! inner x-boundary

  if(mpid%mpirank_3d(1) .eq.0)then

     call bd_frex(0,margin,ro,ix,jx,kx)
     call bd_frex(0,margin,pr,ix,jx,kx)
     call bd_frex(0,margin,vx,ix,jx,kx)
     call bd_frex(0,margin,vy,ix,jx,kx)
     call bd_frex(0,margin,vz,ix,jx,kx)
     
     call bd_frex(0,margin,bx,ix,jx,kx)
     call bd_frex(0,margin,by,ix,jx,kx)
     call bd_frex(0,margin,bz,ix,jx,kx)
     
     call bd_frex(0,margin,phi,ix,jx,kx)
     call bd_frex(0,margin,eta,ix,jx,kx)

  end if
!======================================================================
! outer x-boudary

  if(mpid%mpirank_3d(1) .eq. (mpid%mpisize_3d(1)-1))then

     call bd_frex(1,margin,ro,ix,jx,kx)
     call bd_frex(1,margin,pr,ix,jx,kx)
     
     call bd_frex(1,margin,vx,ix,jx,kx)
     call bd_frex(1,margin,vy,ix,jx,kx)
     call bd_frex(1,margin,vz,ix,jx,kx)

     call bd_frex(1,margin,bx,ix,jx,kx)
     call bd_frex(1,margin,by,ix,jx,kx)
     call bd_frex(1,margin,bz,ix,jx,kx)

!!$     call bd_openbx_cyl(1,margin,bxs,bzs,bx,by,bz &
!!$          ,ix,jx,kx,x,dx,dz)

     call bd_frex(1,margin,phi,ix,jx,kx)
     call bd_frex(1,margin,eta,ix,jx,kx)
  end if

!======================================================================
! inner x-boundary

  if(mpid%mpirank_3d(2) .eq.0)then

     call bd_frey(0,margin,ro,ix,jx,kx)
     call bd_frey(0,margin,pr,ix,jx,kx)
     call bd_frey(0,margin,vx,ix,jx,kx)
     call bd_frey(0,margin,vy,ix,jx,kx)
     call bd_frey(0,margin,vz,ix,jx,kx)
     
     call bd_frey(0,margin,bx,ix,jx,kx)
     call bd_frey(0,margin,by,ix,jx,kx)
     call bd_frey(0,margin,bz,ix,jx,kx)
     
     call bd_frey(0,margin,phi,ix,jx,kx)
     call bd_frey(0,margin,eta,ix,jx,kx)

  end if
!======================================================================
! outer x-boudary

  if(mpid%mpirank_3d(2) .eq. (mpid%mpisize_3d(2)-1))then

     call bd_frey(1,margin,ro,ix,jx,kx)
     call bd_frey(1,margin,pr,ix,jx,kx)
     
     call bd_frey(1,margin,vx,ix,jx,kx)
     call bd_frey(1,margin,vy,ix,jx,kx)
     call bd_frey(1,margin,vz,ix,jx,kx)

     call bd_frey(1,margin,bx,ix,jx,kx)
     call bd_frey(1,margin,by,ix,jx,kx)
     call bd_frey(1,margin,bz,ix,jx,kx)

     call bd_frey(1,margin,phi,ix,jx,kx)
     call bd_frey(1,margin,eta,ix,jx,kx)
  end if

!----------------------------------------------------------------------|
! outer z-boundary
!


  if(mpid%mpirank_3d(3) .eq. (mpid%mpisize_3d(3)-1)) then  
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
! inner z-boundary
!

  if(mpid%mpirank_3d(3) .eq. 0)then

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

  return
end subroutine bnd
