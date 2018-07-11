module bnd

  implicit none
  private

  public :: bnd__exec
  
  
contains


  subroutine bnd__exec(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta)
    
    use mpi_setup, only : mpid, mnull
    use boundary
    
    integer,intent(in) :: margin,ix,jx,kx
    real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr
    real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz
    real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz
    real(8),dimension(ix,jx,kx),intent(inout) :: phi,eta
    integer :: i,j,k

!======================================================================
! inter-process communication by MPI
    call boundary__mpi(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta)

!----------------------------------------------------------------------|
! free boundary in x
    if(mpid%l == mnull) then
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
    if(mpid%r == mnull) then
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
! free boundary in y
    if(mpid%b == mnull) then
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
    if(mpid%f == mnull) then
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
! free boundary in z
    if(mpid%d == mnull) then
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


end module bnd
