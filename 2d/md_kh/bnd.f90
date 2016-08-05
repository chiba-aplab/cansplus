module bnd

  implicit none
  private

  public :: bnd__exec
  
  
contains


  subroutine bnd__exec(margin,ix,jx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta)
    
    use mpi_setup, only : mpid, mnull
    use boundary
    
    integer,intent(in) :: margin,ix,jx
    real(8),dimension(ix,jx),intent(inout) :: ro,pr
    real(8),dimension(ix,jx),intent(inout) :: vx,vy,vz
    real(8),dimension(ix,jx),intent(inout) :: bx,by,bz
    real(8),dimension(ix,jx),intent(inout) :: phi,eta

!======================================================================
! inter-process communication by MPI
! periodic boundary in x
    call boundary__mpi(margin,ix,jx,ro,pr,vx,vy,vz,bx,by,bz,phi)

!----------------------------------------------------------------------|
! free boundary in y
    if(mpid%b == mnull) then
       call bd_frey(0,margin,ro,ix,jx)
       call bd_frey(0,margin,pr,ix,jx)
       call bd_frey(0,margin,vx,ix,jx)
       call bd_frey(0,margin,vy,ix,jx)
       call bd_frey(0,margin,vz,ix,jx)
       call bd_frey(0,margin,bx,ix,jx)
       call bd_frey(0,margin,by,ix,jx)
       call bd_frey(0,margin,bz,ix,jx)
       call bd_frey(0,margin,phi,ix,jx)
       call bd_frey(0,margin,eta,ix,jx)
    end if
    if(mpid%f == mnull) then
       call bd_frey(1,margin,ro,ix,jx)
       call bd_frey(1,margin,pr,ix,jx)
       call bd_frey(1,margin,vx,ix,jx)
       call bd_frey(1,margin,vy,ix,jx)
       call bd_frey(1,margin,vz,ix,jx)
       call bd_frey(1,margin,bx,ix,jx)
       call bd_frey(1,margin,by,ix,jx)
       call bd_frey(1,margin,bz,ix,jx)
       call bd_frey(1,margin,phi,ix,jx)
       call bd_frey(1,margin,eta,ix,jx)
    end if
  
  end subroutine bnd__exec


end module bnd
