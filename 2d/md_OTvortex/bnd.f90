module bnd

  implicit none
  private

  public :: bnd__exec
  
  
contains


  subroutine bnd__exec(margin,ix,jx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta)
    
    use mpi_setup, only : mpid
    use boundary, only : boundary__mpi
    
    integer,intent(in) :: margin,ix,jx
    real(8),dimension(ix,jx),intent(inout) :: ro,pr
    real(8),dimension(ix,jx),intent(inout) :: vx,vy,vz
    real(8),dimension(ix,jx),intent(inout) :: bx,by,bz
    real(8),dimension(ix,jx),intent(inout) :: phi,eta
    integer :: i,j

!======================================================================
! inter-process communication by MPI
! periodic boundary in x & y
    call boundary__mpi(margin,ix,jx,ro,pr,vx,vy,vz,bx,by,bz,phi)
  
  end subroutine bnd__exec


end module bnd
