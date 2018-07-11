module bnd

  implicit none
  private

  public :: bnd__exec
  
  
contains


  subroutine bnd__exec(margin,ix,jx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,x,y &
                      ,roi,pri,vxi,vyi,vzi,bxi,byi,bzi,phii)
    use mpi_setup, only : mpid, mnull
    use boundary
    
    integer,intent(in) :: margin,ix,jx
    real(8),dimension(ix),intent(in) :: x
    real(8),dimension(jx),intent(in) :: y
    real(8),dimension(ix,jx),intent(inout) :: ro,pr
    real(8),dimension(ix,jx),intent(inout) :: vx,vy,vz
    real(8),dimension(ix,jx),intent(inout) :: bx,by,bz
    real(8),dimension(ix,jx),intent(inout) :: phi,eta
    real(8),dimension(ix,jx),intent(in) :: roi,pri
    real(8),dimension(ix,jx),intent(in) :: vxi,vyi,vzi
    real(8),dimension(ix,jx),intent(in) :: bxi,byi,bzi,phii
    integer :: i,j

!======================================================================
! inter-process communication by MPI
! periodic boundary in x
    call boundary__mpi(margin,ix,jx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta)

!----------------------------------------------------------------------|
! boundary in z
    if(mpid%b == mnull)then
     call bd_iniy(0,margin,ro,roi,ix,jx)
     call bd_iniy(0,margin,pr,pri,ix,jx)
     call bd_iniy(0,margin,vx,vxi,ix,jx)
     call bd_iniy(0,margin,vy,vyi,ix,jx)
     call bd_iniy(0,margin,vz,vzi,ix,jx)
     call bd_iniy(0,margin,bx,bxi,ix,jx)
     call bd_iniy(0,margin,by,byi,ix,jx)
     call bd_iniy(0,margin,bz,bzi,ix,jx)
     call bd_iniy(0,margin,phi,phii,ix,jx)
    end if
    if(mpid%f == mnull) then
     call bd_iniy(1,margin,ro,roi,ix,jx)
     call bd_iniy(1,margin,pr,pri,ix,jx)
     call bd_iniy(1,margin,vx,vxi,ix,jx)
     call bd_iniy(1,margin,vy,vyi,ix,jx)
     call bd_iniy(1,margin,vz,vzi,ix,jx)
     call bd_iniy(1,margin,bx,bxi,ix,jx)
     call bd_iniy(1,margin,by,byi,ix,jx)
     call bd_iniy(1,margin,bz,bzi,ix,jx)
     call bd_iniy(1,margin,phi,phii,ix,jx)
    end if

  end subroutine bnd__exec

end module bnd
