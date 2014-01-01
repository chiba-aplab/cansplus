module bnd

  implicit none
  private

  public :: bnd__exec
  
  
contains


  subroutine bnd__exec(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,x,z &
       ,roi,pri,vxi,vyi,vzi,bxi,byi,bzi,phii)
    use mpi_setup, only : mpid, mnull
    use boundary
    
    integer,intent(in) :: margin,ix,jx,kx
    real(8),dimension(ix),intent(in) :: x
    real(8),dimension(kx),intent(in) :: z
    real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr
    real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz
    real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz
    real(8),dimension(ix,jx,kx),intent(inout) :: phi,eta
    real(8),dimension(ix,jx,kx),intent(in) :: roi,pri
    real(8),dimension(ix,jx,kx),intent(in) :: vxi,vyi,vzi
    real(8),dimension(ix,jx,kx),intent(in) :: bxi,byi,bzi,phii
    integer :: i,j,k

!======================================================================
! inter-process communication by MPI
! in this model, x-, y-direction have periodic boundary
    call boundary__mpi(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi)

!----------------------------------------------------------------------|
! outer z-boundary
    if(mpid%t == mnull) then
     call bd_iniz(1,margin,ro,roi,ix,jx,kx)
     call bd_iniz(1,margin,pr,pri,ix,jx,kx)
     call bd_iniz(1,margin,vx,vxi,ix,jx,kx)
     call bd_iniz(1,margin,vy,vyi,ix,jx,kx)
     call bd_iniz(1,margin,vz,vzi,ix,jx,kx)
     call bd_iniz(1,margin,bx,bxi,ix,jx,kx)
     call bd_iniz(1,margin,by,byi,ix,jx,kx)
     call bd_iniz(1,margin,bz,bzi,ix,jx,kx)
     call bd_iniz(1,margin,phi,phii,ix,jx,kx)
    end if

!----------------------------------------------------------------------|
! inner z-boundary
    if(mpid%d == mnull)then
     call bd_iniz(0,margin,ro,roi,ix,jx,kx)
     call bd_iniz(0,margin,pr,pri,ix,jx,kx)
     call bd_iniz(0,margin,vx,vxi,ix,jx,kx)
     call bd_iniz(0,margin,vy,vyi,ix,jx,kx)
     call bd_iniz(0,margin,vz,vzi,ix,jx,kx)
     call bd_iniz(0,margin,bx,bxi,ix,jx,kx)
     call bd_iniz(0,margin,by,byi,ix,jx,kx)
     call bd_iniz(0,margin,bz,bzi,ix,jx,kx)
     call bd_iniz(0,margin,phi,phii,ix,jx,kx)
    end if
  
  end subroutine bnd__exec

end module bnd
