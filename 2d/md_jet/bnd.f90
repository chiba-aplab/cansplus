module bnd

  implicit none
  private

  public :: bnd__exec
  
  
contains


  subroutine bnd__exec(margin,ix,jx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,x)
    
    use mpi_setup, only : mpid, mnull
    use const, only : r_jet,ro_jet,pr_jet,vx_jet,vy_jet,vz_jet,bx_jet,by_jet,bz_jet

    use boundary
    
    integer,intent(in) :: margin,ix,jx
    real(8),dimension(ix),intent(in) :: x
    real(8),dimension(ix,jx),intent(inout) :: ro,pr
    real(8),dimension(ix,jx),intent(inout) :: vx,vy,vz
    real(8),dimension(ix,jx),intent(inout) :: bx,by,bz
    real(8),dimension(ix,jx),intent(inout) :: phi,eta
    real(8) :: pi
    integer :: i,j

!======================================================================
! inter-process communication by MPI
    call boundary__mpi(margin,ix,jx,ro,pr,vx,vy,vz,bx,by,bz,phi)

!======================================================================
! inner x-boundary
    if(mpid%l == mnull)then
       call bd_synpx_car(0,margin,ro,ix,jx)
       call bd_synpx_car(0,margin,pr,ix,jx)
       call bd_synnx_car(0,margin,vx,ix,jx)
       call bd_synnx_car(0,margin,vy,ix,jx)
       call bd_synpx_car(0,margin,vz,ix,jx)
       call bd_synnx_car(0,margin,bx,ix,jx)
       call bd_synnx_car(0,margin,by,ix,jx)
       call bd_synpx_car(0,margin,bz,ix,jx)
       call bd_synpx_car(0,margin,phi,ix,jx)
       call bd_synpx_car(0,margin,eta,ix,jx)
    end if

!======================================================================
! outer x-boudary
    if(mpid%r == mnull)then
       call bd_frex(1,margin,ro,ix,jx)
       call bd_frex(1,margin,pr,ix,jx)
       call bd_frex(1,margin,vx,ix,jx)
       call bd_frex(1,margin,vy,ix,jx)
       call bd_frex(1,margin,vz,ix,jx)
       call bd_frex(1,margin,bx,ix,jx)
       call bd_frex(1,margin,by,ix,jx)
       call bd_frex(1,margin,bz,ix,jx)
       call bd_frex(1,margin,phi,ix,jx)
       call bd_frex(1,margin,eta,ix,jx)
    end if

!----------------------------------------------------------------------|
! free boundary in y
    if(mpid%b == mnull) then
       call bd_synpy_car(0,margin,ro,ix,jx)
       call bd_synpy_car(0,margin,pr,ix,jx)
       call bd_synpy_car(0,margin,vx,ix,jx)
       call bd_synpy_car(0,margin,vy,ix,jx)
       call bd_synny_car(0,margin,vz,ix,jx)
       call bd_synpy_car(0,margin,bx,ix,jx)
       call bd_synpy_car(0,margin,by,ix,jx)
       call bd_synny_car(0,margin,bz,ix,jx)
       call bd_synpy_car(0,margin,phi,ix,jx)
       call bd_synpy_car(0,margin,eta,ix,jx)

       pi = acos(-1d0)
       do j=1,margin
          do i=1,ix
             if (x(i).le.r_jet)then
                ro(i,j) = ro_jet
                pr(i,j) = pr_jet
                vx(i,j) = vx_jet
                vy(i,j) = vy_jet
                vz(i,j) = vz_jet*(sin(x(i)*pi*0.5d0))
                bx(i,j) = bx_jet
                by(i,j) = by_jet
                bz(i,j) = bz_jet
             endif
          enddo
       enddo

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
