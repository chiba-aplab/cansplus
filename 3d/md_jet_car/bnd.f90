module bnd

  implicit none
  private

  public :: bnd__exec
  
  
contains


  subroutine bnd__exec(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,x,y,z)
    
    use mpi_setup, only : mpid, mnull
    use const, only : r_jet,ro_jet,pr_jet,vx_jet,vy_jet,vz_jet,bx_jet,by_jet,bz_jet,pi
    use boundary
    
    integer,intent(in) :: margin,ix,jx,kx
    real(8),dimension(ix),intent(in) :: x
    real(8),dimension(jx),intent(in) :: y
    real(8),dimension(kx),intent(in) :: z
    real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr
    real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz
    real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz
    real(8),dimension(ix,jx,kx),intent(inout) :: phi,eta
    real(8) :: b0s,si,co,r
    integer :: i,j,k

!======================================================================
! inter-process communication by MPI
    call boundary__mpi(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi)

!======================================================================
! inner x-boundary
    if(mpid%l == mnull)then
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
    if(mpid%r == mnull)then
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
! outer z-boundary
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

!----------------------------------------------------------------------|
! inner z-boundary
    if(mpid%d == mnull)then
       call bd_synpz_car(0,margin,ro,ix,jx,kx)
       call bd_synpz_car(0,margin,pr,ix,jx,kx)
       call bd_synpz_car(0,margin,vx,ix,jx,kx)
       call bd_synpz_car(0,margin,vy,ix,jx,kx)
       call bd_synnz_car(0,margin,vz,ix,jx,kx)
       call bd_synpz_car(0,margin,bx,ix,jx,kx)
       call bd_synpz_car(0,margin,by,ix,jx,kx)
       call bd_synnz_car(0,margin,bz,ix,jx,kx)
       call bd_synpz_car(0,margin,phi,ix,jx,kx)
       call bd_synpz_car(0,margin,eta,ix,jx,kx)
       
       do k=1,margin
          do j=1,jx
             do i=1,ix
                r = sqrt(x(i)*x(i)+y(j)*y(j))
                if (r.le.r_jet)then
                   si = y(j)/r
                   co = x(i)/r
                   b0s = by_jet*sin(r*pi)**4d0
                   ro(i,j,k) = ro_jet
                   pr(i,j,k) = pr_jet
                   vx(i,j,k) = vx_jet
                   vy(i,j,k) = vy_jet
                   vz(i,j,k) = vz_jet
                   bx(i,j,k) = -b0s * si
                   by(i,j,k) = b0s * co
                   bz(i,j,k) = bz_jet
                endif
             enddo
          enddo
       enddo
    end if
  
  end subroutine bnd__exec


end module bnd
