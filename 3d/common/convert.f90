module convert

  implicit none
  private

  public :: convert__ptoc, convert__ctop


contains


  subroutine convert__ptoc(ix,jx,kx,gm,ro,pr,vx,vy,vz,bxc,byc,bzc &
                          ,rx,ry,rz,ee)
!======================================================================
! Name :: convert_ptoc
!         convert to conserved
! Input ::
!          ix,jx,kx :: array size
!          gm :: specific heat retio
!          ro,pr,vx,vy,vz,bxc,byc,bzc
!             :: primitice variables
! Output :: 
!          rx,ry,rz,ee :: conserved varialbles
!
!======================================================================
  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: gm
  real(8),dimension(ix,jx,kx),intent(in) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(in) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(in) :: bxc,byc,bzc
  real(8),dimension(ix,jx,kx),intent(out) :: rx,ry,rz,ee

  integer :: i,j,k
  real(8) :: vsq,pb

  do k=1,kx
     do j=1,jx
        do i=1,ix
           rx(i,j,k) = vx(i,j,k)*ro(i,j,k)
           ry(i,j,k) = vy(i,j,k)*ro(i,j,k)
           rz(i,j,k) = vz(i,j,k)*ro(i,j,k)
           
           pb = 0.5d0*(bxc(i,j,k)**2 + byc(i,j,k)**2 + bzc(i,j,k)**2)
           vsq = vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2
           
           ee(i,j,k) = pr(i,j,k)/(gm-1.0d0)+0.5d0*ro(i,j,k)*vsq + pb
        enddo
     enddo
  end do

  end subroutine convert__ptoc


  subroutine convert__ctop(ix,jx,kx,gm,ro,ee,rx,ry,rz,bxc,byc,bzc &
                          ,vx,vy,vz,pr)
!======================================================================
! Name :: convert_ctop
!         convert to primitive
! Input :: 
!          ix,jx,kx :: array size
!          gm :: specific heat retio
!          ro,ee,rx,ry,rz,bxc,byc,bzc
!             :: conserved variables
! Output :: 
!           vx,vy,vz,pr :: primitive variable
!======================================================================
  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: gm
  real(8),dimension(ix,jx,kx),intent(inout) :: ro
  real(8),dimension(ix,jx,kx),intent(in) :: ee
  real(8),dimension(ix,jx,kx),intent(in) :: rx,ry,rz
  real(8),dimension(ix,jx,kx),intent(in) :: bxc,byc,bzc
  real(8),dimension(ix,jx,kx),intent(out) :: vx,vy,vz,pr

  integer :: i,j,k
  real(8) :: rsq,pb,roinverse,vsq
  real(8) :: sign1,sign2,sign3
  real(8) :: temp1,temp2

  do k=1,kx
     do j=1,jx
        do i=1,ix

           vx(i,j,k) = rx(i,j,k)/ro(i,j,k)
           vy(i,j,k) = ry(i,j,k)/ro(i,j,k)
           vz(i,j,k) = rz(i,j,k)/ro(i,j,k)

           pb = 0.5d0*(bxc(i,j,k)**2 + byc(i,j,k)**2 + bzc(i,j,k)**2)
           vsq = vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2
           pr(i,j,k) = (gm-1.0d0)*(ee(i,j,k)-0.5d0*vsq*ro(i,j,k)-pb)

        enddo
     enddo
  end do

  end subroutine convert__ctop


end module convert
