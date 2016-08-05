module convert

  implicit none
  private

  public :: convert__ptoc, convert__ctop


contains


  subroutine convert__ptoc(ix,jx,gm,ro,pr,vx,vy,vz,bxc,byc,bzc &
                          ,rx,ry,rz,ee)
!======================================================================
! Name :: convert_ptoc
!         convert to conserved
! Input ::
!          ix,jx :: array size
!          gm :: specific heat retio
!          ro,pr,vx,vy,vz,bxc,byc,bzc
!             :: primitice variables
! Output :: 
!          rx,ry,rz,ee :: conserved varialbles
!
!======================================================================
  integer,intent(in) :: ix,jx
  real(8),intent(in) :: gm
  real(8),dimension(ix,jx),intent(in) :: ro,pr
  real(8),dimension(ix,jx),intent(in) :: vx,vy,vz
  real(8),dimension(ix,jx),intent(in) :: bxc,byc,bzc
  real(8),dimension(ix,jx),intent(out) :: rx,ry,rz,ee

  integer :: i,j
  real(8) :: vsq,pb

  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j,pb,vsq)
  do j=1,jx
     do i=1,ix
        rx(i,j) = vx(i,j)*ro(i,j)
        ry(i,j) = vy(i,j)*ro(i,j)
        rz(i,j) = vz(i,j)*ro(i,j)
           
        pb = 0.5d0*(bxc(i,j)**2 + byc(i,j)**2 + bzc(i,j)**2)
        vsq = vx(i,j)**2 + vy(i,j)**2 + vz(i,j)**2
           
        ee(i,j) = pr(i,j)/(gm-1.0d0)+0.5d0*ro(i,j)*vsq + pb
     enddo
  enddo
  !$OMP END PARALLEL DO

  end subroutine convert__ptoc


  subroutine convert__ctop(ix,jx,gm,ro,ee,rx,ry,rz,bxc,byc,bzc &
                          ,vx,vy,vz,pr)
!======================================================================
! Name :: convert_ctop
!         convert to primitive
! Input :: 
!          ix,jx :: array size
!          gm :: specific heat retio
!          ro,ee,rx,ry,rz,bxc,byc,bzc
!             :: conserved variables
! Output :: 
!           vx,vy,vz,pr :: primitive variable
!======================================================================
  integer,intent(in) :: ix,jx
  real(8),intent(in) :: gm
  real(8),dimension(ix,jx),intent(inout) :: ro,pr
  real(8),dimension(ix,jx),intent(inout) :: ee
  real(8),dimension(ix,jx),intent(in) :: rx,ry,rz
  real(8),dimension(ix,jx),intent(in) :: bxc,byc,bzc
  real(8),dimension(ix,jx),intent(out) :: vx,vy,vz

  integer :: i,j
  real(8), parameter :: pbeta_min=1d-3
  real(8) :: vsq,pb,roinverse,igm,temppr,signpr,temp1,temp2
  
  igm = 1d0/(gm-1d0)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j,temppr,roinverse,pb,vsq,signpr,temp1,temp2)
  do j=1,jx
     do i=1,ix
        temppr = pr(i,j)
        roinverse = 1d0/ro(i,j)
        vx(i,j) = rx(i,j) * roinverse
        vy(i,j) = ry(i,j) * roinverse
        vz(i,j) = rz(i,j) * roinverse

        pb = 0.5d0*(bxc(i,j)**2 + byc(i,j)**2 + bzc(i,j)**2)
        vsq = vx(i,j)**2 + vy(i,j)**2 + vz(i,j)**2
        pr(i,j) = (gm-1.0d0)*(ee(i,j)-0.5d0*vsq*ro(i,j)-pb)
        signpr = sign(1d0,pr(i,j))
        temp1 = max(0d0,signpr)
        temp2 = -min(0d0,signpr)

        pr(i,j) = max(pbeta_min*pb,temp1*pr(i,j)+temp2*temppr)
        ee(i,j) = pr(i,j)*igm + 0.5d0*vsq*ro(i,j) +pb
     enddo
  enddo
  !$OMP END PARALLEL DO
  end subroutine convert__ctop


end module convert
