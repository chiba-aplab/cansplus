module convert
contains
subroutine convert_ptoc_m(ix,jx,kx,gm,ro,pr,vx,vy,vz,bxc,byc,bzc &
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
  implicit none

!---Input
  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: gm

  real(8),dimension(ix,jx,kx),intent(in) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(in) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(in) :: bxc,byc,bzc

!---Output
  real(8),dimension(ix,jx,kx) :: rx,ry,rz,ee

!---Temp
  real(8) :: vsq,pb
  integer :: i,j,k

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

  return
end subroutine convert_ptoc_m

subroutine convert_ctop_m(ix,jx,kx,gm,ro,ee,rx,ry,rz,bxc,byc,bzc,floor &
     ,vx,vy,vz,pr)
!======================================================================
! Name :: convert_ctop
!         convert to primitive
! Input :: 
!          ix,jx,kx :: array size
!          gm :: specific heat retio
!          ro,ee,rx,ry,rz,bxc,byc,bzc
!             :: conserved variables
!          floor :: minimum number
! Output :: 
!           vx,vy,vz,pr :: primitive variable
!======================================================================
  implicit none

!---Input

  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: gm
  real(8),intent(in) :: floor

  real(8),dimension(ix,jx,kx) :: ro,ee
  real(8),dimension(ix,jx,kx) :: rx,ry,rz
  real(8),dimension(ix,jx,kx) :: bxc,byc,bzc

!---Output
  real(8),dimension(ix,jx,kx) :: vx,vy,vz,pr

!---Temp
  real(8) :: rsq,pb,roinverse,vsq

  integer :: i,j,k

! no if
  real(8) :: sign1,sign2,sign3
  real(8) :: temp1,temp2

  do k=1,kx
     do j=1,jx
        do i=1,ix
! negative density cheak

!!$           sign1 = sign(1.0d0,(ro(i,j,k)-floor))
!!$           
!!$           ro(i,j,k) = max(ro(i,j,k),floor)
!!$           roinverse = 1.0d0/ro(i,j,k)
!!$           vx(i,j,k) = max(sign1,0.0d0)*rx(i,j,k)*roinverse
!!$           vy(i,j,k) = max(sign1,0.0d0)*ry(i,j,k)*roinverse
!!$           vz(i,j,k) = max(sign1,0.0d0)*rz(i,j,k)*roinverse
!!$
           if(ro(i,j,k) < floor)then
              ro(i,j,k) = floor
              vx(i,j,k) = 0.0d0
              vy(i,j,k) = 0.0d0
              vz(i,j,k) = 0.0d0
           else
              vx(i,j,k) = rx(i,j,k)/ro(i,j,k)
              vy(i,j,k) = ry(i,j,k)/ro(i,j,k)
              vz(i,j,k) = rz(i,j,k)/ro(i,j,k)
           endif
           pb = 0.5d0*(bxc(i,j,k)**2 + byc(i,j,k)**2 + bzc(i,j,k)**2)
           vsq = vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2
           pr(i,j,k) = (gm-1.0d0)*(ee(i,j,k)-0.5d0*vsq*ro(i,j,k)-pb)
! negative pressure cheak
!!$           sign1 = sign(1.0d0,(pr(i,j,k)-floor))
!!$
!!$           pr(i,j,k) = max(pr(i,j,k),floor)
!!$
!!$           vx(i,j,k) = max(0.0d0,sign1)*vx(i,j,k)
!!$           vy(i,j,k) = max(0.0d0,sign1)*vy(i,j,k)
!!$           vz(i,j,k) = max(0.0d0,sign1)*vz(i,j,k)

!!$           if(pr(i,j,k) < floor)then
!!$              vx(i,j,k) = 0.0d0
!!$              vz(i,j,k) = 0.0d0
!!$              vsq = vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2
!!$              pr(i,j,k) = (gm-1.0d0)*(ee(i,j,k)-0.5d0*vsq*ro(i,j,k)-pb)
!!$           endif
!!$
!!$           if(pr(i,j,k) < floor)then
!!$              vy(i,j,k) = 0.0d0
!!$              pr(i,j,k) = (gm-1.0d0)*(ee(i,j,k)-0.5d0*vsq*ro(i,j,k)-pb)
!!$           endif

           if(pr(i,j,k) < floor)then
              pr(i,j,k) = floor
              vx(i,j,k) = 0.0d0
              vy(i,j,k) = 0.0d0
              vz(i,j,k) = 0.0d0
           endif
           
        enddo
     enddo
  end do

  return
end subroutine convert_ctop_m

end module convert
