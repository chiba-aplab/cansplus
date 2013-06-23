!
! calcurate magnetic field & divergence B @ cell surface
!
subroutine cal_interface_BP(ix,jx,kx,bxw,phiw &
     ,bx_m,phi_m,ch)
  implicit none

!==========
! INPUT
!========== 
  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: ch

  real(8),dimension(ix,jx,kx,2) :: bxw,phiw

!==========
! OUTPUT
!==========
  real(8),dimension(ix,jx,kx) :: bx_m,phi_m

  integer :: i,j,k

  do k=1,kx
     do j=1,jx
        do i=1,ix
           bx_m(i,j,k) = bxw(i,j,k,1) &
                +(0.5d0*(bxw(i,j,k,2)-bxw(i,j,k,1)) &
                -0.5d0*(phiw(i,j,k,2)-phiw(i,j,k,1))/ch)
           
           phi_m(i,j,k) = phiw(i,j,k,1) &
                +(0.5d0*(phiw(i,j,k,2)-phiw(i,j,k,1)) &
                -0.5d0*ch*(bxw(i,j,k,2)-bxw(i,j,k,1)))
        end do
     end do
  end do
  return
end subroutine cal_interface_BP

  
  
