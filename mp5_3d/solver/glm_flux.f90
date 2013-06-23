!
! calcurate GLM_MHD flux
!
subroutine glm_flux(bx_m,phi_m,ch,fbx,fphi,ix,jx,kx)
  implicit none

  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: ch ! divergence B wave speed

! cell surface Magnetic field and divergence B
  real(8),dimension(ix,jx,kx) :: bx_m,phi_m

! magnetic field flux @ cell surface of normal component
! divergence B flux
  real(8),dimension(ix,jx,kx) :: fbx,fphi

  integer :: i,j,k

  do k=1,kx
     do j=1,jx
        do i=1,ix
           fbx(i,j,k) = phi_m(i,j,k)
           fphi(i,j,k) = bx_m(i,j,k)*ch**2
        end do
     end do
  end do

  return
end subroutine glm_flux

  
