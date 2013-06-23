subroutine phiSourceUpdate(ix,jx,kx,ch,cr,phi,dt)
  implicit none

  Integer,intent(in) :: ix,jx,kx

  real(8),intent(in) :: ch,cr
  real(8),intent(in) :: dt

  real(8),dimension(ix,jx,kx) :: phi
  real(8) :: cp

  integer :: i,j,k

  cp = sqrt(ch*cr)

  do k=1,kx
     do j=1,jx
        do i=1,ix
           phi(i,j,k) = phi(i,j,k)*exp(-dt*ch**2/cp**2)
        end do
     end do
  end do

  return 
end subroutine phiSourceUpdate

