subroutine get_temperature(ix,jx,kx,ro,pr,te,te_factor)

  implicit none

  integer,intent(in) :: ix,jx,kx

  real(8),dimension(ix,jx,kx),intent(in) :: ro,pr

  real(8),dimension(ix,jx,kx) :: te
  real(8),intent(in) :: te_factor

  integer :: i,j,k

  do k=1,kx
     do j=1,jx
        do i=1,ix
		        te(i,j,k) = te_factor*pr(i,j,k)/ro(i,j,k)
        end do
     end do
  end do

end subroutine get_temperature
