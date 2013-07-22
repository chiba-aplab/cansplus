subroutine eta_bnd(ix,jx,kx,x,z,xin,eta,eta0)
  implicit none

  integer,intent(in) :: ix,jx,kx

  real(8),intent(in) :: xin,eta0

  real(8),dimension(ix) :: x
  real(8),dimension(kx) :: z

  real(8),dimension(ix,jx,kx) :: eta

  real(8) :: ss,ssin
  integer :: i,j,k

  ssin = xin + 0.05d0

  do k=1,kx
     do j=1,jx
        do i=1,ix
           ss = sqrt(x(i)**2+z(k)**2)
           if(ss .le. ssin)then
              eta(i,j,k) = eta0
           end if
        end do
     end do
  end do

  return
end subroutine eta_bnd
