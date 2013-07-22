subroutine calTest(margin,ix,jx,kx,x,dx,y,dy,z,dz,qq,qqval)
  implicit none

  integer,intent(in) :: margin,ix,jx,kx

  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: y,dy
  real(8),dimension(kx),intent(in) :: z,dz

  real(8),dimension(ix,jx,kx),intent(in) :: qq

  real(8),intent(inout) :: qqval

  real(8) :: rr,dv

  integer :: i,j,k

  qqval = 0.0d0
  do k=margin+1,kx-margin
     do j=margin+1,jx-margin
        do i=margin+1,ix-margin
           rr = sqrt(x(i)**2)
           dv = rr*dx(i)*dy(j)*dz(k)
           qqval = qqval + qq(i,j,k)*dv
        enddo
     enddo
  enddo

end subroutine calTest
