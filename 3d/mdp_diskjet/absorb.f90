subroutine absorb(ix,jx,kx,x,z,xin,qq,qqi)
  implicit none

  integer,intent(in) :: ix,jx,kx

  real(8),intent(in) :: xin

  real(8),dimension(ix) :: x
  real(8),dimension(kx) :: z

  real(8),dimension(ix,jx,kx) :: qq,qqi

  real(8) :: ai,halfx,ss,dx0

  integer :: i,j,k

  dx0 = 0.01d0
  do k=1,kx
     do j=1,jx
        do i=1,ix
           ss = sqrt(x(i)**2+z(k)**2)
           if(ss .le. xin)then
              ai = 0.1d0*(1.0d0-tanh((ss-xin+5.0d0*dx0)/(2.0d0*dx0)))
              qq(i,j,k) = (1.0d0-ai)*qq(i,j,k)+ai*qqi(i,j,k)
           endif
        enddo
     end do
  end do
  return
end subroutine absorb
