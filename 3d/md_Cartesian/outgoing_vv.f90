subroutine outgoing_vv(ix,jx,kx,x,z,vx,vz,xin)
  implicit none

  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: xin

  real(8),dimension(ix) :: x
  real(8),dimension(kx) :: z

  real(8),dimension(ix,jx,kx) :: vx,vz

  real(8) :: ss

  integer :: i,j,k

  do k=1,kx
     do j=1,jx
        do i=1,ix
           ss = sqrt(x(i)**2+z(k)**2)
           if(ss < xin)then
              vx(i,j,k) = min(0.0d0,vx(i,j,k))
              if(z(k) > 0.0d0)then
                 vz(i,j,k) = min(0.0d0,vz(i,j,k))
              else
                 vz(i,j,k) = max(0.0d0,vz(i,j,k))
              endif
           end if
        end do
     end do
  end do

  return
end subroutine outgoing_vv
