!
! calculate current @ cartesian grid.
! 
subroutine getcurrent(ix,jx,kx,bx,by,bz,dx,dy,dz &
     ,curx,cury,curz)

  implicit none

!==========
! INPUT
!==========
  integer,intent(in) :: ix,jx,kx

  real(8),dimension(ix) :: dx
  real(8),dimension(jx) :: dy
  real(8),dimension(kx) :: dz
  
  real(8),dimension(ix,jx,kx) :: bx,by,bz

!==========
! OUTPUT
!==========

  real(8),dimension(ix,jx,kx) :: curx,cury,curz

  real(8) :: dxm,dym,dzm
  real(8) :: hpi4,pi,inhpi4

  integer :: i,j,k

  real(8) :: sign_cur

  pi = acos(-1.0d0)
  hpi4 = sqrt(4.0d0*pi)
  inhpi4 = 1.0d0/hpi4

! X-component

  do k=2,kx-1
     do j=2,jx-1
        do i=2,ix-1
           dzm = 0.5d0*dz(k-1)+dz(k)+0.5d0*dz(k+1)
           dym = 0.5d0*dy(j-1)+dy(j)+0.5d0*dy(j+1)

           curx(i,j,k) = inhpi4*((bz(i,j+1,k) - bz(i,j-1,k))/dym &
                -(by(i,j,k+1)-by(i,j,k-1))/dzm)

        end do
     end do
  end do

! Y-component

  do k=2,kx-1
     do j=2,jx-1
        do i=2,ix-1
           dxm = 0.5d0*dx(i-1)+dx(i)+0.5d0*dx(i+1)
           dzm = 0.5d0*dz(k-1)+dz(k)+0.5d0*dz(k+1)

           cury(i,j,k) = inhpi4*((bx(i,j,k+1)-bx(i,j,k-1))/dzm &
                - (bz(i+1,j,k)-bz(i-1,j,k))/dxm)
        end do
     end do
  end do

! Z-component

  do k=2,kx-1
     do j=2,jx-1
        do i=2,ix-1
           dxm = 0.5d0*dx(i-1)+dx(i)+0.5d0*dx(i+1)
           dym = 0.5d0*dy(j-1)+dy(j)+0.5d0*dy(j+1)

           curz(i,j,k) = inhpi4*((by(i+1,j,k)-by(i-1,j,k))/dxm &
                - (bx(i,j+1,k)-bx(i,j-1,k))/dym)
        end do
     end do
  end do

  return
end subroutine getcurrent

