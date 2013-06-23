subroutine bd_openbx_cyl(mbnd,margin,bxi,byi,bzi,bxc,byc,bzc &
     ,ix,jx,kx,x,dx,dy,dz)

  implicit none

  integer,intent(in) :: mbnd
  integer,intent(in) :: margin
  integer,intent(in) :: ix,jx,kx

  real(8),dimension(ix) :: x,dx
  real(8),dimension(jx) :: dy
  real(8),dimension(kx) :: dz

  real(8),dimension(ix,jx,kx) :: bxi,byi,bzi
  real(8),dimension(ix,jx,kx) :: bxc,byc,bzc

  real(8) :: bxi0

  real(8) :: dxodx,dxody,dxodz

  integer :: ibnd

  integer :: i2,i1
  integer :: i,j,k

  if(mbnd .eq. 0)then
     ibnd = 1+margin
     do k=2,kx
        do j=2,jx
           do i=1,margin
              i2 = ibnd-i+1
              i1 = ibnd-i
              dxodx = (x(i2)+0.5d0*dx(i2))/(x(i2)-0.5d0*dx(i2))
              dxody = dx(i2)/(dy(j)*(x(i2)-0.5d0*dx(i2)))
              dxodz = (dx(i2)/dz(k))*(x(i2)/(x(i2)-0.5d0*dx(i2)))
              bxi(i1,j,k) = dxodx*bxi(i2,j,k) &
                   + dxody*(byi(i2,j,k)-byi(i2,j-1,k)) &
                   + dxodz*(bzi(i2,j,k)-bzi(i2,j,k-1))
              bxc(i2,j,k) = 0.5d0*(bxi(i2,j,k)+bxi(i1,j,k))
           enddo
           dxodx = (x(1)+0.5d0*dx(1))/(x(1)-0.5d0*dx(1))
           dxody = dx(1)/(dy(j)*(x(1)-0.5d0*dx(1)))
           dxodz = (dx(1)/dz(k))*(x(1)/(x(1)-0.5d0*dx(1)))
           bxi0 = dxodx*bxi(1,j,k) &
                + dxody*(byi(1,j,k)-byi(1,j-1,k)) &
                + dxodz*(bzi(1,j,k)-bzi(1,j,k-1))
           bxc(1,j,k) = 0.5d0*(bxi(1,j,k)+bxi0)
        enddo
     enddo
  else
     ibnd = ix-margin
     do k=2,kx
        do j=2,jx
           do i=1,margin
              i1 = ibnd+i
              i2 = ibnd+i-1
              dxodx = (x(i1)-0.5d0*dx(i1))/(x(i1)+0.5d0*dx(i1))
              dxody = dx(i1)/(dy(j)*(x(i1)+0.5d0*dx(i1)))
              dxodz = (dx(i1)/dz(k))*(x(i1)/(x(i1)+0.5d0*dx(i1)))
              bxi(i1,j,k) = dxodx*bxi(i2,j,k) &
                   - dxody*(byi(i1,j,k)-byi(i1,j-1,k)) &
                   - dxodz*(bzi(i1,j,k)-bzi(i1,j,k-1))
              bxc(i1,j,k) = 0.5d0*(bxi(i1,j,k)+bxi(i2,j,k))
           enddo
        enddo
     enddo
  endif
  return
end subroutine bd_openbx_cyl
