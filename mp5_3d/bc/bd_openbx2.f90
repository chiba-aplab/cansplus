subroutine bd_openbx2(mbnd,margin,bxi,byi,bzi,bxc,byc,bzc &
     ,ix,jx,kx,dx,dy,dz)

  implicit none

  integer,intent(in) :: mbnd
  integer,intent(in) :: margin
  integer,intent(in) :: ix,jx,kx

  real(8),dimension(ix) :: dx
  real(8),dimension(jx) :: dy
  real(8),dimension(kx) :: dz

  real(8),dimension(ix,jx,kx) :: bxi,byi,bzi
  real(8),dimension(ix,jx,kx) :: bxc,byc,bzc

  real(8) :: bxi0

  real(8) :: dxody,dxodz

  integer :: ibnd

  integer :: i,j,k

  if(mbnd .eq. 0)then
     call bd_frex(0,margin,byi,ix,jx,kx)
     call bd_frex(0,margin,byc,ix,jx,kx)
     call bd_frex(0,margin,bzi,ix,jx,kx)
     call bd_frex(0,margin,bzc,ix,jx,kx)
     call bd_frex(0,margin,bxi,ix,jx,kx)
     call bd_frex(0,margin,bxc,ix,jx,kx)
     ibnd = 1 + margin
     do i=1,margin
        do k=2,kx
           do j=2,jx
              dxody = dx(ibnd-i+1)/dy(j)
              dxodz = dx(ibnd-i+1)/dz(k)
              bxi(ibnd-i,j,k) = bxi(ibnd-i+1,j,k) &
                   + dxody*(byi(ibnd-i+1,j,k)-byi(ibnd-i+1,j-1,k)) &
                   + dxodz*(bzi(ibnd-i+1,j,k)-bzi(ibnd-i+1,j,k-1))
              bxc(ibnd-i+1,j,k) = 0.5d0*(bxi(ibnd-i+1,j,k)+bxi(ibnd-i,j,k))
           enddo
        enddo
     enddo
  else
     call bd_frex(1,margin,byi,ix,jx,kx)
     call bd_frex(1,margin,byc,ix,jx,kx)
     call bd_frex(1,margin,bzi,ix,jx,kx)
     call bd_frex(1,margin,bzc,ix,jx,kx)
     call bd_frex(1,margin,bxi,ix,jx,kx)
     call bd_frex(1,margin,bxc,ix,jx,kx)
     ibnd=ix-margin
     do i=1,margin
        do k=2,kx
           do j=2,jx
              dxody = dx(ibnd+i)/dy(j)
              dxodz = dx(ibnd+i)/dz(k)
              bxi(ibnd+i,j,k) = bxi(ibnd+i-1,j,k) &
                   - dxody*(byi(ibnd+i,j,k)-byi(ibnd+i,j-1,k)) &
                   - dxodz*(bzi(ibnd+i,j,k)-bzi(ibnd+i,j,k-1))
              bxc(ibnd+i,j,k) = 0.5d0*(bxi(ibnd+i,j,k)+bxi(ibnd+i-1,j,k))
           enddo
        enddo
     enddo
  endif
  return
end subroutine bd_openbx2
