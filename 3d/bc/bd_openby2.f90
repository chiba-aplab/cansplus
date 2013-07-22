subroutine bd_openby2(mbnd,margin,bxi,byi,bzi,bxc,byc,bzc &
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

  real(8) :: byi0

  real(8) :: dyodx, dyodz

  integer :: jbnd

  integer :: i,j,k

  if (mbnd .eq. 0)then
     call bd_frey(0,margin,bxi,ix,jx,kx)
     call bd_frey(0,margin,bxc,ix,jx,kx)
     call bd_frey(0,margin,byi,ix,jx,kx)
     call bd_frey(0,margin,byc,ix,jx,kx)
     call bd_frey(0,margin,bzi,ix,jx,kx)
     call bd_frey(0,margin,bzc,ix,jx,kx)
     jbnd = 1 + margin
     do j=1,margin
        do k=2,kx
           do i=2,ix
              dyodx = dy(jbnd-j+1)/dx(i)
              dyodz = dy(jbnd-j+1)/dz(k)
              byi(i,jbnd-j,k) = byi(i,jbnd-j+1,k) &
                   + dyodx*(bxi(i,jbnd-j+1,k)-bxi(i-1,jbnd-j+1,k)) &
                   + dyodz*(bzi(i,jbnd-j+1,k)-bzi(i,jbnd-j+1,k-1))
              byc(i,jbnd-j+1,k) = 0.5d0*(byi(i,jbnd-j+1,k)+byi(i,jbnd-j,k))
           enddo
        enddo
     enddo
  else
     call bd_frey(1,margin,bxi,ix,jx,kx)
     call bd_frey(1,margin,bxc,ix,jx,kx)
     call bd_frey(1,margin,byi,ix,jx,kx)
     call bd_frey(1,margin,byc,ix,jx,kx)
     call bd_frey(1,margin,bzi,ix,jx,kx)
     call bd_frey(1,margin,bzc,ix,jx,kx)
     jbnd = jx-margin
     do j=1,margin     
        do k=2,kx
           do i=2,ix
              dyodx = dy(jbnd+j)/dx(i)
              dyodz = dy(jbnd+j)/dz(k)
              byi(i,jbnd+j,k) = byi(i,jbnd+j-1,k) &
                   - dyodx*(bxi(i,jbnd+j,k)-bxi(i-1,jbnd+j,k)) &
                   - dyodz*(bzi(i,jbnd+j,k)-bzi(i,jbnd+j,k-1))
              byc(i,jbnd+j,k) = 0.5d0*(byi(i,jbnd+j,k)+byi(i,jbnd+j-1,k))
           enddo
        enddo
     enddo
  endif
  return
end subroutine bd_openby2
