subroutine bd_openbz2(mbnd,margin,bxi,byi,bzi,bxc,byc,bzc &
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

  real(8) :: bzi0

  real(8) :: dzodx,dzody

  integer :: kbnd

  integer :: i,j,k

  if ( mbnd .eq. 0)then
     call bd_frez(0,margin,byi,ix,jx,kx)
     call bd_frez(0,margin,byc,ix,jx,kx)
     call bd_frez(0,margin,bxi,ix,jx,kx)
     call bd_frez(0,margin,bxc,ix,jx,kx)
     call bd_frez(0,margin,bzi,ix,jx,kx)
     call bd_frez(0,margin,bzc,ix,jx,kx)
     kbnd = 1+margin
     do k=1,margin
        do j=2,jx
           do i=2,ix
              dzodx = dz(kbnd-k+1)/dx(i)
              dzody = dz(kbnd-k+1)/dy(j)
              bzi(i,j,kbnd-k) = bzi(i,j,kbnd-k+1) &
                   + dzodx*(bxi(i,j,kbnd-k+1)-bxi(i-1,j,kbnd-k+1)) &
                   + dzody*(byi(i,j,kbnd-k+1)-byi(i,j-1,kbnd-k+1))
              bzc(i,j,kbnd-k+1) = 0.5d0*(bzi(i,j,kbnd-k+1)+bzi(i,j,kbnd-k))
           enddo
        enddo
     enddo
  else
     call bd_frez(1,margin,byi,ix,jx,kx)
     call bd_frez(1,margin,byc,ix,jx,kx)
     call bd_frez(1,margin,bxi,ix,jx,kx)
     call bd_frez(1,margin,bxc,ix,jx,kx)
     call bd_frez(1,margin,bzi,ix,jx,kx)
     call bd_frez(1,margin,bzc,ix,jx,kx)
     kbnd = kx-margin
     do k=1,margin
        do j=2,jx
           do i=2,ix
              dzodx = dz(kbnd+k)/dx(i)
              dzody = dz(kbnd+k)/dy(j)
              bzi(i,j,kbnd+k) = bzi(i,j,kbnd+k-1) &
                   - dzodx*(bxi(i,j,kbnd+k)-bxi(i-1,j,kbnd+k)) &
                   - dzody*(byi(i,j,kbnd+k)-byi(i,j-1,kbnd+k))
              bzc(i,j,kbnd+k) = 0.5d0*(bzi(i,j,kbnd+k)+bzi(i,j,kbnd+k-1))
           enddo
        enddo
     enddo
  endif
  return
end subroutine bd_openbz2
