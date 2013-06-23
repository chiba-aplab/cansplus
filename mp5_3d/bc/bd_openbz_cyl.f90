subroutine bd_openbz_cyl(mbnd,margin,bxi,byi,bzi,bxc,byc,bzc &
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

  real(8) :: bzi0

  real(8) :: dzodx,dzody
  real(8) :: xright,xleft

  integer :: kbnd
  
  integer :: k1,k2
  integer :: i,j,k
  
  if(mbnd.eq.0)then
     kbnd = 1+margin
     do k=1,margin
        do j=2,jx
           do i=2,ix
              k2 = kbnd-k+1
              k1 = kbnd-k
              xright = (x(i)+0.5d0*dx(i))/x(i)
              xleft = (x(i)-0.5d0*dx(i))/x(i)
              dzodx = dz(k2)/dx(i)
              dzody = dz(k2)/(x(i)*dy(j))
              bzi(i,j,k1) = bzi(i,j,k2) &
                   + dzodx*(xright*bxi(i,j,k2)-xleft*bxi(i,j,k2)) &
                   + dzody*(byi(i,j,k2)-byi(i,j-1,k2)) 
              bzc(i,j,k2) = 0.5d0*(bzi(i,j,k1)+bzi(i,j,k2))
           enddo
        enddo
     enddo
     do j=2,jx
        do i=2,ix
           xright = (x(i)+0.5d0*dx(i))/x(i)
           xleft = (x(i)-0.5d0*dx(i))/x(i)
           dzodx = dz(1)/dx(i)
           dzody = dz(1)/(x(i)*dy(j))
           bzi0 = bzi(i,j,1) &
                + dzodx*(xright*bxi(i,j,1)-xleft*bxi(i,j,1)) &
                + dzody*(byi(i,j,1)-byi(i,j-1,1))
           bzc(i,j,1) = 0.5d0*(bzi(i,j,1)+bzi0)
        enddo
     enddo
  else
     kbnd = kx-margin
     do k=1,margin
        do j=2,jx
           do i=2,ix
              k2 = kbnd+k
              k1 = kbnd+k-1
              xright = (x(i)+0.5d0*dx(i))/x(i)
              xleft = (x(i)-0.5d0*dx(i))/x(i)
              dzodx = dz(k2)/dx(i)
              dzody = dz(k2)/(x(i)*dy(j))
              bzi(i,j,k2) = bzi(i,j,k1) &
                   - dzodx*(xright*bxi(i,j,k2)-xleft*bxi(i,j,k2)) &
                   - dzody*(byi(i,j,k2)-byi(i,j-1,k2)) 
              bzc(i,j,k2) = 0.5d0*(bzi(i,j,k1)+bzi(i,j,k2))
           enddo
        enddo
     enddo
  endif
  return
end subroutine bd_openbz_cyl
