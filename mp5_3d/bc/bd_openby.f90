subroutine bd_openby(mbnd,margin,bxi,byi,bzi,bxc,byc,bzc &
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

  real(8) :: dyodx,dyodz

  integer :: jbnd
  integer :: i,j,k

  if (mbnd .eq. 0)then
     call bd_frey(0,margin,bxi,ix,jx,kx)
     call bd_frey(0,margin,bzi,ix,jx,kx)
     call bd_frey(0,margin,bxc,ix,jx,kx)
     call bd_frey(0,margin,bzc,ix,jx,kx)
     jbnd = 1+margin
     do k=2,kx
        do j=1,margin
           do i=2,ix
              dyodx = dy(jbnd-j+1)/dx(i)
              dyodz = dz(jbnd-j+1)/dz(k)
              byi(i,jbnd-j,k) = byi(i,jbnd-j+1,k) &
                   +dyodx*(bxi(i,jbnd-j+1,k)-bxi(i-1,jbnd-j+1,k)) &
                   +dyodz*(bzi(i,jbnd-j+1,k)-bzi(i,jbnd-j+1,k-1))

              byc(i,jbnd-j+1,k) = 0.5d0*(byi(i,jbnd-j+1,k)+byi(i,jbnd-j,k))
           enddo
        enddo
     enddo
     do k=2,kx
        do i=2,ix
           dyodx = dy(1)/dx(i)
           dyodz = dy(1)/dz(k)
           byi0 = byi(i,1,k) &
                +dyodx*(bxi(i,1,k)-bxi(i-1,1,k)) &
                +dyodz*(bzi(i,1,k)-bzi(i,1,k-1))
           byc(i,1,k) = 0.5d0*(byi0 + byi(i,1,k))
        enddo
     enddo
! i=1,boundary
     do k=2,kx
        do j=1,margin
           dyodx = dy(jbnd-j+1)/dx(1)
           dyodz = dy(jbnd-j+1)/dz(k)
           byi(1,jbnd-j,k) = byi(1,jbnd-j+1,k) !&
!                +dyodx*(bxi(2,jbnd-j+1,k)-bxi(1,jbnd-j+1,k)) &
!                +dyodz*(bzi(1,jbnd-j+1,k)-bzi(1,jbnd-j+1,k-1))
           byc(1,jbnd-j+1,k) = 0.5d0*(byi(1,jbnd-j+1,k)+byi(1,jbnd-j,k))
        enddo
        dyodx = dy(1)/dx(1)
        dyodz = dy(1)/dz(k)
        byi0 = byi(1,1,k) !&
!             +dyodx*(bxi(2,1,k)-bxi(1,1,k)) &
!             +dyodz*(bzi(1,1,k)-bzi(1,1,k-1))
        byc(1,1,k) = 0.5d0*(byi0 + byi(1,1,k))
     enddo
! k=1,boundary
     do i=2,ix
        do j=1,margin
           dyodx = dy(jbnd-j+1)/dx(i)
           dyodz = dy(jbnd-j+1)/dz(1)
           byi(i,jbnd-j,1) = byi(i,jbnd-j+1,1) !&
!                + dyodx*(bxi(i,jbnd-j+1,1)-bxi(i-1,jbnd-j+1,1)) &
!                + dyodz*(bzi(i,jbnd-j+1,2)-bzi(i,jbnd-j+1,1))
           byc(i,jbnd-j+1,1) = 0.5d0*(byi(i,jbnd-j+1,1)+byi(i,jbnd-j,1))
        enddo
        dyodx = dy(1)/dx(i)
        dyodz = dy(1)/dz(1)
        byi0 = byi(i,1,1) !&
!             +dyodx*(bxi(i,1,1)-bxi(i-1,1,1)) &
!             +dyodz*(bzi(i,1,2)-bzi(i,1,1))
        byc(i,1,1) = 0.5d0*(byi0 + byi(i,1,1))
     enddo

! i=k=1,boundary
     do j=1,margin
        dyodx = dy(jbnd-j+1)/dx(1)
        dyodz = dy(jbnd-j+1)/dz(1)
        byi(1,jbnd-j,1) = byi(1,jbnd-j+1,1) !&
!             +dyodx*(bxi(2,jbnd-j+1,1)-bxi(1,jbnd-j+1,1)) &
!             +dyodz*(bzi(1,jbnd-j+1,2)-bzi(1,jbnd-j+1,1))
        byc(1,jbnd-j+1,1) = 0.5d0*(byc(1,jbnd-j+1,1)+byc(1,jbnd-j,1))
     enddo
     dyodx = dy(1)/dx(1)
     dyodz = dy(1)/dz(1)
     byi0 = byi(1,1,1) !&
!          +dyodx*(bxi(2,1,1)-bxi(1,1,1)) &
!          +dyodz*(bzi(1,1,2)-bzi(1,1,1))
     byc(1,1,1) = 0.5d0*(byi0 + byi(1,1,1))

  else
     call bd_frey(1,margin,bxi,ix,jx,kx)
     call bd_frey(1,margin,bzi,ix,jx,kx)
     call bd_frey(1,margin,bxc,ix,jx,kx)
     call bd_frey(1,margin,bzc,ix,jx,kx)
     jbnd = jx-margin
     do k=2,kx
        do j=1,margin
           do i=2,ix
              dyodz = dy(jbnd+j)/dz(k)
              dyodx = dy(jbnd+j)/dx(i)
              byi(i,jbnd+j,k) = byi(i,jbnd+j-1,k) &
                   -dyodx*(bxi(i,jbnd+j,k) - bxi(i-1,jbnd+j,k)) &
                   -dyodz*(bzi(i,jbnd+j,k) - bzi(i,jbnd+j,k-1))
              byc(i,jbnd+j,k) = 0.5d0*(byi(i,jbnd+j,k) + byi(i,jbnd+j-1,k))
           enddo
        enddo
     enddo
! k=1,boundary
     do j=1,margin
        do i=2,ix
           dyodz = dy(jbnd+j)/dz(1)
           dyodx = dy(jbnd+j)/dx(i)
           byi(i,jbnd+j,1) = byi(i,jbnd+j-1,1)! &
!                -dyodx*(bxi(i,jbnd+j,1) - bxi(i-1,jbnd+j,1)) &
!                -dyodz*(bzi(i,jbnd+j,2) - bzi(i,jbnd+j,1))
           byc(i,jbnd+j,1) = 0.5d0*(byi(i,jbnd+j,1) + byi(i,jbnd+j-1,1))
        enddo
     enddo
! i=1,boundary
     do k=2,kx
        do j=1,margin
           dyodz = dy(jbnd+j)/dz(k)
           dyodz = dy(jbnd+j)/dx(1)
           byi(1,jbnd+j,k) = byi(1,jbnd+j-1,k) !&
!                -dyodx*(bxi(2,jbnd+j,k) - bxi(1,jbnd+j,k)) &
!                -dyodz*(bzi(1,jbnd+j,k) - bzi(1,jbnd+j,k-1))
           byc(1,jbnd+j,k) = 0.5d0*(byi(1,jbnd+j,k) + byi(1,jbnd+j-1,k))
        enddo
     enddo
! i=k=1,boundary
     do j=1,margin
        dyodz = dy(jbnd+j)/dz(1)
        dyodx = dy(jbnd+j)/dx(1)
        byi(1,jbnd+j,1) = byi(1,jbnd+j-1,1) !&
!             +dyodx*(bxi(2,jbnd+j,1) - bxi(1,jbnd+j,1)) &
!             +dyodz*(bzi(1,jbnd+j,2) - bzi(1,jbnd+j,1))
        byc(1,jbnd+j,1) = 0.5d0*(byi(1,jbnd+j,1) + byi(1,jbnd+j-1,1))
     enddo
  endif

  return
end subroutine bd_openby
