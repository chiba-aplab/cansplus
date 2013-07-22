subroutine bd_openbz(mbnd,margin,bxi,byi,bzi,bxc,byc,bzc &
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

  if (mbnd .eq. 0)then
     call bd_frez(0,margin,bxi,ix,jx,kx)
     call bd_frez(0,margin,byi,ix,jx,kx)
     call bd_frez(0,margin,bxc,ix,jx,kx)
     call bd_frez(0,margin,byc,ix,jx,kx)
     kbnd = 1+margin
     do k=1,margin
        do j=2,jx
           do i=2,ix
              dzodx = dz(kbnd-k+1)/dx(i)
              dzody = dz(kbnd-k+1)/dy(j)
              bzi(i,j,kbnd-k) = bzi(i,j,kbnd-k+1) &
                   +dzodx*(bxi(i,j,kbnd-k+1)-bxi(i-1,j,kbnd-k+1)) &
                   +dzody*(byi(i,j,kbnd-k+1)-byi(i,j-1,kbnd-k+1))

              bzc(i,j,kbnd-k+1) = 0.5d0*(bzi(i,j,kbnd-k+1)+bzi(i,j,kbnd-k))
           enddo
        enddo
     enddo
     do j=2,jx
        do i=2,ix
           dzodx = dz(1)/dx(i)
           dzody = dz(1)/dy(j)
           bzi0 = bzi(i,j,1) &
                +dzodx*(bxi(i,j,1)-bxi(i-1,j,1)) &
                +dzody*(byi(i,j,1)-byi(i,j-1,1))
           bzc(i,j,1) = 0.5d0*(bzi0 + bzi(i,j,1))
        enddo
     enddo
! i=1,boundary
     do j=2,jx
        do k=1,margin
           dzodx = dz(kbnd-k+1)/dx(1)
           dzody = dz(kbnd-k+1)/dy(k)
           bzi(1,j,kbnd-k) = bzi(1,j,kbnd-k+1) !&
!                +dzodx*(bxi(2,j,kbnd-k+1)-bxi(1,j,kbnd-k+1)) &
!                +dzody*(byi(1,j,kbnd-k+1)-byi(1,j-1,kbnd-k+1))
           bzc(1,j,kbnd-k+1) = 0.5d0*(bzi(1,j,kbnd-k+1)+bzi(1,j,kbnd-k))
        enddo
        dzodx = dz(1)/dx(1)
        dzody = dz(1)/dy(j)
        bzi0 = bzi(1,j,1) !&
!             +dzodx*(bxi(2,j,1)-bxi(1,j,1)) &
!             +dzody*(byi(1,j,1)-byi(1,j-1,1))
        bzc(1,j,1) = 0.5d0*(bzi0 + bzi(1,j,1))
     enddo
! j=1,boundary
     do i=2,ix
        do k=1,margin
           dzodx = dz(kbnd-k+1)/dx(i)
           dzody = dz(kbnd-k+1)/dy(1)
           bzi(i,1,kbnd-k) = bzi(i,1,kbnd-k+1) !&
!                +dzodx*(bxi(i,1,kbnd-k+1)-bxi(i-1,1,kbnd-k+1))  &
!                +dzody*(byi(i,2,kbnd-k+1)-byi(i,1,kbnd-k+1))
           bzc(i,1,kbnd-k+1) = 0.5d0*(bzi(i,1,kbnd-k+1)+bzi(i,1,kbnd-k))
        enddo
        dzodx = dz(1)/dx(i)
        dzody = dz(1)/dy(1)
        bzi0 = bzi(i,1,1) !&
!             +dzodx*(bxi(i,1,1)-bxi(i-1,1,1)) &
!             +dzody*(byi(i,2,1)-byi(i,1,1))
        bzc(i,1,1) = 0.5d0*(bzi0 + bzi(i,1,1))
     enddo
! i=j=1,boundary
     do k=1,margin
        dzodx = dz(kbnd-k+1)/dx(1)
        dzody = dz(kbnd-k+1)/dy(1)
        bzi(1,1,kbnd-k) = bzi(1,1,kbnd-k+1) !&
!             +dzodx*(bxi(2,1,kbnd-k+1)-bxi(1,1,kbnd-k+1)) &
!             +dzody*(byi(1,2,kbnd-k+1)-byi(1,1,kbnd-k+1))
        bzc(1,1,kbnd-k+1) = 0.5d0*(bzi(1,1,kbnd-k+1)+bzi(1,1,kbnd-k))
     enddo
     dzodx = dz(1)/dx(1)
     dzody = dz(1)/dy(1)
     bzi0 = bzi(1,1,1) !&
!          +dzodx*(bxi(2,1,1)-bxi(1,1,1)) &
!          +dzody*(byi(1,2,1)-byi(1,1,1))
     bzc(1,1,1) = 0.5d0*(bzi0 + bzi(1,1,1))
  else
     call bd_frez(1,margin,bxi,ix,jx,kx)
     call bd_frez(1,margin,byi,ix,jx,kx)
     call bd_frez(1,margin,bxc,ix,jx,kx)
     call bd_frez(1,margin,byc,ix,jx,kx)
     kbnd = kx-margin
     do k=1,margin
        do j=2,jx
           do i=2,ix
              dzodx = dz(kbnd+k)/dx(i)
              dzody = dz(kbnd+k)/dy(j)
              bzi(i,j,kbnd+k) = bzi(i,j,kbnd+k-1) &
                   -dzodx*(bxi(i,j,kbnd+k)-bxi(i-1,j,kbnd+k)) &
                   -dzody*(byi(i,j,kbnd+k)-byi(i,j-1,kbnd+k))
              bzc(i,j,kbnd+k) = 0.5d0*(bzi(i,j,kbnd+k) + bzi(i,j,kbnd+k-1))
           enddo
        enddo
     enddo
! i=1,boundary
     do k=1,margin
        do j=2,jx
           dzodx = dz(kbnd+k)/dx(1)
           dzody = dz(kbnd+k)/dy(j)
           bzi(1,j,kbnd+k) = bzi(1,j,kbnd+k-1) !&
!                -dzodx*(bxi(2,j,kbnd+k) - bxi(1,j,kbnd+k)) &
!                -dzody*(byi(1,j,kbnd+k) - byi(1,j-1,kbnd+k))
           bzc(1,j,kbnd+k) = 0.5d0*(bzi(1,j,kbnd+k) + bzi(1,j,kbnd+k-1))
        enddo
     enddo
! j=1,boundary
     do k=1,margin
        do i=2,ix
           dzodx = dz(kbnd+k)/dx(i)
           dzody = dz(kbnd+k)/dy(1)
           bzi(i,1,kbnd+k) = bzi(i,1,kbnd+k-1) !&
!                -dzodx*(bxi(i,1,kbnd+k)-bxi(i-1,1,kbnd+k)) &
!                -dzody*(byi(i,2,kbnd+k)-byi(i,1,kbnd+k))
           bzc(i,1,kbnd+k) = 0.5d0*(bzi(i,1,kbnd+k) + bzi(i,1,kbnd+k-1))
        enddo
     enddo
! i=j=1,boundary
     do k=1,margin
        dzodx = dz(kbnd+k)/dx(1)
        dzody = dz(kbnd+k)/dy(1)
        bzi(1,1,kbnd+k) = bzi(1,1,kbnd+k-1) !&
!             +dzodx*(bxi(2,1,kbnd+k) - bxi(1,1,kbnd+k)) &
!             +dzody*(byi(1,2,kbnd+k) - byi(1,1,kbnd+k))
        bzc(1,1,kbnd+k) = 0.5d0*(bzi(1,1,kbnd+k) + bzi(1,1,kbnd+k-1))
     enddo
  endif

  return
end subroutine bd_openbz
