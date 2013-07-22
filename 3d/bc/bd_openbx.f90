subroutine bd_openbx(mbnd,margin,bxi,byi,bzi,bxc,byc,bzc &
     ,ix,jx,kx,dx,dy,dz)
!=====================================================================
! all magnetic-field component boundary
!=====================================================================
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


  if (mbnd .eq. 0)then
     call bd_frex(0,margin,byi,ix,jx,kx)
     call bd_frex(0,margin,bzi,ix,jx,kx)
     call bd_frex(0,margin,byc,ix,jx,kx)
     call bd_frex(0,margin,bzc,ix,jx,kx)
     ibnd = 1+margin
     do k=2,kx
        do j=2,jx
           do i=1,margin
              dxodz = dx(ibnd-i+1)/dz(k)
              dxody = dx(ibnd-i+1)/dy(j)
              bxi(ibnd-i,j,k) = bxi(ibnd-i+1,j,k) &
                  +dxody*(byi(ibnd-i+1,j,k)-byi(ibnd-i+1,j-1,k)) &
                  +dxodz*(bzi(ibnd-i+1,j,k)-bzi(ibnd-i+1,j,k-1))
              bxc(ibnd-i+1,j,k) = 0.5d0*(bxi(ibnd-i+1,j,k)+bxi(ibnd-i,j,k))
           enddo
           dxodz = dx(1)/dz(k)
           dxody = dx(1)/dy(j)
           bxi0 = bxi(1,j,k) &
                +dxody*(byi(1,j,k)-byi(1,j-1,k)) &
                +dxodz*(bzi(1,j,k)-bzi(1,j,k-1))
           bxc(1,j,k) = 0.5d0*(bxi0 + bxi(1,j,k))
        enddo
     enddo
! k=1,boundary
     do j=2,jx
        do i=1,margin
           dxodz = dx(ibnd-i+1)/dz(1)
           dxody = dx(ibnd-i+1)/dy(j)
           bxi(ibnd-i,j,1) = bxi(ibnd-i+1,j,1) !&
!                +dxody*(byi(ibnd-i+1,j,1)-byi(ibnd-i+1,j-1,1)) &
!                +dxodz*(bzi(ibnd-i+1,j,2)-bzi(ibnd-i+1,j,1))
           bxc(ibnd-i+1,j,k) = 0.5d0*(bxi(ibnd-i+1,j,1)+bxi(ibnd-i,j,1))
        enddo
        dxodz = dx(1)/dz(1)
        dxody = dx(1)/dy(j)
        bxi0 = bxi(1,j,1) !&
!             +dxody*(byi(1,j,1)-byi(1,j-1,1)) &
!             +dxodz*(bzi(1,j,2)-bzi(1,j,1))
        bxc(1,j,1) = 0.5d0*(bxi0 + bxi(1,j,1))
     enddo
! j=1,boundary
     do k=2,kx
        do i=1,margin
           dxodz = dx(ibnd-i+1)/dz(k)
           dxody = dx(ibnd-i+1)/dy(1)
           bxi(ibnd-i,1,k) = bxi(ibnd-i+1,1,k) !&
!                +dxody*(byi(ibnd-i+1,2,k)-byi(ibnd-i+1,1,k)) &
!                +dxodz*(bzi(ibnd-i+1,1,k)-bzi(ibnd-i+1,1,k-1))
           bxc(ibnd-i+1,j,k) = 0.5d0*(bxi(ibnd-i+1,1,k)+bxi(ibnd-i,1,k))
        enddo
        dxodz = dx(1)/dz(k)
        dxody = dx(1)/dy(1)
        bxi0 = bxi(1,1,k) !&
!             +dxody*(byi(1,2,k)-byi(1,1,k)) &
!             +dxodz*(bzi(1,1,k)-bzi(1,1,k-1))
        bxc(1,1,k) = 0.5d0*(bxi0 + bxi(1,1,k))
     enddo
! j=k=1,boundary
     do i=1,margin
        dxodz = dx(ibnd-i+1)/dz(1)
        dxody = dx(ibnd-i+1)/dy(1)
        bxi(ibnd-i,1,1) = bxi(ibnd-i+1,1,1) !&
!             +dxody*(byi(ibnd-i+1,2,1)-byi(ibnd-i+1,1,1)) &
!             +dxodz*(bzi(ibnd-i+1,1,2)-bzi(ibnd-i+1,1,1))
        bxc(ibnd-i+1,j,k) = 0.5d0*(bxi(ibnd-i+1,1,1) + bxi(ibnd-i,1,1))
     enddo
     dxodz = dx(1)/dz(1)
     dxody = dy(1)/dy(1)
     bxi0 = bxi(1,1,1) !&
!          +dxody*(byi(1,2,1)-byi(1,1,1)) &
!          +dxodz*(bzi(1,1,2)-bzi(1,1,1))
     bxc(1,1,1) = 0.5d0*(bxi0 + bxi(1,1,1))

  else
     call bd_frex(1,margin,byi,ix,jx,kx)
     call bd_frex(1,margin,bzi,ix,jx,kx)
     call bd_frex(1,margin,byc,ix,jx,kx)
     call bd_frex(1,margin,bzc,ix,jx,kx)
     ibnd = ix-margin
     do k=2,kx
        do j=2,jx
           do i=1,margin
              dxodz = dx(ibnd+i)/dz(k)
              dxody = dx(ibnd+i)/dy(j)
              bxi(ibnd+i,j,k) = bxi(ibnd+i-1,j,k) &
                   - dxody*(byi(ibnd+i,j,k)-byi(ibnd+i,j-1,k)) &
                   - dxodz*(bzi(ibnd+i,j,k)-bzi(ibnd+i,j,k-1))
              bxc(ibnd+i,j,k) = 0.5d0*(bxi(ibnd+i,j,k)+bxi(ibnd+i-1,j,k))
           enddo
        enddo
     enddo
! k=1,boundary
     do j=2,jx
        do i=1,margin
           dxodz = dx(ibnd+i)/dz(k)
           dxody = dx(ibnd+i)/dy(j)
           bxi(ibnd+i,j,1) = bxi(ibnd+i-1,j,1) !&
!                +dxody*(byi(ibnd-i+1,j,1)-byi(ibnd-i+1,j-1,1)) &
!                +dxodz*(bzi(ibnd-i+1,j,2)-bzi(ibnd-i+1,j,1))
           bxc(ibnd+i,j,k) = 0.5d0*(bxi(ibnd+i-1,j,1)+bxi(ibnd+i,j,1))
        enddo
     enddo
! j=1,boundary
     do k=2,kx
        do i=1,margin
           dxodz = dx(ibnd+i)/dz(k)
           dxody = dx(ibnd+i)/dy(1)
           bxi(ibnd+i,1,k) = bxi(ibnd+i-1,1,k) !&
!                +dxody*(byi(ibnd-i+1,2,k)-byi(ibnd-i+1,1,k)) &
!                +dxodz*(bzi(ibnd-i+1,1,k)-bzi(ibnd-i+1,1,k-1))
           bxc(ibnd+i,j,k) = 0.5d0*(bxi(ibnd+i-1,1,k)+bxi(ibnd+i,1,k))
        enddo
     enddo
! j=k=1,boundary
     do i=1,margin
        dxodz = dx(ibnd+i)/dz(1)
        dxody = dx(ibnd+i)/dy(1)
        bxi(ibnd+i,1,1) = bxi(ibnd+i-1,1,1) !&
!             +dxody*(byi(ibnd-i+1,2,1)-byi(ibnd-i+1,1,1)) &
!             +dxodz*(bzi(ibnd-i+1,1,2)-bzi(ibnd-i+1,1,1))
        bxc(ibnd+i,j,k) = 0.5d0*(bxi(ibnd+i-1,1,1) + bxi(ibnd+i,1,1))
     enddo
  endif

  return
end subroutine bd_openbx
