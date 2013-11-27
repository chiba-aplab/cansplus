module getEta

  implicit none
  private

  public :: getEta__anomalous


contains


  subroutine getEta__anomalous(ix,jx,kx,ro,bx,by,bz,dx,dy,dz,eta0,vc,eta,curx,cury,curz)

  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: eta0,vc
  real(8),dimension(ix,jx,kx),intent(in) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(in) :: ro
  real(8),dimension(ix),intent(in) :: dx
  real(8),dimension(jx),intent(in) :: dy
  real(8),dimension(kx),intent(in) :: dz
  real(8),dimension(ix,jx,kx),intent(out) :: eta,curx,cury,curz

  integer :: i,j,k
  real(8) :: cur_abs,vd,etamax,flag


! In Cartesian Coordinate
  call getcurrent__cart(bx,by,bz,ix,jx,kx,dx,dy,dz &
                       ,curx,cury,curz)
  
  etamax=0.01d0
  do k=2,kx-1
     do j=2,jx-1
        do i=2,ix-1
           cur_abs = sqrt(+curx(i,j,k)*curx(i,j,k)+cury(i,j,k)*cury(i,j,k) &
                          +curz(i,j,k)*curz(i,j,k))
           vd = cur_abs/ro(i,j,k)

           flag = max(sign(1.0d0,ro(i,j,k)-1.0d-2),0.0d0)
           eta(i,j,k) = flag*min(etamax,eta0*(max((vd/vc-1.0d0),0.0d0)**2))
        end do
     end do
  end do
  
  end subroutine getEta__anomalous


  subroutine getcurrent__cart(bx,by,bz,ix,jx,kx,dx,dy,dz &
                             ,curx,cury,curz)

  integer,intent(in) :: ix,jx,kx
  real(8),dimension(ix),intent(in) :: dx
  real(8),dimension(jx),intent(in) :: dy
  real(8),dimension(kx),intent(in) :: dz
  real(8),dimension(ix,jx,kx),intent(in) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(out):: curx,cury,curz

  integer :: i,j,k
  real(8) :: dxm,dym,dzm
  real(8) :: hpi4,pi,inhpi4
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

  end subroutine getcurrent__cart


end module getEta
