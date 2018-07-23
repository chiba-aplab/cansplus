module getEta

  implicit none
  private

  public :: getEta__anomalous,getcurrent__cart


contains


  subroutine getEta__anomalous(ix,jx,ro,bx,by,bz,dx,dy,eta0,vc,eta,curx,cury,curz)

  integer,intent(in) :: ix,jx
  real(8),intent(in) :: eta0,vc
  real(8),dimension(ix,jx),intent(in) :: bx,by,bz
  real(8),dimension(ix,jx),intent(in) :: ro
  real(8),dimension(ix),intent(in) :: dx
  real(8),dimension(jx),intent(in) :: dy
  real(8),dimension(ix,jx),intent(out) :: eta,curx,cury,curz

  integer :: i,j
  real(8) :: cur_abs,vd,etamax,flag


! In Cartesian Coordinate
  call getcurrent__cart(bx,by,bz,ix,jx,dx,dy,curx,cury,curz)
  
  etamax=0.01d0
  do j=2,jx-1
     do i=2,ix-1
        cur_abs = sqrt(+curx(i,j)*curx(i,j)+cury(i,j)*cury(i,j) &
                       +curz(i,j)*curz(i,j))
        vd = cur_abs/ro(i,j)

        flag = max(sign(1.0d0,ro(i,j)-1.0d-2),0.0d0)
        eta(i,j) = flag*min(etamax,eta0*(max((vd/vc-1.0d0),0.0d0)**2))
     end do
  end do
  
  end subroutine getEta__anomalous


  subroutine getcurrent__cart(bx,by,bz,ix,jx,dx,dy,curx,cury,curz)

  integer,intent(in) :: ix,jx
  real(8),dimension(ix),intent(in) :: dx
  real(8),dimension(jx),intent(in) :: dy
  real(8),dimension(ix,jx),intent(in) :: bx,by,bz
  real(8),dimension(ix,jx),intent(out):: curx,cury,curz

  integer :: i,j
  real(8) :: dxm,dym
  real(8) :: sign_cur

! X-component
  do j=2,jx-1
     do i=2,ix-1
        dym = 0.5d0*dy(j-1)+dy(j)+0.5d0*dy(j+1)
        curx(i,j) = (bz(i,j+1) - bz(i,j-1))/dym
     end do
  end do

! Y-component
  do j=2,jx-1
     do i=2,ix-1
        dxm = 0.5d0*dx(i-1)+dx(i)+0.5d0*dx(i+1)
        cury(i,j) = - (bz(i+1,j)-bz(i-1,j))/dxm
     end do
  end do

! Z-component
  do j=2,jx-1
     do i=2,ix-1
        dxm = 0.5d0*dx(i-1)+dx(i)+0.5d0*dx(i+1)
        dym = 0.5d0*dy(j-1)+dy(j)+0.5d0*dy(j+1)
        curz(i,j) = ( (by(i+1,j)-by(i-1,j))/dxm &
                     -(bx(i,j+1)-bx(i,j-1))/dym)
     end do
  end do

  end subroutine getcurrent__cart


end module getEta
