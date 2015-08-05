module getEta

  implicit none
  private

  public :: getEta__anomalous


contains


  subroutine getEta__anomalous(ix,jx,kx,ro,bx,by,bz,x,dx,dy,dz &
                               ,eta0,vc,eta,curx,cury,curz)
             
  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: eta0,vc
  real(8),dimension(ix,jx,kx),intent(in) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(in) :: ro
  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: dy
  real(8),dimension(kx),intent(in) :: dz
  real(8),dimension(ix,jx,kx),intent(inout) :: eta,curx,cury,curz

  integer :: i,j,k
  real(8) :: cur_abs,vd,etamax,flag


! In Cylindrical coordinate
  call getcurrent(bx,by,bz,ix,jx,kx,x,dx,dy,dz &
                      ,curx,cury,curz)
  
  etamax=0.01d0
  !$OMP PARALLEL DO
  !$OMP PRIVATE(i,j,cur_abs,vd)
  do k=1,kx
     do j=1,jx
        do i=1,ix
           cur_abs = sqrt(+curx(i,j,k)*curx(i,j,k)+cury(i,j,k)*cury(i,j,k) &
                          +curz(i,j,k)*curz(i,j,k))
           vd = cur_abs/ro(i,j,k)
           eta(i,j,k) = min(etamax,eta0*(max((vd/vc-1.0d0),0.0d0)**2))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  end subroutine getEta__anomalous


  subroutine getcurrent(bx,by,bz,ix,jx,kx,x,dx,dy,dz &
                            ,curx,cury,curz)

  integer,intent(in) :: ix,jx,kx
  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: dy
  real(8),dimension(kx),intent(in) :: dz
  real(8),dimension(ix,jx,kx),intent(in) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(out) :: curx,cury,curz

  integer :: i,j,k
  real(8) :: ddx,ddy,ddz
  real(8) :: pi,hpi4,inhpi4
  real(8) :: line1,line2

  pi = acos(-1.0d0)
  hpi4 = sqrt(4.0d0*pi)
  inhpi4 = 1.0d0/hpi4

!-- start OpenMP----------------------------------------------------
  !$OMP PARALLEL
 
! x-component
  !$OMP DO PRIVATE(i,j,ddy,ddz)
  do k=2,kx-1
     do j=2,jx-1
        do i=2,ix-1
           ddy = 0.5d0*dy(j-1)+dy(j)+0.5d0*dy(j+1)
           ddz = 0.5d0*dz(k-1)+dz(k)+0.5d0*dz(k+1)

           curx(i,j,k) = (bz(i,j+1,k)-bz(i,j-1,k))/ddy &
                - (by(i,j,k+1)-by(i,j,k-1))/ddz
        enddo
     enddo
  enddo
  !$OMP END DO NOWAIT

! y-component
  !$OMP DO PRIVATE(i,j,ddx,ddz)
  do k=2,kx-1
     do j=2,jx-1
        do i=2,ix-1
           ddz = 0.5d0*dz(k-1)+dz(k)+0.5d0*dz(k+1)
           ddx = 0.5d0*dx(i-1)+dx(i)+0.5d0*dx(i+1)
           
           cury(i,j,k) = (bx(i,j,k+1)-bx(i,j,k-1))/ddz &
                - (bz(i+1,j,k)-bz(i-1,j,k))/ddx
        enddo
     enddo
  enddo
  !$OMP END DO NOWAIT
 
! z-component
  !$OMP DO PRIVATE(i,j,ddx,ddy)
  do k=2,kx-1
     do j=2,jx-1
        do i=2,ix-1
           ddx = 0.5d0*dx(i-1)+dx(i)+0.5d0*dx(i+1)
           ddy = 0.5d0*dy(j-1)+dy(j)+0.5d0*dy(j+1)

           curz(i,j,k) = (by(i+1,j,k)-by(i-1,j,k))/ddx &
                - (bx(i,j+1,k)-bx(i,j-1,k))/ddy
        enddo
     enddo
  enddo
  !$OMP END DO NOWAIT

!-- end OpenMP
  !$OMP END PARALLEL DO

 
 
 end subroutine getcurrent


end module getEta
