module getEta

  implicit none
  private

  public :: getEta__anomalous


contains


  subroutine getEta__anomalous(ix,jx,ro,bx,by,bz,x,dx,dy &
                               ,eta0,vc,eta,curx,cury,curz)
             
  integer,intent(in) :: ix,jx
  real(8),intent(in) :: eta0,vc
  real(8),dimension(ix,jx),intent(in) :: bx,by,bz
  real(8),dimension(ix,jx),intent(in) :: ro
  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: dy
  real(8),dimension(ix,jx),intent(inout) :: eta,curx,cury,curz

  integer :: i,j
  real(8) :: cur_abs,vd,etamax,flag


! In Cylindrical coordinate
  call getcurrent(bx,by,bz,ix,jx,kx,x,dx,dy &
                      ,curx,cury,curz)
  
  etamax=0.01d0
  !$OMP PARALLEL DO
  !$OMP PRIVATE(i,j,cur_abs,vd)
     do j=1,jx
        do i=1,ix
           cur_abs = sqrt(+curx(i,j)*curx(i,j)+cury(i,j)*cury(i,j) &
                          +curz(i,j)*curz(i,j))
           vd = cur_abs/ro(i,j)
           eta(i,j) = min(etamax,eta0*(max((vd/vc-1.0d0),0.0d0)**2))
        end do
     end do
  !$OMP END PARALLEL DO
  
  end subroutine getEta__anomalous


  subroutine getcurrent(bx,by,bz,ix,jx,kx,x,dx,dy &
                            ,curx,cury,curz)

  integer,intent(in) :: ix,jx,kx
  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: dy
  real(8),dimension(ix,jx),intent(in) :: bx,by,bz
  real(8),dimension(ix,jx),intent(out) :: curx,cury,curz

  integer :: i,j
  real(8) :: ddx,ddy
  real(8) :: pi,hpi4,inhpi4
  real(8) :: line1,line2

  pi = acos(-1.0d0)
  hpi4 = sqrt(4.0d0*pi)
  inhpi4 = 1.0d0/hpi4

!-- start OpenMP----------------------------------------------------
  !$OMP PARALLEL
 
! x-component
  !$OMP DO PRIVATE(i,j,ddy)
  do j=2,jx-1
     do i=2,ix-1
        ddy = 0.5d0*dy(j-1)+dy(j)+0.5d0*dy(j+1)
        curx(i,j) = inhpi4*(bz(i,j+1)-bz(i,j-1))/ddy 
     enddo
  enddo
  !$OMP END DO NOWAIT

! y-component
  !$OMP DO PRIVATE(i,j,ddx)
  do j=2,jx-1
     do i=2,ix-1
        ddx = 0.5d0*dx(i-1)+dx(i)+0.5d0*dx(i+1)
        cury(i,j) = inhpi4*(-(bz(i+1,j)-bz(i-1,j))/ddx
     enddo
  enddo
  !$OMP END DO NOWAIT
 
! z-component
  !$OMP DO PRIVATE(i,j,ddx,ddy)
  do j=2,jx-1
     do i=2,ix-1
        ddx = 0.5d0*dx(i-1)+dx(i)+0.5d0*dx(i+1)
        ddy = 0.5d0*dy(j-1)+dy(j)+0.5d0*dy(j+1)
        curz(i,j) = inhpi4*(+(by(i+1,j)-by(i-1,j))/ddx &
                            -(bx(i,j+1)-bx(i,j-1))/ddy )
     enddo
  enddo
  !$OMP END DO NOWAIT

!-- end OpenMP
  !$OMP END PARALLEL DO

 
 
 end subroutine getcurrent


end module getEta
