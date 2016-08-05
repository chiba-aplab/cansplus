module model

  use const
  use mpi_setup, only : mpid

  implicit none
  private

  public :: model_setup


contains

  
  subroutine model_setup(ro,pr,vx,vy,vz,bx,by,bz,phi &
                        ,x,dx,xm,y,dy,ym,eta,min_dx)

!---Input & Output
  real(8),dimension(ix)   ,intent(out) :: x,dx
  real(8),dimension(0:ix) ,intent(out) :: xm  
  real(8),dimension(jx)   ,intent(out) :: y,dy
  real(8),dimension(0:jx) ,intent(out) :: ym  
  real(8),dimension(ix,jx),intent(out) :: ro,pr
  real(8),dimension(ix,jx),intent(out) :: vx,vy,vz
  real(8),dimension(ix,jx),intent(out) :: bx,by,bz
  real(8),dimension(ix,jx),intent(out) :: eta,phi
  real(8),                 intent(out) :: min_dx

  integer :: i,j
  integer :: ig,jg
  integer :: izero,jzero
  real(8),dimension(igx) :: xg,dxg
  real(8),dimension(jgx) :: yg,dyg
  real(8),dimension(0:igx) :: xmg
  real(8),dimension(0:jgx) :: ymg

! set global grid 
  do i=1,igx
     dxg(i) = dxg0
  enddo
  do j=1,jgx
     dyg(j) = dyg0
  enddo
! XYZ position
  izero = margin+1
  xmg(izero) = xmin+dxg(izero)
  do i=izero,igx-1
     xmg(i+1) = xmg(i)+dxg(i+1)
  enddo
  do i=izero-1,0,-1
     xmg(i) = xmg(i+1)-dxg(i+1)
  enddo
  do i=1,igx
     xg(i) = 0.5d0*(xmg(i)+xmg(i-1))
  enddo
  jzero = margin+1
  ymg(jzero) = ymin+dyg(jzero)
  do j=jzero,jgx-1
     ymg(j+1) = ymg(j)+dyg(j+1)
  enddo
  do j=jzero-1,0,-1
     ymg(j) = ymg(j+1)-dyg(j+1)
  enddo
  do j=1,jgx
     yg(j) = 0.5d0*(ymg(j)+ymg(j-1))
  enddo

! set individual grid 
  do i=1,ix
     ig = mpid%mpirank_2d(1)*(ix-2*margin)+i
     x(i)=xg(ig)
     dx(i) = dxg(ig)
  enddo
  do j=1,jx
     jg = mpid%mpirank_2d(2)*(jx-2*margin)+j
     y(j) = yg(jg)
     dy(j) = dyg(jg)
  enddo
  do i=0,ix
     ig = mpid%mpirank_2d(1)*(ix-2*margin)+i
     xm(i)=xmg(ig)
  enddo
  do j=0,jx
     jg = mpid%mpirank_2d(2)*(jx-2*margin)+j
     ym(j) = ymg(jg)
  enddo

! calculate min_dx
  min_dx = min(minval(dxg),minval(dyg))
  
!----------------------------------------------------------------------|
! set initial model
  do j=1,jx
     do i=1,ix
        ro(i,j) = gm*gm
        pr(i,j) = gm
        vx(i,j) = -sin(y(j))
        vy(i,j) = sin(x(i))
        vz(i,j) = 0.D0
        bx(i,j) = -sin(y(j))
        by(i,j) = sin(2.*x(i))
        bz(i,j) = 0.d0
        phi(i,j) = 0.d0
        eta(i,j) = 0.d0
     enddo
  enddo

end subroutine model_setup


end module model
