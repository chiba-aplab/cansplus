module model

  use const
  use mpi_setup, only : mpid

  implicit none
  private

  public :: model_setup


contains

  
  subroutine model_setup(ro,pr,vx,vy,vz,bx,by,bz,phi         &
                      ,roi,pri,vxi,vyi,vzi,bxi,byi,bzi &
                      ,x,dx,xm,y,dy,ym,z,dz,zm &
                      ,gx,gy,gz,eta,min_dx)

  implicit none

!---Input & Output
  real(8),dimension(ix)      ,intent(out) :: x,dx
  real(8),dimension(0:ix)    ,intent(out) :: xm  
  real(8),dimension(jx)      ,intent(out) :: y,dy
  real(8),dimension(0:jx)    ,intent(out) :: ym  
  real(8),dimension(kx)      ,intent(out) :: z,dz
  real(8),dimension(0:kx)    ,intent(out) :: zm  
  real(8),dimension(ix,jx,kx),intent(out) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(out) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(out) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(out) :: eta,phi
  real(8),dimension(ix,jx,kx),intent(out) :: gx,gy,gz
  real(8),dimension(ix,jx,kx),intent(out) :: roi,pri
  real(8),dimension(ix,jx,kx),intent(out) :: vxi,vyi,vzi
  real(8),dimension(ix,jx,kx),intent(out) :: bxi,byi,bzi
  real(8),                    intent(out) :: min_dx

  integer :: i,j,k
  integer :: ig,jg,kg
  integer :: izero,jzero,kzero
  real(8) :: ss
  real(8) :: dzmax
  real(8),dimension(igx) :: xg,dxg
  real(8),dimension(0:igx) :: xmg
  real(8),dimension(jgx) :: yg,dyg
  real(8),dimension(0:jgx) :: ymg
  real(8),dimension(kgx) :: zg,dzg
  real(8),dimension(0:kgx) :: zmg
  real(8),dimension(ix,jx,kx) :: gpot
  real(8),dimension(igx,jgx,kgx) :: gxg,gyg,gzg,gpotg

!---Step 1a.-------------------------------------------------------------|
! set global x-grid 
!

  do i=1,igx
     dxg(i) = dxg0
  enddo

  izero = margin+1

! X origin
  xmg(izero) = xmin + dxg(izero)
  do i=izero,igx-1
     xmg(i+1) = xmg(i)+dxg(i+1)
  enddo
  do i=izero-1,0,-1
     xmg(i) = xmg(i+1)-dxg(i+1)
  enddo
  do i=1,igx
     xg(i) = 0.5d0*(xmg(i)+xmg(i-1))
  enddo

!---Step 1b.-------------------------------------------------------------|
! set global y-grid
  do j=1,jgx
     dyg(j) = dyg0
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

!---Step 1c.-------------------------------------------------------------|
! set global z-grid
  dzmax = 10.0d0*dzg0
  do,k=1,kgx
     dzg(k) = dzg0
  enddo
  do k=margin+1,kgx
     dzg(k) = dzg(k-1)*ratio_z
     if(dzg(k).gt.dzmax) dzg(k)=dzmax
  enddo
  do k=kgx-margin,kgx
     dzg(k)=dzg(kgx-margin)
  enddo
  do k=margin,1,-1
     dzg(k) = dzg(k+1)
  enddo

! Z origin
  kzero = margin+1
  zmg(kzero) = zmin + dzg(kzero)

  do k=kzero,kgx-1
     zmg(k+1) = zmg(k)+dzg(k+1)
  enddo
  do k=kzero-1,0,-1
     zmg(k) = zmg(k+1)-dzg(k+1)
  enddo
  do k=2,kgx
     zg(k) = 0.5d0*(zmg(k)+zmg(k-1))
  enddo
  zg(1) = 0.5d0*(zmg(1) + (zmg(1)-dzg0))

!---Step 2a.-------------------------------------------------------------|
! set individual x-grid 
  do i=1,ix
     ig = mpid%mpirank_3d(1)*(ix-2*margin)+i
     x(i)=xg(ig)
     dx(i) = dxg(ig)
  enddo
  do i=0,ix
     ig = mpid%mpirank_3d(1)*(ix-2*margin)+i
     xm(i) = xmg(ig)
  end do

!---Step 2b.-------------------------------------------------------------|
! set individual y-grid 
  do j=1,jx
     jg = mpid%mpirank_3d(2)*(jx-2*margin)+j
     y(j) = yg(jg)
     dy(j) = dyg(jg)
  enddo
  do j=0,jx
     jg = mpid%mpirank_3d(2)*(jx-2*margin)+j
     ym(j) = ymg(jg)
  end do

!---Step 2c.-------------------------------------------------------------|
! set individual z-grid 
  do k=1,kx
     kg = mpid%mpirank_3d(3)*(kx-2*margin)+k
     z(k) = zg(kg)
     dz(k) = dzg(kg)
  enddo
  do k=0,kx
     kg = mpid%mpirank_3d(3)*(kx-2*margin)+k
     zm(k) = zmg(kg)
  end do

! calculate min_dx
  min_dx = min(minval(dxg),minval(dyg),minval(dzg))

!---Step 3a.----------------------------------------------------------|
! set gravitation
  !$OMP PARALLEL PRIVATE(i,j,ss,ig,jg,kg)
  !$OMP DO
  do k=1,kgx
     do j=1,jgx
        do i=1,igx
           gpotg(i,j,k) = 0.0d0
           gxg(i,j,k) = 0.0d0
           gyg(i,j,k) = 0.0d0
           gzg(i,j,k) = 0.0d0
        enddo
     enddo
  enddo
  !$OMP END DO 

!----------------------------------------------------------------------|
  !$OMP DO
  do k=1,kgx
     do j=1,jgx
        do i=1,igx
           ss = sqrt(xg(i)**2+zg(k)**2)
           if( ss .gt. sseps) then
              gpotg(i,j,k) = -grav/(ss-ssg)
              gxg(i,j,k) = -(grav/((ss-ssg)**2))*(xg(i)/ss)
              gzg(i,j,k) = -(grav/((ss-ssg)**2))*(zg(k)/ss)
           else 
              if( ss .gt. (sseps*0.5d0))then
                 gpotg(i,j,k) = -grav*sseps/((sseps-ssg)**2)*(2.0d0-ss/sseps-ssg/sseps)
                 gxg(i,j,k) = -grav/((sseps-ssg)**2)*(xg(i)/ss)
                 gzg(i,j,k) = -grav/((sseps-ssg)**2)*(zg(k)/ss)
              else
                 gpotg(i,j,k) = -grav/((sseps-ssg)**2)*(1.5d0*sseps-ssg)
                 gxg(i,j,k) = 0.0d0
                 gyg(i,j,k) = 0.0d0
                 gzg(i,j,k) = 0.0d0
              endif
           endif
        end do
     end do
  end do
  !$OMP END DO

!---Step 3a.----------------------------------------------------------|
! set individual gravitation
  !$OMP DO
  do k=1,kx
     do j=1,jx
        do i=1,ix
           ig = mpid%mpirank_3d(1)*(ix-2*margin)+i
           jg = mpid%mpirank_3d(2)*(jx-2*margin)+j
           kg = mpid%mpirank_3d(3)*(kx-2*margin)+k
           gpot(i,j,k) = gpotg(ig,jg,kg)
           gx(i,j,k) = gxg(ig,jg,kg)
           gy(i,j,k) = gyg(ig,jg,kg)
           gz(i,j,k) = gzg(ig,jg,kg)
        enddo
     enddo
  enddo
  !$OMP END DO

!----------------------------------------------------------------------|
! set initial model
  !$OMP DO
  do k=1,kx
     do j=1,jx
        do i=1,ix
           ro(i,j,k) = ro_am
           pr(i,j,k) = pr_am
           vx(i,j,k) = vx_am
           vy(i,j,k) = vy_am
           vz(i,j,k) = vz_am
           bx(i,j,k) = bx_am
           by(i,j,k) = by_am
           bz(i,j,k) = bz_am
           phi(i,j,k) = 0d0
           eta(i,j,k) = 0d0
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  roi(1:ix,1:jx,1:kx) = ro(1:ix,1:jx,1:kx)
  pri(1:ix,1:jx,1:kx) = pr(1:ix,1:jx,1:kx)
  vxi(1:ix,1:jx,1:kx) = vx(1:ix,1:jx,1:kx)
  vyi(1:ix,1:jx,1:kx) = vy(1:ix,1:jx,1:kx)
  vzi(1:ix,1:jx,1:kx) = vz(1:ix,1:jx,1:kx)
  bxi(1:ix,1:jx,1:kx) = bx(1:ix,1:jx,1:kx)
  byi(1:ix,1:jx,1:kx) = by(1:ix,1:jx,1:kx)
  bzi(1:ix,1:jx,1:kx) = bz(1:ix,1:jx,1:kx)

  return
end subroutine model_setup

end module model
