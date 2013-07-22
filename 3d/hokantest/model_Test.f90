subroutine model_Test(idf,mpid,igx,jgx,kgx,ix,jx,kx,margin &
       ,ro,pr,vx,vy,vz,bx,by,bz,phi,gm &
       ,x,dx,xm,dxm,y,dy,ym,dym,z,dz,zm,dzm &
       ,gx,gz,mf_params,dtout,tend,cr,eta0,vc,eta,xin)
  use convertCoordinate
  use mpi_domain_xz

  implicit none
  include 'mpif.h'

  integer,intent(in) :: ix,jx,kx
  integer,intent(in) :: mf_params

  integer,intent(in) :: margin
  integer,intent(in) :: igx,jgx,kgx
  
  type(mpidomain) :: mpid

  integer :: idf
  real(8) :: tend,dtout

  real(8) :: gm,cr,eta0,vc,xin
  real(8) :: pi,pi2,pi4,hpi4

!---Input & Output
  real(8),dimension(ix) :: x,dx,xm,dxm
  real(8),dimension(jx) :: y,dy,ym,dym
  real(8),dimension(kx) :: z,dz,zm,dzm

  real(8),dimension(igx) :: xg,dxg,xmg,dxmg
  real(8),dimension(jgx) :: yg,dyg,ymg,dymg
  real(8),dimension(kgx) :: zg,dzg,zmg,dzmg

  real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz

  real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(inout) :: eta,phi

  real(8),dimension(ix,jx,kx),intent(inout) :: gx,gz
  real(8),dimension(ix,jx,kx) :: gpot
  real(8),dimension(igx,jgx,kgx) :: gxg,gzg,gpotg

!---Temp coordinate
  real(8),dimension(ix,jx,kx) :: curx,cury,curz

  real(8) :: dxg0,dyg0,dzg0

  real(8) :: xmin,ymin,zmin
  real(8) :: xmax,ymax,zmax

  real(8) :: dxmax,dymax,dzmax

  integer :: izero,jzero,kzero

  real(8) :: ratio

!-------Temp phys
  real(8) :: t0,b0,p0

  real(8) :: temp,v2,pb
  real(8) :: grav, sseps, ss, ssg,sstmp

  integer :: i,j,k
  integer :: ig,jg,kg

  real(8) :: max,min !test

! disk jet params
  real(8) :: beta

!----test1012
  real(8) :: x1,x2,x3

  real(8) :: roCar, prCar
  real(8) :: v1Car, v2Car, v3Car
  real(8) :: b1Car, b2Car, b3Car

  real(8) :: ro0,pr0
  real(8) :: v10,v20,v30
  real(8) :: b10,b20,b30

  real(8) :: ro1,pr1
  real(8) :: v11,v21,v31
  real(8) :: b11,b21,b31
  
  real(8) :: wtr,smoothFactor
!----test1012

!---Step 0.--------------------------------------------------------------|
! set parameter
!

  pi = acos(-1.0d0)
  gm=5.0d0/3.0d0
  cr = 0.18d0
  pi2 = 2.0d0*pi
  pi4 = 4.0d0*pi
  hpi4 = sqrt(pi4)

  beta = 1.0d0

!----test1012
  wtr = 0.2d0
!!$!----left
!!$  ro0 = 1.08d0
!!$  pr0 = 0.95d0
!!$  v10 = 1.2d0
!!$  v20 = 0.01d0
!!$  v30 = 0.5d0
!!$  b10 = 2.0d0/hpi4
!!$  b20 = 3.6d0/hpi4
!!$  b30 = 2.0d0/hpi4
!!$!-----right
!!$  ro1 = 1.0d0
!!$  pr1 = 1.0d0
!!$  v11 = 0.0d0
!!$  v21 = 0.0d0
!!$  v31 = 0.0d0
!!$  b11 = 2.0d0/hpi4
!!$  b21 = 4.0d0/hpi4
!!$  b31 = 2.0d0/hpi4
!----Model2
  ro1 = 1.0d0
  pr1 = 1.0d0/gm
  v11 = -1.0d0
  v21 = 0.0d0
  v31 = 0.0d0
  b11 = 2.0d0*sqrt(pr1/beta)
  b21 = 0.0d0
  b31 = 0.0d0

  ro0 = 0.1d0
  pr0 = pr1
  v10 = v11
  v20 = 0.0d0
  v30 = 0.0d0
  b10 = b11
  b20 = 0.0d0
  b30 = 0.0d0



!---Step 1.--------------------------------------------------------------|
! set grid
!
  dxg0 = 0.01d0
  dzg0 = 0.01d0

  xmin = 0.0d0
!  xmin = 1.6d0
  ymin = 0.0d0
  ymax = pi2
  zmin = 0.0d0

!---Step 1a.-------------------------------------------------------------|
! set global x-grid 
!
  dxmax = 10.0d0*dxg0
  ratio = 1.05d0
  do i=margin,igx
     dxg(i) = dxg0
  enddo

  do i=100+margin,igx-margin
     dxg(i) = dxg(i-1)*ratio
     if(dxg(i).gt.dxmax) dxg(i)=dxmax
  enddo
  do i=igx-margin+1,igx
     dxg(i)=dxg(igx-margin)
  enddo

  do i=1,margin-1
     dxg(margin-i) = dxg(margin-i+1)
  enddo

  izero = margin+1

! X origin

  dxg(izero) = 4.0d0*dxg0
  dxg(izero-1) = 4.0d0*dxg0
  xg(izero) = xmin+0.5d0*dxg(izero)
  xmg(izero) = xmin

  do i=1,igx-1
     dxmg(i) = 0.5d0*(dxg(i)+dxg(i+1))
  enddo
  dxmg(igx)=dxmg(igx-1)

  do i=izero+1,igx
     xmg(i) = xmg(i-1)+dxg(i-1)
  enddo
  do i=izero-1,1,-1
     xmg(i) = xmg(i+1)-dxg(i)
  enddo

  do i=izero+1,igx
     xg(i) = xg(i-1)+dxmg(i-1)
  enddo
  do i=izero-1,1,-1
     xg(i) = xg(i+1)-dxmg(i)
  enddo

!---Step 1b.-------------------------------------------------------------|
! set global y-grid
!

  dyg0 = (ymax-ymin)/real(jgx-margin*2)

  do j=1,jgx
     dyg(j) = dyg0
  enddo

  do j=1,jgx-1
     dymg(j) = dyg(j+1)
  enddo
  dymg(jgx)=dymg(jgx-1)
  jzero = margin+1


!!$  yg(jzero) = ymin
!!$  ymg(jzero) = ymin+0.5d0*dyg(jzero)
!----test1012
  yg(jzero) = ymin+0.5d0*dyg(jzero)
  ymg(jzero) = ymin+dyg(jzero)
!----test1012

  do j=jzero+1,jgx
     ymg(j) = ymg(j-1)+dyg(j-1)
  enddo
  do j=jzero-1,1,-1
     ymg(j) = ymg(j+1)-dyg(j)
  enddo

  do j=jzero+1,jgx
     yg(j) = yg(j-1)+dymg(j-1)
  enddo
  do j=jzero-1,1,-1
     yg(j) = yg(j+1)-dymg(j)
  enddo


!---Step 1c.-------------------------------------------------------------|
! set global z-grid
!

  dzmax = 10.0d0*dzg0
  ratio = 1.05d0

  do k=1,kgx
     dzg(k) = dzg0
  enddo

  do k=kgx/2+50,kgx-margin
     dzg(k) = dzg(k-1)*ratio
     if(dzg(k).gt.dzmax) dzg(k)=dzmax
  enddo
  do k=kgx-margin,kgx
     dzg(k)=dzg(kx-margin)
  enddo

  do k=kgx/2-49,margin,-1
     dzg(k) = dzg(k+1)*ratio
     if(dzg(k).gt.dzmax) dzg(k)=dzmax
  end do
  do k=margin,1,-1
     dzg(k) = dzg(k+1)
  enddo

  do k=1,kgx-1
     dzmg(k) = 0.5d0*(dzg(k+1)+dzg(k))
  enddo

  dzmg(kgx)=dzmg(kgx-1)

! Z origin

  kzero = kgx/2

  zg(kzero) = zmin-0.5d0*dzg(kzero)
  zmg(kzero) = zmin

  do k=kzero+1,kgx
     zmg(k) = zmg(k-1)+dzg(k-1)
  enddo
  do k=kzero-1,1,-1
     zmg(k) = zmg(k+1)-dzg(k)
  enddo

  do k=kzero+1,kgx
     zg(k) = zg(k-1)+dzmg(k-1)
  enddo
  do k=kzero-1,1,-1
     zg(k) = zg(k+1)-dzmg(k)
  enddo


!---Step 2a.-------------------------------------------------------------|
! set individual x-grid 
!
  do i=1,ix
     ig = mpid%mpirank_2d(1)*(ix-2*margin)+i
     x(i)=xg(ig)
     dx(i) = dxg(ig)
     dxm(i) = dxmg(ig)
     xm(i) = xmg(ig)
  enddo

!---Step 2b.-------------------------------------------------------------|
! set individual y-grid 
!
  do j=1,jx
     jg = j
     y(j) = yg(jg)
     dy(j) = dyg(jg)
     dym(j) = dymg(jg)
     ym(j) = ymg(jg)
  enddo

!---Step 2c.-------------------------------------------------------------|
! set individual z-grid 
!
  do k=1,kx
     kg=mpid%mpirank_2d(2)*(kx-2*margin)+k
     z(k) = zg(kg)
     dz(k) = dzg(kg)
     dzm(k) = dzmg(kg)
     zm(k) = zmg(kg)
  enddo

!---Step 3a.----------------------------------------------------------|
! set gravitation

  do k=1,kgx
     do j=1,jgx
        do i=1,igx
           gpotg(i,j,k) = 0.0d0
           gxg(i,j,k) = 0.0d0
           gzg(i,j,k) = 0.0d0
        enddo
     enddo
  enddo

!---Step 3a.----------------------------------------------------------|
! set individual gravitation

  do k=1,kx
     do j=1,jx
        do i=1,ix
           ig = mpid%mpirank_2d(1)*(ix-2*margin)+i
           jg = j
           kg = mpid%mpirank_2d(2)*(kx-2*margin)+k
           gpot(i,j,k) = gpotg(ig,jg,kg)
           gx(i,j,k) = gxg(ig,jg,kg)
           gz(i,j,k) = gzg(ig,jg,kg)
        enddo
     enddo
  enddo

!----------------------------------------------------------------------|
! set initial mocel
! 
  do k=1,kx
     do j=1,jx
        do i=1,ix

!           smoothFactor = 0.5d0*(1.0d0+tanh((x(i)*cos(y(j)))/wtr))

!!$           roCar = ro0 + (ro1-ro0)*smoothFactor
!!$           prCar = pr0 + (pr1-pr0)*smoothFactor
!!$           v1Car = v10 + (v11-v10)*smoothFactor
!!$           v2Car = v20 + (v21-v20)*smoothFactor
!!$           v3Car = v30 + (v31-v30)*smoothFactor
!!$           b1Car = b10 + (b11-b10)*smoothFactor
!!$           b2Car = b20 + (b21-b20)*smoothFactor
!!$           b3Car = b30 + (b31-b30)*smoothFactor
           call CylinderToCartesian(x(i),y(j),z(k),x1,x2,x3)
           if(sqrt((x1-0.2d0)**2+x2**2) < 0.1d0)then
              roCar = ro1
              prCar = pr1
              v1Car = v11
              v2Car = v21
              v3Car = v31
              b1Car = b11
              b2Car = b21
              b3Car = b31
           else
              roCar = ro0
              prCar = pr0
              v1Car = v10
              v2Car = v20
              v3Car = v30
              b1Car = b10
              b2Car = b20
              b3Car = b30
           end if
           ro(i,j,k) = roCar
           pr(i,j,k) = prCar
           vx(i,j,k) = v1Car*cos(y(j)) + v2Car*sin(y(j))
           vy(i,j,k) = -v1Car*sin(y(j)) + v2Car*cos(y(j))
           vz(i,j,k) = v3Car
           bx(i,j,k) = b1Car*cos(y(j)) + b2Car*sin(y(j))
           by(i,j,k) = -b1Car*sin(y(j)) + b2Car*cos(y(j))
           bz(i,j,k) = b3Car
           
           phi(i,j,k) = 0.0d0
        enddo
     enddo
  enddo

end subroutine model_Test
