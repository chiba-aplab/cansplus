subroutine model_diskjet(idf,mpid,igx,jgx,kgx,ix,jx,kx,margin &
     ,ro,pr,vx,vy,vz,bx,by,bz,gm &
     ,x,dx,xm,dxm,y,dy,ym,dym,z,dz,zm,dzm &
     ,gx,gz,mf_params,dtout,tend,cr,eta0,vc,eta,xin)
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
  real(8),dimension(ix,jx,kx),intent(inout) :: eta

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
  real(8) :: grav, sseps, ss, ssg

  integer :: i,j,k
  integer :: ig,jg,kg

  real(8) :: max,min !test

! disk jet params
  real(8) :: aa,rn,eth,tec0,roc0,tec00,emg,psi0
  
  real(8) :: roc,prc,rod,prd,vyd,te
  real(8) :: xleft1,xleft2
  real(8) :: zleft1,zleft2

!---Step 0.--------------------------------------------------------------|
! set parameter
!

  pi = acos(-1.0d0)
  gm=5.0d0/3.0d0
  cr = 0.18d0
  pi2 = 2.0d0*pi
  pi4 = 4.0d0*pi
  hpi4 = sqrt(pi4)

  aa=0.d0
  rn=3.d0
  eth=0.05d0
  tec0=1.d0
  roc0=1.d-3
  emg=2.d-3
  eta0 = 0.0d0
  vc = 10.0d0 

!---Step 1.--------------------------------------------------------------|
! set grid
!
  dxg0 = 0.01d0
  dzg0 = 0.01d0

!  xmin = 0.0d0+0.0000001d0
  xmin = 0.0d0
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

  do i=igx/2+margin,igx-margin
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

  dxg(izero) = 3.0d0*dxg0
  dxg(izero-1) = 3.0d0*dxg0
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

  yg(jzero) = ymin
  ymg(jzero) = ymin+0.5d0*dyg(jzero)

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

  do k=kgx/2+100,kgx-margin
     dzg(k) = dzg(k-1)*ratio
     if(dzg(k).gt.dzmax) dzg(k)=dzmax
  enddo
  do k=kgx-margin+1,kgx
     dzg(k)=dzg(kx-margin)
  enddo

  do k=kgx/2-99,margin,-1
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
  grav = 1.0d0
  sseps = 0.2d0
  ssg = 0.0d0

  xin = sseps

  do k=1,kgx
     do j=1,jgx
        do i=1,igx
           ss = sqrt(xg(i)**2+zg(k)**2)
           if( ss .gt. sseps) then
              gpotg(i,j,k) = -grav/(ss-ssg)
              gxg(i,j,k) = -(grav/(ss-ssg)**2)*(xg(i)/ss)
              gzg(i,j,k) = -(grav/(ss-ssg)**2)*(zg(k)/ss)
           else
              if (ss.gt.0.5d0*sseps)then
                 gpotg(i,j,k) = -(grav/sseps-(ss-sseps)/sseps**2)
                 gxg(i,j,k) = -(grav/(sseps-ssg)**2)*(xg(i)/sseps)
                 gzg(i,j,k) = -(grav/(sseps-ssg)**2)*(zg(k)/sseps)
              else
                 gpotg(i,j,k) = -1.5d0/sseps
                 gxg(i,j,k) = 0.0d0
                 gzg(i,j,k) = 0.0d0
              end if
           endif
        enddo
     enddo
  enddo

!--------------------
! gravPot_boundary

!--------------------
! grav_boundary

!----------
! X

!----------
! Z

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

  tec00=tec0*gm
  
  psi0=-1.0d0+1.0d0/2.0d0/(1.0d0-aa)+(rn+1.0d0)*eth
  b0=sqrt(emg)
  
  do k=1,kx
     do j=1,jx
        do i=1,ix
           ss = sqrt(x(i)**2+z(k)**2)
           
           roc=roc0*exp(-gm/tec00*(gpot(i,j,k)+1.0d0))
           prc=roc*tec00/gm
           
           rod = 0.0d0
           prd = 0.0d0
           vyd = 0.0d0
           if (x(i).gt.sseps) then
              te= ( psi0+1.0d0/ss-1.0d0/2.0d0 &
                   / (1.0d0-aa)*x(i)**(2.0d0*aa-2.0d0) ) &
                   / (rn+1.0d0)*gm
              if (te.gt.0.0d0) then
                 rod = (te/gm/eth)**rn
                 prd = te * rod /gm
                 vyd = x(i)**(aa-1.0d0)
              endif
           endif

           ro(i,j,k) = roc+rod
           pr(i,j,k) = prc+prd
           
           vx(i,j,k) = 0.0d0
           vy(i,j,k) = vyd
           vz(i,j,k) = 0.0d0
           
           bx(i,j,k) = 0.0d0
           by(i,j,k) = 0.0d0
           bz(i,j,k) = b0
        enddo
     enddo
  enddo

!---
! compute current and cal eta
!

  call perturb(2,ix,jx,kx,vx,vy,vz,x,dx,dy,dz)

  write(*,*) 'model1'
  call getcurrent_cyl(bx,by,bz,ix,jx,kx,x,dx,dy,dz &
       ,curx,cury,curz)
  write(*,*) 'model2'
  call getEta(ix,jx,kx,ro,curx,cury,curz,eta0,vc,eta)

  return
end subroutine model_diskjet
