subroutine model_AD(idf,mpid,igx,jgx,kgx,ix,jx,kx,margin &
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
  real(8) :: grav, sseps, ss, ssg,sstmp

  integer :: i,j,k
  integer :: ig,jg,kg

  real(8) :: max,min !test

! disk jet params
  real(8) :: betab,k0,h0,r0
  real(8) :: l0,rohalo,cshalo,v0

  real(8) :: psi0,pr0,cs0d,ca0d
  real(8) :: roc,prc,pot0
  real(8) :: rod,prd,tmp
  real(8) :: vyd,byd
!---Step 0.--------------------------------------------------------------|
! set parameter
!

  pi = acos(-1.0d0)
  gm=5.0d0/3.0d0
  cr = 0.18d0
  pi2 = 2.0d0*pi
  pi4 = 4.0d0*pi
  hpi4 = sqrt(pi4)

  r0 = 10.0d0
  k0 = 0.0005d0
  betab = 100.0d0
  rohalo = 0.001d0
  cshalo = 3.0d0/5.0d0

  eta0 = 0.0005d0
  vc = 8.0d0 

!---Step 1.--------------------------------------------------------------|
! set grid
!
  dxg0 = 0.1d0
  dzg0 = 0.1d0

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

  grav = 0.5d0
  sseps = 3.0d0
  ssg = 1.0d0

  xin = 3.0d0
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
                 gpotg(i,j,k) = -grav*(1.0d0/(sseps-ssg) &
                      -(ss-sseps)/(sseps-ssg)**2)
!                 gpotg(i,j,k) = -grav*(ss**2-3.0d0*sseps**2+2.0d0*sseps)/(2.0d0*sseps*(sseps-ssg)**2)
                 gxg(i,j,k) = -(grav/(sseps-ssg)**2)*(xg(i)/ss)
                 gzg(i,j,k) = -(grav/(sseps-ssg)**2)*(zg(k)/ss)
              else
                 ss = 0.5d0*sseps
!                 gpotg(i,j,k) = grav*(ss**2-3.0d0*sseps**2+2.0d0*sseps)/(2.0d0*sseps*(sseps-ssg)**2)
                 gpotg(i,j,k) = -grav*(1.0d0/(sseps-ssg) &
                      -(ss-sseps)/(sseps-ssg)**2)
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

  l0 = sqrt(0.5d0*r0)*r0/(r0-1.0d0)
  v0 = l0/r0
  h0 = 2.0d0*k0/betab/r0**(2.0d0*(gm-1.0d0))
  pr0 = k0
  cs0d = gm*pr0
  ca0d = h0*(r0**2)**(gm-1.0d0)

  pot0 = -0.5d0/(r0-1.0d0)
  psi0 = -0.5d0/(r0-1.0d0)+0.5d0*l0**2/r0**2+cs0d/(gm-1.0d0) &
       +0.5d0*gm*ca0d/(gm-1.0d0)
  
  do k=1,kx
     do j=1,jx
        do i=1,ix
           ss = sqrt(x(i)**2+z(k)**2)
           roc = rohalo*exp(-(gm/cshalo**2)*(gpot(i,j,k)-pot0))
           prc = (cshalo**2/gm)*roc

           rod = 0.0d0
           prd = 0.0d0
           vyd = 0.0d0
           byd = 0.0d0
           if (x(i).gt.sseps) then
              tmp = psi0+0.5d0/(ss-1.0d0)-0.5d0*l0**2/x(i)**2
              if(tmp > 0.0d0)then
                 rod = (tmp/(k0*(gm/(gm-1.0d0)) &
                      *(1.0d0+x(i)**(2.0d0*(gm-1.0d0)) &
                      /(betab*r0**(2.0d0*(gm-1.0d0))))))**(1.0d0/(gm-1.0d0))
                 prd = k0*rod**gm
                 vyd = l0/x(i)
                 byd = sqrt(rod*h0*(rod*x(i)**2)**(gm-1.0d0))
              endif
           endif
           ro(i,j,k) = roc+rod
           pr(i,j,k) = prc+prd
           
           vx(i,j,k) = 0.0d0
           vy(i,j,k) = vyd
           vz(i,j,k) = 0.0d0
           
           bx(i,j,k) = 0.0d0
           by(i,j,k) = byd
           bz(i,j,k) = 0.0d0
        enddo
     enddo
  enddo

!---
! compute current and cal eta
!

  call perturb(2,ix,jx,kx,vx,vy,vz,x,dx,dy,dz,v0)

  call getcurrent_cyl(bx,by,bz,ix,jx,kx,x,dx,dy,dz &
       ,curx,cury,curz)
  call getEta(ix,jx,kx,ro,curx,cury,curz,eta0,vc,eta)

  return
end subroutine model_AD
