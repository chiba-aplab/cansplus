subroutine model_asahina(idf,mpid,igx,jgx,kgx,ix,jx,kx,margin &
       ,ro,pr,vx,vy,vz,bx,by,bz,phi,gm &
       ,x,dx,y,dy,z,dz &
       ,gx,gy,gz,mf_params,dtout,tend,cr,eta0,vc,eta,xin &
       ,ccx,ccy,ccz)
  use mpi_domain

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
  real(8),dimension(ix) :: x,dx,dxm
  real(8),dimension(0:ix) :: xm  
  real(8),dimension(jx) :: y,dy,dym
  real(8),dimension(0:jx) :: ym  
  real(8),dimension(kx) :: z,dz,dzm
  real(8),dimension(0:kx) :: zm  
  real(8),dimension(5,2,ix) :: ccx
  real(8),dimension(5,2,jx) :: ccy
  real(8),dimension(5,2,kx) :: ccz


  real(8),dimension(igx) :: xg,dxg,dxmg
  real(8),dimension(0:igx) :: xmg
  real(8),dimension(jgx) :: yg,dyg,dymg
  real(8),dimension(0:jgx) :: ymg
  real(8),dimension(kgx) :: zg,dzg,dzmg
  real(8),dimension(0:kgx) :: zmg

  real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz

  real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(inout) :: eta,phi

  real(8),dimension(ix,jx,kx),intent(inout) :: gx,gy,gz
  real(8),dimension(ix,jx,kx) :: gpot
  real(8),dimension(igx,jgx,kgx) :: gxg,gyg,gzg,gpotg

!---Temp coordinate
  real(8),dimension(ix,jx,kx) :: curx,cury,curz

  real(8) :: ratio
  real(8) :: dxg0,dyg0,dzg0

  real(8) :: xmin,ymin,zmin
  real(8) :: xmax,ymax,zmax

  real(8) :: dxmax,dymax,dzmax

  integer :: izero,jzero,kzero,lan,nnx,nc,mdx


!-------Temp phys
  real(8) :: t0,b0,p0

  real(8) :: temp,v2,pb
  real(8) :: grav, sseps, ss, ssg,sstmp

  integer :: i,j,k
  integer :: ig,jg,kg

! asahina
  real(8) :: ro_am,pr_am,te_am
  real(8) :: ro_am1

  real(8) :: rmax,xc
!---Step 0.--------------------------------------------------------------|
! set parameter
!

  pi = acos(-1.0d0)
  gm=5.0d0/3.0d0
  cr = 0.18d0
  pi2 = 2.0d0*pi
  pi4 = 4.0d0*pi
  hpi4 = sqrt(pi4)

  ro_am = 1.6d-24*0.1478709138d0
  ro_am1 = 1.6d-24*6.881726147d0
  te_am = 1.0d4
  pr_am = 6.881726147d0*1.38d-16*2d2

!---Step 1.--------------------------------------------------------------|
! set grid
!
!  xmin = 0.0d0+0.0000001d0
  xmax = 15.0d0
  xmin = 0.0d0
  ymin = -xmax*3d18
  ymax = xmax
  zmin = 0.0d0
  zmax = 60.0d0

  ratio = 1.1d0
  rmax = 5.0d0
  nc = log(rmax)/log(ratio)
  xc =5.0d0
  lan = ix/2-margin
  mdx = ix/2+1

  nnx = (ratio*((1.0d0-rmax)/(1.0d0-ratio))+rmax*real(lan-nc))&
       *xc/(xmax+xc*(rmax-1.0d0))

  dxg0 = 3d18*xc/real(nnx)
  dyg0 = dxg0
  dzg0 = 3d18*(zmax-zmin)/real(kgx-margin*2)

!---Step 1a.-------------------------------------------------------------|
! set global x-grid 
!

  dxmax = rmax*dxg0
 
  do i=mdx,mdx+nnx
     dxg(i) = dxg0
  enddo

  do i=mdx+nnx+1,igx
     dxg(i) = min(dxg(i-1)*ratio,dxmax)
  enddo

  do i=1,mdx-1
     dxg(i)=dxg(igx-i+1)
  enddo

!  do i=igx-margin+1,igx
!     dxg(i)=dxg(igx-margin)
!  enddo

!  do i=1,margin
!     dxg(margin-i+1) = dxg(margin+1)
!  enddo

!  izero = margin+1
  izero = mdx

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
!

!  dyg0 = (ymax-ymin)/real(jgx-margin*2)

  do j=1,jgx
     dyg(j) = dyg0
  enddo

  do j=1,jgx-1
     dymg(j) = dyg(j+1)
  enddo
  dymg(jgx)=dymg(jgx-1)

  jzero = margin+1
  ymg(jzero) = ymin+0.5d0*dyg(jzero)

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
!

  do k=1,kgx
     dzg(k) = dzg0
  enddo

  do k=kgx-margin,kgx
     dzg(k)=dzg(kx-margin)
  enddo

  do k=margin,1,-1
     dzg(k) = dzg(k+1)
  enddo

  do k=1,kgx-1
     dzmg(k) = 0.5d0*(dzg(k+1)+dzg(k))
  enddo

  dzmg(kgx)=dzmg(kgx-1)

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
!
  do i=1,ix
     ig = mpid%mpirank_3d(1)*(ix-2*margin)+i
     x(i)=xg(ig)
     dx(i) = dxg(ig)
     dxm(i) = dxmg(ig)
  enddo
  do i=0,ix
     ig = mpid%mpirank_3d(1)*(ix-2*margin)+i
     xm(i) = xmg(ig)
  end do

!---Step 2b.-------------------------------------------------------------|
! set individual y-grid 
!
  do j=1,jx
     jg = mpid%mpirank_3d(2)*(jx-2*margin)+j
     y(j) = yg(jg)
     dy(j) = dyg(jg)
     dym(j) = dymg(jg)
  enddo

  do j=0,jx
     jg = mpid%mpirank_3d(2)*(jx-2*margin)+j
     ym(j) = ymg(jg)
  end do
!---Step 2c.-------------------------------------------------------------|
! set individual z-grid 
!
  do k=1,kx
     kg=mpid%mpirank_3d(3)*(kx-2*margin)+k
     z(k) = zg(kg)
     dz(k) = dzg(kg)
     dzm(k) = dzmg(kg)
  enddo
  do k=0,kx
     kg=mpid%mpirank_3d(3)*(kx-2*margin)+k
     zm(k) = zmg(kg)
  end do

  call reconstructionConstant(margin,ix,x,xm,dx,ccx)  
  call reconstructionConstant(margin,jx,y,ym,dy,ccy)  
  call reconstructionConstant(margin,kx,z,zm,dz,ccz)  

!---Step 3a.----------------------------------------------------------|
! set gravitation

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

  grav = 1.0d0
  sseps = 0.2d0
  ssg = 1.0d0

  xin = sseps

!---Step 3a.----------------------------------------------------------|
! set individual gravitation

  do k=1,kx
     do j=1,jx
        do i=1,ix
           ig = mpid%mpirank_3d(1)*(ix-2*margin)+i
           jg = mpid%mpirank_3d(1)*(jx-2*margin)+j
           kg = mpid%mpirank_3d(2)*(kx-2*margin)+k
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
           if(z(k) .le. 30.0d0*3d18)then
              pr(i,j,k) = pr_am
              ro(i,j,k) = ro_am
           else
              pr(i,j,k) = pr_am
              ro(i,j,k) = ro_am1
           end if
!!$        if(z(k) > 30.0d0)then
!!$           ro(i,j,k) = 10.0d0*ro_am
!!$        end if
           vx(i,j,k) = 0.0d0
           vy(i,j,k) = 0.0d0
           vz(i,j,k) = 0.0d0
           
!        bx(i,j,k) = 0.0d0
           by(i,j,k) = 0.0d0
!        bz(i,j,k) = 0.0d0
!        r=sqrt(x(i)*x(i) + z(k)*z(k)+ep**2)
           bx(i,j,k)=0d0
           bz(i,j,k)=1d-30
           
           phi(i,j,k) = 0.0d0
           eta(i,j,k) = 0.0d0
        enddo
     enddo
  enddo


!!$  call getcurrent_cyl(bx,by,bz,ix,jx,kx,x,dx,dy,dz &
!!$       ,curx,cury,curz)
!!$  call getEta(ix,jx,kx,ro,curx,cury,curz,eta0,vc,eta)

  return
end subroutine model_asahina
