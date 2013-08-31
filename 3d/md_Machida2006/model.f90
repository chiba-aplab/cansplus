module  model

  implicit none
  private

  public :: model_machida

contains

subroutine  model_machida(igx,jgx,kgx&
       ,ro,pr,vx,vy,vz,bx,by,bz,phi &
       ,roi,pri,vxi,vyi,vzi,bxi,byi,bzi &
       ,x,dx,y,dy,z,dz &
       ,gx,gz,mf_params,eta &
       ,ccx,ccy,ccz )


  use mpi_domain_xz
  use const
  use boundary

  implicit none

! const
!  integer,intent(in) :: ix,jx,kx
  integer,intent(in) :: mf_params

! const
!  integer,intent(in) :: margin
  integer,intent(in) :: igx,jgx,kgx

  type(mpidomain) :: mpid
  integer :: merr 

!  here
!  integer :: idf
!  real(8) :: tend,dtout

! const
!  real(8) :: gm,cr,eta0,vc,xin
! const  
!real(8) :: pi,pi2,pi4,hpi4

!---Input & Output
  real(8),dimension(ix) :: x,dx,dxm,xx
  real(8),dimension(0:ix) :: xm
  real(8),dimension(jx) :: y,dy,dym
  real(8),dimension(0:jx) :: ym
  real(8),dimension(kx) :: z,dz,dzm
  real(8),dimension(0:kx) :: zm
  real(8),dimension(5,2,ix) :: ccx
  real(8),dimension(5,2,jx) :: ccy
  real(8),dimension(5,2,kx) :: ccz

!  real(8),intent(in) :: rg_nrmlx,rohalo,te_factor
!  real(8),intent(in) :: nrmlv,nrmlte,boltzmann_const,Navo,mmw

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

  real(8),dimension(ix,jx,kx),intent(inout) :: gx,gz
  real(8),dimension(ix,jx,kx) :: gpot
  real(8),dimension(igx,jgx,kgx) :: gxg,gzg,gpotg

  real(8),dimension(ix,jx,kx) :: roi,pri
  real(8),dimension(ix,jx,kx) :: vxi,vyi,vzi
  real(8),dimension(ix,jx,kx) :: bxi,byi,bzi
!---Temp coordinate
  real(8),dimension(ix,jx,kx) :: curx,cury,curz

! const  here
  real(8) :: dyg0

! const
!  real(8) :: xmin,ymin,zmin
!  real(8) :: xmax,ymax,zmax

!  const
!  real(8) :: dxmax,dymax,dzmax
!  integer :: ugrid_xmax, ugrid_zmax

  integer :: izero,jzero,kzero

!  const
!  real(8) :: ratio_x, ratio_z

!-------Temp phys
  real(8) :: ss!grav,sseps, ssg   const

  integer :: i,j,k
  integer :: ig,jg,kg

! machida
  real(8) :: tmp !kk,beta0    const

  real(8) :: psi0, pot0 !, tec0  const

  real(8) :: roc,prc,vyc
  real(8) :: rod,prd,vyd,byd

  real(8) :: bbAbsMax

!add
!  cons
!  real(8) :: factorc
!  const
!  real(8) :: aa

!---Step 0.--------------------------------------------------------------|
! set parameter
!

!  const
!  factorc=boltzmann_const*Navo/mmw*nrmlte/nrmlv/nrmlv

!  const
!  vc=0.9d0*2.998d10/nrmlv
!  vc=0.2d0*2.998d10/nrmlv
!  vc=0.5d0*2.998d10/nrmlv

!  const
!  pi = acos(-1.0d0)
!  gm=5.0d0/3.0d0
!  pi2 = 2.0d0*pi
!  cr = 0.18d0
!  pi4 = 4.0d0*pi
!  hpi4 = sqrt(pi4)

! const
!  eta0=4.0d0*pi*0.0001d0

!   const
!  tec0 = 5.0d-1
!  const
!  beta0 = 100.0d0
!  kk = 0.05d0
!  aa=0.2d0

!---Step 1.--------------------------------------------------------------|
! set grid

! const
!  dxg0 = 0.01d0
!  dzg0 = 0.01d0

!  xmin = 0.0d0
!  ymin = 0.0d0
!  ymax = pi2
!  zmin = 0.0d0

!---Step 1a.-------------------------------------------------------------|
! set global x-grid 
!

!  const
!  dxmax = 10.0d0*dxg0
!  ratio_x = 1.05d0
!  ugrid_max = 96

!  do i=margin+1,2*margin
!       dxg(i) = 4.0d0*dxg0
!  enddo
!
!  do i=2*margin+1,2*margin+19
!    dxg(i) = dxg0
!  enddo


!  here
  dxg(margin+1)=4.0d0*dxg0

!  do i=margin+2,margin+96
  do i=margin+2,margin+ugrid_xmax
     dxg(i) = dxg0
  enddo

!  do i=margin+97,igx-margin
  do i=margin+1+ugrid_xmax,igx-margin
     dxg(i) = dxg(i-1)*ratio_x
     if(dxg(i).gt.dxmax) dxg(i)=dxmax
  enddo
  do i=igx-margin+1,igx
     dxg(i)=dxg(igx-margin)
  enddo

  do i=0,margin-1
     dxg(margin-i) = dxg(margin-i+1)
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

!  const
!  dzmax = 10.0d0*dzg0
!  ratio_z = 1.07d0
!  ugrid_zmax=45

  do,k=1,kgx
     dzg = dzg0
  enddo

!  do k=int(kgx/2)+45,kgx
  do k=int(kgx/2)+ugrid_zmax,kgx
     dzg(k) = dzg(k-1)*ratio_z
     if(dzg(k).gt.dzmax) dzg(k)=dzmax
  enddo
  do k=kgx-margin,kgx
     dzg(k)=dzg(kx-margin)
  enddo

!  do k=int(kgx/2)-44,margin,-1
  do k=int(kgx/2)-1-ugrid_zmax,margin,-1
     dzg(k) = dzg(k+1)*ratio_z
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
     ig = mpid%mpirank_2d(1)*(ix-2*margin)+i
     x(i)=xg(ig)
     dx(i) = dxg(ig)
     dxm(i) = dxmg(ig)
  enddo
  do i=0,ix
     ig = mpid%mpirank_2d(1)*(ix-2*margin)+i
     xm(i) = xmg(ig)
  end do

!---Step 2b.-------------------------------------------------------------|
! set individual y-grid 
!
  do j=1,jx
     jg = j
     y(j) = yg(jg)
     dy(j) = dyg(jg)
     dym(j) = dymg(jg)
  enddo

  do j=0,jx
     jg = j
     ym(j) = ymg(jg)
  end do
!---Step 2c.-------------------------------------------------------------|
! set individual z-grid 
!
  do k=1,kx
     kg=mpid%mpirank_2d(2)*(kx-2*margin)+k
     z(k) = zg(kg)
     dz(k) = dzg(kg)
     dzm(k) = dzmg(kg)

  enddo
  do k=0,kx
     kg=mpid%mpirank_2d(2)*(kx-2*margin)+k
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
           gzg(i,j,k) = 0.0d0
        enddo
     enddo
  enddo

!----------------------------------------------------------------------|
! const 
!  ssg=rg_nrmlx
!  sseps = 0.2d0
!  grav = (1.0d0-ssg)**2

!  const
!  xin = sseps
!  xin = 0.2d0   !2rg

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
                 gpotg(i,j,k) =-grav*sseps/((sseps-ssg)**2)*(2.0d0-ss/sseps-ssg/sseps)
                 gxg(i,j,k) = -grav/((sseps-ssg)**2)*(xg(i)/ss)
                 gzg(i,j,k) = -grav/((sseps-ssg)**2)*(zg(k)/ss)
              else
                 gpotg(i,j,k) = -grav/((sseps-ssg)**2)*(1.5d0*sseps-ssg)
                 gxg(i,j,k) = 0.0d0
                 gzg(i,j,k) = 0.0d0
              endif
           endif
        end do
     end do
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---g y
  call bd_pery(margin,gpotg,igx,jgx,kgx)
  call bd_pery(margin,gxg,igx,jgx,kgx)
  call bd_pery(margin,gzg,igx,jgx,kgx)

!---g x
  call bd_synpx(0,margin,gpotg,igx,jgx,kgx)
  call bd_synnx(0,margin,gxg,igx,jgx,kgx)
  call bd_synpx(0,margin,gzg,igx,jgx,kgx)

  call bd_frex(1,margin,gpotg,igx,jgx,kgx)
  call bd_consx(1,margin,0.0d0,gxg,igx,jgx,kgx)
  call bd_frex(1,margin,gzg,igx,jgx,kgx)

!---g z
  call bd_frez(0,margin,gpotg,igx,jgx,kgx)
  call bd_frez(0,margin,gxg,igx,jgx,kgx)
  call bd_consz(0,margin,0.0d0,gzg,igx,jgx,kgx)

  call bd_frez(1,margin,gpotg,igx,jgx,kgx)
  call bd_frez(1,margin,gxg,igx,jgx,kgx)
  call bd_consz(1,margin,0.0d0,gzg,igx,jgx,kgx)

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

  bbAbsMax = 0.0d0
  pot0=-grav/(1.0d0-ssg)
  psi0 = pot0+1.0d0/(2.0d0-2.0d0*aa)+gm*kk/(gm-1.0d0)*(1+1.d0/beta0)
  do k=1,kx
     do j=1,jx
        do i=1,ix
           !---- set corona ---------------------
           ss = sqrt(x(i)**2+z(k)**2)
!           roc = rohalo*exp(-(gpot(i,j,k)-pot0)/tec0)
           roc = rohalo*exp(-(gpot(i,j,k)-pot0)/(factorc*tec0))
           prc = tec0*roc/te_factor

           !---- set torus ----------------------
           rod = 0.0d0
           prd = 0.0d0
           vyd = 0.0d0
           byd = 0.0d0

           if (ss.gt.sseps) then
              tmp = psi0-gpot(i,j,k)-x(i)**(2.0d0*aa-2.0d0)/(2.0d0-2.0d0*aa)
              if(tmp > 0.0d0)then
                 rod = (tmp/(kk*(gm/(gm-1.0d0)) &
                      *(1.0d0+(x(i)**(2.0d0*(gm-1.0d0)))/beta0)))**(1.0d0/(gm-1.0d0))
                 prd = kk*rod**gm
                 vyd = x(i)**(aa-1.0d0)
                 byd = sqrt(rod*2*kk/beta0*(rod*x(i)**2)**(gm-1.0d0))
                 if(dabs(byd) > bbAbsMax)then
                    bbAbsMax = dabs(byd)
                 endif
              endif
           endif


           ro(i,j,k) = rod+roc
           pr(i,j,k) = prd+prc
           vx(i,j,k) = 0.0d0
           vy(i,j,k) = vyd
           vz(i,j,k) = 0.0d0
           bx(i,j,k) = 0.0d0
           by(i,j,k) = byd
           bz(i,j,k) = 0.0d0
           phi(i,j,k) = 0.0d0

           eta(i,j,k) = 0.0d0
        enddo
     enddo
  enddo
           roi = ro
           pri = pr
           vxi = 0.0d0
           vyi = 0.0d0
           vzi = 0.0d0
           bxi = 0.0d0
           byi = by
           bzi = 0.0d0

  call perturb(1,mpid,bx,by,bz,x,dx,dy,dz,bbAbsMax)
       
end subroutine model_machida

!
!----------------------------------------------------------------
!

subroutine perturb(iperturb,mpid,vx,vy,vz,x,dx,dy,dz,v0)
! perturbation
  use mpi_domain_xz
  use const
  implicit none

!  const
!  integer,intent(in) :: ix,jx,kx
  integer,intent(in) :: iperturb

  type(mpidomain) :: mpid

  real(8),intent(in) :: v0

  real(8),dimension(ix) :: x,dx
  real(8),dimension(jx) :: dy
  real(8),dimension(kx) :: dz

  real(8),dimension(ix,jx,kx) :: vx,vy,vz

  real(8),dimension(ix,jx,kx) :: dax,day,daz
  real(8),dimension(ix,jx,kx) :: dvx,dvy,dvz
  real(8),dimension(ix,jx,kx) :: dvxc,dvyc,dvzc

  integer :: n
  integer,dimension(:),allocatable :: seed

  real(8) :: amp
  real(8) :: dtody,dtodz,dtodx
  real(8) :: line1,line2
  real(8) :: temp
  real(8) :: dvxmax,dvymax,dvzmax
  real(8) :: vmax

  integer :: i,j,k

  amp = 0.01d0*v0

  call random_seed(size = n)
  allocate(seed(n))
  seed = 123+mpid%mpirank
  call random_seed(PUT = seed)
  call random_number(dax)
  call random_number(day)
  call random_number(daz)

  if(iperturb .eq. 1)then
     do k=1,kx
        do j=1,jx
           do i=1,ix
              dvx(i,j,k) = 0.0d0
              dvy(i,j,k) = 0.0d0
              dvz(i,j,k) = 0.0d0
           enddo
        enddo
     enddo


     dvxmax = 0.0d0
     do k=2,kx-1
        do j=1,jx
           do i=1,ix
              if(vy(i,j,k) .ne. 0.0d0)then
                 dtodz = 1.0d0/dz(k)
                 dvx(i,j,k) = dvx(i,j,k)-dtodz*(day(i,j,k-1)-day(i,j,k))
                 temp = dabs(dvx(i,j,k))
                 if(temp > dvxmax)then
                    dvxmax = temp
                 endif
              endif
           enddo
        enddo
     enddo

     dvymax = 0.0d0
     do k=2,kx-1
        do j=1,jx
           do i=2,ix-1
              if(vy(i,j,k) .ne. 0.0d0)then
                 dtodz = 1.0d0/dz(k)
                 dtodx = 1.0d0/dx(i)
                 dvy(i,j,k) = dvy(i,j,k)-(dtodz)*(dax(i,j,k)-dax(i,j,k-1)) &
                      -(dtodx)*(-daz(i,j,k)+daz(i-1,j,k))
                 temp = dabs(dvy(i,j,k))
                 if(temp > dvymax)then
                    dvymax = temp
                 endif
              endif
           enddo
        enddo
     enddo

     dvzmax = 0.0d0
     do k=1,kx
        do j=1,jx
           do i=2,ix-1
              if(vy(i,j,k) .ne. 0.0d0)then
                 line1 = dabs(x(i)+0.5d0*dx(i))
                 line2 = dabs(x(i)-0.5d0*dx(i))
                 dtodx = 1.0d0/(x(i)*dx(i))
                 dvz(i,j,k) = dvz(i,j,k)-(dtodx)*(line1*day(i,j,k)-line2*day(i-1,j,k))
                 temp = dabs(dvz(i,j,k))
                 if(temp > dvzmax)then
                    dvzmax = temp
                 endif
              endif
           enddo
        enddo
     enddo

     do k=1,kx
        do j=1,jx
           do i=1,ix
              if(dvxmax .ne. 0.0d0)then
                 dvx(i,j,k) = amp*(dvx(i,j,k)/dvxmax)
              else
                 dvx(i,j,k) = 0.0d0
              endif
              if(dvymax .ne. 0.0d0)then
                 dvy(i,j,k) = amp*(dvy(i,j,k)/dvymax)
              else
                 dvy(i,j,k) = 0.0d0
              endif
              if(dvzmax .ne. 0.0d0)then
                 dvz(i,j,k) = amp*(dvz(i,j,k)/dvzmax)
              else
                 dvz(i,j,k) = 0.0d0
              endif
           enddo
        enddo
     enddo

     do k=1,kx-1
        do j=1,jx-1
           do i=1,ix-1
              dvxc(i,j,k) = 0.5d0*(dvx(i,j,k)+dvx(i-1,j,k))
              dvyc(i,j,k) = dvy(i,j,k)
              dvzc(i,j,k) = 0.5d0*(dvz(i,j,k)+dvz(i,j,k-1))
           enddo
        enddo
     enddo

     do k=1,kx-1
        do j=1,jx-1
           do i=1,ix-1
              vx(i,j,k) = vx(i,j,k) + dvxc(i,j,k)
              vy(i,j,k) = vy(i,j,k) + dvyc(i,j,k)
              vz(i,j,k) = vz(i,j,k) + dvzc(i,j,k)
           enddo
        enddo
     enddo
  else if(iperturb .eq. 2)then
     do k=1,kx
        do j=1,jx
           do i=1,ix
              if(dabs(vy(i,j,k)) > 0.0d0)then
                 vx(i,j,k) = vx(i,j,k) + 2.0d0*amp*(dax(i,j,k)-0.5d0)
                 vy(i,j,k) = vy(i,j,k) + 2.0d0*amp*(day(i,j,k)-0.5d0)
                 vz(i,j,k) = vz(i,j,k) + 2.0d0*amp*(daz(i,j,k)-0.5d0)
              endif
           enddo
        enddo
     enddo
  else
     do k=1,kx
        do j=1,jx
           do i=1,ix
              dvx(i,j,k) = 0.0d0
              dvy(i,j,k) = 0.0d0
              dvz(i,j,k) = 0.0d0
           enddo
        enddo
     enddo
  endif

  return
end subroutine perturb

end module model
