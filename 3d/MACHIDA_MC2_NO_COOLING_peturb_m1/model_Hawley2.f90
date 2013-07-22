subroutine model_Hawley2(idf,mpid,igx,jgx,kgx,ix,jx,kx,margin &
     ,ro,pr,vx,vy,vz,bx,by,bz,phi,gm &
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
  real(8),dimension(ix,jx,kx),intent(inout) :: eta,phi
  real(8),dimension(ix,jx,kx) :: bxs,bys,bzs

  real(8),dimension(ix,jx,kx),intent(inout) :: gx,gz
  real(8),dimension(ix,jx,kx) :: gpot,gpot_edge
  real(8),dimension(igx,jgx,kgx) :: gxg,gzg,gpotg,gpotg_edge

!---Temp coordinate
  real(8) :: ro_ed,ro_edge_chk,ssr,ssz
  real(8) :: line1,line2
  real(8) :: inbeta,inbetamax,alpha,inbetamaxg
  real(8) :: romin,dtodz,dtodx
  integer :: merr

  real(8),dimension(ix,jx,kx) :: ro_edge
  real(8),dimension(ix,jx,kx) :: curx,cury,curz
  real(8),dimension(ix,jx,kx) :: aphi
  real(8),dimension(ix,jx,kx) :: bxs_temp,bzs_temp

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
  real(8) :: betab,k0,h0,rkep,rin
  real(8) :: l0,rohalo,cshalo,v0

  real(8) :: psi0,pr0,cs0d,ca0d
  real(8) :: roc,prc,pot0
  real(8) :: rod,prd,tmp,tmp2,rod0
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

  romin = 0.005d0

  rkep = 6.5d0
  rin = 4.5d0
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
  xmin = 1.6d0
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

  grav = 0.5d0
  sseps = 3.0d0
  ssg = 1.0d0

  xin = 2.0d0
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
           ss = sqrt((xg(i)+0.5d0*dxg(i))**2+(zg(k)+0.5d0*dzg(k))**2)
           if( ss .gt. sseps) then
              gpotg_edge(i,j,k) = -grav/(ss-ssg)
           else
              if (ss.gt.0.5d0*sseps)then
                 gpotg_edge(i,j,k) = -grav*(1.0d0/(sseps-ssg) &
                      -(ss-sseps)/(sseps-ssg)**2)
              else
                 ss = 0.5d0*sseps
                 gpotg_edge(i,j,k) = -grav*(1.0d0/(sseps-ssg) &
                      -(ss-sseps)/(sseps-ssg)**2)
              end if
           end if
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
           gpot_edge(i,j,k)=gpotg_edge(ig,jg,kg)
           gpot(i,j,k) = gpotg(ig,jg,kg)
           gx(i,j,k) = gxg(ig,jg,kg)
           gz(i,j,k) = gzg(ig,jg,kg)
        enddo
     enddo
  enddo

!----------------------------------------------------------------------|
! set initial mocel
! 


  l0 = sqrt(0.5d0*rkep)*rkep/(rkep-1.0d0)
  v0 = l0/rkep
  h0 = 2.0d0*k0/betab/rkep**(2.0d0*(gm-1.0d0))
  pr0 = k0
  cs0d = gm*pr0
  ca0d = h0*(rkep**2)**(gm-1.0d0)

  pot0 = -0.5d0/(rkep-1.0d0)
  psi0 = -0.5d0/(rin-1.0d0)+0.5d0*l0**2/rin**2
  
  do k=1,kx
     do j=1,jx
        do i=1,ix
           ss = sqrt(x(i)**2+z(k)**2)
           roc = rohalo*exp(-(gm/cshalo**2)*(gpot(i,j,k)-pot0))
           prc = (cshalo**2/gm)*roc
           ro_ed = rohalo*exp(-(gm/cshalo**2)*(gpot(i,j,k)-pot0))

           rod = 0.0d0
           prd = 0.0d0
           vyd = 0.0d0
           byd = 0.0d0
           if (x(i).gt.sseps) then
              tmp = psi0+0.5d0/(ss-1.0d0)-0.5d0*l0**2/x(i)**2
              tmp2 = psi0+0.5d0/(rkep-1.0d0)-0.5d0*l0**2/rkep**2
              rod0 = (tmp2/(k0*(gm/(gm-1.0d0))))**(1.0d0/(gm-1.0d0))
              if(tmp > 0.0d0)then
!!$                 rod = (tmp/(k0*(gm/(gm-1.0d0)) &
!!$                      *(1.0d0+x(i)**(2.0d0*(gm-1.0d0)) &
!!$                      /(betab*rkep**(2.0d0*(gm-1.0d0))))))**(1.0d0/(gm-1.0d0))
                 rod = ((tmp/(k0*(gm/(gm-1.0d0))))**(1.0d0/(gm-1.0d0)))!/rod0
                 prd = k0*rod**gm
                 vyd = l0/x(i)
!                 byd = sqrt(rod*h0*(rod*x(i)**2)**(gm-1.0d0))
!                 byd = sqrt(2.0d0*prd/betab)
              endif
           endif
           ro(i,j,k) = roc+rod
           pr(i,j,k) = prc+prd
           
           vx(i,j,k) = 0.0d0
           vy(i,j,k) = vyd
           vz(i,j,k) = 0.0d0
           
           bx(i,j,k) = 0.0d0
           by(i,j,k) = 0.0d0
           bz(i,j,k) = 0.0d0

           rod = 0.0d0
           ss = sqrt((x(i)+0.5d0*dx(i))**2+(z(k)+0.5d0*dz(k))**2)
           if (x(i).gt.sseps) then
              tmp = psi0+0.5d0/(ss-1.0d0)-0.5d0*l0**2/(x(i)+0.5d0*dx(i))**2
              tmp2 = psi0+0.5d0/(rkep-1.0d0)-0.5d0*l0**2/rkep**2
              rod0 = (tmp2/(k0*(gm/(gm-1.0d0))))**(1.0d0/(gm-1.0d0))
              if(tmp > 0.0d0)then
                 rod = ((tmp/(k0*(gm/(gm-1.0d0))))**(1.0d0/(gm-1.0d0)))!/rod0
              endif
           end if
           ro_edge(i,j,k) = ro_ed+rod
        enddo
     enddo
  enddo

  do k=1,kx-1
     do j=1,jx-1
        do i=1,ix-1
           aphi(i,j,k) = 0.0d0
           bxs(i,j,k) = 0.0d0
           bys(i,j,k) = 0.0d0
           bzs(i,j,k) = 0.0d0
           bxs_temp(i,j,k) = 0.0d0
           bzs_temp(i,j,k) = 0.0d0

           bx(i,j,k) = 0.0d0
           by(i,j,k) = 0.0d0
           bz(i,j,k) = 0.0d0
           phi(i,j,k) = 0.0d0
           eta(i,j,k) = 0.0d0
           ssr = x(i)+0.5d0*dx(i)
           ssz = z(k)+0.5d0*dz(k)
           ss = sqrt(ssr**2 + ssz**2)
           ro_edge_chk = min(ro(i,j,k),ro(i+1,j,k),ro(i,j,k+1),ro(i+1,j,k+1))
           if(ssr .gt. rin)then
              if(ro_edge_chk .gt. romin)then
                 aphi(i,j,k) = ro_edge(i,j,k)
              end if
           end if
        end do
     end do
  end do

  do k=2,kx-1
     do j=2,jx-1
        do i=2,ix-1
           dtodz = 1.0d0/dz(k)
           bxs_temp(i,j,k) = -dtodz*(aphi(i,j,k)-aphi(i,j,k-1))
           line1 = x(i)+0.5d0*dx(i)
           line2 = x(i)-0.5d0*dx(i)
           dtodx = 1.0d0/(x(i)*dx(i))
           bzs_temp(i,j,k) = dtodx*(line1*aphi(i,j,k)-line2*aphi(i-1,j,k))
        end do
     end do
  end do

  call exchangeMpixz_res(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bxs_temp,bys,bzs_temp &
       ,phi,eta,merr)

  do k=1,kx
     do j=1,jx
        do i=2,ix
           bx(i,j,k) = 0.5d0*(bxs_temp(i,j,k)+bxs_temp(i-1,j,k))
        enddo
     enddo
  enddo

  do k=1,kx
     do j=2,jx
        do i=1,ix
           by(i,j,k) = 0.5d0*(bys(i,j,k)+bys(i,j-1,k))
        enddo
     enddo
  enddo

  do k=2,kx
     do j=1,jx
        do i=1,ix
           bz(i,j,k) = 0.5d0*(bzs_temp(i,j,k)+bzs_temp(i,j,k-1))
        enddo
     enddo
  enddo

  inbetamax = -1.0d0
  do k=1,kx
     do j=1,jx
        do i=1,ix
           pb = 0.5d0*(bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2)
           inbeta = pb/pr(i,j,k)
           if(inbeta > inbetamax)then
              inbetamax = inbeta
           endif
        end do
     end do
  end do

  call mpi_allreduce(inbetamax,inbetamaxg,1,mpi_double_precision,mpi_max &
       ,mpi_comm_world,merr)

  inbeta = inbetamaxg
  alpha = sqrt(1.0d0/(inbeta*betab))
  do k=1,kx
     do j=1,jx
        do i=1,ix
           bxs(i,j,k) = alpha*bxs_temp(i,j,k)
           bzs(i,j,k) = alpha*bzs_temp(i,j,k)
        end do
     end do
  end do

  do k=1,kx
     do j=1,jx
        do i=2,ix
           bx(i,j,k) = 0.5d0*(bxs(i,j,k)+bxs(i-1,j,k))
        enddo
     enddo
  enddo

  do k=1,kx
     do j=2,jx
        do i=1,ix
           by(i,j,k) = 0.5d0*(bys(i,j,k)+bys(i,j-1,k))
        enddo
     enddo
  enddo

  do k=2,kx
     do j=1,jx
        do i=1,ix
           bz(i,j,k) = 0.5d0*(bzs(i,j,k)+bzs(i,j,k-1))
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
end subroutine model_Hawley2
