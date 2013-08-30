subroutine model( )
  use mpi_domain_xz
  use const
  use init
  use boundary

  implicit none

  integer :: idf

!---Input & Output
  real(8),dimension(ix) :: dxm,xx
  real(8),dimension(0:ix) :: xm  
  real(8),dimension(jx) :: dym
  real(8),dimension(0:jx) :: ym  
  real(8),dimension(kx) :: dzm
  real(8),dimension(0:kx) :: zm  

  real(8),dimension(igx) :: xg,dxg,dxmg
  real(8),dimension(0:igx) :: xmg
  real(8),dimension(jgx) :: yg,dyg,dymg
  real(8),dimension(0:jgx) :: ymg
  real(8),dimension(kgx) :: zg,dzg,dzmg
  real(8),dimension(0:kgx) :: zmg

  real(8),dimension(ix,jx,kx) :: gpot
  real(8),dimension(igx,jgx,kgx) :: gxg,gzg,gpotg

!---Temp coordinate
  real(8),dimension(ix,jx,kx) :: curx,cury,curz

  real(8) :: dxmax,dymax,dzmax

  integer :: izero,jzero,kzero

  real(8) :: ss

!-------Temp phys
  integer :: ig,jg,kg

!---Step 1a.-------------------------------------------------------------|
! set global x-grid 
!

  dxmax = 10.0d0*dxg0

  dxg(margin+1)=4.0d0*dxg0

  do i=margin+2,margin+96
     dxg(i) = dxg0
  enddo

  do i=margin+97,igx-margin
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

  dzmax = 10.0d0*dzg0

  do,k=1,kgx
     dzg = dzg0
  enddo

  do k=margin+1,kgx
     dzg(k) = dzg(k-1)*ratio_z
     if(dzg(k).gt.dzmax) dzg(k)=dzmax
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
                 gzg(i,j,k) = 0.0d0
              endif
           endif
        end do
     end do
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        enddo
     enddo
  enddo

  roi = ro
  pri = pr
  vxi = vx
  vyi = vy
  vzi = vz
  bxi = bx
  byi = by
  bzi = bz

  return
end subroutine model
