subroutine model_kh(idf,mpid,igx,jgx,kgx,ix,jx,kx,margin &
       ,ro,pr,vx,vy,vz,bx,by,bz,phi,gm &
       ,x,dx,y,dy,z,dz &
       ,gx,gz,mf_params,dtout,tend,cr,eta0,vc,eta,xin &
       ,ccx,ccy,ccz)
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
  real(8),dimension(ix) :: x,dx,dxm,xx
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

  integer :: i,j,k,temp_i
  integer :: ig,jg,kg

  real(8) :: max,min !test

! machida
  real(8) :: l0, kk, beta0, hh,tmp
  real(8) :: rohalo, ro0, vk0, va0

  real(8) :: psi0, cshalo, pot0, tec0, tec00

  real(8) :: roc,prc,vyc
  real(8) :: rod,prd,vyd,byd

  real(8) :: bbAbsMax

  real(8) :: ro0,ro1,u0,u1,pr0,pr1,wtr
  real(8) :: rlambda,amp,mdum,alpha1,alpha2,wkx,omegar,omegai,vre,vim,rannum
  real(8) :: ma,mm,mq,mr,rnml,kk,mseed



!---Step 0.--------------------------------------------------------------|
! set parameter
!

  pi = acos(-1.0d0)
  gm=5.0d0/3.0d0
  cr = 0.18d0
  pi2 = 2.0d0*pi
  pi4 = 4.0d0*pi
  hpi4 = sqrt(pi4)

  tec0 = 1.0d0
  tec00 = tec0*gm

  beta0 = 100.0d0
  ro0 = 1.0d0
  vk0 = 1.0d0
  kk = 0.05d0
  hh = 2.0d0*kk/beta0
  l0 = 1.0d0
  rohalo = 0.001d0
  cshalo = 1.0d0
  va0 = hh
  
  psi0 = -1.0d0+0.5d0+gm*kk/(gm-1.0d0)+0.5d0*hh*gm/(gm-1.0)
  pot0 = -1.0d0
!---Step 1.--------------------------------------------------------------|
! set grid
!
  dxg0 = 1.0d0/real(ix-2*margin)
  dzg0 = 1.0d0/real(kx*32-2*margin)

  xmin = 1.0d0
  ymin = 0.0d0
  ymax = pi2/32.0d0
  zmin = 0.0d0

!---Step 1a.-------------------------------------------------------------|
! set global x-grid 
!
  dxmax = 10.0d0*dxg0
  do i=1,igx
     dxg(i) = dxg0
  enddo

  izero = margin+1

! X origin

  xmg(izero) = xmin + dxg(izero)

  do i=izero,igx-1
     xmg(i+1) = xmg(i)+dxg(i+1)
  enddo

  do i=izero-1,1,-1
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

  do j=jzero-1,1,-1
     ymg(j) = ymg(j+1)-dyg(j+1)
  enddo

  do j=1,jgx
     yg(j) = 0.5d0*(ymg(j)+ymg(j-1))
  enddo

!---Step 1c.-------------------------------------------------------------|
! set global z-grid
!

  dzmax = 10.0d0*dzg0

  do k=1,kgx
     dzg(k) = dzg0
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

  do k=kzero-1,1,-1
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


!---g y
!  call bd_pery(margin,gpotg,igx,jgx,kgx)
!  call bd_pery(margin,gxg,igx,jgx,kgx)
!  call bd_pery(margin,gzg,igx,jgx,kgx)

!---g x
!  call bd_synpx(0,margin,gpotg,igx,jgx,kgx)
!  call bd_synnx(0,margin,gxg,igx,jgx,kgx)
!  call bd_synpx(0,margin,gzg,igx,jgx,kgx)

!  call bd_frex(1,margin,gpotg,igx,jgx,kgx)
!  call bd_consx(1,margin,0.0d0,gxg,igx,jgx,kgx)  
!  call bd_frex(1,margin,gzg,igx,jgx,kgx)

!---g z
!  call bd_frez(0,margin,gpotg,igx,jgx,kgx)
!  call bd_frez(0,margin,gxg,igx,jgx,kgx)
!  call bd_consz(0,margin,0.0d0,gzg,igx,jgx,kgx)  

!  call bd_frez(1,margin,gpotg,igx,jgx,kgx)
!  call bd_frez(1,margin,gxg,igx,jgx,kgx)
!  call bd_consz(1,margin,0.0d0,gzg,igx,jgx,kgx)  

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

  ro0=1.0d0
  ro1=0.8d0
  u0=0.0d0
  u1=0.5d0
  pr0=1.0d0/gm
  pr1=1.0d0/gm
  wtr=0.01d0

  do k=1,kx
     do j=1,jx
        do i=1,ix

           ro(i,j,k) = ro0+(ro1-ro0)*(1+tanh(z(k)/wtr))/2
           pr(i,j,k) = pr0+(pr1-pr0)*(1+tanh(z(k)/wtr))/2
           vx(i,j,k) = u0+(u1-u0)*(1+tanh(z(k)/wtr))/2
           vz(i,j,k) = 0.0d0


           vy(i,j,k) = 0.0d0
           bx(i,j,k) = 0.0d0
           by(i,j,k) = 0.0d0
           bz(i,j,k) = 0.0d0
           phi(i,j,k) = 0.0d0

           eta(i,j,k) = 0.0d0
        enddo
     enddo
  enddo


!  call perturb2(1,mpid,ix,jx,kx,bx,by,bz,x,dx,dy,dz,bbAbsMax)
!!$
!!$  call getcurrent_cyl(bx,by,bz,ix,jx,kx,x,dx,dy,dz &
!!$       ,curx,cury,curz)
!!$  call getEta(ix,jx,kx,ro,curx,cury,curz,eta0,vc,eta)

!----------------------------------------------------------------------|
! perturbation
! 

      rlambda=1.d0/2.d0
      amp=1.d-2

      if (.false.) then
      mdum=-1
      do k=1,kx
      do j=1,jx
      do i=1,ix
      ma=16807
      mm=2147483647
      mq=127773
      mr=2836
      rnml=4.656612875245796924105750827168d-10

      kk=mdum/mq
      mseed=ma*(mdum-kk*mq)-mr*kk
      if (mseed.lt.0) mseed=mdum+mm
      rannum=2.0d0*rnml*mseed

!         vz(i,j,k) = amp*rannum
      enddo
      enddo
      enddo

      else

      alpha1=ro0/(ro0+ro1)
      alpha2=ro1/(ro0+ro1)

      wkx=2*pi/rlambda
      omegar=wkx*(alpha1*u0+alpha2*u1)
      omegai=sqrt(wkx*wkx*alpha1*alpha2*(u0-u1)*(u0-u1))

      do k=1,kx
      do j=1,jx
      do i=1,ix
       if (z(k).lt.0) then
         vre=amp*exp(wkx*z(k))*((wkx*u0-omegar)*cos(wkx*x(i))+ &
             omegai*sin(wkx*x(i)))
         vim=amp*exp(wkx*z(k))*((wkx*u0-omegar)*sin(wkx*x(i))- &
             omegai*cos(wkx*x(i)))
!         vz(i,j,k)=vz(i,j,k)+vre
!         vx(i,j,k)=vx(i,j,k)-vim
!         pr(i,j,k)=pr(i,j,k) &
!              -ro(i,j,k)/wkx*(omegai*vre-(wkx*u0-omegar)*vim)
       else
         vre=amp*exp(-wkx*z(k))*((wkx*u1-omegar)*cos(wkx*x(i)) &
             +omegai*sin(wkx*x(i)))
         vim=amp*exp(-wkx*z(k))*((wkx*u1-omegar)*sin(wkx*x(i)) &
             -omegai*cos(wkx*x(i)))
!         vz(i,j,k)=vz(i,j,k)+vre
!         vx(i,j,k)=vx(i,j,k)+vim
!         pr(i,j,k)=pr(i,j,k) &
!             +ro(i,j,k)/wkx*(omegai*vre-(wkx*u1-omegar)*vim)
        endif
      enddo
      enddo
      enddo

      endif







  return
end subroutine model_kh
