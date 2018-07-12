module  model

  use const
  use mpi_setup, only : mpid

  implicit none
  private

  public :: model_setup


contains


  subroutine  model_setup(ro,pr,vx,vy,vz,bx,by,bz,phi &
                         ,roi,pri,vxi,vyi,vzi,bxi,byi,bzi,phii &
                         ,x,dx,xm,y,dy,ym,z,dz,zm &
                         ,gx,gz,eta,min_dx)

!---Input & Output
  real(8),dimension(ix),      intent(out) :: x,dx
  real(8),dimension(0:ix),    intent(out) :: xm
  real(8),dimension(jx),      intent(out) :: y,dy
  real(8),dimension(0:jx),    intent(out) :: ym
  real(8),dimension(kx),      intent(out) :: z,dz
  real(8),dimension(0:kx),    intent(out) :: zm
  real(8),dimension(ix,jx,kx),intent(out) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(out) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(out) :: bx,by,bz


  real(8),dimension(ix,jx,kx),intent(out) :: eta,phi
  real(8),dimension(ix,jx,kx),intent(out) :: gx,gz
  real(8),dimension(ix,jx,kx),intent(out) :: roi,pri
  real(8),dimension(ix,jx,kx),intent(out) :: vxi,vyi,vzi
  real(8),dimension(ix,jx,kx),intent(out) :: bxi,byi,bzi
  real(8),dimension(ix,jx,kx),intent(out) :: phii
  real(8),                    intent(out) :: min_dx



  integer :: i,j,k
  integer :: ig,jg,kg
  integer :: izero,jzero,kzero
  real(8) :: ss
  real(8) :: tmp 
  real(8) :: psi0, pot0 
  real(8) :: roc,prc,vyc
  real(8) :: rod,prd,vyd,bxd,byd,bzd
  real(8) :: bbAbsMax
  real(8),dimension(ix,jx,kx) :: curx,cury,curz
  real(8),dimension(igx) :: xg,dxg
  real(8),dimension(0:igx) :: xmg
  real(8),dimension(jgx) :: yg,dyg
  real(8),dimension(0:jgx) :: ymg
  real(8),dimension(kgx) :: zg,dzg
  real(8),dimension(0:kgx) :: zmg
  real(8),dimension(ix,jx,kx) :: gpot
  real(8),dimension(igx,jgx,kgx) :: gxg,gzg,gpotg

  real(8),dimension(ix,jx,kx) :: pr_per, flag_torus
  real(8),parameter :: x0=1.d0
  real(8) :: Lnrml,L_r
!---Step 1a.-------------------------------------------------------------|
! set global x-grid 
!
  dxg(margin+1)=4.0d0*dxg0

!$OMP PARALLEL
!$OMP DO
  do i=margin+2,margin+ugrid_xmax
     dxg(i) = dxg0
  enddo
!$OMP END DO
!$OMP END PARALLEL

  do i=margin+1+ugrid_xmax,igx-margin
     dxg(i) = dxg(i-1)*ratio_x
     if(dxg(i) > dxmax) dxg(i)=dxmax
  enddo

  do i=igx-margin+1,igx
     dxg(i)=dxg(igx-margin)
  enddo
  do i=0,margin-1
     dxg(margin-i) = dxg(margin+i+1)
  enddo

! X origin
  izero = margin+1
  xmg(izero) = xmin + dxg(izero)

  do i=izero,igx-1
     xmg(i+1) = xmg(i)+dxg(i+1)
  enddo
  do i=izero-1,0,-1
     xmg(i) = xmg(i+1)-dxg(i+1)
  enddo

!$OMP PARALLEL DO
  do i=1,igx
     xg(i) = 0.5d0*(xmg(i)+xmg(i-1))
  enddo
!$OMP END PARALLEL DO

!---Step 1b.-------------------------------------------------------------|
! set global y-grid

!$OMP PARALLEL DO
  do j=1,jgx
     dyg(j) = dyg0
  enddo
!$OMP END PARALLEL DO

! Y origin
  jzero = margin+1
  ymg(jzero) = ymin+0.5d0*dyg(jzero)
  do j=jzero,jgx-1
     ymg(j+1) = ymg(j)+dyg(j+1)
  enddo
  do j=jzero-1,0,-1
     ymg(j) = ymg(j+1)-dyg(j+1)
  enddo

!$OMP PARALLEL DO
  do j=1,jgx
     yg(j) = 0.5d0*(ymg(j)+ymg(j-1))
  enddo
!$OMP END PARALLEL DO

!---Step 1c.-------------------------------------------------------------|
! set global z-grid
!$OMP PARALLEL DO
  do,k=1,kgx
     dzg(k) = dzg0
  enddo
!$OMP END PARALLEL DO

  do k=int(kgx/2)+1+ugrid_zmax,kgx
     dzg(k) = dzg(k-1)*ratio_z
     if(dzg(k).gt.dzmax) dzg(k)=dzmax
  enddo
  do k=kgx-margin,kgx
     dzg(k)=dzg(kgx-margin)
  enddo
  do k=int(kgx/2)-ugrid_zmax,margin,-1
     dzg(k) = dzg(k+1)*ratio_z
     if(dzg(k).gt.dzmax) dzg(k)=dzmax
  end do
  do k=margin,1,-1
     dzg(k) = dzg(k+1)
  enddo

! Z origin
  kzero = kgx/2+1
  zmg(kzero) = zmin + dzg(kzero)
  do k=kzero,kgx-1
     zmg(k+1) = zmg(k)+dzg(k+1)
  enddo
  do k=kzero-1,0,-1
     zmg(k) = zmg(k+1)-dzg(k+1)
  enddo
!$OMP PARALLEL DO
  do k=1,kgx
     zg(k) = 0.5d0*(zmg(k)+zmg(k-1))
  enddo
!$OMP END PARALLEL DO

!---Step 2a.-------------------------------------------------------------|
! set individual x-grid 
!$OMP PARALLEL
!$OMP DO PRIVATE(ig)
  do i=1,ix
     ig = mpid%mpirank_3d(1)*(ix-2*margin)+i
     x(i)=xg(ig)
     dx(i) = dxg(ig)
  enddo
!$OMP END DO
!$OMP DO PRIVATE(ig)
  do i=0,ix
     ig = mpid%mpirank_3d(1)*(ix-2*margin)+i
     xm(i) = xmg(ig)
  end do
!$OMP END DO

!---Step 2b.-------------------------------------------------------------|
! set individual y-grid 
!$OMP DO PRIVATE(jg)
  do j=1,jx
     jg = mpid%mpirank_3d(2)*(jx-2*margin)+j
     y(j) = yg(jg)
     dy(j) = dyg(jg)
  enddo
!$OMP END DO
!$OMP DO PRIVATE(jg)
  do j=0,jx
     jg = mpid%mpirank_3d(2)*(jx-2*margin)+j
     ym(j) = ymg(jg)
  end do
!$OMP END DO

!---Step 2c.-------------------------------------------------------------|
! set individual z-grid 
!$OMP DO PRIVATE(kg)
  do k=1,kx
     kg=mpid%mpirank_3d(3)*(kx-2*margin)+k
     z(k) = zg(kg)
     dz(k) = dzg(kg)
  enddo
!$OMP END DO
!$OMP DO PRIVATE(kg)
  do k=0,kx
     kg=mpid%mpirank_3d(3)*(kx-2*margin)+k
     zm(k) = zmg(kg)
  end do
!$OMP END DO
!$OMP END PARALLEL

! calculate min_dx
  min_dx = min(minval(dxg),minval(dzg))
!$OMP PARALLEL DO REDUCTION(min:min_dx)
  do j=1,jgx
     do i=margin+1,igx
        min_dx=min(min_dx,xg(i)*dyg(j))
     enddo
  enddo
!$OMP END PARALLEL DO

!---Step 3a.----------------------------------------------------------|
  ! set gravitation

  !$OMP PARALLEL
  !$OMP DO &
  !$OMP PRIVATE(i,j,ss)
  do k=1,kgx
     do j=1,jgx
        do i=1,igx
           gpotg(i,j,k) = 0.0d0
           gxg(i,j,k) = 0.0d0
           gzg(i,j,k) = 0.0d0
        enddo
     enddo
  enddo
  !$OMP END DO
  
  !$OMP DO &
  !$OMP PRIVATE(i,j,ss)
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
  !$OMP END DO

!---Step 3b.----------------------------------------------------------|
  ! set individual gravitation
  !$OMP DO &
  !$OMP PRIVATE(i,j,ig,jg,kg)
  do k=1,kx
     do j=1,jx
        do i=1,ix
           ig = mpid%mpirank_3d(1)*(ix-2*margin)+i
           jg = mpid%mpirank_3d(2)*(jx-2*margin)+j
           kg = mpid%mpirank_3d(3)*(kx-2*margin)+k
           gpot(i,j,k) = gpotg(ig,jg,kg)
           gx(i,j,k) = gxg(ig,jg,kg)
           gz(i,j,k) = gzg(ig,jg,kg)
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
!----------------------------------------------------------------------|
! set initial model
!  bbAbsMax = 0.0d0
  pot0=-grav/(1.0d0-ssg)
  psi0 = pot0+1.0d0/(2.0d0-2.0d0*aa)+gm*kk/(gm-1.0d0)*(1+1.d0/beta0)

!---2D model with poloidal magnetic field -------------------------
  if(pol_2D)then
     !$OMP PARALLEL
     !$OMP DO &
     !$OMP PRIVATE(i,j) &
     !$OMP PRIVATE(ss,roc,prc,rod,prd,vyd,bxd,byd,bzd,Lnrml,L_r)
     do k=1,kx
        do j=1,jx
           do i=1,ix
              !---- set corona ---------------------
              ss = sqrt(x(i)**2+z(k)**2)
              roc = rohalo*exp(-(gpot(i,j,k)-pot0)/tec0)
              prc = tec0*roc/gm
              
              !---- set torus ----------------------
              rod = 0.0d0
              prd = 0.0d0
              vyd = 0.0d0
              bxd = 0.0d0
              byd = 0.0d0
              bzd = 0.0d0
              Lnrml = x0**0.5
              flag_torus(i,j,k) = 0.d0
              if (ss.gt.sseps) then
                 L_r = Lnrml*x(i)**aa
                 rod = (x0/ss - 0.5d0*x0**2/x(i)**2 - 0.5d0/dst) &
                      /(n_poly+1.d0)/a_poly/x0
                 if (rod.le.0.d0) then
                    rod=0.d0
                    prd=0.d0
                    vyd=0.d0
                    bxd=0.d0
                    byd=0.d0
                    bzd=0.d0
                    goto 100
                 endif
                 rod = rod**n_poly ! i.e., rod =  rod**(1.d0/(gm-1.d0))
                 prd = a_poly*rod**gm
                 
                 bxd =2.d0/beta0/a_poly*n_poly/(n_poly +1.d0)*z(k)/ss**3 &
                      *rod**(-gm + 3.d0)
                 byd = 0.d0
                 bzd = rod**2/beta0/x(i) + 2.d0/beta0/a_poly*n_poly &
                      /(n_poly+1.d0)*(-x(i)/ss**3 + x0/x(i)**3) &
                      *rod**(-gm+3.d0)
                 vyd = L_r/x(i)                 
                 flag_torus(i,j,k) = 1
                 !endif
              endif
100           continue

              ro(i,j,k) = rod+roc
              pr(i,j,k) = prd+prc
              pr_per(i,j,k) = prd+prc 
              vx(i,j,k) = 0.0d0
              vy(i,j,k) = vyd
              vz(i,j,k) = 0.0d0
              bx(i,j,k) = bxd
              by(i,j,k) = byd
              bz(i,j,k) = bzd
              phi(i,j,k) = 0.0d0
              eta(i,j,k) = 0.0d0
           enddo
        enddo
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
     roi = ro
     pri = pr
     vxi = 0.0d0
     vyi = 0.0d0
     vzi = 0.0d0
     bxi = bx
     byi = 0.0d0
     bzi = bz
     phii = 0.0d0

     call perturb2(2,pr_per)
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j) 
     do k=1,kx
        do j=1,jx
           do i=1,ix
              if (flag_torus(i,j,k) == 1) then
                 pr(i,j,k) = pr_per(i,j,k)
              endif
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
endif
!---3D model with toroldal magnetic field -------------------------
if(tor_3D)then
   !$OMP PARALLEL DO &
   !$OMP PRIVATE(i,j) &
   !$OMP PRIVATE(ss,roc,prc,rod,prd,vyd,byd,tmp) &
   !$OMP REDUCTION(max:bbAbsMax)
   do k=1,kx
      do j=1,jx
         do i=1,ix
            !---- set corona ---------------------
            ss = sqrt(x(i)**2+z(k)**2)
            roc = rohalo*exp(-(gpot(i,j,k)-pot0)/tec0)
            !           prc = tec0*roc
            prc = tec0*roc/gm

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
   !$OMP END PARALLEL DO

  call perturb(1,bx,by,bz,x,dx,dy,dz,bbAbsMax)

endif
       
end subroutine model_setup

!----------------------------------------------------------------
! perturbation
!----------------------------------------------------------------
  subroutine perturb(iperturb,vx,vy,vz,x,dx,dy,dz,v0)

  integer,intent(in) :: iperturb
  real(8),intent(in) :: v0
  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: dy
  real(8),dimension(kx),intent(in) :: dz
  real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz

  integer :: i,j,k
  integer :: n
  integer,dimension(:),allocatable :: seed
  real(8) :: amp
  real(8) :: dtody,dtodz,dtodx
  real(8) :: line1,line2
  real(8) :: temp
  real(8) :: dvxmax,dvymax,dvzmax
  real(8) :: vmax
  real(8),dimension(ix,jx,kx) :: dax,day,daz
  real(8),dimension(ix,jx,kx) :: dvx,dvy,dvz
  real(8),dimension(ix,jx,kx) :: dvxc,dvyc,dvzc

  amp = 0.01d0*v0

  call random_seed(size = n)
  allocate(seed(n))
  seed = 123+mpid%mpirank
  call random_seed(PUT = seed)
  call random_number(dax)
  call random_number(day)
  call random_number(daz)

  if(iperturb .eq. 1)then
     !$OMP PARALLEL
     !$OMP DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,jx
           do i=1,ix
              dvx(i,j,k) = 0.0d0
              dvy(i,j,k) = 0.0d0
              dvz(i,j,k) = 0.0d0
           enddo
        enddo
     enddo
     !$OMP END DO

     dvxmax = 0.0d0
     !$OMP DO &
     !$OMP PRIVATE(i,j,dtodz,temp) &
     !$OMP REDUCTION(max:dvxmax)
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
     !$OMP END DO 
     
     dvymax = 0.0d0
     
     !$OMP DO &
     !$OMP PRIVATE(i,j,dtodz,dtodx,temp) &
     !$OMP REDUCTION(max:dvymax)
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
     !$OMP END DO 

     dvzmax = 0.0d0
     !$OMP DO &
     !$OMP PRIVATE(i,j,line1,line2,dtodx,temp) &
     !$OMP REDUCTION(max:dvzmax)
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
     !$OMP END DO 

     !$OMP DO &
     !$OMP PRIVATE(i,j) 
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
     !$OMP END DO 

     !$OMP DO &
     !$OMP PRIVATE(i,j)
     do k=2,kx-1
        do j=1,jx-1
           do i=2,ix-1
              dvxc(i,j,k) = 0.5d0*(dvx(i,j,k)+dvx(i-1,j,k))
              dvyc(i,j,k) = dvy(i,j,k)
              dvzc(i,j,k) = 0.5d0*(dvz(i,j,k)+dvz(i,j,k-1))
           enddo
        enddo
     enddo
     !$OMP END DO 

     !$OMP DO &
     !$OMP PRIVATE(i,j)
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
              vx(i,j,k) = vx(i,j,k) + dvxc(i,j,k)
              vy(i,j,k) = vy(i,j,k) + dvyc(i,j,k)
              vz(i,j,k) = vz(i,j,k) + dvzc(i,j,k)
           enddo
        enddo
     enddo
     !$OMP END DO 
     !$OMP END PARALLEL
  else if(iperturb .eq. 2)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
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
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,jx
           do i=1,ix
              dvx(i,j,k) = 0.0d0
              dvy(i,j,k) = 0.0d0
              dvz(i,j,k) = 0.0d0
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

  end subroutine perturb

  subroutine perturb2(iperturb,vx)

  integer,intent(in) :: iperturb
!  real(8),intent(in) :: v0
!  real(8),dimension(ix),intent(in) :: x,dx
!  real(8),dimension(jx),intent(in) :: dy
!  real(8),dimension(kx),intent(in) :: dz
!  real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(inout) :: vx

  integer :: i,j,k
  integer :: n
  integer,dimension(:),allocatable :: seed
  real(8) :: amp
  real(8) :: dtody,dtodz,dtodx
  real(8) :: line1,line2
  real(8) :: temp
  real(8) :: dvxmax,dvymax,dvzmax
  real(8) :: vmax
  real(8),dimension(ix,jx,kx) :: dax
!  real(8),dimension(ix,jx,kx) :: dvx,dvy,dvz
!  real(8),dimension(ix,jx,kx) :: dvxc,dvyc,dvzc

  amp = 0.01d0

  call random_seed(size = n)
  allocate(seed(n))
  seed = 123+mpid%mpirank
  call random_seed(PUT = seed)
  call random_number(dax)
!  call random_number(day)
!  call random_number(daz)

  if(iperturb .eq. 2)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,jx
           do i=1,ix
!              if(dabs(vy(i,j,k)) > 0.0d0)then
                 vx(i,j,k) = vx(i,j,k)*(1.d0 + 2.0d0*amp*(dax(i,j,k)-0.5d0))
!                 vy(i,j,k) = vy(i,j,k) + 2.0d0*amp*(day(i,j,k)-0.5d0)
!                 vz(i,j,k) = vz(i,j,k) + 2.0d0*amp*(daz(i,j,k)-0.5d0)
!              endif
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
!   else
!      do k=1,kx
!         do j=1,jx
!            do i=1,ix
!               dvx(i,j,k) = 0.0d0
!  !             dvy(i,j,k) = 0.0d0
!  !             dvz(i,j,k) = 0.0d0
!            enddo
!         enddo
!      enddo
  endif

  end subroutine perturb2


end module model
