!======================================================================|
subroutine model_mpi_cart3d &
     & (igx, jgx, kgx, ix, jx, kx, margin &
     & , ro, pr, vx, vy, vz, bx, by, bz, phi &
     & , gm &
     & , x, dx, y, dy, z, dz &
     & , gx, gz, cr, eta0, vc, eta, xin &
     & , ccx, ccy, ccz, rg_nrmlx, rohalo, te_limit, te_factor, nrmlv &
     & , nrmlte, boltzmann_const, Navo, mmw &
     & , rank, topo_coords)
!======================================================================!
! NAME : model_mpi_cart3d
!
! PURPOSE: Set cordinates and initial condition
!
! INPUTS
!
! INOUTPUTS
!    
! OUTPUTS
!
! INTERNAL VARIABLES
!
! UPDATE LOG:
! 2012/12/30 initial coding (by H. ODA)
!======================================================================!
!  use mpi_domain_xz
  use cans_type
  use cans_mpi

  implicit none
  include 'mpif.h'

  integer,intent(in) :: ix,jx,kx

  integer,intent(in) :: margin
  integer,intent(in) :: igx,jgx,kgx

  type(cans_mpi_rank),intent(in) :: rank
  integer,intent(in) :: topo_coords(n_mpi_dims)

  real(8),intent(out) :: gm,cr,eta0,vc,xin

!---Input & Output
  real(8),dimension(ix),intent(out) :: x,dx
  real(8),dimension(jx),intent(out) :: y,dy
  real(8),dimension(kx),intent(out) :: z,dz

  real(8),dimension(5,2,ix),intent(out) :: ccx
  real(8),dimension(5,2,jx),intent(out) :: ccy
  real(8),dimension(5,2,kx),intent(out) :: ccz

  real(8),intent(in) :: rg_nrmlx,rohalo,te_limit,te_factor
  real(8),intent(in) :: nrmlv,nrmlte,boltzmann_const,Navo,mmw

  real(8),dimension(ix) :: dxm
  real(8),dimension(jx) :: dym
  real(8),dimension(kx) :: dzm
  real(8),dimension(0:ix) :: xm  
  real(8),dimension(0:jx) :: ym  
  real(8),dimension(0:kx) :: zm  

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
  real(8) :: l0, kk, ll, beta0, hh,tmp,tmp2
  real(8) :: ro0, vk0, va0

  real(8) :: psi0, cshalo, pot0, tec0, tec00

  real(8) :: roc,prc,vyc
  real(8) :: rod,prd,vyd,byd

  real(8) :: bbAbsMax

!add
  real(8) :: factorc
  real(8) :: aa
  
  integer :: perturb_type, i_random_seed

!----------------------------------------------------------------------|

!---Step 0.--------------------------------------------------------------|
! set parameter
!

  factorc=boltzmann_const*Navo/mmw*nrmlte/nrmlv/nrmlv

  !pi = acos(-1.0d0)
  gm=5.0d0/3.0d0
  !pi2 = 2.0d0*pi
  cr = 0.18d0

  tec0 = 5.0d-1

  beta0 = 100.0d0
  kk = 0.05d0
  ll = 1.0d0 
  aa=0.2d0 
  hh = 2.0d0*kk/beta0

!---Step 1.--------------------------------------------------------------|
! set grid
!
  dxg0 = 0.01d0
  dzg0 = 0.01d0

  xmin = 0.0d0
  ymin = 0.0d0
  ymax = pi2
  zmin = 0.0d0

!---Step 1a.-------------------------------------------------------------|
! set global x-grid 
!

  dxmax = 10.0d0*dxg0
  ratio = 1.05d0

  dxg(margin+1)=4.0d0*dxg0

  do i=margin+2,margin+96
	dxg(i) = dxg0
  enddo

  do i=margin+97,igx-margin
     dxg(i) = dxg(i-1)*ratio
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

  dzmax = 10.0d0*dzg0
  ratio = 1.07d0

  kzero = int(kgx/2)+1

  do,k=1,kgx
     !dzg = dzg0
     dzg(k) = dzg0
  enddo

  !do k=int(kgx/2)+45,kgx
  do k = kzero + 45, kgx - margin
     !dzg(k) = dzg(k-1) * ratio
     !if(dzg(k).gt.dzmax) dzg(k)=dzmax
     dzg(k) = min(dzg(k-1) * ratio, dzmax)
  enddo

  do k=kgx-margin,kgx
     !dzg(k)=dzg(kx-margin)
     dzg(k)=dzg(kgx-margin)
  enddo

  !do k=int(kgx/2)-44,margin,-1
  do k = kzero - 46, margin, -1
     !dzg(k) = dzg(k+1)*ratio
     !if(dzg(k).gt.dzmax) dzg(k)=dzmax
     dzg(k) = min(dzg(k+1) * ratio, dzmax)
  end do
  do k=margin,1,-1
     dzg(k) = dzg(k+1)
  enddo

  do k=1,kgx-1
     dzmg(k) = 0.5d0*(dzg(k+1)+dzg(k))
  enddo

  dzmg(kgx)=dzmg(kgx-1)

! Z origin

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
  do i = 1, ix
     !ig = mpid%mpirank_2d(1)*(ix-2*margin)+i
     ig     = topo_coords(1) * (ix - 2 * margin) + i
     x(i)   = xg(ig)
     dx(i)  = dxg(ig)
     dxm(i) = dxmg(ig)
  enddo

  do i = 0, ix
     !ig = mpid%mpirank_2d(1)*(ix-2*margin)+i
     ig    = topo_coords(1) * (ix - 2 * margin) + i
     xm(i) = xmg(ig)
  end do

!---Step 2b.-------------------------------------------------------------|
! set individual y-grid 
!
  do j = 1, jx
     !jg = j
     jg     = topo_coords(2) * (jx - 2 * margin) + j
     y(j)   = yg(jg)
     dy(j)  = dyg(jg)
     dym(j) = dymg(jg)
  enddo

  do j = 0, jx
     !jg = j
     jg    = topo_coords(2) * (jx - 2 * margin) + j
     ym(j) = ymg(jg)
  end do

!---Step 2c.-------------------------------------------------------------|
! set individual z-grid 
!
  do k = 1, kx
     !kg=mpid%mpirank_2d(2)*(kx-2*margin)+k
     kg     = topo_coords(3) * (kx - 2 * margin) + k
     z(k)   = zg(kg)
     dz(k)  = dzg(kg)
     dzm(k) = dzmg(kg)

  enddo

  do k = 0, kx
     !kg=mpid%mpirank_2d(2)*(kx-2*margin)+k
     kg    = topo_coords(3) * (kx - 2 * margin) + k
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











!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ssg=rg_nrmlx
  sseps = 0.2d0
  grav = (1.0d0-ssg)**2
  xin = 0.2d0	!2rg

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

  do k = 1, kx
     do j = 1, jx
        do i = 1, ix
           !ig = mpid%mpirank_2d(1)*(ix-2*margin)+i
           !jg = j
           !kg = mpid%mpirank_2d(2)*(kx-2*margin)+k
           ig = topo_coords(1) * (ix - 2 * margin) + i
           jg = topo_coords(2) * (jx - 2 * margin) + j
           kg = topo_coords(3) * (kx - 2 * margin) + k
           gpot(i,j,k) = gpotg(ig, jg, kg)
           gx(i,j,k) = gxg(ig, jg, kg)
           gz(i,j,k) = gzg(ig, jg, kg)
        enddo
     enddo
  enddo

  
!----------------------------------------------------------------------|
! set initial mocel
! 

  bbAbsMax = 0.0d0
  pot0=-grav/(1.0d0-ssg)
  psi0 = pot0+1.0d0/(2.0d0-2.0d0*aa)+gm*kk/(gm-1.0d0)+0.5d0*hh*gm/(gm-1.0)
  do k=1,kx
     do j=1,jx
        do i=1,ix
           ss = sqrt(x(i)**2+z(k)**2)
           roc = rohalo*exp(-(gpot(i,j,k)-pot0)/(factorc*tec0))
           prc = tec0*roc/te_factor

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
                 byd = sqrt(rod*hh*(rod*x(i)**2)**(gm-1.0d0))
                 if(dabs(byd) > bbAbsMax)then
                    bbAbsMax = dabs(byd)
                 endif
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
           phi(i,j,k) = 0.0d0

           eta(i,j,k) = 0.0d0
        enddo
     enddo
  enddo

  perturb_type  = 1
  i_random_seed = rank%nnn
  call perturb_vector(perturb_type,i_random_seed,ix,jx,kx,bx,by,bz,x,dx,dy,dz,bbAbsMax)

end subroutine model_mpi_cart3d
