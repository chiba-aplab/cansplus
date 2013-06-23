subroutine gas_fall(idf,mpid,igx,jgx,kgx,ix,jx,kx,margin &
       ,ro,pr,vx,vy,vz,bx,by,bz,phi,gm &
       ,x,dx,y,dy,z,dz &
       ,gx,gz,mf_params,cr,eta0,vc,eta,xin &
       ,ccx,ccy,ccz,rg_nrmlx,rohalo,te_limit,te_factor,nrmlv &
       ,nrmlte,boltzmann_const,Navo,mmw)
  use mpi_domain_xz

  implicit none
  include 'mpif.h'

  integer,intent(in) :: ix,jx,kx
  integer,intent(in) :: mf_params

  integer,intent(in) :: margin
  integer,intent(in) :: igx,jgx,kgx
  
  type(mpidomain) :: mpid

  integer :: idf

  real(8),intent(in) :: gm,cr,eta0,vc,xin
  real(8) :: pi,pi2,pi4,hpi4

!---Input & Output
  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: y,dy
  real(8),dimension(kx),intent(in) :: z,dz
  real(8),dimension(5,2,ix),intent(in) :: ccx
  real(8),dimension(5,2,jx),intent(in) :: ccy
  real(8),dimension(5,2,kx),intent(in) :: ccz

  real(8),intent(in) :: rg_nrmlx,rohalo,te_limit,te_factor
  real(8),intent(in) :: nrmlv,nrmlte,boltzmann_const,Navo,mmw

  real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz

  real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(inout) :: eta,phi

  real(8),dimension(ix,jx,kx),intent(inout) :: gx,gz

!---Temp coordinate
  real(8),dimension(ix,jx,kx):: curx,cury,curz



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

  real(8) :: rogas,prgas,vxgas,vygas,vzgas,bxgas,bygas,bzgas

!add
  real(8) :: factorc
  real(8) :: aa

  real(8):: dv1,dv2,dv3,dv4,dv5,dv6

!---Step 0.--------------------------------------------------------------|
! set parameter
!

  factorc=boltzmann_const*Navo/mmw*nrmlte/nrmlv/nrmlv

  pi = acos(-1.0d0)
  pi2 = 2.0d0*pi

  tec0 = 5.0d-1
!  tec0 = 1.0d-2
!  tec0 = 1.0d-3
!  cshalo = 3.0d0/5.0d0

  beta0 = 100.0d0
!  beta0 = 1.0d10
!  ro0 = 1.0d0
!  vk0 = 1.0d0
!  kk = 0.05d0
  kk = 0.05d0
!  kk = /nrmlv/nrmlv
  ll = 1.0d0 
  aa=0.2d0 
  hh = 2.0d0*kk/beta0
     
  ssg=rg_nrmlx
  sseps = 0.2d0
  grav = (1.0d0-ssg)**2
  
!----------------------------------------------------------------------|
! set gas
! 
!  te_limit=0.0d0
!  pot0=-1.0d0
  pot0=-grav/(1.0d0-ssg)
!  psi0 = pot0+0.5d0*ll**2+gm*kk/(gm-1.0d0)+0.5d0*hh*gm/(gm-1.0)
  psi0 = pot0+1.0d0/(2.0d0-2.0d0*aa)+gm*kk/(gm-1.0d0)+0.5d0*hh*gm/(gm-1.0)
  do k=2,kx-1
     do j=2,jx-1
        do i=2,ix-1
           ss = sqrt((x(i)*cos(y(j))-3.5d0)**2+(x(i)*sin(y(j))-3.5d0)**2+(z(k))**2)
           rogas=0.0d0
           prgas=0.0d0
           vxgas=0.0d0
           vygas=0.0d0
           vzgas=0.0d0
           bxgas=0.0d0
           bygas=0.0d0
           bzgas=0.0d0

           if (ss.le.0.25d0) then
            rogas=30.0d0
            dv1=dx(i-1)*x(i-1)*dy(j)*dz(k)
            dv2=dx(i+1)*x(i+1)*dy(j)*dz(k)
            dv3=dx(i)*x(i)*dy(j-1)*dz(k)
            dv4=dx(i)*x(i)*dy(j+1)*dz(k)
            dv5=dx(i)*x(i)*dy(j)*dz(k-1)
            dv6=dx(i)*x(i)*dy(j)*dz(k+1)
            prgas=pr(i-1,j,k)*dv1+pr(i+1,j,k)*dv2+pr(i,j-1,k)*dv3+pr(i,j+1,k)*dv4+pr(i,j,k-1)*dv5+pr(i,j,k+1)*dv6
            prgas=prgas/(dv1+dv2+dv3+dv4+dv5+dv6)
!            vzgas=-1.0d0/(6.28d0*5.0d0)
            vzgas=0.0d0
            vygas=0.5d0*1.0d0/(3.5d0*1.414d0)
!              tmp = psi0-gpot(i,j,k)-x(i)**(2.0d0*aa-2.0d0)/(2.0d0-2.0d0*aa)
!              if(tmp > 0.0d0)then
!                 rogas = (tmp/(kk*(gm/(gm-1.0d0))))**(1.0d0/(gm-1.0d0))
!                 prgas = kk*rogas**gm
!                 vygas = x(i)**(aa-1.0d0)
!              endif
           endif


           ro(i,j,k) = ro(i,j,k)+rogas
           pr(i,j,k) = pr(i,j,k)+prgas
           vx(i,j,k) = vx(i,j,k)+vxgas
           vy(i,j,k) = vy(i,j,k)+vygas
           vz(i,j,k) = vz(i,j,k)+vzgas
           bx(i,j,k) = bx(i,j,k)+bxgas
           by(i,j,k) = by(i,j,k)+bygas
           bz(i,j,k) = bz(i,j,k)+bzgas
        enddo
     enddo
  enddo


  return
end subroutine gas_fall
