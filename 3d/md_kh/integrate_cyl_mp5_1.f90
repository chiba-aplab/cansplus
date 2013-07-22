subroutine integrate_cyl_mp5_1(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,dt &
     ,gx,gz,floor,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
     ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
     ,eta0,vc,eta,ccx,ccy,ccz)
  use convert
  implicit none

!--Input
  integer,intent(in) :: ix,jx,kx,margin

  real(8),intent(in) :: ch,cr
  real(8) :: cp

  real(8),intent(in) :: dt,gm,eta0,vc
  real(8) :: dts

  real(8),intent(in) :: floor

  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: y,dy
  real(8),dimension(kx),intent(in) :: z,dz
  real(8),dimension(5,2,ix),intent(in) :: ccx
  real(8),dimension(5,2,jx),intent(in) :: ccy
  real(8),dimension(5,2,kx),intent(in) :: ccz

  real(8),dimension(ix,jx,kx),intent(in) :: gx,gz

!-- using flux
  real(8),dimension(ix,jx,kx) :: ro1,pr1,vx1,vy1,vz1
  real(8),dimension(ix,jx,kx) :: bx1,by1,bz1

  real(8),dimension(ix,jx,kx) :: phi1

!-- Input & output
  real(8),dimension(ix,jx,kx),intent(in) :: ro,pr,vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(in) :: bx,by,bz

  real(8),dimension(ix,jx,kx),intent(in) :: phi
  real(8),dimension(ix,jx,kx) :: eta
!--Temporary vari
!  real(8),dimension(ix,jx,kx) :: curx,cury,curz

!-conserved variable
  real(8),dimension(ix,jx,kx) :: rx,ry,rz,ee
  real(8),dimension(ix,jx,kx) :: rx1,ry1,rz1,ee1

!-source
  real(8),dimension(ix,jx,kx) :: sro,srx,sry,srz
  real(8),dimension(ix,jx,kx) :: see,sphi,sbz

!-surface variables
  real(8),dimension(ix,jx,kx,2) :: row,prw,vxw,vyw,vzw
  real(8),dimension(ix,jx,kx,2) :: bxw,byw,bzw,phiw
  real(8),dimension(ix,jx,kx) :: bx_m,by_m,bz_m,phi_m
  
!-Numerical flux
!x-component
  real(8),dimension(ix,jx,kx) :: frox,feex,frxx,fryx,frzx
  real(8),dimension(ix,jx,kx) :: fbyx,fbzx,fbxx,fphix
!y-component
  real(8),dimension(ix,jx,kx) :: froy,feey,frxy,fryy,frzy
  real(8),dimension(ix,jx,kx) :: fbyy,fbzy,fbxy,fphiy
!z-component
  real(8),dimension(ix,jx,kx) :: froz,feez,frxz,fryz,frzz
  real(8),dimension(ix,jx,kx) :: fbxz,fbyz,fbzz,fphiz

!-other temporary variables

  real(8) :: dtodx,dtody,dtodz

  integer :: mdir
  integer :: i,j,k

  real(8) :: inversex

  real(8) :: pt

  real(8) :: pi,hpi4,inhpi4
!-----Step 0.----------------------------------------------------------|
! primitive to conserve
!

  cp = sqrt(ch*cr)
  pi = acos(-1.0d0)
  hpi4 = sqrt(4.0d0*pi)
  inhpi4 = 1.0d0/hpi4

  call convert_ptoc_m(ix,jx,kx,gm,ro,pr,vx,vy,vz,bx,by,bz &
       ,rx,ry,rz,ee)

!!$  call getcurrent_cyl(bxf,byf,bzf,ix,jx,kx,x,dx,dy,dz &
!!$       ,curx,cury,curz)
!!$
!!$  call getEta(ix,jx,kx,rof,curx,cury,curz,eta0,vc,eta)

!-----Step 1a.---------------------------------------------------------|
! Compute flux in x-direction
! set L/R state at x-direction
!
  mdir = 1

  call MP5_reconstruction_charGlmMhd(mdir,ix,jx,kx,ro,pr &
       ,vx,vy,vz,bx,by,bz,phi &
       ,ch,gm,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw,ccx,ccy,ccz)

  call cal_interface_BP(ix,jx,kx,bxw,phiw &
       ,bx_m,phi_m,ch)
  call glm_flux(bx_m,phi_m,ch,fbxx,fphix,ix,jx,kx)
  
  call hlld_flux(row,prw,vxw,vyw,vzw,bx_m,byw,bzw,gm,ix,jx,kx,floor &
       ,frox,feex,frxx,fryx,frzx,fbyx,fbzx)

!!$  call cal_resflux(mdir,ix,jx,kx,fbyx,curz,eta,-1.0d0 &
!!$       ,fbyxr)
!!$  call cal_resflux(mdir,ix,jx,kx,fbzx,cury,eta,+1.0d0 &
!!$       ,fbzxr)

!-----Step 1b.---------------------------------------------------------|
! compute flux at y-direction
! set L/R state at y-direction
!
  mdir = 2

  call MP5_reconstruction_charGlmMhd(mdir,ix,jx,kx,ro,pr &
       ,vy,vz,vx,by,bz,bx,phi &
       ,ch,gm,row,prw,vyw,vzw,vxw,byw,bzw,bxw,phiw,ccx,ccy,ccz)

  call cal_interface_BP(ix,jx,kx,byw,phiw &
       ,by_m,phi_m,ch)

  call glm_flux(by_m,phi_m,ch,fbyy,fphiy,ix,jx,kx)

  call hlld_flux(row,prw,vyw,vzw,vxw,by_m,bzw,bxw,gm,ix,jx,kx,floor &
       ,froy,feey,fryy,frzy,frxy,fbzy,fbxy)

!!$  call cal_resflux(mdir,ix,jx,kx,fbzy,curz,eta,-1.0d0 &
!!$       ,fbzyr)
!!$  call cal_resflux(mdir,ix,jx,kx,fbxy,curx,eta,+1.0d0 &
!!$       ,fbxyr)

!-----Step 1c.---------------------------------------------------------|
! compute flux at z-direction
! set L/R state at z-direction
!
  mdir = 3

  call MP5_reconstruction_charGlmMhd(mdir,ix,jx,kx,ro,pr &
       ,vz,vx,vy,bz,bx,by,phi &
       ,ch,gm,row,prw,vzw,vxw,vyw,bzw,bxw,byw,phiw,ccx,ccy,ccz)

  call cal_interface_BP(ix,jx,kx,bzw,phiw &
       ,bz_m,phi_m,ch)

  call glm_flux(bz_m,phi_m,ch,fbzz,fphiz,ix,jx,kx)

  call hlld_flux(row,prw,vzw,vxw,vyw,bz_m,bxw,byw,gm,ix,jx,kx,floor &
       ,froz,feez,frzz,frxz,fryz,fbxz,fbyz)

!!$  call cal_resflux(mdir,ix,jx,kx,fbxz,cury,eta,-1.0d0 &
!!$       ,fbxzr)
!!$  call cal_resflux(mdir,ix,jx,kx,fbyz,curx,eta,+1.0d0 &
!!$       ,fbyzr)

!-----Step 2a.---------------------------------------------------------|
! half time step update cell center variables using x-flux
!
  do k=3,kx-2
     do j=3,jx-2
        do i=3,ix-2
           dtodx = dt/dx(i)
           ro1(i,j,k) = ro(i,j,k)+dtodx*(frox(i-1,j,k)-frox(i,j,k))
           ee1(i,j,k) = ee(i,j,k)+dtodx*(feex(i-1,j,k)-feex(i,j,k))
           rx1(i,j,k) = rx(i,j,k)+dtodx*(frxx(i-1,j,k)-frxx(i,j,k))
           ry1(i,j,k) = ry(i,j,k)+dtodx*(fryx(i-1,j,k)-fryx(i,j,k))
           rz1(i,j,k) = rz(i,j,k)+dtodx*(frzx(i-1,j,k)-frzx(i,j,k))
           bx1(i,j,k) = bx(i,j,k)+dtodx*(fbxx(i-1,j,k)-fbxx(i,j,k))
           by1(i,j,k) = by(i,j,k)+dtodx*(fbyx(i-1,j,k)-fbyx(i,j,k))
           bz1(i,j,k) = bz(i,j,k)+dtodx*(fbzx(i-1,j,k)-fbzx(i,j,k))
           phi1(i,j,k) = phi(i,j,k)+dtodx*(fphix(i-1,j,k)-fphix(i,j,k))
        enddo
     enddo
  enddo

!-----Step 2b.---------------------------------------------------------|
! half time step update cell center variables using z-flux
!

  do k=3,kx-2
     do j=3,jx-2
        do i=3,ix-2
           dtody = dt/(x(i)*dy(j))
           ro1(i,j,k) = ro1(i,j,k)+dtody*(froy(i,j-1,k)-froy(i,j,k))
           ee1(i,j,k) = ee1(i,j,k)+dtody*(feey(i,j-1,k)-feey(i,j,k))
           rx1(i,j,k) = rx1(i,j,k)+dtody*(frxy(i,j-1,k)-frxy(i,j,k))
           ry1(i,j,k) = ry1(i,j,k)+dtody*(fryy(i,j-1,k)-fryy(i,j,k))
           rz1(i,j,k) = rz1(i,j,k)+dtody*(frzy(i,j-1,k)-frzy(i,j,k))
           bx1(i,j,k) = bx1(i,j,k)+dtody*(fbxy(i,j-1,k)-fbxy(i,j,k))
           bz1(i,j,k) = bz1(i,j,k)+dtody*(fbzy(i,j-1,k)-fbzy(i,j,k))
           by1(i,j,k) = by1(i,j,k)+dtody*(fbyy(i,j-1,k)-fbyy(i,j,k))
           phi1(i,j,k) = phi1(i,j,k)+dtody*(fphiy(i,j-1,k)-fphiy(i,j,k))
        enddo
     enddo
  enddo

!-----Step 2c.---------------------------------------------------------|
! half time step update cell center variables using z-flux
!

  do k=3,kx-2
     do j=3,jx-2
        do i=3,ix-2
           dtodz = dt/dz(k)
           ro1(i,j,k) = ro1(i,j,k)+dtodz*(froz(i,j,k-1)-froz(i,j,k))
           ee1(i,j,k) = ee1(i,j,k)+dtodz*(feez(i,j,k-1)-feez(i,j,k))
           rx1(i,j,k) = rx1(i,j,k)+dtodz*(frxz(i,j,k-1)-frxz(i,j,k))
           ry1(i,j,k) = ry1(i,j,k)+dtodz*(fryz(i,j,k-1)-fryz(i,j,k))
           rz1(i,j,k) = rz1(i,j,k)+dtodz*(frzz(i,j,k-1)-frzz(i,j,k))
           bx1(i,j,k) = bx1(i,j,k)+dtodz*(fbxz(i,j,k-1)-fbxz(i,j,k))
           by1(i,j,k) = by1(i,j,k)+dtodz*(fbyz(i,j,k-1)-fbyz(i,j,k))
           bz1(i,j,k) = bz1(i,j,k)+dtodz*(fbzz(i,j,k-1)-fbzz(i,j,k))
           phi1(i,j,k) = phi1(i,j,k)+dtodz*(fphiz(i,j,k-1)-fphiz(i,j,k))
        enddo
     enddo
  enddo
!-----Step 3.----------------------------------------------------------|
! add source term
! 


  do k=3,kx-2
     do j=3,jx-2
        do i=3,ix-2
           inversex = 1.0d0/x(i)
! density
           sro(i,j,k) = -0.5d0*(frox(i,j,k)+frox(i-1,j,k))*inversex
! x-momentum
           srx(i,j,k) = -(ro(i,j,k)*(vx(i,j,k)**2-vy(i,j,k)**2) &
                +(by(i,j,k)**2-bx(i,j,k)**2))*inversex &
                +ro(i,j,k)*gx(i,j,k)
! y-momentum
           sry(i,j,k) = -(fryx(i,j,k)+fryx(i-1,j,k))*inversex
! z-momentum
           srz(i,j,k) = -0.5d0*(frzx(i,j,k)+frzx(i-1,j,k))*inversex &
                + ro(i,j,k)*gz(i,j,k)
! z-magnetic
           sbz(i,j,k) = -0.5d0*(fbzx(i-1,j,k)+fbzx(i,j,k))*inversex
! energy
           see(i,j,k) = -0.5d0*(feex(i,j,k)+feex(i-1,j,k))*inversex &
                +ro(i,j,k)*(vx(i,j,k)*gx(i,j,k)+vz(i,j,k)*gz(i,j,k))
! phi
           sphi(i,j,k) = -0.5d0*(fphix(i,j,k)+fphix(i-1,j,k))*inversex
        end do
     end do
  end do
  !------
  ! Step 4.
  ! add source term
  !------
  do k=3,kx-2
     do j=3,jx-2
        do i=3,ix-2
           dts = dt
           ro1(i,j,k) = ro1(i,j,k)+dts*sro(i,j,k)
           rx1(i,j,k) = rx1(i,j,k)+dts*srx(i,j,k)
           ry1(i,j,k) = ry1(i,j,k)+dts*sry(i,j,k)
           rz1(i,j,k) = rz1(i,j,k)+dts*srz(i,j,k)
           bz1(i,j,k) = bz1(i,j,k)+dts*sbz(i,j,k)
           ee1(i,j,k) = ee1(i,j,k)+dts*see(i,j,k)
           phi1(i,j,k) = phi1(i,j,k)+dts*sphi(i,j,k)
           phi1(i,j,k) = phi1(i,j,k)*exp(-dts*ch**2/cp**2)
        end do
     end do
  end do

!-----Step 5.----------------------------------------------------------|
! conserved to primitive
!
  call convert_ctop_m(ix,jx,kx,gm,ro1,ee1,rx1,ry1,rz1,bx1,by1,bz1,floor &
       ,vx1,vy1,vz1,pr1)

  return
end subroutine integrate_cyl_mp5_1
