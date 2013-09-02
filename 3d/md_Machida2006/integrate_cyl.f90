subroutine integrate_cyl(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,dt &
     ,gx,gz,floor,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
     ,roi,pri,vxi,vyi,vzi,bxi,byi,bzi &
     ,eta0,vc,eta,ccx,ccy,ccz,RadCool,te_factor,time,rohalo,swtch_t,xin)

  use convert
  use mpi_domain_xz
  implicit none

!--Input
  integer,intent(in) :: ix,jx,kx,margin
  integer :: merr
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

  real(8),intent(in) :: RadCool,te_factor,time,rohalo,swtch_t,xin
  real(8),dimension(ix,jx,kx),intent(in) :: roi,pri,vxi,vyi,vzi
  real(8),dimension(ix,jx,kx),intent(in) :: bxi,byi,bzi
  real(8),dimension(ix,jx,kx),intent(in) :: gx,gz

!-- using flux
  real(8),dimension(ix,jx,kx) :: ro1,pr1,vx1,vy1,vz1
  real(8),dimension(ix,jx,kx) :: bx1,by1,bz1

  real(8),dimension(ix,jx,kx) :: phi1

!-- Input & output
  real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr,vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz

  real(8),dimension(ix,jx,kx),intent(inout) :: phi
  real(8),dimension(ix,jx,kx) :: eta
!--Temporary vari
  real(8),dimension(ix,jx,kx) :: curx,cury,curz

!-conserved variable
  real(8),dimension(ix,jx,kx) :: rx,ry,rz,ee
  real(8),dimension(ix,jx,kx) :: rx1,ry1,rz1,ee1

!-source
  real(8) :: sro,srx,sry,srz
  real(8) :: see,sphi,sbz

!-surface variables
  real(8),dimension(ix,jx,kx,2) :: row,prw,vxw,vyw,vzw
  real(8),dimension(ix,jx,kx,2) :: bxw,byw,bzw,phiw
  real(8),dimension(ix,jx,kx) :: bx_m,by_m,bz_m,phi_m
!-Numerical flux
!x-component
  real(8),dimension(ix,jx,kx) :: frox,feex,frxx,fryx,frzx
  real(8),dimension(ix,jx,kx) :: fbyx,fbzx,fbxx,fphix
  real(8),dimension(ix,jx,kx) :: fbyxr,fbzxr,fbxxr,feexr
!y-component
  real(8),dimension(ix,jx,kx) :: froy,feey,frxy,fryy,frzy
  real(8),dimension(ix,jx,kx) :: fbyy,fbzy,fbxy,fphiy
  real(8),dimension(ix,jx,kx) :: fbyyr,fbzyr,fbxyr,feeyr
!z-component
  real(8),dimension(ix,jx,kx) :: froz,feez,frxz,fryz,frzz
  real(8),dimension(ix,jx,kx) :: fbxz,fbyz,fbzz,fphiz
  real(8),dimension(ix,jx,kx) :: fbxzr,fbyzr,fbzzr,feezr

!-other temporary variables

  real(8) :: dtodx,dtody,dtodz,hdt

  integer :: mdir
  integer :: i,j,k,n

  real(8) :: inversex             !1/x

  real(8) :: pi,hpi4,inhpi4

  real(8) ::te

  real(8) :: ratio,limit

  ratio=10000.0d0
  limit=0.0d0

!-----Step 0.----------------------------------------------------------|
! primitive to conserve
!

  cp = sqrt(ch*cr)
  pi = acos(-1.0d0)
  hpi4 = sqrt(4.0d0*pi)
  inhpi4 = 1.0d0/hpi4

  call convert_ptoc_m(ix,jx,kx,gm,ro,pr,vx,vy,vz,bx,by,bz &
       ,rx,ry,rz,ee)

  ro1=ro
  rx1=rx
  ry1=ry
  rz1=rz
  bx1=bx
  by1=by
  bz1=bz
  ee1=ee
  phi1=phi

do n=1,2
  call getcurrent_cyl(bx,by,bz,ix,jx,kx,x,dx,dy,dz &
       ,curx,cury,curz)

  call getEta2(ix,jx,kx,ro,curx,cury,curz,eta0,vc,eta,rohalo)

!-----Step 1a.---------------------------------------------------------|
! Compute flux in x-direction
! set L/R state at x-direction
!
  mdir = 1

  call lr_state_MP5_2(mdir,ix,jx,kx,ro,pr &
       ,vx,vy,vz,bx,by,bz,phi &
       ,ch,gm,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw,ccx,ccy,ccz)

  call MP5toMC2(ix,jx,kx,x,dx,y,dy,z,dz &
       ,ro,ro,pr,vx,vy,vz,bx,by,bz,phi &
       ,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw &
       ,mdir,floor,ratio,xin)
 
  call cal_interface_BP(ix,jx,kx,bxw,phiw &
       ,bx_m,phi_m,ch)
  call glm_flux(bx_m,phi_m,ch,fbxx,fphix,ix,jx,kx)
  
  call hlld_flux(row,prw,vxw,vyw,vzw,bx_m,byw,bzw,gm,ix,jx,kx,floor &
       ,frox,feex,frxx,fryx,frzx,fbyx,fbzx)

  call cal_resflux(mdir,ix,jx,kx,fbyx,curz,eta,-1.0d0 &
       ,fbyxr)
  call cal_resflux(mdir,ix,jx,kx,fbzx,cury,eta,+1.0d0 &
       ,fbzxr)
  call cal_resflux_fee(mdir,ix,jx,kx,feex,curx,cury,curz,bx,by,bz,eta &
       ,feexr)
!-----Step 1b.---------------------------------------------------------|
! compute flux at y-direction
! set L/R state at y-direction
!
  mdir = 2

  call lr_state_MP5_2(mdir,ix,jx,kx,ro,pr &
       ,vy,vz,vx,by,bz,bx,phi &
       ,ch,gm,row,prw,vyw,vzw,vxw,byw,bzw,bxw,phiw,ccx,ccy,ccz)

  call MP5toMC2(ix,jx,kx,x,dx,y,dy,z,dz &
       ,ro,ro,pr,vx,vy,vz,bx,by,bz,phi &
       ,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw &
       ,mdir,floor,ratio,xin)

  call cal_interface_BP(ix,jx,kx,byw,phiw &
       ,by_m,phi_m,ch)

  call glm_flux(by_m,phi_m,ch,fbyy,fphiy,ix,jx,kx)

  call hlld_flux(row,prw,vyw,vzw,vxw,by_m,bzw,bxw,gm,ix,jx,kx,floor &
       ,froy,feey,fryy,frzy,frxy,fbzy,fbxy)

  call cal_resflux(mdir,ix,jx,kx,fbzy,curz,eta,-1.0d0 &
       ,fbzyr)
  call cal_resflux(mdir,ix,jx,kx,fbxy,curx,eta,+1.0d0 &
       ,fbxyr)
  call cal_resflux_fee(mdir,ix,jx,kx,feey,curx,cury,curz,bx,by,bz,eta &
       ,feeyr)
!-----Step 1c.---------------------------------------------------------|
! compute flux at z-direction
! set L/R state at z-direction
!
  mdir = 3
  
  call lr_state_MP5_2(mdir,ix,jx,kx,ro,pr &
       ,vz,vx,vy,bz,bx,by,phi &
       ,ch,gm,row,prw,vzw,vxw,vyw,bzw,bxw,byw,phiw,ccx,ccy,ccz)

  call MP5toMC2(ix,jx,kx,x,dx,y,dy,z,dz &
       ,ro,ro,pr,vx,vy,vz,bx,by,bz,phi &
       ,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw &
       ,mdir,floor,ratio,xin)

  call cal_interface_BP(ix,jx,kx,bzw,phiw &
       ,bz_m,phi_m,ch)

  call glm_flux(bz_m,phi_m,ch,fbzz,fphiz,ix,jx,kx)

  call hlld_flux(row,prw,vzw,vxw,vyw,bz_m,bxw,byw,gm,ix,jx,kx,floor &
       ,froz,feez,frzz,frxz,fryz,fbxz,fbyz)

  call cal_resflux(mdir,ix,jx,kx,fbxz,cury,eta,-1.0d0 &
       ,fbxzr)
  call cal_resflux(mdir,ix,jx,kx,fbyz,curx,eta,+1.0d0 &
       ,fbyzr)
  call cal_resflux_fee(mdir,ix,jx,kx,feez,curx,cury,curz,bx,by,bz,eta &
       ,feezr)

!-----Step 2.---------------------------------------------------------|
! half time step update cell center variables using flux
!
  hdt=dt*0.5d0*n
  do k=margin+1,kx-margin
     do j=margin+1,jx-margin
        do i=margin+1,ix-margin

           inversex = 1.0d0/x(i)
! source term
! density
           sro = -0.5d0*(frox(i,j,k)+frox(i-1,j,k))*inversex
! x-momentum
           srx = -(ro(i,j,k)*(vx(i,j,k)**2-vy(i,j,k)**2) &
                +(by(i,j,k)**2-bx(i,j,k)**2))*inversex &
                +ro(i,j,k)*gx(i,j,k)
! y-momentum
           sry = -(fryx(i,j,k)+fryx(i-1,j,k))*inversex
! z-momentum
           srz = -0.5d0*(frzx(i,j,k)+frzx(i-1,j,k))*inversex &
                + ro(i,j,k)*gz(i,j,k)
! z-magnetic
           sbz = -0.5d0*(fbzxr(i-1,j,k)+fbzxr(i,j,k))*inversex
! energy
           see = -0.5d0*(feexr(i,j,k)+feexr(i-1,j,k))*inversex &
                +ro(i,j,k)*(vx(i,j,k)*gx(i,j,k)+vz(i,j,k)*gz(i,j,k))
!           if (time .gt. swtch_t) then
           if (ro(i,j,k) .gt. rohalo) then
           te = te_factor*pr(i,j,k)/ro(i,j,k)
           see = see - RadCool*(ro(i,j,k)**2)*sqrt(te)
           endif
!           endif
! phi
           sphi = -0.5d0*(fphix(i,j,k)+fphix(i-1,j,k))*inversex

! update
           dtodx = hdt/dx(i)
           dtody = hdt/(x(i)*dy(j))
           dtodz = hdt/dz(k)

           ro(i,j,k) = ro1(i,j,k)+dtodx*(frox(i-1,j,k)-frox(i,j,k)) &
                +dtody*(froy(i,j-1,k)-froy(i,j,k)) &
                +dtodz*(froz(i,j,k-1)-froz(i,j,k)) &
                +hdt*sro
           ee(i,j,k) = ee1(i,j,k)+dtodx*(feexr(i-1,j,k)-feexr(i,j,k)) &
                +dtody*(feeyr(i,j-1,k)-feeyr(i,j,k)) &
                +dtodz*(feezr(i,j,k-1)-feezr(i,j,k)) &
                +hdt*see
           rx(i,j,k) = rx1(i,j,k)+dtodx*(frxx(i-1,j,k)-frxx(i,j,k)) &
                +dtody*(frxy(i,j-1,k)-frxy(i,j,k)) &
                +dtodz*(frxz(i,j,k-1)-frxz(i,j,k)) &
                +hdt*srx
           ry(i,j,k) = ry1(i,j,k)+dtodx*(fryx(i-1,j,k)-fryx(i,j,k)) &
                +dtody*(fryy(i,j-1,k)-fryy(i,j,k)) &
                +dtodz*(fryz(i,j,k-1)-fryz(i,j,k)) &
                +hdt*sry
           rz(i,j,k) = rz1(i,j,k)+dtodx*(frzx(i-1,j,k)-frzx(i,j,k)) &
                +dtody*(frzy(i,j-1,k)-frzy(i,j,k)) &
                +dtodz*(frzz(i,j,k-1)-frzz(i,j,k)) &
                +hdt*srz
           bx(i,j,k) = bx1(i,j,k)+dtodx*(fbxx(i-1,j,k)-fbxx(i,j,k)) &
                +dtody*(fbxyr(i,j-1,k)-fbxyr(i,j,k)) &
                +dtodz*(fbxzr(i,j,k-1)-fbxzr(i,j,k))
           by(i,j,k) = by1(i,j,k)+dtodx*(fbyxr(i-1,j,k)-fbyxr(i,j,k)) &
                +dtody*(fbyy(i,j-1,k)-fbyy(i,j,k)) &
                +dtodz*(fbyzr(i,j,k-1)-fbyzr(i,j,k))
           bz(i,j,k) = bz1(i,j,k)+dtodx*(fbzxr(i-1,j,k)-fbzxr(i,j,k)) &
                +dtody*(fbzyr(i,j-1,k)-fbzyr(i,j,k)) &
                +dtodz*(fbzz(i,j,k-1)-fbzz(i,j,k)) &
                +hdt*sbz
           phi(i,j,k) = (phi1(i,j,k)+dtodx*(fphix(i-1,j,k)-fphix(i,j,k)) &
                +dtody*(fphiy(i,j-1,k)-fphiy(i,j,k)) &
                +dtodz*(fphiz(i,j,k-1)-fphiz(i,j,k)) &
                +hdt*sphi)*exp(-hdt*ch**2/cp**2)
        enddo
     enddo
  enddo

!-----Step 3.----------------------------------------------------------|
! conserved to primitive
!
  call convert_ctop_m(ix,jx,kx,gm,ro,ee,rx,ry,rz,bx,by,bz,floor &
       ,vx,vy,vz,pr)

  call exchangeMpixz(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz &
       ,phi,merr)

  call bnd(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,x,z &
             ,xin,roi,pri,vxi,vyi,vzi,bxi,byi,bzi)
enddo

return
end subroutine integrate_cyl
