module integrate

  implicit none
  private

  public :: integrate__TVDRK3


contains


  subroutine integrate__TVDRK3(margin,ix,jx,gm,dx,dy,dt          &
                              ,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cp &
                              ,eta,ccx,ccy)

  use convert
  use lr_state, only : lr_state__MP5, lr_state__MSCL2
  use flux_calc
  use bnd

  integer,intent(in) :: ix,jx,margin
  real(8),intent(in) :: ch,cp
  real(8),intent(in) :: dt,gm
  real(8),dimension(ix),intent(in) :: dx
  real(8),dimension(jx),intent(in) :: dy
  real(8),dimension(5,2,ix),intent(in) :: ccx
  real(8),dimension(5,2,jx),intent(in) :: ccy
  real(8),dimension(ix,jx),intent(inout) :: ro,pr,vx,vy,vz
  real(8),dimension(ix,jx),intent(inout) :: bx,by,bz
  real(8),dimension(ix,jx),intent(inout) :: phi
  real(8),dimension(ix,jx),intent(inout) :: eta

!-- using flux
  real(8),dimension(ix,jx) :: ro1,pr1,vx1,vy1,vz1
  real(8),dimension(ix,jx) :: bx1,by1,bz1
  real(8),dimension(ix,jx) :: phi1
!-conserved variable
  real(8),dimension(ix,jx) :: rx,ry,rz,ee
  real(8),dimension(ix,jx) :: rx1,ry1,rz1,ee1
!-surface variables
  real(8),dimension(ix,jx,2) :: row,prw,vxw,vyw,vzw
  real(8),dimension(ix,jx,2) :: bxw,byw,bzw,phiw
  real(8),dimension(ix,jx) :: bx_m,by_m,bz_m,phi_m
!-Numerical flux
!x-component
  real(8),dimension(ix,jx) :: frox,feex,frxx,fryx,frzx
  real(8),dimension(ix,jx) :: fbyx,fbzx,fbxx,fphix
!y-component
  real(8),dimension(ix,jx) :: froy,feey,frxy,fryy,frzy
  real(8),dimension(ix,jx) :: fbyy,fbzy,fbxy,fphiy
!z-component
  real(8),dimension(ix,jx) :: froz,feez,frxz,fryz,frzz
  real(8),dimension(ix,jx) :: fbxz,fbyz,fbzz,fphiz

!-other temporary variables
  integer :: mdir
  integer :: i,j,n
  real(8), parameter :: fac=1.D0/12.D0
  real(8) :: dtodx,dtody,k1,k2

!-----Step 0.----------------------------------------------------------|
! primitive to conserve
  call convert__ptoc(ix,jx,gm,ro,pr,vx,vy,vz,bx,by,bz &
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

  do n=1,3

!-----Step 1a.---------------------------------------------------------|
! Compute flux in x-direction
! set L/R state at x-direction
!
  mdir = 1

  call lr_state__MP5(mdir,ix,jx,ro,pr,vx,vy,vz,bx,by,bz,phi &
       ,ch,gm,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw,ccx,ccy)
!  call lr_state__MSCL2(mdir,ix,jx,ro,pr &
!       ,vx,vy,vz,bx,by,bz,phi &
!       ,ch,gm,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw,dx,dy)

  call flux_calc__bp(ix,jx,bxw,phiw,bx_m,phi_m,ch)
  call flux_calc__glm(bx_m,phi_m,ch,fbxx,fphix,ix,jx)
  
  call flux_calc__hlld(row,prw,vxw,vyw,vzw,bx_m,byw,bzw,gm,margin,ix,jx &
                      ,frox,feex,frxx,fryx,frzx,fbyx,fbzx)

!-----Step 1b.---------------------------------------------------------|
! compute flux at y-direction
! set L/R state at y-direction
!
  mdir = 2

  call lr_state__MP5(mdir,ix,jx,ro,pr,vy,vz,vx,by,bz,bx,phi &
       ,ch,gm,row,prw,vyw,vzw,vxw,byw,bzw,bxw,phiw,ccx,ccy)
!  call lr_state__MSCL2(mdir,ix,jx,ro,pr &
!       ,vy,vz,vx,by,bz,bx,phi &
!       ,ch,gm,row,prw,vyw,vzw,vxw,byw,bzw,bxw,phiw,dx,dy)

  call flux_calc__bp(ix,jx,byw,phiw,by_m,phi_m,ch)

  call flux_calc__glm(by_m,phi_m,ch,fbyy,fphiy,ix,jx)

  call flux_calc__hlld(row,prw,vyw,vzw,vxw,by_m,bzw,bxw,gm,margin,ix,jx &
                      ,froy,feey,fryy,frzy,frxy,fbzy,fbxy)


!-----Step 2.---------------------------------------------------------|
! TVDRK substep
  k1 = fac*(-7.D0*n*n+30.D0*n-23.D0)
  k2 = fac*(+7.D0*n*n-30.D0*n+35.D0)
  do j=margin+1,jx-margin
     do i=margin+1,ix-margin

        dtodx = dt/dx(i)
        dtody = dt/dy(j)

        ro(i,j) = k1*ro1(i,j)+k2*(+ro(i,j) &
             +dtodx*(frox(i-1,j)-frox(i,j))  &
             +dtody*(froy(i,j-1)-froy(i,j)) )
        ee(i,j) = k1*ee1(i,j)+k2*(+ee(i,j)  &
             +dtodx*(feex(i-1,j)-feex(i,j)) &
             +dtody*(feey(i,j-1)-feey(i,j)) )
        rx(i,j) = k1*rx1(i,j)+k2*(+rx(i,j) &
             +dtodx*(frxx(i-1,j)-frxx(i,j))  &
             +dtody*(frxy(i,j-1)-frxy(i,j)) )
        ry(i,j) = k1*ry1(i,j)+k2*(+ry(i,j) &
             +dtodx*(fryx(i-1,j)-fryx(i,j))  &
             +dtody*(fryy(i,j-1)-fryy(i,j)) )
        rz(i,j) = k1*rz1(i,j)+k2*(+rz(i,j) &
             +dtodx*(frzx(i-1,j)-frzx(i,j))  &
             +dtody*(frzy(i,j-1)-frzy(i,j)) )
        bx(i,j) = k1*bx1(i,j)+k2*(+bx(i,j)  &
             +dtodx*(fbxx(i-1,j)-fbxx(i,j))   &
             +dtody*(fbxy(i,j-1)-fbxy(i,j)) )
        by(i,j) = k1*by1(i,j)+k2*(+by(i,j)  &
             +dtodx*(fbyx(i-1,j)-fbyx(i,j)) &
             +dtody*(fbyy(i,j-1)-fbyy(i,j)) )
        bz(i,j) = k1*bz1(i,j)+k2*(+bz(i,j)  &
             +dtodx*(fbzx(i-1,j)-fbzx(i,j)) &
             +dtody*(fbzy(i,j-1)-fbzy(i,j)) )
        phi(i,j) = k1*phi1(i,j)+k2*(+phi(i,j) &
             +dtodx*(fphix(i-1,j)-fphix(i,j))   &
             +dtody*(fphiy(i,j-1)-fphiy(i,j))   &
             )*exp(-dt*ch**2/cp**2)
     enddo
  enddo

!-----Step 3.----------------------------------------------------------|
! conserved to primitive
!
  call convert__ctop(ix,jx,gm,ro,ee,rx,ry,rz,bx,by,bz &
                    ,vx,vy,vz,pr)
  call bnd__exec(margin,ix,jx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta)

  enddo

  end subroutine integrate__TVDRK3


end module integrate
