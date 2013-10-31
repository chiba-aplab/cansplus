module flux_calc

  implicit none
  private

  public :: flux_calc__hll, flux_calc__hlld, flux_calc__glm, flux_calc__bp, &
            flux_calc__fbres, flux_calc__feres


contains


  subroutine flux_calc__hll(row,prw,vxw,vyw,vzw,bx,byw,bzw,gm,ix,jx,kx &
                           ,fro,fee,frx,fry,frz,fby,fbz)

  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: gm
! primitive variables :: 1 = left state , 2 = right state
  real(8),dimension(ix,jx,kx,2),intent(in) :: row,prw,vxw,vyw,vzw
  real(8),dimension(ix,jx,kx,2),intent(in) :: byw,bzw
  real(8),dimension(ix,jx,kx),intent(in) :: bx ! magnetic field 
  real(8),dimension(ix,jx,kx),intent(out) :: fro,fee,frx,fry,frz ! numerical flux
  real(8),dimension(ix,jx,kx),intent(out) :: fby,fbz

!----- U -----
! qql :: left state
! qqr :: right state
  real(8) :: rol,vxl,vyl,vzl,byl,bzl,ptl,eel 
  real(8) :: ror,vxr,vyr,vzr,byr,bzr,ptr,eer 
  real(8) :: rxl,ryl,rzl 
  real(8) :: rxr,ryr,rzr
  real(8) :: bxs,bxsq
  real(8) :: pbl,pbr,ptst,prl,prr
  real(8) :: gmpl,gmpr,gpbl,gpbr

!----- U* ----
! qqlst :: left state
! qqrst :: right state
  real(8) :: rost,vxst,vyst,vzst,byst,bzst,eest
  real(8) :: rxst,ryst,rzst
  real(8) :: sl,sr
!----- flux ---
! fqql :: left physical flux
! fqqr :: right physical flux
! fluxqq :: intemediate HLLD flux (OUTPUT)
  real(8) :: frol,frxl,fryl,frzl
  real(8) :: fbyl,fbzl,feel
  real(8) :: fror,frxr,fryr,frzr
  real(8) :: fbyr,fbzr,feer
  real(8) :: bsql,bsqr
  real(8) :: cfl,cfr
  integer :: i,j,k

  do k=2,kx-2
     do j=2,jx-2
        do i=2,ix-2

!----- Step 0. ----------------------------------------------------------|
! set L/R-state
!
           bxs = bx(i,j,k)
           bxsq = bxs**2
!---- Left state
        
           rol = row(i,j,k,1)
           vxl = vxw(i,j,k,1)
           vyl = vyw(i,j,k,1)
           vzl = vzw(i,j,k,1)
           byl = byw(i,j,k,1)
           bzl = bzw(i,j,k,1)
           
           bsql = bxs**2+byl**2+bzl**2
           pbl = 0.5d0*(bxs**2 + byl**2 + bzl**2)
           prl = prw(i,j,k,1)
           ptl = prl + pbl
           
           eel = prl/(gm-1.0d0) &
                + 0.5d0*rol*(vxl**2+vyl**2+vzl**2) &
                + pbl
           rxl = rol*vxl
           ryl = rol*vyl
           rzl = rol*vzl
           
!---- Right state

           ror = row(i,j,k,2)
           vxr = vxw(i,j,k,2)
           vyr = vyw(i,j,k,2)
           vzr = vzw(i,j,k,2)
           byr = byw(i,j,k,2)
           bzr = bzw(i,j,k,2)
           
           bsqr = bxs**2+byr**2+bzr**2
           pbr = 0.5d0*(bxs**2 + byr**2 + bzr**2)
           prr = prw(i,j,k,2)
           
           ptr = prr + pbr
           eer = prr/(gm-1.0d0) &
                + 0.5d0*(ror*(vxr**2+vyr**2+vzr**2)) &
                + pbr
           
           rxr = ror*vxr
           ryr = ror*vyr
           rzr = ror*vzr
           
!----- Step 1. ----------------------------------------------------------|
! Compute wave left & right wave speed
!
           
           gmpl = gm*prl
           gmpr = gm*prr
           
           gpbl = gmpl+2.0d0*pbl
           gpbr = gmpr+2.0d0*pbr
           
           cfl = sqrt((gpbl + sqrt(gpbl**2-4.0d0*gmpl*bxsq))/(2.0d0*rol))
           cfr = sqrt((gpbr + sqrt(gpbr**2-4.0d0*gmpr*bxsq))/(2.0d0*ror))
           
           sl = min(vxl,vxr)-max(cfl,cfr)
           sr = max(vxl,vxr)+max(cfl,cfr)
!----- Step 2. ----------------------------------------------------------|
! compute L/R fluxs
!
! Left value
           frol = rxl
           frxl = rxl*vxl + ptl -bxsq
           fryl = rol*vxl*vyl - bxs*byl
           frzl = rol*vxl*vzl - bxs*bzl
           feel = vxl*(eel + ptl -bxsq) - bxs*(vyl*byl + vzl*bzl)
           fbyl = byl*vxl - bxs*vyl
           fbzl = bzl*vxl - bxs*vzl
           
! Right value
           fror = rxr
           frxr = rxr*vxr + ptr -bxsq
           fryr = ror*vxr*vyr - bxs*byr
           frzr = ror*vxr*vzr - bxs*bzr
           feer = vxr*(eer + ptr -bxsq) - bxs*(vyr*byr + vzr*bzr)
           fbyr = byr*vxr - bxs*vyr
           fbzr = bzr*vxr - bxs*vzr

!----- Step 3. ----------------------------------------------------------|
! return upwind flux
!

           if (sl >= 0.0d0) then
              fro(i,j,k) = frol
              frx(i,j,k) = frxl
              fry(i,j,k) = fryl
              frz(i,j,k) = frzl
              fee(i,j,k) = feel
              fby(i,j,k) = fbyl
              fbz(i,j,k) = fbzl
              cycle
           endif
           
           if (sr <= 0.0d0) then
              fro(i,j,k) = fror
              frx(i,j,k) = frxr
              fry(i,j,k) = fryr
              frz(i,j,k) = frzr
              fee(i,j,k) = feer
              fby(i,j,k) = fbyr
              fbz(i,j,k) = fbzr
              cycle
           endif
!----- Step 4. ----------------------------------------------------------|
! compute HLL flux
!
           
           fro(i,j,k) = (sr*frol-sl*fror+sr*sl*(ror-rol))/(sr-sl)
           fee(i,j,k) = (sr*feel-sl*feer+sr*sl*(eer-eel))/(sr-sl)
           
           frx(i,j,k) = (sr*frxl-sl*frxr+sr*sl*(rxr-rxl))/(sr-sl)
           fry(i,j,k) = (sr*fryl-sl*fryr+sr*sl*(ryr-ryl))/(sr-sl)
           frz(i,j,k) = (sr*frzl-sl*frzr+sr*sl*(rzr-rzl))/(sr-sl)
           
           fby(i,j,k) = (sr*fbyl-sl*fbyr+sr*sl*(byr-byl))/(sr-sl)
           fbz(i,j,k) = (sr*fbzl-sl*fbzr+sr*sl*(bzr-bzl))/(sr-sl)
        end do
     end do
  end do
  
  end subroutine flux_calc__hll


  subroutine flux_calc__hlld(row,prw,vxw,vyw,vzw,bx,byw,bzw,gm,ix,jx,kx,zfloor &
                            ,fro,fee,frx,fry,frz,fby,fbz)
!=====================================================================
!
! Purpose
! Calculation of Numerical Flux by HLLD method
!
! Reference
!    T.Miyoshi, K.Kusano, JCP, 208, 315 (2005)
!
! Input
!   ix,jx,kx :: array size
!   row,prw,vxw,vyw,vzw,byw,bzw
!            :: primitive value array(ix,jx,kx,2)
!               qqw(i,j,k,1) = qq(i,j,k)'s left state
!               qqw(i,j,k,2) = qq(i,j,k)'s right state
!   bx       :: magnetic field normal to cell surface
!   gm       :: specific heat rate
! Output
!   f**      :: numerical flux
!=====================================================================

  integer,intent(in) :: ix,jx,kx ! array size
  real(8),intent(in) :: gm       ! specific heat rate
  real(8),intent(in) :: zfloor 
! primitive variables :: 1 = left state , 2 = right state
  real(8),dimension(ix,jx,kx,2),intent(in) :: row,prw,vxw,vyw,vzw
  real(8),dimension(ix,jx,kx,2),intent(in) :: byw,bzw
  real(8),dimension(ix,jx,kx),intent(in) :: bx ! magnetic field 
  real(8),dimension(ix,jx,kx),intent(out) :: fro,fee,frx,fry,frz ! numerical flux
  real(8),dimension(ix,jx,kx),intent(out) :: fby,fbz

!----- U -----
! qql :: left state
! qqr :: right state
  real(8) :: rol,vxl,vyl,vzl,byl,bzl,ptl,eel 
  real(8) :: ror,vxr,vyr,vzr,byr,bzr,ptr,eer 
  real(8) :: rxl,ryl,rzl 
  real(8) :: rxr,ryr,rzr
  real(8) :: bxs,bxsq
  real(8) :: pbl,pbr,ptst,prl,prr
!----- U* ----
! qqlst :: left state
! qqrst :: right state
  real(8) :: rolst,vxlst,vylst,vzlst,bylst,bzlst,eelst
  real(8) :: rorst,vxrst,vyrst,vzrst,byrst,bzrst,eerst
  real(8) :: rxlst,rylst,rzlst
  real(8) :: rxrst,ryrst,rzrst
!----- U** ---
! qqldst :: left state
! qqrdst :: right state
  real(8) :: roldst,vxldst,vyldst,vzldst,byldst,bzldst,eeldst
  real(8) :: rordst,vxrdst,vyrdst,vzrdst,byrdst,bzrdst,eerdst
  real(8) :: rxldst,ryldst,rzldst
  real(8) :: rxrdst,ryrdst,rzrdst

!----- flux ---
! fqql :: left physical flux
! fqqr :: right physical flux
! fluxqq :: intemediate HLLD flux (OUTPUT)
  real(8) :: frol,frxl,fryl,frzl
  real(8) :: fbyl,fbzl,feel
  real(8) :: fror,frxr,fryr,frzr
  real(8) :: fbyr,fbzr,feer

!----- wave speed ---
! sl :: left-going fastest signal velocity
! sr :: right-going fastest signal velocity
! sm :: contact discontinuity velocity
! slst :: left-going alfven velocity
! srst :: right-going alfven velocity
  real(8) :: sm,sl,sr,slst,srst

! cfl :: left-state Fast wave velocity
! cfr :: right-sate Fast wave velocity
! csl :: left-state Slow wave velocity
! csr :: right-state Slow wave velocity
! cal :: left-state Alfven wave velocity
! car :: right-state Alfven wave velocity
! cacl :: left-state Sound wave velocity
! cacr :: right-state Sound wave velocity
  real(8) :: cfl,cfr,csl,csr,cal,car,cacl,cacr 

!--------------------
! temporary variables
  real(8) :: gmpl,gmpr
  real(8) :: gpbl,gpbr
  real(8) :: cfmax,maxspd
  real(8) :: sdl,sdr,sdml,sdmr
  real(8) :: vdbstl,vdbstr
  real(8) :: sqrtrol,sqrtror
  real(8) :: invsumro
  real(8) :: signbx,inverse_sqrtro,msqrtro
  real(8) :: bsq,sqrt_bsq
  real(8) :: bsqr,bsql
  real(8) :: temp,tmp
  real(8) :: signBy,signBz,sign_fst
  real(8) :: temp_fst,sqrt_fst,temp_fa,sqrt_fa
  real(8) :: alpha_1,alpha_2,alpha_3  
  
  integer :: i,j,k  
  integer :: nanflag
! no if
  real(8) :: sign1,sign2,sign3
  real(8) :: temp1,temp2
  real(8) :: msl,mslst,msm,msrst,msr

  do k=2,kx-2
     do j=2,jx-2
        do i=2,ix-2

!----- Step 0. ----------------------------------------------------------|
! set L/R-state
!
           bxs = bx(i,j,k)
           bxsq = bxs**2
!---- Left state
        
           rol = row(i,j,k,1)
           vxl = vxw(i,j,k,1)
           vyl = vyw(i,j,k,1)
           vzl = vzw(i,j,k,1)
           byl = byw(i,j,k,1)
           bzl = bzw(i,j,k,1)
           
           bsql = bxs**2+byl**2+bzl**2
           pbl = 0.5d0*(bxs**2 + byl**2 + bzl**2)
           prl = prw(i,j,k,1)
           ptl = prl + pbl
           
           eel = prl/(gm-1.0d0) &
                + 0.5d0*rol*(vxl**2+vyl**2+vzl**2) &
                + pbl
           rxl = rol*vxl
           ryl = rol*vyl
           rzl = rol*vzl
           
!---- Right state

           ror = row(i,j,k,2)
           vxr = vxw(i,j,k,2)
           vyr = vyw(i,j,k,2)
           vzr = vzw(i,j,k,2)
           byr = byw(i,j,k,2)
           bzr = bzw(i,j,k,2)
           
           bsqr = bxs**2+byr**2+bzr**2
           pbr = 0.5d0*(bxs**2 + byr**2 + bzr**2)
           prr = prw(i,j,k,2)
           
           ptr = prr + pbr
           eer = prr/(gm-1.0d0) &
                + 0.5d0*(ror*(vxr**2+vyr**2+vzr**2)) &
                + pbr
           
           rxr = ror*vxr
           ryr = ror*vyr
           rzr = ror*vzr
           
!----- Step 1. ----------------------------------------------------------|
! Compute wave left & right wave speed
!
           
           gmpl = gm*prl
           gmpr = gm*prr
           
           gpbl = gmpl+2.0d0*pbl
           gpbr = gmpr+2.0d0*pbr
           
           cfl = sqrt((gpbl + sqrt(gpbl**2-4.0d0*gmpl*bxsq))/(2.0d0*rol))
           cfr = sqrt((gpbr + sqrt(gpbr**2-4.0d0*gmpr*bxsq))/(2.0d0*ror))
           
           sl = min(vxl,vxr)-max(cfl,cfr)
           sr = max(vxl,vxr)+max(cfl,cfr)

!----- Step 2. ----------------------------------------------------------|
! compute L/R fluxs
!

! Left value
           frol = rxl
           frxl = rxl*vxl + ptl -bxsq
           fryl = rol*vxl*vyl - bxs*byl
           frzl = rol*vxl*vzl - bxs*bzl
           feel = vxl*(eel + ptl -bxsq) - bxs*(vyl*byl + vzl*bzl)
           fbyl = byl*vxl - bxs*vyl
           fbzl = bzl*vxl - bxs*vzl

! Right value
           fror = rxr
           frxr = rxr*vxr + ptr -bxsq
           fryr = ror*vxr*vyr - bxs*byr
           frzr = ror*vxr*vzr - bxs*bzr
           feer = vxr*(eer + ptr -bxsq) - bxs*(vyr*byr + vzr*bzr)
           fbyr = byr*vxr - bxs*vyr
           fbzr = bzr*vxr - bxs*vzr

!----- Step 3. ----------------------------------------------------------|
! return upwind flux
!
           
           if (sl >= 0.0d0) then
              fro(i,j,k) = frol
              frx(i,j,k) = frxl
              fry(i,j,k) = fryl
              frz(i,j,k) = frzl
              fee(i,j,k) = feel
              fby(i,j,k) = fbyl
              fbz(i,j,k) = fbzl
              cycle
           endif
        
           if (sr <= 0.0d0) then
              fro(i,j,k) = fror
              frx(i,j,k) = frxr
              fry(i,j,k) = fryr
              frz(i,j,k) = frzr
              fee(i,j,k) = feer
              fby(i,j,k) = fbyr
              fbz(i,j,k) = fbzr
              cycle
           endif

!----- Step 4. ----------------------------------------------------------|
! compute middle and alfven wave
!
           sdl = sl - vxl
           sdr = sr - vxr
        
           sm = (sdr*ror*vxr - sdl*rol*vxl - ptr + ptl) &
                /(sdr*ror - sdl*rol)
           
           sdml = sl - sm
           sdmr = sr - sm
        
!----- Step 5. ----------------------------------------------------------|
! compute intermediate states
!

!              ptst = ptl + rol*sdl*(sdl-sdml)
           ptst = (sdr*ror*ptl-sdl*rol*ptr+rol*ror*sdr*sdl*(vxr-vxl)) &
                /(sdr*ror-sdl*rol)

!----- Step 5A. ----------------------------------------------------------|
! compute Ul*
!

           temp_fst = (rol*sdl*sdml - bxsq)
           sign1 = sign(1.0d0,temp_fst+zfloor)
           sign2 = sign(1.0d0,temp_fst-zfloor)
           sign3 = sign1*sign2
           
           tmp = bxs*(sdl-sdml)/(temp_fst+min(0.0d0,sign3))
           rolst = max(0.0d0,sign3)*(rol*sdl/sdml) &
                - min(0.0d0,sign3)*rol
           vxlst = max(0.0d0,sign3)*(sm) &
                - min(0.0d0,sign3)*vxl
           rxlst = rolst*vxlst
           
           vylst = max(0.0d0,sign3)*(vyl-byl*tmp) &
                - min(0.0d0,sign3)*vyl
           rylst = rolst*vylst
           vzlst = max(0.0d0,sign3)*(vzl-bzl*tmp) &
                - min(0.0d0,sign3)*vzl
           rzlst = rolst*vzlst
           
           tmp = (rol*sdl**2-bxsq)/(temp_fst+min(0.0d0,sign3))
           bylst = max(0.0d0,sign3)*(byl*tmp) &
                - min(0.0d0,sign3)*byl
           bzlst = max(0.0d0,sign3)*(bzl*tmp) &
                - min(0.0d0,sign3)*bzl
           
           vdbstl = (vxlst*bxs+vylst*bylst+vzlst*bzlst)
           eelst = max(0.0d0,sign3)*((sdl*eel - ptl*vxl + ptst*sm + bxs*(vxl*bxs+vyl*byl+vzl*bzl-vdbstl))/sdml) &
                - min(0.0d0,sign3)*eel
           
!----- Step 5B. ----------------------------------------------------------|
! compute Ur*
!

           temp_fst = (ror*sdr*sdmr - bxsq)
           sign1 = sign(1.0d0,temp_fst+zfloor)
           sign2 = sign(1.0d0,temp_fst-zfloor)
           sign3 = sign1*sign2
           
           tmp = bxs*(sdr-sdmr)/(temp_fst+min(0.0d0,sign3))
           rorst = max(0.0d0,sign3)*(ror*sdr/sdmr) &
                - min(0.0d0,sign3)*ror
           vxrst = max(0.0d0,sign3)*(sm) &
                - min(0.0d0,sign3)*vxr
           rxrst = rorst*vxrst
           
           vyrst = max(0.0d0,sign3)*(vyr-byr*tmp) &
                - min(0.0d0,sign3)*vyr
           ryrst = rorst*vyrst
           vzrst = max(0.0d0,sign3)*(vzr-bzr*tmp) &
                - min(0.0d0,sign3)*vzr
           rzrst = rorst*vzrst
           
           tmp = (ror*sdr**2-bxsq)/(temp_fst+min(0.0d0,sign3))
           byrst = max(0.0d0,sign3)*(byr*tmp) &
                - min(0.0d0,sign3)*byr
           bzrst = max(0.0d0,sign3)*(bzr*tmp) &
                - min(0.0d0,sign3)*bzr
           
           vdbstr = (vxrst*bxs+vyrst*byrst+vzrst*bzrst)
           eerst = max(0.0d0,sign3)*((sdr*eer - ptr*vxr  &
                + ptst*sm + bxs*(vxr*bxs+vyr*byr+vzr*bzr-vdbstr))/sdmr) &
                - min(0.0d0,sign3)*eer
           
!----- Step 5C. ----------------------------------------------------------|
! compute Ul** and Ur**
!
           sqrtrol = sqrt(rolst)
           sqrtror = sqrt(rorst)
           
           slst = sm - dabs(bxs)/sqrtrol
           srst = sm + dabs(bxs)/sqrtror
           
           if (bxs == 0.0d0)then
              roldst = rolst
              rxldst = rxlst
              ryldst = rylst
              rzldst = rzlst
              vxldst = vxlst
              vyldst = vylst
              vzldst = vzlst
              eeldst = eelst
              byldst = bylst
              bzldst = bzlst
              
              rordst = rorst
              rxrdst = rxrst
              ryrdst = ryrst
              rzrdst = rzrst
              vxrdst = vxrst
              vyrdst = vyrst
              vzrdst = vzrst
              eerdst = eerst
              byrdst = byrst
              bzrdst = bzrst
           else
              invsumro = 1.0d0/(sqrtrol + sqrtror)
              signbx = sign(1.0d0,bxs)
              
              roldst = rolst
              rordst = rorst
              
              rxldst = rxlst
              rxrdst = rxrst
              vxldst = vxlst
              vxrdst = vxrst
              
              temp = invsumro*(sqrtrol*vylst + sqrtror*vyrst + signbx*(byrst-bylst))
              vyldst = temp
              vyrdst = temp
              ryldst = roldst * temp
              ryrdst = rordst * temp
              
              temp = invsumro*(sqrtrol*vzlst + sqrtror*vzrst + signbx*(bzrst-bzlst))
              vzldst = temp
              vzrdst = temp
              rzldst = roldst * temp
              rzrdst = rordst * temp
              
              temp = invsumro*(sqrtrol*byrst + sqrtror*bylst  &
                   + signbx*sqrtrol*sqrtror*(vyrst-vylst))
              byldst = temp
              byrdst = temp
              
              temp = invsumro*(sqrtrol*bzrst + sqrtror*bzlst &
                   + signbx*sqrtrol*sqrtror*(vzrst-vzlst))
              
              bzldst = temp
              bzrdst = temp
              
              temp = sm*bxs + vyldst*byldst + vzldst*bzldst
              eeldst = eelst - sqrtrol*signbx*(vdbstl - temp)
              eerdst = eerst + sqrtror*signbx*(vdbstr - temp)
              
           endif
!----- Step 6. ----------------------------------------------------------|
! compute flux
           if (slst >= 0.0d0)then
              fro(i,j,k) = frol + sl*(rolst - rol)
              frx(i,j,k) = frxl + sl*(rxlst - rxl)
              fry(i,j,k) = fryl + sl*(rylst - ryl)
              frz(i,j,k) = frzl + sl*(rzlst - rzl)
              fee(i,j,k) = feel + sl*(eelst - eel)
              fby(i,j,k) = fbyl + sl*(bylst - byl)
              fbz(i,j,k) = fbzl + sl*(bzlst - bzl)
           else if(sm >= 0.0d0) then
              temp = slst - sl
              fro(i,j,k) = frol - sl*rol - temp*rolst + slst*roldst
              frx(i,j,k) = frxl - sl*rxl - temp*rxlst + slst*rxldst
              fry(i,j,k) = fryl - sl*ryl - temp*rylst + slst*ryldst
              frz(i,j,k) = frzl - sl*rzl - temp*rzlst + slst*rzldst
              fee(i,j,k) = feel - sl*eel - temp*eelst + slst*eeldst
              fby(i,j,k) = fbyl - sl*byl - temp*bylst + slst*byldst
              fbz(i,j,k) = fbzl - sl*bzl - temp*bzlst + slst*bzldst
           else if (srst > 0.0d0) then
              temp = srst - sr
              fro(i,j,k) = fror - sr*ror - temp*rorst + srst*rordst
              frx(i,j,k) = frxr - sr*rxr - temp*rxrst + srst*rxrdst
              fry(i,j,k) = fryr - sr*ryr - temp*ryrst + srst*ryrdst
              frz(i,j,k) = frzr - sr*rzr - temp*rzrst + srst*rzrdst
              fee(i,j,k) = feer - sr*eer - temp*eerst + srst*eerdst
              fby(i,j,k) = fbyr - sr*byr - temp*byrst + srst*byrdst
              fbz(i,j,k) = fbzr - sr*bzr - temp*bzrst + srst*bzrdst
           else
              fro(i,j,k) = fror + sr*(rorst - ror)
              frx(i,j,k) = frxr + sr*(rxrst - rxr)
              fry(i,j,k) = fryr + sr*(ryrst - ryr)
              frz(i,j,k) = frzr + sr*(rzrst - rzr)
              fee(i,j,k) = feer + sr*(eerst - eer)
              fby(i,j,k) = fbyr + sr*(byrst - byr)
              fbz(i,j,k) = fbzr + sr*(bzrst - bzr)
           endif
        end do
     end do
  enddo

  end subroutine flux_calc__hlld


  subroutine flux_calc__glm(bx_m,phi_m,ch,fbx,fphi,ix,jx,kx)

  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: ch ! divergence B wave speed
! cell surface Magnetic field and divergence B
  real(8),dimension(ix,jx,kx), intent(in) :: bx_m,phi_m
! magnetic field flux @ cell surface of normal component
! divergence B flux
  real(8),dimension(ix,jx,kx), intent(out) :: fbx,fphi

  integer :: i,j,k

  do k=1,kx
     do j=1,jx
        do i=1,ix
           fbx(i,j,k) = phi_m(i,j,k)
           fphi(i,j,k) = bx_m(i,j,k)*ch**2
        end do
     end do
  end do

  end subroutine flux_calc__glm


  subroutine flux_calc__bp(ix,jx,kx,bxw,phiw &
                          ,bx_m,phi_m,ch)

  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: ch
  real(8),dimension(ix,jx,kx,2),intent(in) :: bxw,phiw
  real(8),dimension(ix,jx,kx),intent(out) :: bx_m,phi_m

  integer :: i,j,k

! calcurate magnetic field & divergence B @ cell surface
  do k=1,kx
     do j=1,jx
        do i=1,ix
           bx_m(i,j,k) = bxw(i,j,k,1) &
                +(0.5d0*(bxw(i,j,k,2)-bxw(i,j,k,1)) &
                -0.5d0*(phiw(i,j,k,2)-phiw(i,j,k,1))/ch)
           
           phi_m(i,j,k) = phiw(i,j,k,1) &
                +(0.5d0*(phiw(i,j,k,2)-phiw(i,j,k,1)) &
                -0.5d0*ch*(bxw(i,j,k,2)-bxw(i,j,k,1)))
        end do
     end do
  end do

  end subroutine flux_calc__bp


  subroutine flux_calc__fbres(mdir,ix,jx,kx,fbx,curx,eta,pm &
                             ,fbx_res)

  integer,intent(in) :: ix,jx,kx,mdir
! +1 or -1, consistency between Electric fields and numerical flux
  real(8),intent(in) :: pm
  real(8),dimension(ix,jx,kx),intent(in) :: fbx,curx,eta
  real(8),dimension(ix,jx,kx),intent(out) :: fbx_res

  integer :: i,j,k
  real(8) :: fres

! Average eta*current at cell curface
! Flux at i+1/2
  if(mdir == 1)then
     do k=2,kx-2
        do j=2,jx-2
           do i=2,ix-2
              fres = 0.5d0*pm*(eta(i,j,k)*curx(i,j,k) &
                   +eta(i+1,j,k)*curx(i+1,j,k))

              fbx_res(i,j,k) = fbx(i,j,k)+fres
           end do
        end do
     end do
! Flux at j+1/2
  else if(mdir == 2)then
     do k=2,kx-2
        do j=2,jx-2
           do i=2,ix-2
              fres = 0.5d0*pm*(eta(i,j,k)*curx(i,j,k) &
                   +eta(i,j+1,k)*curx(i,j+1,k))
              fbx_res(i,j,k) = fbx(i,j,k)+fres
           end do
        end do
     end do
! Flux at k+1/2
  else
     do k=2,kx-2
        do j=2,jx-2
           do i=2,ix-2
              fres = 0.5d0*pm*(eta(i,j,k)*curx(i,j,k) &
                   +eta(i,j,k+1)*curx(i,j,k+1))
              fbx_res(i,j,k) = fbx(i,j,k)+fres
           end do
        end do
     end do
  end if

  end subroutine flux_calc__fbres


  subroutine flux_calc__feres(mdir,ix,jx,kx,fee,curx,cury,curz,bx,by,bz,eta &
                            ,fee_res)

  integer,intent(in) :: ix,jx,kx,mdir
  real(8),dimension(ix,jx,kx),intent(in) :: fee,curx,cury,curz,bx,by,bz,eta
  real(8),dimension(ix,jx,kx),intent(out) :: fee_res

  integer :: i,j,k
  real(8) :: fres

! Average eta*current at cell curface
! Flux at i+1/2
  if(mdir == 1)then
     do k=2,kx-2
        do j=2,jx-2
           do i=2,ix-2
              fres =0.5d0*(eta(i,j,k)*(cury(i,j,k)*bz(i,j,k)-curz(i,j,k)*by(i,j,k)) &
                   +eta(i+1,j,k)*(cury(i+1,j,k)*bz(i+1,j,k)-curz(i+1,j,k)*by(i+1,j,k)))
              fee_res(i,j,k) = fee(i,j,k)+fres
           end do
        end do
     end do
! Flux at j+1/2
  else if(mdir == 2)then
     do k=2,kx-2
        do j=2,jx-2
           do i=2,ix-2
              fres =0.5d0*(eta(i,j,k)*(curz(i,j,k)*bx(i,j,k)-curx(i,j,k)*bz(i,j,k)) &
                   +eta(i,j+1,k)*(curz(i,j+1,k)*bx(i,j+1,k)-curx(i,j+1,k)*bz(i,j+1,k)))
              fee_res(i,j,k) = fee(i,j,k)+fres
           end do
        end do
     end do
! Flux at k+1/2
  else
     do k=2,kx-2
        do j=2,jx-2
           do i=2,ix-2
              fres =0.5d0*(eta(i,j,k)*(curx(i,j,k)*by(i,j,k)-cury(i,j,k)*bx(i,j,k)) &
                   +eta(i,j,k+1)*(curx(i,j,k+1)*by(i,j,k+1)-cury(i,j,k+1)*bx(i,j,k+1)))
              fee_res(i,j,k) = fee(i,j,k)+fres
           end do
        end do
     end do
  end if

  end subroutine flux_calc__feres


end module flux_calc
