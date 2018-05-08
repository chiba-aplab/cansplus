module flux_calc

  implicit none
  private

  public :: flux_calc__hll, flux_calc__hlld, flux_calc__glm, flux_calc__bp, &
            flux_calc__fbres, flux_calc__feres

  real(8), parameter :: eps=1d-40


contains


  subroutine flux_calc__hll(row,prw,vxw,vyw,vzw,bx,byw,bzw,gm,margin,ix,jx,kx &
                           ,fro,fee,frx,fry,frz,fby,fbz)

  integer,intent(in) :: ix,jx,kx,margin
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
  real(8) :: pbl,pbr,prl,prr
  real(8) :: gmpl,gmpr,gpbl,gpbr

!----- U* ----
! qqlst :: left state
! qqrst :: right state
  real(8) :: sl,sr
!----- flux ---
! fqql :: left physical flux
! fqqr :: right physical flux
! fluxqq :: intermediate HLLD flux (OUTPUT)
  real(8) :: frol,frxl,fryl,frzl
  real(8) :: fbyl,fbzl,feel
  real(8) :: fror,frxr,fryr,frzr
  real(8) :: fbyr,fbzr,feer
  real(8) :: bsql,bsqr
  real(8) :: cfl,cfr
  integer :: i,j,k

  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j) &
  !$OMP PRIVATE(bxs,bxsq) &
  ! --- l state ---
  !$OMP PRIVATE(rol,vxl,vyl,vzl,byl,bzl,bsql,pbl,prl,ptl,eel,rxl,ryl,rzl) &
  ! --- r state ---
  !$OMP PRIVATE(ror,vxr,vyr,vzr,byr,bzr,bsqr,pbr,prr,ptr,eer,rxr,ryr,rzr) &
  !--- step 1----
  !$OMP PRIVATE(gmpl,gmpr,gpbl,gpbr,cfl,cfr,sl,sr) &
  !--- step 2 ---
  ! - left
  !$OMP PRIVATE(frol,frxl,fryl,frzl,feel,fbyl,fbzl) &
  ! - right
  !$OMP PRIVATE(fror,frxr,fryr,frzr,feer,fbyr,fbzr)
  do k=margin,kx-margin
     do j=margin,jx-margin
        do i=margin,ix-margin
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
           
           cfl = sqrt((gpbl + sqrt((gmpl-2.0d0*pbl)**2+4.0d0*gmpl*(byl**2+bzl**2)))/(2.0d0*rol))
           cfr = sqrt((gpbr + sqrt((gmpr-2.0d0*pbr)**2+4.0d0*gmpr*(byr**2+bzr**2)))/(2.0d0*ror))
           
           sl = min(vxl,vxr)-max(cfl,cfr)
           sr = max(vxl,vxr)+max(cfl,cfr)
!----- Step 2. ----------------------------------------------------------|
! compute L/R fluxes
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
  !$OMP END PARALLEL DO
  
  end subroutine flux_calc__hll


  subroutine flux_calc__hlld(row,prw,vxw,vyw,vzw,bx,byw,bzw,gm,margin,ix,jx,kx &
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

  integer,intent(in) :: ix,jx,kx,margin ! array size
  real(8),intent(in) :: gm       ! specific heat rate
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
!  real(8) :: csl,csr,cal,car,cacl,cacr 
  real(8) :: cfl,cfr

!--------------------
! temporary variables
  real(8) :: gmpl,gmpr,gpbl,gpbr
  real(8) :: sdl,sdr,sdml,sdmr,isdml,isdmr,rosdl,rosdr
  real(8) :: vdbstl,vdbstr
  real(8) :: sqrtrol,sqrtror,invsumro
  real(8) :: signbx
  real(8) :: temp,temp1
  real(8) :: temp_fst
  
  integer :: i,j,k  
! no if
  real(8) :: sign1,igm,itf,maxs1,mins1,abbx
  real(8) :: msl,mslst,msrst,msr

  igm = 1.0d0/(gm-1.0d0)


  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j) &
  !$OMP PRIVATE(bxs,bxsq) &
  ! --- l state ---
  !$OMP PRIVATE(rol,vxl,vyl,vzl,byl,bzl,pbl,prl,ptl,rxl,ryl,rzl,eel) &
  ! --- r state ---
  !$OMP PRIVATE(ror,vxr,vyr,vzr,byr,bzr,pbr,prr,ptr,rxr,ryr,rzr,eer) &
  !--- step 1----
  !$OMP PRIVATE(gmpl,gmpr,gpbl,gpbr,cfl,cfr,sl,sr) &
  !--- step 2 ---
  ! - left
  !$OMP PRIVATE(frol,frxl,fryl,frzl,feel,fbyl,fbzl) &
  ! - right
  !$OMP PRIVATE(fror,frxr,fryr,frzr,feer,fbyr,fbzr) &
  !--- step 4 ---
  !$OMP PRIVATE(sdl,sdr,rosdl,rosdr,temp,sm,sdml,sdmr) &
  !--- step 5 ---
  !$OMP PRIVATE(ptst)   &
  !--- step 5a ---
  !$OMP PRIVATE(temp_fst,sign1,maxs1,mins1,itf,isdml,rolst,vxlst,rxlst) &
  !$OMP PRIVATE(vylst,rylst,vzlst,rzlst,bylst,bzlst,vdbstl,eelst) &
  !--- step 5b ---
  !$OMP PRIVATE(isdmr,rorst,vxrst,rxrst) &
  !$OMP PRIVATE(vyrst,ryrst,vzrst,rzrst,byrst,bzrst,vdbstr,eerst) &
  !--- step 5c ---
  !$OMP PRIVATE(sqrtrol,sqrtror,abbx,slst,srst,signbx,invsumro,roldst,rordst) &
  !$OMP PRIVATE(rxldst,rxrdst,vxldst,vxrdst,vyldst,vyrdst,ryldst,ryrdst) &
  !$OMP PRIVATE(vzldst,vzrdst,rzldst,rzrdst,byldst,byrdst,bzldst,bzrdst) &
  !$OMP PRIVATE(eeldst,eerdst) &
  !--- step 6 ---
  !$OMP PRIVATE(msl,mslst,msrst,msr,temp1) 
  do k=margin,kx-margin
     do j=margin,jx-margin
        do i=margin,ix-margin

!----- Step 0. ----------------------------------------------------------|
! set L/R-state

           bxs = bx(i,j,k)
           bxsq = bxs*bxs

!---- Left state
           rol = row(i,j,k,1)
           vxl = vxw(i,j,k,1)
           vyl = vyw(i,j,k,1)
           vzl = vzw(i,j,k,1)
           byl = byw(i,j,k,1)
           bzl = bzw(i,j,k,1)
           
           pbl = 0.5d0*(bxsq+byl*byl+bzl*bzl)
           prl = prw(i,j,k,1)
           ptl = prl + pbl

           rxl = rol*vxl
           ryl = rol*vyl
           rzl = rol*vzl
           
           eel = prl*igm &
                + 0.5d0*(rxl*vxl+ryl*vyl+rzl*vzl) &
                + pbl
           
!---- Right state
           ror = row(i,j,k,2)
           vxr = vxw(i,j,k,2)
           vyr = vyw(i,j,k,2)
           vzr = vzw(i,j,k,2)
           byr = byw(i,j,k,2)
           bzr = bzw(i,j,k,2)
           
           pbr = 0.5d0*(bxsq+byr*byr+bzr*bzr)
           prr = prw(i,j,k,2)
           ptr = prr + pbr

           rxr = ror*vxr
           ryr = ror*vyr
           rzr = ror*vzr

           eer = prr*igm &
                + 0.5d0*(rxr*vxr+ryr*vyr+rzr*vzr) &
                + pbr
           
!----- Step 1. ----------------------------------------------------------|
! Compute wave left & right wave speed
!
           gmpl = gm*prl
           gmpr = gm*prr
           
           gpbl = gmpl+2.0d0*pbl
           gpbr = gmpr+2.0d0*pbr
           
           cfl = sqrt((gpbl + sqrt((gmpl-2.0d0*pbl)**2+4.0d0*gmpl*(byl**2+bzl**2)))*0.5d0/rol)
           cfr = sqrt((gpbr + sqrt((gmpr-2.0d0*pbr)**2+4.0d0*gmpr*(byr**2+bzr**2)))*0.5d0/ror)
           
           sl = min(vxl,vxr)-max(cfl,cfr)
           sr = max(vxl,vxr)+max(cfl,cfr)

!----- Step 2. ----------------------------------------------------------|
! compute L/R fluxes
!
! Left value
           frol = rxl
           frxl = rxl*vxl + ptl -bxsq
           fryl = rxl*vyl - bxs*byl
           frzl = rxl*vzl - bxs*bzl
           feel = vxl*(eel + ptl -bxsq) - bxs*(vyl*byl + vzl*bzl)
           fbyl = byl*vxl - bxs*vyl
           fbzl = bzl*vxl - bxs*vzl

! Right value
           fror = rxr
           frxr = rxr*vxr + ptr -bxsq
           fryr = rxr*vyr - bxs*byr
           frzr = rxr*vzr - bxs*bzr
           feer = vxr*(eer + ptr -bxsq) - bxs*(vyr*byr + vzr*bzr)
           fbyr = byr*vxr - bxs*vyr
           fbzr = bzr*vxr - bxs*vzr

!----- Step 3. ----------------------------------------------------------|
! return upwind flux
!

!----- Step 4. ----------------------------------------------------------|
! compute middle and alfven wave
!
           sdl = sl - vxl
           sdr = sr - vxr
           rosdl = rol*sdl
           rosdr = ror*sdr
           temp = 1.0d0/(rosdr - rosdl)
           sm = (rosdr*vxr - rosdl*vxl - ptr + ptl)*temp
           
           sdml = sl - sm
           sdmr = sr - sm
        
!----- Step 5. ----------------------------------------------------------|
! compute intermediate states
!
           ptst = (rosdr*ptl-rosdl*ptr+rosdl*rosdr*(vxr-vxl))*temp

!----- Step 5A. ----------------------------------------------------------|
! compute Ul*
!
           temp_fst = rosdl*sdml - bxsq
           sign1 = sign(1.0d0,abs(temp_fst)-eps)
           maxs1 = max(0.0d0,sign1)
           mins1 = min(0.0d0,sign1)
           itf = 1.0d0/(temp_fst+mins1)
           isdml = 1.0d0/sdml

           temp = bxs*(sdl-sdml)*itf
           rolst = maxs1*(rosdl*isdml) &
                - mins1*rol
           vxlst = maxs1*sm &
                - mins1*vxl
           rxlst = rolst*vxlst
           
           vylst = maxs1*(vyl-byl*temp) &
                - mins1*vyl
           rylst = rolst*vylst
           vzlst = maxs1*(vzl-bzl*temp) &
                - mins1*vzl
           rzlst = rolst*vzlst
           
           temp = (rosdl*sdl-bxsq)*itf
           bylst = maxs1*(byl*temp) &
                - mins1*byl
           bzlst = maxs1*(bzl*temp) &
                - mins1*bzl
           
           vdbstl = vxlst*bxs+vylst*bylst+vzlst*bzlst
           eelst = maxs1*((sdl*eel - ptl*vxl &
                + ptst*sm + bxs*(vxl*bxs+vyl*byl+vzl*bzl-vdbstl))*isdml) &
                - mins1*eel
           
!----- Step 5B. ----------------------------------------------------------|
! compute Ur*
!
           temp_fst = rosdr*sdmr - bxsq
           sign1 = sign(1.0d0,abs(temp_fst)-eps)
           maxs1 = max(0.0d0,sign1)
           mins1 = min(0.0d0,sign1)
           itf = 1.0d0/(temp_fst+mins1)
           isdmr = 1.0d0/sdmr
           
           temp = bxs*(sdr-sdmr)*itf
           rorst = maxs1*(rosdr*isdmr) &
                - mins1*ror
           vxrst = maxs1*sm &
                - mins1*vxr
           rxrst = rorst*vxrst
           
           vyrst = maxs1*(vyr-byr*temp) &
                - mins1*vyr
           ryrst = rorst*vyrst
           vzrst = maxs1*(vzr-bzr*temp) &
                - mins1*vzr
           rzrst = rorst*vzrst
           
           temp = (rosdr*sdr-bxsq)*itf
           byrst = maxs1*(byr*temp) &
                - mins1*byr
           bzrst = maxs1*(bzr*temp) &
                - mins1*bzr
           
           vdbstr = vxrst*bxs+vyrst*byrst+vzrst*bzrst
           eerst = maxs1*((sdr*eer - ptr*vxr  &
                + ptst*sm + bxs*(vxr*bxs+vyr*byr+vzr*bzr-vdbstr))*isdmr) &
                - mins1*eer
           
!----- Step 5C. ----------------------------------------------------------|
! compute Ul** and Ur**
!
           sqrtrol = sqrt(rolst)
           sqrtror = sqrt(rorst)

           abbx = abs(bxs)
           slst = sm - abbx/sqrtrol
           srst = sm + abbx/sqrtror

           signbx = sign(1.0d0,bxs)           
           sign1 = sign(1.0d0,abbx-eps)
           maxs1 = max(0d0,sign1)
           mins1 = -min(0d0,sign1)
           invsumro = maxs1/(sqrtrol + sqrtror)

           roldst = rolst
           rordst = rorst

           rxldst = rxlst
           rxrdst = rxrst
           vxldst = vxlst
           vxrdst = vxrst

           temp = invsumro*(sqrtrol*vylst + sqrtror*vyrst + signbx*(byrst-bylst))
           vyldst = vylst*mins1 + temp
           vyrdst = vyrst*mins1 + temp
           ryldst = rylst*mins1 + roldst * temp
           ryrdst = ryrst*mins1 + rordst * temp

           temp = invsumro*(sqrtrol*vzlst + sqrtror*vzrst + signbx*(bzrst-bzlst))
           vzldst = vzlst*mins1 + temp
           vzrdst = vzrst*mins1 + temp
           rzldst = rzlst*mins1 + roldst * temp
           rzrdst = rzrst*mins1 + rordst * temp

           temp = invsumro*(sqrtrol*byrst + sqrtror*bylst  &
                + signbx*sqrtrol*sqrtror*(vyrst-vylst))
           byldst = bylst*mins1 + temp
           byrdst = byrst*mins1 + temp
              
           temp = invsumro*(sqrtrol*bzrst + sqrtror*bzlst &
                + signbx*sqrtrol*sqrtror*(vzrst-vzlst))
           bzldst = bzlst*mins1 + temp
           bzrdst = bzrst*mins1 + temp
              
           temp = sm*bxs + vyldst*byldst + vzldst*bzldst
           eeldst = eelst - sqrtrol*signbx*(vdbstl - temp)*maxs1
           eerdst = eerst + sqrtror*signbx*(vdbstr - temp)*maxs1
              
!----- Step 6. ----------------------------------------------------------|
! compute flux
              
           sign1 = sign(1.0d0,sm)
           maxs1 = max(0.0d0,sign1)
           mins1 = -min(0.0d0,sign1)

           msl = min(sl,0.0d0)
           mslst = min(slst,0.0d0)
           msrst = max(srst,0.0d0)
           msr = max(sr,0.0d0)
           temp = mslst-msl
           temp1 = msrst-msr

           fro(i,j,k) = (frol-msl*rol-rolst*temp+roldst*mslst)*maxs1 &
                +(fror-msr*ror-rorst*temp1+rordst*msrst)*mins1
           frx(i,j,k) = (frxl-msl*rxl-rxlst*temp+rxldst*mslst)*maxs1 &
                +(frxr-msr*rxr-rxrst*temp1+rxrdst*msrst)*mins1
           fry(i,j,k) = (fryl-msl*ryl-rylst*temp+ryldst*mslst)*maxs1 &
                +(fryr-msr*ryr-ryrst*temp1+ryrdst*msrst)*mins1
           frz(i,j,k) = (frzl-msl*rzl-rzlst*temp+rzldst*mslst)*maxs1 &
                +(frzr-msr*rzr-rzrst*temp1+rzrdst*msrst)*mins1
           fee(i,j,k) = (feel-msl*eel-eelst*temp+eeldst*mslst)*maxs1 &
                +(feer-msr*eer-eerst*temp1+eerdst*msrst)*mins1
           fby(i,j,k) = (fbyl-msl*byl-bylst*temp+byldst*mslst)*maxs1 &
                +(fbyr-msr*byr-byrst*temp1+byrdst*msrst)*mins1
           fbz(i,j,k) = (fbzl-msl*bzl-bzlst*temp+bzldst*mslst)*maxs1 &
                +(fbzr-msr*bzr-bzrst*temp1+bzrdst*msrst)*mins1
        end do
     end do
  enddo
  !$OMP END PARALLEL DO
  

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


  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j)
  do k=1,kx
     do j=1,jx
        do i=1,ix
           fbx(i,j,k) = phi_m(i,j,k)
           fphi(i,j,k) = bx_m(i,j,k)*ch**2
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  end subroutine flux_calc__glm


  subroutine flux_calc__bp(ix,jx,kx,bxw,phiw &
                          ,bx_m,phi_m,ch)

  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: ch
  real(8),dimension(ix,jx,kx,2),intent(in) :: bxw,phiw
  real(8),dimension(ix,jx,kx),intent(out) :: bx_m,phi_m

  integer :: i,j,k

  ! calcurate magnetic field & divergence B @ cell surface
  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j)
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
  !$OMP END PARALLEL DO

  end subroutine flux_calc__bp


  subroutine flux_calc__fbres(mdir,margin,ix,jx,kx,fbx,curx,eta,pm &
                             ,fbx_res)

  integer,intent(in) :: ix,jx,kx,mdir,margin
! +1 or -1, consistency between Electric fields and numerical flux
  real(8),intent(in) :: pm
  real(8),dimension(ix,jx,kx),intent(in) :: fbx,curx,eta
  real(8),dimension(ix,jx,kx),intent(out) :: fbx_res

  integer :: i,j,k
  real(8) :: fres

! Average eta*current at cell curface
! Flux at i+1/2
  if(mdir == 1)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,fres)
     do k=margin,kx-margin
        do j=margin,jx-margin
           do i=margin,ix-margin
              fres = 0.5d0*pm*(eta(i,j,k)*curx(i,j,k) &
                   +eta(i+1,j,k)*curx(i+1,j,k))

              fbx_res(i,j,k) = fbx(i,j,k)+fres
           end do
        end do
     end do
     !$OMP END PARALLEL DO
! Flux at j+1/2
  else if(mdir == 2)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,fres)
     do k=margin,kx-margin
        do j=margin,jx-margin
           do i=margin,ix-margin
              fres = 0.5d0*pm*(eta(i,j,k)*curx(i,j,k) &
                   +eta(i,j+1,k)*curx(i,j+1,k))
              fbx_res(i,j,k) = fbx(i,j,k)+fres
           end do
        end do
     end do
     !$OMP END PARALLEL DO
! Flux at k+1/2
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,fres)
     do k=margin,kx-margin
        do j=margin,jx-margin
           do i=margin,ix-margin
              fres = 0.5d0*pm*(eta(i,j,k)*curx(i,j,k) &
                   +eta(i,j,k+1)*curx(i,j,k+1))
              fbx_res(i,j,k) = fbx(i,j,k)+fres
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  end subroutine flux_calc__fbres


  subroutine flux_calc__feres(mdir,margin,ix,jx,kx,fee,curx,cury,curz,bx,by,bz,eta &
                            ,fee_res)

  integer,intent(in) :: ix,jx,kx,mdir,margin
  real(8),dimension(ix,jx,kx),intent(in) :: fee,curx,cury,curz,bx,by,bz,eta
  real(8),dimension(ix,jx,kx),intent(out) :: fee_res

  integer :: i,j,k
  real(8) :: fres

! Average eta*current at cell curface
! Flux at i+1/2
  if(mdir == 1)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,fres)
     do k=margin,kx-margin
        do j=margin,jx-margin
           do i=margin,ix-margin
              fres =0.5d0*(eta(i,j,k)*(cury(i,j,k)*bz(i,j,k)-curz(i,j,k)*by(i,j,k)) &
                   +eta(i+1,j,k)*(cury(i+1,j,k)*bz(i+1,j,k)-curz(i+1,j,k)*by(i+1,j,k)))
              fee_res(i,j,k) = fee(i,j,k)+fres
           end do
        end do
     end do
     !$OMP END PARALLEL DO
! Flux at j+1/2
  else if(mdir == 2)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,fres)
     do k=margin,kx-margin
        do j=margin,jx-margin
           do i=margin,ix-margin
              fres =0.5d0*(eta(i,j,k)*(curz(i,j,k)*bx(i,j,k)-curx(i,j,k)*bz(i,j,k)) &
                   +eta(i,j+1,k)*(curz(i,j+1,k)*bx(i,j+1,k)-curx(i,j+1,k)*bz(i,j+1,k)))
              fee_res(i,j,k) = fee(i,j,k)+fres
           end do
        end do
     end do
     !$OMP END PARALLEL DO
! Flux at k+1/2
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,fres)
     do k=margin,kx-margin
        do j=margin,jx-margin
           do i=margin,ix-margin
              fres =0.5d0*(eta(i,j,k)*(curx(i,j,k)*by(i,j,k)-cury(i,j,k)*bx(i,j,k)) &
                   +eta(i,j,k+1)*(curx(i,j,k+1)*by(i,j,k+1)-cury(i,j,k+1)*bx(i,j,k+1)))
              fee_res(i,j,k) = fee(i,j,k)+fres
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  end subroutine flux_calc__feres


end module flux_calc
