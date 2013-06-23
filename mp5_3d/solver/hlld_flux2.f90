!=====================================================================
subroutine hlld_flux2(row,prw,vxw,vyw,vzw,bx,byw,bzw,gm,ix,jx,kx,zfloor &
     ,fro,fee,frx,fry,frz,fby,fbz)
!=====================================================================
!
! Name :: hlld_flux
!
! Purpose
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

  implicit none

!----> INPUT
  integer,intent(in) :: ix,jx,kx ! array size
  real(8),intent(in) :: gm       ! specific heat rate
  real(8),intent(in) :: zfloor 

! primitive variables :: 1 = left state , 2 = right state
! 
  real(8),dimension(ix,jx,kx,2) :: row,prw,vxw,vyw,vzw
  real(8),dimension(ix,jx,kx,2) :: byw,bzw

  real(8),dimension(ix,jx,kx) :: bx ! magnetic field 
!----< INPUT

!----> OUTPUT

  real(8),dimension(ix,jx,kx) :: fro,fee,frx,fry,frz ! numerical flux
  real(8),dimension(ix,jx,kx) :: fby,fbz

!----< OUTPUT

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
  real(8) :: temp_fst,sqrt_fst,temp_fa,sqrt_fa,temp_alf
  real(8) :: alpha_1,alpha_2,alpha_3  
  
  integer :: i,j,k  

  integer :: nanflag
! no if

  real(8) :: sign1,sign2,sign3,sign4
  real(8) :: temp1,temp2

  real(8) :: msl,mslst,msm,msrst,msr

  integer :: flag
  integer :: tflag1,tflag2,tflag3
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
           rolst = rol*sdl/sdml
              
           vxlst = sm
           rxlst = rolst*vxlst

           if( dabs(rol*sdl*sdml/bxsq-1.0d0) < zfloor)then
              vylst = vyl
              rylst = rolst*vylst
              vzlst = vzl
              rzlst = rolst*vzlst

              bylst = byl
              bzlst = bzl
           else
              temp = bxs*(sdl-sdml)/(rol*sdl*sdml-bxsq)
              vylst = vyl-byl*temp
              rylst = rolst*vylst
              vzlst = vzl-bzl*temp
              rzlst = rolst*vzlst
              temp = (rol*sdl**2-bxsq)/(rol*sdl*sdml-bxsq)
              bylst = byl*temp
              bzlst = bzl*temp
           endif
           vdbstl = (vxlst*bxs+vylst*bylst+vzlst*bzlst)
           eelst = (sdl*eel-ptl*vxl+ptst*sm+bxs*(vxl*bxs+vyl*byl+vzl*bzl-vdbstl))/sdml
!----- Step 5B. ----------------------------------------------------------|
! compute Ur*
!
           rorst = ror*sdr/sdmr
              
           vxrst = sm
           rxrst = rorst*vxrst

           if( dabs(ror*sdr*sdmr/bxsq-1.0d0) < zfloor)then
              vyrst = vyr
              ryrst = rorst*vyrst
              vzrst = vzr
              rzrst = rorst*vzrst

              byrst = byr
              bzrst = bzr
           else
              temp = bxs*(sdr-sdmr)/(ror*sdr*sdmr-bxsq)
              vyrst = vyr-byr*temp
              ryrst = rorst*vyrst
              vzrst = vzr-bzr*temp
              rzrst = rorst*vzrst
              temp = (ror*sdr**2-bxsq)/(ror*sdr*sdmr-bxsq)
              byrst = byr*temp
              bzrst = bzr*temp
           endif
           vdbstr = (vxrst*bxs+vyrst*byrst+vzrst*bzrst)
           eerst = (sdr*eer-ptr*vxr+ptst*sm+bxs*(vxr*bxs+vyr*byr+vzr*bzr-vdbstr))/sdmr
           
!----- Step 5C. ----------------------------------------------------------|
! compute Ul** and Ur**
!
           sqrtrol = sqrt(rolst)
           sqrtror = sqrt(rorst)
           
           slst = sm - dabs(bxs)/sqrtrol
           srst = sm + dabs(bxs)/sqrtror

           roldst = rolst
           rordst = rorst

           vxldst = vxlst
           vxrdst = vxrst

           rxldst = rxlst
           rxrdst = rxrst
           
           if (dabs(bxs) < zfloor)then
              vyldst = vylst
              vzldst = vzlst
              ryldst = rylst
              rzldst = rzlst
              eeldst = eelst
              byldst = bylst
              bzldst = bzlst
              
              vyrdst = vyrst
              vzrdst = vzrst
              ryrdst = ryrst
              rzrdst = rzrst
              eerdst = eerst
              byrdst = byrst
              bzrdst = bzrst
           else
              invsumro = 1.0d0/(sqrtrol + sqrtror)
              signbx = sign(1.0d0,bxs)
              
              vxldst = vxlst
              vxrdst = vxrst
              
              temp = invsumro*(sqrtrol*vylst+sqrtror*vyrst+signbx*(byrst-bylst))
              vyldst = temp
              vyrdst = temp
              ryldst = roldst*temp
              ryrdst = rordst*temp
              
              temp = invsumro*(sqrtrol*vzlst+sqrtror*vzrst+signbx*(bzrst-bzlst))
              vzldst = temp
              vzrdst = temp
              rzldst = roldst*temp
              rzrdst = rordst*temp
              
              temp = invsumro*(sqrtrol*byrst+sqrtror*bylst  &
                   +signbx*sqrtrol*sqrtror*(vyrst-vylst))
              byldst = temp
              byrdst = temp
              
              temp = invsumro*(sqrtrol*bzrst+sqrtror*bzlst &
                   +signbx*sqrtrol*sqrtror*(vzrst-vzlst))
              bzldst = temp
              bzrdst = temp
              
              temp =sm*bxs+vyldst*byldst+vzldst*bzldst
              eeldst = eelst-sqrtrol*signbx*(vdbstl-temp)
              eerdst = eerst+sqrtror*signbx*(vdbstr-temp)
              
           endif
!----- Step 6. ----------------------------------------------------------|
! compute flux
!
!!$        msl = min(sl,0.0d0)
!!$        mslst = min(slst,0.0d0)
!!$        msm = min(sm,0.0d0)
!!$        msrst = min(srst,0.0d0)
!!$        msr = min(sr,0.0d0)
!!$
!!$        fro(i,j,k) = frol-msl*rol+ msr*ror &
!!$             -(rolst*(mslst-msl)+roldst*(msm-mslst) &
!!$             +rordst*(msrst-msm)+rorst*(msr-msrst)) 
!!$        frx(i,j,k) = frxl-msl*rxl+ msr*rxr &
!!$             -(rxlst*(mslst-msl)+rxldst*(msm-mslst) &
!!$             +rxrdst*(msrst-msm)+rxrst*(msr-msrst)) 
!!$        fry(i,j,k) = fryl-msl*ryl+ msr*ryr &
!!$             -(rylst*(mslst-msl)+ryldst*(msm-mslst) &
!!$             +ryrdst*(msrst-msm)+ryrst*(msr-msrst)) 
!!$        frz(i,j,k) = frzl-msl*rzl+ msr*rzr &
!!$             -(rzlst*(mslst-msl)+rzldst*(msm-mslst) &
!!$             +rzrdst*(msrst-msm)+rzrst*(msr-msrst)) 
!!$        fee(i,j,k) = feel-msl*eel+ msr*eer &
!!$             -(eelst*(mslst-msl)+eeldst*(msm-mslst) &
!!$             +eerdst*(msrst-msm)+eerst*(msr-msrst)) 
!!$        fby(i,j,k) = fbyl-msl*byl+ msr*byr &
!!$             -(bylst*(mslst-msl)+byldst*(msm-mslst) &
!!$             +byrdst*(msrst-msm)+byrst*(msr-msrst)) 
!!$        fbz(i,j,k) = fbzl-msl*bzl+ msr*bzr &
!!$             -(bzlst*(mslst-msl)+bzldst*(msm-mslst) &
!!$             +bzrdst*(msrst-msm)+bzrst*(msr-msrst)) 
           
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
end subroutine hlld_flux2
