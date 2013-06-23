subroutine hll_flux(row,prw,vxw,vyw,vzw,bx,byw,bzw,gm,ix,jx,kx,zfloor &
     ,fro,fee,frx,fry,frz,fby,fbz)

  implicit none
  
  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: gm

  real(8),intent(in) :: zfloor

! primitive variables :: 1 = left state , 2 = right state
! 
  real(8),dimension(ix,jx,kx,2),intent(in) :: row,prw,vxw,vyw,vzw
  real(8),dimension(ix,jx,kx,2),intent(in) :: byw,bzw

  real(8),dimension(ix,jx,kx),intent(in) :: bx ! magnetic field 

!----> OUTPUT

  real(8),dimension(ix,jx,kx) :: fro,fee,frx,fry,frz ! numerical flux
  real(8),dimension(ix,jx,kx) :: fby,fbz

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
  return
end subroutine hll_flux
