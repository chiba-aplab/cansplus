! Log
! 2013.07.28 make-up, unimplimentation flag (MC2 or MP5)  
!
subroutine integrate(FLAG_accuracy,mpid,margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,dt,hdt&
 ,gx,gz,floor,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
 ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
 ,roi,pri,vxi,vyi,vzi,bxi,byi,bzi &
 ,eta0,vc,eta,ccx,ccy,ccz,RadCool,te_factor,time,rohalo,swtch_t,xin)

  use mpi_domain_xz
  
  implicit none

  include 'mpif.h'

  integer,intent(in) :: ix,jx,kx,margin
  real(8),intent(in) :: floor

  type(mpidomain),intent(in) :: mpid
  integer :: merr

  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: y,dy
  real(8),dimension(kx),intent(in) :: z,dz
  real(8),dimension(5,2,ix),intent(in) :: ccx
  real(8),dimension(5,2,jx),intent(in) :: ccy
  real(8),dimension(5,2,kx),intent(in) :: ccz

  real(8),intent(in) :: RadCool,te_factor,time,rohalo,swtch_t,xin
  real(8),intent(in) :: ch,cr
  real(8),intent(in) :: dt,hdt,gm,eta0,vc

  real(8),dimension(ix,jx,kx),intent(in) :: gx,gz

!-- Initial state
  real(8),dimension(ix,jx,kx),intent(in) :: roi,pri
  real(8),dimension(ix,jx,kx),intent(in) :: vxi,vyi,vzi
  real(8),dimension(ix,jx,kx),intent(in) :: bxi,byi,bzi

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
  real(8),dimension(ix,jx,kx) :: curx,cury,curz

!-other temporary variables
  integer :: i,j,k
  integer :: FLAG_accuracy

!------2stage version-----------------------------------------------

  if(FLAG_accuracy == 0) then
  call integrate_cyl_mp5_1st(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,hdt &
       ,gx,gz,floor,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
       ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
       ,eta0,vc,eta,ccx,ccy,ccz,RadCool,te_factor,time,rohalo,swtch_t,xin)

  call exchangeMpixz2(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1 &
       ,phi1,merr)

  call  bnd(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1,eta,x)

  call  bnd_absoub(ix,jx,kx,x,z,xin,roi,pri,vxi,vyi,vzi,bxi,byi,bzi  &
        ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1)

  call integrate_cyl_mp5_2nd(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,dt  &
       ,gx,gz,floor,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1  &
       ,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr  &
       ,eta0,vc,eta,ccx,ccy,ccz,RadCool,te_factor,time,rohalo,swtch_t,xin)

  call  exchangeMpixz2(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz  &
        ,phi,merr)

  call  bnd(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,x)

  call  bnd_absoub(ix,jx,kx,x,z,xin,roi,pri,vxi,vyi,vzi,bxi,byi,bzi  &
        ,ro,pr,vx,vy,vz,bx,by,bz)

!  else
  endif

!!$!---1st
!!$  call integrate_cyl_mp5_1(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,dt &
!!$       ,gx,gz,floor,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
!!$       ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
!!$       ,eta0,vc,eta,ccx,ccy,ccz)
!!$
!!$  call exchangeMpixz2(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1 &
!!$       ,phi1,merr)

!!$  call bnd(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1,eta,x)
!!$
!!$  call bnd_absoub(ix,jx,kx,x,z,xin,roi,pri,vxi,vyi,vzi,bxi,byi,bzi &
!!$       ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1)
!!$
!!$!---2nd
!!$  call integrate_cyl_mp5_2(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,dt &
!!$       ,gx,gz,floor,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
!!$       ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
!!$       ,eta0,vc,eta,ccx,ccy,ccz)
!!$
!!$  call exchangeMpixz2(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1 &
!!$       ,phi1,merr)
!!$
!!$  call bnd(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1,eta,x)
!!$
!!$  call bnd_absoub(ix,jx,kx,x,z,xin,roi,pri,vxi,vyi,vzi,bxi,byi,bzi &
!!$       ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1)
!!$
!!$!----3rd
!!$  call integrate_cyl_mp5_3(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,dt &
!!$       ,gx,gz,floor,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
!!$       ,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
!!$       ,eta0,vc,eta,ccx,ccy,ccz)

!----------------------------------------------------------------------|
return
end subroutine integrate
