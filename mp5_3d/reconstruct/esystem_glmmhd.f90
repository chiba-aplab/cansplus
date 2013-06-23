!=====================================================================
subroutine esystem_glmmhd(lem,rem,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,gm)
!=====================================================================
! version 1.0
!
!=====================================================================
  implicit none

  real(8),intent(in) :: gm,ch
  real(8),intent(in) :: ro,pr,vx,vy,vz
  real(8),intent(in) :: bx,by,bz,phi

  real(8),dimension(9,9) :: lem,rem
  real(8),dimension(9) :: eigenvalue

  real(8) :: roi,btsq,vaxsq,asq
  real(8) :: ct2,tsum,tdif,cf2_cs2
  real(8) :: cfsq,cf,cssq,cs
  real(8) :: bt,bet2,bet3
  real(8) :: alpha_f,alpha_s

  real(8) :: na,qf,qs,af_prm,as_prm

  real(8) :: sqrtro,s,a,af,as
  real(8) :: vax

  real(8) :: temp,ich

  integer :: flag
  roi = 1.0d0/ro
  btsq = by**2+bz**2
  vaxsq = (bx**2)*roi
  asq = gm*pr*roi

! fast and slow speed
  ct2 = btsq*roi
  tsum = vaxsq + ct2 + asq
  tdif = vaxsq + ct2 - asq
  cf2_cs2 = sqrt(tdif**2 + 4.0d0*asq*ct2)

  cfsq = 0.5d0*(tsum + cf2_cs2)
  cf = sqrt(cfsq)

  cssq = asq*vaxsq/cfsq
  cs = sqrt(cssq)

! compute beta

  bt = sqrt(btsq)
  if ( bt .eq. 0.0d0) then
     bet2 = 1.0d0/sqrt(2.0d0)
     bet3 = 1.0d0/sqrt(2.0d0)
  else
     bet2 = by/bt
     bet3 = bz/bt
  endif

! compute alpha

  if(cf2_cs2 .eq. 0.0d0) then
     alpha_f = 1.0d0
     alpha_s = 0.0d0
  else if ((asq -cssq) <= 0.0d0)then
     alpha_f = 0.0d0
     alpha_s = 1.0d0
  else if ((cfsq - asq) <= 0.0d0)then
     alpha_f = 1.0d0
     alpha_s = 0.0d0
  else
     alpha_f = sqrt((asq -cssq)/cf2_cs2)
     alpha_s = sqrt((cfsq -asq)/cf2_cs2)
  endif

! compute Q(s),A(s)

  sqrtro = sqrt(ro)
  s = sign(1.0d0,bx)
  a = sqrt(asq)
  qf = cf*alpha_f*s
  qs = cs*alpha_s*s
  af = a*alpha_f*sqrtro
  as = a*alpha_s*sqrtro

! compute eigenvalue

  vax = sqrt(vaxsq)
  eigenvalue(1) = -ch
  eigenvalue(2) = vx - cf
  eigenvalue(3) = vx - vax
  eigenvalue(4) = vx - cs
  eigenvalue(5) = vx
  eigenvalue(6) = vx + cs
  eigenvalue(7) = vx + vax
  eigenvalue(8) = vx + cf
  eigenvalue(9) = ch

! Right-eigenvector
!  ich = 1.0d0/(ch)

  rem(1,1) = ro*alpha_f
  rem(1,2) = 0.0d0
  rem(1,3) = ro*alpha_s
  rem(1,4) = 1.0d0
  rem(1,5) = rem(1,3)
  rem(1,6) = 0.0d0
  rem(1,7) = 0.0d0
  rem(1,8) = rem(1,1)
  rem(1,9) = 0.0d0

  rem(2,1) = -cf*alpha_f
  rem(2,2) = 0.0d0
  rem(2,3) = -cs*alpha_s
  rem(2,4) = 0.0d0
  rem(2,5) = -rem(2,3)
  rem(2,6) = 0.0d0
  rem(2,7) = 0.0d0
  rem(2,8) = -rem(2,1)
  rem(2,9) = 0.0d0

  rem(3,1) = qs*bet2
  rem(3,2) = -bet3
  rem(3,3) = -qf*bet2
  rem(3,4) = 0.0d0
  rem(3,5) = -rem(3,3)
  rem(3,6) = 0.0d0
  rem(3,7) = bet3
  rem(3,8) = -rem(3,1)
  rem(3,9) = 0.0d0

  rem(4,1) = qs*bet3
  rem(4,2) = bet2
  rem(4,3) = -qf*bet3
  rem(4,4) = 0.0d0
  rem(4,5) = -rem(4,3)
  rem(4,6) = 0.0d0
  rem(4,7) = -bet2
  rem(4,8) = -rem(4,1)
  rem(4,9) = 0.0d0

  rem(5,1) = ro*asq*alpha_f
  rem(5,2) = 0.0d0
  rem(5,3) = ro*asq*alpha_s
  rem(5,4) = 0.0d0
  rem(5,5) = rem(5,3) 
  rem(5,6) = 0.0d0
  rem(5,7) = 0.0d0
  rem(5,8) = rem(5,1)
  rem(5,9) = 0.0d0

  rem(6,1) = 0.0d0
  rem(6,2) = 0.0d0
  rem(6,3) = 0.0d0
  rem(6,4) = 0.0d0
  rem(6,5) = 0.0d0
  rem(6,6) = 1.0d0
  rem(6,7) = 0.0d0
  rem(6,8) = 0.0d0
  rem(6,9) = 1.0d0


  rem(7,1) = as*bet2
  rem(7,2) = -bet3*s*sqrtro
  rem(7,3) = -af*bet2
  rem(7,4) = 0.0d0
  rem(7,5) = rem(7,3)
  rem(7,6) = 0.0d0
  rem(7,7) = rem(7,2) 
  rem(7,8) = rem(7,1)
  rem(7,9) = 0.0d0

  rem(8,1) = as*bet3
  rem(8,2) = bet2*s*sqrtro
  rem(8,3) = -af*bet3
  rem(8,4) = 0.0d0
  rem(8,5) = rem(8,3)
  rem(8,6) = 0.0d0
  rem(8,7) = rem(8,2)
  rem(8,8) = rem(8,1)
  rem(8,9) = 0.0d0

  rem(9,1) = 0.0d0
  rem(9,2) = 0.0d0
  rem(9,3) = 0.0d0
  rem(9,4) = 0.0d0
  rem(9,5) = 0.0d0
  rem(9,6) = -ch
  rem(9,7) = 0.0d0
  rem(9,8) = 0.0d0
  rem(9,9) = ch

! Left eigenvector

  na = 0.5d0/asq
  qf = na*qf
  qs = na*qs
  af_prm = na*af*roi
  as_prm = na*as*roi

  lem(1,1) = 0.0d0
  lem(1,2) = -na*cf*alpha_f
  lem(1,3) = qs*bet2
  lem(1,4) = qs*bet3
  lem(1,5) = na*alpha_f*roi
  lem(1,6) = 0.0d0
  lem(1,7) = as_prm*bet2
  lem(1,8) = as_prm*bet3
  lem(1,9) = 0.0d0

  lem(2,1) = 0.0d0
  lem(2,2) = 0.0d0
  lem(2,3) = -0.5d0*bet3
  lem(2,4) = 0.5d0*bet2
  lem(2,5) = 0.0d0
  lem(2,6) = 0.0d0
  lem(2,7) = -0.5d0*bet3*s/sqrtro
  lem(2,8) = 0.5d0*bet2*s/sqrtro
  lem(2,9) = 0.0d0

  lem(3,1) = 0.0d0
  lem(3,2) = -na*cs*alpha_s
  lem(3,3) = -qf*bet2
  lem(3,4) = -qf*bet3
  lem(3,5) = na*alpha_s*roi
  lem(3,6) = 0.0d0
  lem(3,7) = -af_prm*bet2
  lem(3,8) = -af_prm*bet3
  lem(3,9) = 0.0d0

  lem(4,1) = 1.0d0
  lem(4,2) = 0.0d0
  lem(4,3) = 0.0d0
  lem(4,4) = 0.0d0
  lem(4,5) = -1.0d0/asq
  lem(4,6) = 0.0d0
  lem(4,7) = 0.0d0
  lem(4,8) = 0.0d0
  lem(4,9) = 0.0d0

  lem(5,1) = 0.0d0
  lem(5,2) = -lem(3,2)
  lem(5,3) = -lem(3,3)
  lem(5,4) = -lem(3,4)
  lem(5,5) = lem(3,5) 
  lem(5,6) = 0.0d0
  lem(5,7) = lem(3,7)
  lem(5,8) = lem(3,8)
  lem(5,9) = 0.0d0

  lem(6,1) = 0.0d0
  lem(6,2) = 0.0d0
  lem(6,3) = 0.0d0
  lem(6,4) = 0.0d0
  lem(6,5) = 0.0d0
  lem(6,6) = 0.5d0
  lem(6,7) = 0.0d0
  lem(6,8) = 0.0d0
  lem(6,9) = -0.5d0/ch

  lem(7,1) = 0.0d0
  lem(7,2) = 0.0d0
  lem(7,3) = -lem(2,3)
  lem(7,4) = -lem(2,4)
  lem(7,5) = 0.0d0
  lem(7,6) = 0.0d0
  lem(7,7) = lem(2,7)
  lem(7,8) = lem(2,8)
  lem(7,9) = 0.0d0

  lem(8,1) = 0.0d0
  lem(8,2) = -lem(1,2)
  lem(8,3) = -lem(1,3)
  lem(8,4) = -lem(1,4)
  lem(8,5) = lem(1,5)
  lem(8,6) = 0.0d0
  lem(8,7) = lem(1,7)
  lem(8,8) = lem(1,8)
  lem(8,9) = 0.0d0

  lem(9,1) = 0.0d0
  lem(9,2) = 0.0d0
  lem(9,3) = 0.0d0
  lem(9,4) = 0.0d0
  lem(9,5) = 0.0d0
  lem(9,6) = 0.5d0
  lem(9,7) = 0.0d0
  lem(9,8) = 0.0d0
  lem(9,9) = 0.5d0/ch

end subroutine esystem_glmmhd
