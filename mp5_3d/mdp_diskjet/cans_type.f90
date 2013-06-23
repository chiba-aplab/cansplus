module cans_type
!======================================================================!
! NAME: cans_type
!
! PURPOSE: Set types of variables
!
! CONTENTS: kinds
!!           coord 
!!           primitive_variable
!!           conserved_variable
!!           constant
!
! UPDATE LOG:
!======================================================================!

  implicit none 
  
!-- Symbolic names for kind types of sigle and double precision real --!
  integer,parameter :: sp = selected_real_kind(6,37)
  integer,parameter :: dp = selected_real_kind(15,307)

!-- Symbolic names for kind types of 
!                                single and double precision complex --!
  integer,parameter :: spc = kind((1.0,1.0))
  integer,parameter :: dpc = kind((1.0d0,1.0d0))

!-- Symbolic name for kind type of default logical --!
  integer,parameter :: lgt = kind(.true.)
  
!-- Set user's choice of kind type (default: double precision) --!
!   use sp for single precision, dp for double precision
  integer,parameter :: real_kind = dp
  integer,parameter :: cmplx_kind = dpc

!-- Frequently used mathematical constants --!
  real(real_kind),parameter :: pi   = acos(-1.0d0)
  real(real_kind),parameter :: pi2  = 1.0d0 * pi
  real(real_kind),parameter :: pi4i = 1.0d0 / pi / 4.0d0
  real(real_kind),parameter :: pi8i = 1.0d0 / pi / 8.0d0



!-- Set structure type for cell --!
!      ct: position of cell center
!      if: position of cell interface
!      iv: cell interval (the interval between cell centers)
!      sz: cell size (the interval between cell interfaces)
!  type coord
!     real(real_kind) :: ct,if,iv,sz
!  end type coord


!-- Set structure type for primitive variables in 3D MHD --!
!      ro : density
!      vx : x-velocity 
!      vy : y-velocity
!      bx : x-magnetic field
!      by : y-magnetic field
!      pr : pressure
!      er : radiative energy density
  type primitive_variable
     real(real_kind) :: ro,vx,vy,bx,by,pr,er
  end type primitive_variable


!-- Set structure type for conserved variables in 1D Rad.-MHD --!
!      ro : density
!      mx : x-momentum
!      my : y-momnetum
!      bx : x-magnetic field
!      by : y-magnetic field
!      ee : internal + kinetic + magnetic energy
!      et : internal + kinetic + magnetic + radiation energy
!      er : radiative energy density
!  type conserved_variable
!     real(real_kind) :: ro,mx,my,bx,by,ee,et,er
!  end type conserved_variable


!-- Set structure type for physical constants in 1D Rad.-MHD --!
!      gm : polytoropic index (specific heat ratio), gamma
!      mu : dimensionless mean molecular weight, mu
!      ke : opacity for the electron scattering, kappa_e
!      lu : unit of length scale, l
!      vu : unit of velocity, v_u
!      kb : Boltzmann constant, k_b
!      ms : mass of particle, m
!      ar : radiation constant, a_R 
!  type constant
!     real(real_kind) :: gm,mu,ke,lu,vu,kb,ms,ar
!  end type constant


end module cans_type
