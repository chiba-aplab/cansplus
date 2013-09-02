module init

  use const
  use mpi_domain_xz

  implicit none
  private

  public :: initialize

  real(8),public,dimension(ix) :: x,dx
  real(8),public,dimension(jx) :: y,dy
  real(8),public,dimension(kx) :: z,dz
  real(8),public,dimension(5,2,ix) :: ccx
  real(8),public,dimension(5,2,jx) :: ccy
  real(8),public,dimension(5,2,kx) :: ccz

! cell center variables
  real(8),public,dimension(ix,jx,kx) :: ro,pr
  real(8),public,dimension(ix,jx,kx) :: vx,vy,vz
  real(8),public,dimension(ix,jx,kx) :: bx,by,bz,phi
  real(8),public,dimension(ix,jx,kx) :: eta

  real(8),public,dimension(ix,jx,kx) :: ro1,pr1
  real(8),public,dimension(ix,jx,kx) :: vx1,vy1,vz1
  real(8),public,dimension(ix,jx,kx) :: bx1,by1,bz1,phi1

! gravitation
  real(8),public,dimension(ix,jx,kx) :: gx,gy,gz
! initial variables
  real(8),public,dimension(ix,jx,kx) :: roi,pri
  real(8),public,dimension(ix,jx,kx) :: vxi,vyi,vzi
  real(8),public,dimension(ix,jx,kx) :: bxi,byi,bzi

  real(8), public :: ch

!----------------------------------------------------------------------|
!  initialize counters
  integer, public :: nd,ns,merr,ns1,ns2,nscount
  real(8), public :: time,timep

!----------------------------------------------------------------------|
!   time control parameters
!     nstop : number of total time steps for the run
  integer, public :: mwflag,mw,nt1,nt2

  real(8), public :: dt,hdt
  real(8), public :: dtg
!-----------------------------------------------------------  

!----------------------------------------------------------------------|
! MPI variable
  type(mpidomain), public :: mpid


contains


subroutine initialize

  use model, only : model_setup

  implicit none

!----------------------------------------------------------------------|
!   for MPI
!  
  call mpi_setup(mpid,mpisize_x,mpisize_z)

!----------------------------------------------------------------------|
!   setup numerical model (grid, initial conditions, etc.)
!
  call model_setup(ro,pr,vx,vy,vz,bx,by,bz,phi &
       ,roi,pri,vxi,vyi,vzi,bxi,byi,bzi &
       ,x,dx,y,dy,z,dz &
       ,gx,gz,eta,ccx,ccy,ccz ) 

  call exchangeMpixz(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz &
         ,phi,merr)
  call bnd(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,x,z &
               ,xin,roi,pri,vxi,vyi,vzi,bxi,byi,bzi)

!----------------------------------------------------------------------|
!  initialize counters

  nd=1
  time  = 0.0d0
  timep = 0.0d0

  ns    = 0
  merr  = 0
  nscount=0

  dt = tend
               
end subroutine initialize


end module init

