module init
  use const
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

! gravitation
  real(8),public,dimension(ix,jx,kx) :: gx,gy,gz

  real(8),public,dimension(ix,jx,kx) :: roi,pri
  real(8),public,dimension(ix,jx,kx) :: vxi,vyi,vzi
  real(8),public,dimension(ix,jx,kx) :: bxi,byi,bzi

  real(8),public :: ch,cp,min_dx

!----------------------------------------------------------------------|
!  initialize counters

  integer,public :: nd,ns,ns1,ns2
  real(8),public :: time,timep

!----------------------------------------------------------------------|
!   time control parameters
  integer,public :: mwflag,mw,nt1,nt2
  real(8),public :: dt,dtg


contains


  subroutine initialize
  use lr_state, only : reconstructionConstant
  use model, only : model_setup
  use mpi_setup, only : mpi_setup__init, mpi_setup__init_cyl
  use bnd

  real(8),dimension(0:ix) :: xm
  real(8),dimension(0:jx) :: ym
  real(8),dimension(0:kx) :: zm

!----------------------------------------------------------------------|
!   for MPI
  call mpi_setup__init(mpisize_x,mpisize_y,mpisize_z,pbcheck)
!----------------------------------------------------------------------|

!----------------------------------------------------------------------|
!   setup numerical model (grid, initial conditions, etc.)
!----------------------------------------------------------------------|
  call model_setup(ro,pr,vx,vy,vz,bx,by,bz,phi &
       ,roi,pri,vxi,vyi,vzi,bxi,byi,bzi &
       ,x,dx,xm,y,dy,ym,z,dz,zm,gx,gy,gz,eta,min_dx) 

  call bnd__exec(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,x,y,z)

!----------------------------------------------------------------------|
!  cal reconstruction constant for MP5
!----------------------------------------------------------------------|
  call reconstructionConstant(margin,ix,xm,dx,ccx)
  call reconstructionConstant(margin,jx,ym,dy,ccy)
  call reconstructionConstant(margin,kx,zm,dz,ccz)
!----------------------------------------------------------------------|

!----------------------------------------------------------------------|
!  initialize counters
!----------------------------------------------------------------------|
  nd=1
  time  = 0.0d0
  timep = 0.0d0
  ns    = 0
  dt = tend

  end subroutine initialize

end module init
