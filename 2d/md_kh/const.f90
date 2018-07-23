module const

  implicit none

! physical constants
  real(8),parameter :: pi = acos(-1.0d0), pi2 = 2d0*pi
  real(8),parameter :: gm = 5d0/3d0 ! specific heat retio

! time control parameters
  logical,parameter :: restart = .false. ! if .true. then start from restart data
  integer,parameter :: nstop = 100000
  real(8),parameter :: tend = 100.0d0, dtout = 5
  real(8),parameter :: dtmin = 1d-10! minimum time step
  real(8),parameter :: safety = 0.3d0 ! CFL number
  character(*),parameter :: input_dir = "./data/", output_dir = "./data/"

! Cell & MPI
  integer,parameter :: margin = 3 ! for 5th order interpolation
  integer,parameter :: ix = 63+2*margin,jx=80+2*margin
  integer,parameter :: mpisize_x = 2, mpisize_y = 2
  integer,parameter :: igx = ix*mpisize_x-2*margin*(mpisize_x-1)
  integer,parameter :: jgx = jx*mpisize_y-2*margin*(mpisize_y-1)
  real(8),parameter :: xmin = 0.d0, ymin = -10.d0
  real(8),parameter :: xmax = 2.d0*pi/0.4d0, ymax = +10.d0
  real(8),parameter :: dxg0 = (xmax-xmin)/real(igx-margin*2)
  real(8),parameter :: dyg0 = (ymax-ymin)/real(jgx-margin*2)
  !TRUE if periodic boundary condition is applied. (1:x, 2:y, 3:z)
  logical,parameter :: pbcheck(2) = (/.true., .false./)

! Parameters for the Kelvin-Helmholtz instability
  real(8),parameter :: ro0   = 1.d0  !at Y=ymax
  real(8),parameter :: b0    = 1.d0  !at Y=ymax
  real(8),parameter :: beta  = 2.d0 !at Y=ymax, Pressure in B0^2/4pi = beta/2
  real(8),parameter :: v0    = 1.d0*sqrt(1.d0+0.5*gm*beta)  !Velocity difference in VA0
  real(8),parameter :: rr    = 1.d0 !rho ratio (ro(y=ymin)/ro(y=ymax))
  real(8),parameter :: br    = 1.d0  !field strength ratio
  real(8),parameter :: theta = pi/2. !field elevation angle
  real(8),parameter :: lmd   = 1.d0  !velocity shear width
  real(8),parameter :: eta0  = 0.0d0
  real(8),parameter :: vc    = 0.0d0

end module const
