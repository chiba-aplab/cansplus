module const

  implicit none

! physical constants
  real(8),parameter :: pi = acos(-1.0d0), pi2 = 2d0*pi
  real(8),parameter :: gm = 5d0/3d0 ! specific heat retio

! time control parameters
  logical,parameter :: restart = .false. ! if .true. then start from restart data
  integer,parameter :: nstop = 10000
  real(8),parameter :: tend = 250.0d0, dtout = 10.0d0 
  real(8),parameter :: dtmin = 1d-10! minimum time step
  real(8),parameter :: safety = 0.3d0 ! CFL number
  character(*),parameter :: input_dir = "./data/", output_dir = "./data/"

! Cell & MPI
  integer,parameter :: margin = 3 ! for 5th order interpolation
  integer,parameter :: ix = 75+2*margin,jx=75+2*margin,kx=1+2*margin
  integer,parameter :: mpisize_x = 8, mpisize_y = 1, mpisize_z = 1
  integer,parameter :: igx = ix*mpisize_x-2*margin*(mpisize_x-1)
  integer,parameter :: jgx = jx*mpisize_y-2*margin*(mpisize_y-1)
  integer,parameter :: kgx = kx*mpisize_z-2*margin*(mpisize_z-1)
  real(8),parameter :: xmin = 0.d0, ymin = 0.d0, zmin =0.d0
  real(8),parameter :: xmax = 120.d0, ymax = 15.d0, zmax = 1.d0
  real(8),parameter :: dxg0 = (xmax-xmin)/real(igx-margin*2)
  real(8),parameter :: dyg0 = (ymax-ymin)/real(jgx-margin*2)
  real(8),parameter :: dzg0 = (zmax-zmin)/real(kgx-margin*2)
  !TRUE if periodic boundary condition is applied. (1:x, 2:y, 3:z)
  logical,parameter :: pbcheck(3) = (/.false., .false., .false./)

! Parameters for magnetic reconnection
  real(8),parameter :: ro0   = 1.d0  !at Y=ymax
  real(8),parameter :: b0    = 1.d0  !at Y=ymax
  real(8),parameter :: beta  = 0.2d0 !at Y=ymax, Pressure in B0^2/4pi = beta/2
  real(8),parameter :: lmd   = 1.d0  !current sheet width
  real(8),parameter :: eta0  = 0.001d0          !resistivity 
  real(8),parameter :: eta1  = 1.0d0/60.0d0     !resistivity 
  real(8),parameter :: amp   = 0.03d0     !perturbation 

end module const
