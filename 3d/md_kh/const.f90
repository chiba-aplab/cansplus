module const

  implicit none

! physical constants
  real(8),parameter :: pi = acos(-1.0d0), pi2 = 2d0*pi
  real(8),parameter :: gm = 5d0/3d0 ! specific heat retio
  real(8),parameter :: safety = 0.3d0 ! CFL number

! time control parameters
  logical           :: restart = .false. ! if .true. then start from restart data
  real(8),parameter :: tend = 200.0d0, dtout = tend/400
  integer,parameter :: nstop = 100000
  real(8),parameter :: dtmin = 1d-10! minimum time step
  character(*),parameter :: input_dir = "./data/", output_dir = "./data/"
! grid & MPI
  integer,parameter :: margin = 3 ! for 5th order interpolation
  integer,parameter :: ix = 64+2*margin,jx=40+2*margin,kx=64+2*margin
  integer,parameter :: mpisize_x = 2, mpisize_y = 4,mpisize_z = 2
  integer,parameter :: igx = ix*mpisize_x-2*margin*(mpisize_x-1)
  integer,parameter :: jgx = jx*mpisize_y-2*margin*(mpisize_y-1)
  integer,parameter :: kgx = kx*mpisize_z-2*margin*(mpisize_z-1)
  real(8),parameter :: xmin = 0.d0, ymin = -10.d0, zmin =0.d0
  real(8),parameter :: xmax = 16.d0, ymax = +10.d0, zmax = 16.d0
  real(8),parameter :: dxg0 = (xmax-xmin)/real(igx-margin*2)
  real(8),parameter :: dyg0 = (ymax-ymin)/real(jgx-margin*2)
  real(8),parameter :: dzg0 = (zmax-zmin)/real(kgx-margin*2)
  !TRUE if periodic boundary condition is applied. (1:x, 2:y, 3:z)
  logical,parameter :: pbcheck(3) = (/.true., .false., .true./)

! Parameters for the Kelvin-Helmholtz instability
  real(8),parameter :: ro0   = 1.d0  !at Y=ymax
  real(8),parameter :: b0    = 1.d0  !at Y=ymax
  real(8),parameter :: beta  = 100.d0 !at Y=ymax, Pressure in B0^2/4pi = beta/2
  real(8),parameter :: v0    = sqrt(1.d0+0.5*gm*beta)  !Velocity difference in VA0
  real(8),parameter :: rr    = 1.d0 !rho ratio (ro(y=ymin)/ro(y=ymax))
  real(8),parameter :: br    = 1.d0  !field strength ratio
  real(8),parameter :: theta = pi/2. !field elevation angle
  real(8),parameter :: lmd   = 1.d0  !velocity shear width

end module const
