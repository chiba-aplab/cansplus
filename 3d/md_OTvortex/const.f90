module const

  implicit none

! physical constants
  real(8),parameter :: pi = acos(-1.0d0), pi2 = 2d0*pi
  real(8),parameter :: gm = 5d0/3d0 ! specific heat retio
  real(8),parameter :: safety = 0.3d0 ! CFL number

! time control parameters
  logical           :: restart = .false. ! if .true. then start from restart data
  real(8),parameter :: tend = pi, dtout = tend/10
  integer,parameter :: nstop = 100000
  real(8),parameter :: dtmin = 1d-10! minimum time step
  character(*),parameter :: input_dir = "./data/", output_dir = "./data/"
! grid & MPI
  integer,parameter :: margin = 3 ! for 5th order interpolation
  integer,parameter :: ix = 50+2*margin,jx=100+2*margin,kx=1+2*margin
  integer,parameter :: mpisize_x = 4, mpisize_y = 2,mpisize_z = 1
  integer,parameter :: igx = ix*mpisize_x-2*margin*(mpisize_x-1)
  integer,parameter :: jgx = jx*mpisize_y-2*margin*(mpisize_y-1)
  integer,parameter :: kgx = kx*mpisize_z-2*margin*(mpisize_z-1)
  real(8),parameter :: xmin = 0.d0, ymin = 0.d0, zmin =0.d0
  real(8),parameter :: xmax = pi2, ymax = pi2, zmax = 1.D0
  real(8),parameter :: dxg0 = (xmax-xmin)/real(igx-margin*2)
  real(8),parameter :: dyg0 = (ymax-ymin)/real(jgx-margin*2)
  real(8),parameter :: dzg0 = (zmax-zmin)/real(kgx-margin*2)
  !TRUE if periodic boundary condition is applied. (1:x, 2:y, 3:z)
  logical,parameter :: pbcheck(3) = (/.true., .true., .false./)

end module const
