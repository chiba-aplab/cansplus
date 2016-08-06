module const

  implicit none
 
! physical constants
  real(8),parameter :: pi= acos(-1.0d0) 
  real(8),parameter :: pi2=2.0d0*pi

!   time control parameters
  logical,parameter :: restart = .false.     ! if ".true." then start from restart data
  integer,parameter :: nstop = int(2.d6)       !number of total time steps for the run
  real(8),parameter :: tend  = 40.d0
  real(8),parameter :: dtout = 1.d0
  real(8),parameter :: safety=0.2d0
  real(8),parameter :: dtmin=1.d-10
  character(*),parameter :: input_dir = "./data/", output_dir = "./data/"

! CELL & MPI
  integer,parameter :: margin=3
  integer,parameter :: ix = 40+2*margin,jx=100+2*margin
  integer,parameter :: mpisize_x=2, mpisize_y=2
  integer,parameter :: igx = ix*mpisize_x-2*margin*(mpisize_x-1)
  integer,parameter :: jgx = jx*mpisize_y-2*margin*(mpisize_y-1)
  !TRUE if periodic boundary condition is applied. (1:x, 2:y, 3:z)
  logical,parameter :: pbcheck(2) = (/.true., .false./)
  ! size (uniform grid)
  real(8),parameter :: xmin = -20.d0, ymin = -25.0d0
  real(8),parameter :: xmax = -xmin,  ymax = -ymin  
  real(8),parameter :: dxg0 = (xmax-xmin)/real(igx-margin*2)
  real(8),parameter :: dyg0 = (ymax-ymin)/real(jgx-margin*2)

! model parameter
  ! temperature
  real(8),parameter :: ytr=6.0d0
  real(8),parameter :: wtr=0.5d0
  real(8),parameter :: tcor=25.d0
  ! pressure ratio
  real(8),parameter :: alpha=1.0d0   ! inverse beta , pressure ratio magnetic/gas 
  real(8),parameter :: yf1=-2.5d0,yf2=2.5d0,wf0=0.5d0
  real(8),parameter :: gm=1.05d0
  ! gravity
  real(8),parameter :: g0 = 1.d0/gm
  real(8),parameter :: wg = 0.5d0
  ! maximum physical value @y=0.d0
  real(8),parameter :: ro0=1.0d0, pr0=1.d0/gm
  ! perturbation
  real(8),parameter ::  amp=0.05d0
  real(8),parameter ::  xlamd=20.0d0,ylamd=0.d0
  real(8),parameter ::  yptb1=-4.d0,yptb2=-1.d0,yptb3=1.d0,yptb4=4.d0
  real(8),parameter ::  wptb1=0.5d0,wptb2=0.5d0
  ! resistivity
  real(8),parameter :: eta0 = 0.0d0  ! upper limit of resistivity
  real(8),parameter :: vc=0.0d0      ! threshold 

end module const 
