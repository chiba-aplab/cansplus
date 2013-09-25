module const

 implicit none
 
  integer,parameter :: margin=3
  integer,parameter :: ix = 256+2*margin,jx=64+2*margin,kx=16+2*margin

  real(8),parameter :: cr=0.18d0

! physical parameter
  real(8),parameter :: pi= acos(-1.0d0)            ! circular constant
  real(8),parameter :: pi2=2.0d0*pi
! real(8),parameter :: pi4=4.0d0*pi, hpi4=sqrt(pi4)
  real(8),parameter :: gm=5.d0/3.d0                ! specific heat retio
  real(8),parameter :: beta0=100.d0                ! pressure ratio or plasma beta

  real(8),parameter :: floor= 1.d-6                 ! minimum value

!--------------------------------------------------------------------------
!  MPI
  integer,parameter :: mpisize_x=1
  integer,parameter :: mpisize_y=4                  ! unimplement for mpi
  integer,parameter :: mpisize_z=4
  integer,parameter :: igx = ix*mpisize_x-2*margin*(mpisize_x-1)
  integer,parameter :: jgx = jx*mpisize_y-2*margin*(mpisize_y-1)
  integer,parameter :: kgx = kx*mpisize_z-2*margin*(mpisize_z-1)

!--------------------------------------------------------------------------
!   time control parameters

  logical :: restart = .false.     ! if ".true." then start from restart data
  
  real(8),parameter :: tend  = 0.01d0
  real(8),parameter :: dtout = tend/(5.0d0)
  integer,parameter :: nstop = 20            !number of total time steps for the run

!  cooling switching     here モデル依存?
  real(8),parameter :: swtch_t=62.8d0

  real(8),parameter :: safety=0.2d0
  real(8),parameter :: dtmin=1.d-10
!--------------------------------------------------------------------------
! real(8),parameter :: ~~~



! model parameter

  !  size 
  real(8),parameter :: xmin = 0.0d0
!  real(8),parameter :: xmax
  real(8),parameter :: ymin = 0.0d0
  real(8),parameter :: ymax = pi2
  real(8),parameter :: zmin = 0.0d0  
!  real(8),parameter :: zmax = 0.0d0  

  ! set uniform grid
  real(8),parameter :: dxg0 = 0.01d0
  real(8),parameter :: dyg0 = (ymax-ymin)/real(jgx-margin*2)
  real(8),parameter :: dzg0 = 0.01d0

  ! set global non-uniform grid  
  real(8),parameter :: dxmax = 10.0d0*dxg0
  real(8),parameter :: ratio_x = 1.05d0     ! in case of uniform, ratio_x=1.0d0
  integer,parameter :: ugrid_xmax=96        ! in case of uniform, ugrid_xmax=0  
  real(8),parameter :: dzmax = 10.0d0*dxg0
  real(8),parameter :: ratio_z = 1.05d0     ! in case of uniform, ratio_z=1.0d0
  integer,parameter :: ugrid_zmax=96        ! in case of uniform, ugrid_zmax=0 


!--cooling
  real(8),parameter :: velocity_c = 2.9979d0*(10.0d0**10)         !cgs
  real(8),parameter :: mass_bh = 10.0d0                           !per mass_solar
  real(8),parameter :: mass_proton = 1.6726d0*(10.0d0**(-24))     !cgs
  real(8),parameter :: boltzmann_const = 1.3807d0*(10.0d0**(-16)) !cgs
  real(8),parameter :: Navo = 6.0221d0*(10.0d0**23)               !Avogadro's constant
  real(8),parameter :: mmw = 1.0d0                                !mean molecular weight
  real(8),parameter :: rg = 3.0d0*(10.0d0**6)*(mass_bh/10.0d0)    !cgs Schwarzschild radius
  real(8),parameter :: rg_nrmlx = 1.0d0/10.0d0                    !per nrmlx

  real(8),parameter :: nrmlro = 0.29d0*8.3d0*(10.0d0**(-7))/(mass_bh/10.0d0)  !normalized ro
  real(8),parameter :: nrmlte = mass_proton*(velocity_c**2)/boltzmann_const   !normalized te
  real(8),parameter :: nrmlx = rg/rg_nrmlx                                    !normalized x
  real(8),parameter :: nrmlv =sqrt(0.5d0*rg_nrmlx/((1.0d0-rg_nrmlx)**2))*velocity_c    !normalized v
  real(8),parameter :: nrmlt = nrmlx/nrmlv                                    !normalized t
  real(8),parameter :: nrmlee = nrmlro*nrmlv*nrmlv                            !normalized ee

  real(8),parameter :: RadCool =nrmlro*sqrt(nrmlte)*nrmlx/(nrmlv**3)*(6.2d0*(10.0d0**20)) !cooling factor
  real(8),parameter :: te_factor =mmw*(nrmlv**2)/Navo/mass_proton/(velocity_c**2)         !factor for te


  real(8),parameter :: xin = 0.2d0   ! inner boundary
  ! resistivity
  real(8),parameter :: eta0 = 4.0d0*pi*1.0d-4      ! upper limit of resistivity
  real(8),parameter :: vc=0.2d0*2.998d10/nrmlv!here! threshold 
  ! gravity
  real(8),parameter :: ssg=rg_nrmlx    ! 
  real(8),parameter :: sseps = 0.2d0   ! switch gravity model for Black hole center (or innerboundary)
  real(8),parameter :: grav = (1.0d0-ssg)**2

  ! corona (or halo)
  real(8),parameter :: factorc=boltzmann_const*Navo/mmw*nrmlte/nrmlv/nrmlv
  real(8),parameter :: tec0 = 5.0d-1
  real(8),parameter :: rohalo = 1.0d-4  ! density ratio of halo to torus
  ! torus 
  real(8),parameter :: aa=0.2d0   ! dependence of raius on angular momentum
  real(8),parameter :: kk=0.05d0  ! ratio of thermal energy to bulk energy     


end module const 
