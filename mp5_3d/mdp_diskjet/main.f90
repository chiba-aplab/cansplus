program main
!======================================================================|
! PURPOSE
!    
!
! VARIABLES
!    margin: [integer] margin outside the boundary
!    ix, jx, kx : [integer] grid numbers
!
!    nd      [integer] : counter for the output
!    ns      [integer] : counter for the time step
!    nstop   [integer] : ending step number
!    mpi_err [integer] : error counter
!    mw      [integer] : output flag
!    mwflag  [integer] : error flag
!    mf_**   [integer] : unit specifier for **
!    nt1     [integer] :
!    nt2     [integer] :
!
!    dtout [real(8)] : interval of the output
!    tend  [real(8)] : ending time
!    dtmin [real(8)] : minimum value of dt
!
!    time  [real(8)] :
!    timep [real(8)] : variable to store the time at previous step
!    
!    safety: [real(8)] courant number
!    qav: [real(8)] artificial visosity
!    
!    dt [real(8)] : time interval
!    x  (ix) [real(8)] : x-coordinate at the grid point, i
!    xm (ix) [real(8)] : x-coordinate at the cell surface, i+1/2
!    dx (ix) [real(8)] : surface interval at i
!    dxm(ix) [real(8)] : grid interval at i+1/2
!    y  (jx) [real(8)] : y-coordinate at the grid point, j
!    ym (jx) [real(8)] : y-coordinate at the cell surface, j+1/2
!    dy (jx) [real(8)] : surface interval at j
!    dym(jx) [real(8)] : grid interval at j+1/2
!    z  (kx) [real(8)] : z-coordinate at the grid point, k
!    zm (kx) [real(8)] : z-coordinate at the cell surface, k+1/2
!    dz (kx) [real(8)] : surface interval at k
!    dzm(kx) [real(8)] : grid interval at k+1/2
!
!    ro (ix, jx, kx) [real(8)] : density
!    pr (ix, jx, kx) [real(8)] : pressure
!    vx (ix, jx, kx) [real(8)] : velocity
!    vy (ix, jx, kx) [real(8)] : velocisy
!    vz (ix, jx, kx) [real(8)] : velocisy
!    bx (ix, jx, kx) [real(8)] : magnetic field
!    by (ix, jx, kx) [real(8)] : magnetic field
!    bz (ix, jx, kx) [real(8)] : magnetic field
!    phi(ix, jx, kx) [real(8)] : 
!
!    gm [real(8)] : polytoropic index gamma (specific heat ratio)
!
!>>> for MPI >>>!
!    n_mpi_dims  [integer] : number of dimension of domain decomposition
!    n_tiles     [integer] : number of elements of each dimension
!    periods     [lgt]     : specifier wheter the grid is 
!                              periodic (.true.) or not (.false.) 
!                              in each dimension
!    reorder     [lgt]     : ranking may reordered (.true.) 
!                              or not (.false.)
!<<< for MPI <<<!
!    
!     
! UPDATE LOG: 
! 2012/12/28 implement 3D MPI parallelization in 3D Cartesian topology (by H. ODA)
! 
!======================================================================|
!     array definitions
!======================================================================|

!>>> for MPI >>>!
  use cans_mpi
!<<< for MPI <<<!

  use cans_type

  use openfile


  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer,parameter :: FLAG_ModAverage=0 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
  integer,parameter :: margin = 3

  integer,parameter :: ix = 64+2*margin,jx=64+2*margin,kx=64+2*margin

  integer :: igx,jgx,kgx
  integer :: ix0,jx0,kx0,nx0

  real(8),dimension(ix) :: x,dx
  real(8),dimension(jx) :: y,dy
  real(8),dimension(kx) :: z,dz
  real(8),dimension(5,2,ix) :: ccx
  real(8),dimension(5,2,jx) :: ccy
  real(8),dimension(5,2,kx) :: ccz

! cell center variables
  real(8),dimension(ix,jx,kx) :: ro,pr
  real(8),dimension(ix,jx,kx) :: vx,vy,vz
  real(8),dimension(ix,jx,kx) :: bx,by,bz,phi
  real(8),dimension(ix,jx,kx) :: eta

  real(8),dimension(ix,jx,kx) :: ro1,pr1
  real(8),dimension(ix,jx,kx) :: vx1,vy1,vz1
  real(8),dimension(ix,jx,kx) :: bx1,by1,bz1,phi1
  
! gravitation
  real(8),dimension(ix,jx,kx) :: gx,gy,gz

  real(8),dimension(ix,jx,kx) :: roi,pri
  real(8),dimension(ix,jx,kx) :: vxi,vyi,vzi
  real(8),dimension(ix,jx,kx) :: bxi,byi,bzi

  real(8) :: ch,cr

  real(8) :: gm ! specific heat retio
  real(8) :: eta0,vc

  real(8) :: floor ! minimum value


! test1028
  integer,parameter :: iCalNum = 5
  integer,parameter :: iInterNum = 6
  integer,parameter :: mf_calPhysValue = 401

  real(8),dimension(ix,jx,kx) :: pb,ee,am,rx,ry,rz
  real(8),dimension(iCalNum) :: avVal
  real(8),dimensioN(iInterNum) :: avInterVal

!--av
  real(8) :: baseAV
  real(8) :: pbAV,prAV
  real(8) :: ebAV,roAv,eeAv,amAv

  real(8),parameter :: xmin_MAV=0.7d0
  real(8),parameter :: xmax_MAV=1.3d0
  real(8),parameter :: zmin_MAV=0.0d0
  real(8),parameter :: zmax_MAV=0.3d0
!--tcr
  real(8),dimension(kx) :: baseTCR

  real(8),dimension(kx) :: bxTCR,byTCR,bzTCR
  real(8),dimension(kx) :: rxTCR,ryTCR,rzTCR
  real(8),dimension(kx) :: roTCR,pbTCR,prTCR

  real(8),parameter :: xObs1 = 0.3d0
  real(8),parameter :: xObs2 = 0.5d0
  real(8),parameter :: xObs3 = 0.75d0
  real(8),parameter :: xObs4 = 1.0d0
  real(8),parameter :: xObs5 = 1.5d0

!--cooling
  real(8),parameter :: velocity_c = 2.9979d0*(10.0d0**10)                                       !cgs
  real(8),parameter :: mass_bh = 10.0d0                                                         !per mass_solar
  real(8),parameter :: mass_proton = 1.6726d0*(10.0d0**(-24))                                   !cgs
  real(8),parameter :: boltzmann_const = 1.3807d0*(10.0d0**(-16))                               !cgs
  real(8),parameter :: Navo = 6.0221d0*(10.0d0**23)                                             !Avogadro's constant
!  real(8),parameter :: mmw = 0.62d0									!mean molecular weight
  real(8),parameter :: mmw = 1.0d0                                                              !mean molecular weight
  real(8),parameter :: rg = 3.0d0*(10.0d0**6)*(mass_bh/10.0d0)                                  !cgs Schwarzschild radius
!  real(8),parameter :: rg_nrmlx = 1.0d0/50.0d0					         	!per nrmlx
  real(8),parameter :: rg_nrmlx = 1.0d0/10.0d0                                                  !per nrmlx
!  real(8),parameter :: rg_nrmlx = 1.0d0/1000.0d0                                                  !per nrmlx


  real(8),parameter :: nrmlro = 0.29d0*8.3d0*(10.0d0**(-7))/(mass_bh/10.0d0)            !normalized ro
  real(8),parameter :: nrmlte = mass_proton*(velocity_c**2)/boltzmann_const             !normalized te
  real(8),parameter :: nrmlx = rg/rg_nrmlx                                              !normalized x
  real(8),parameter :: nrmlv = sqrt(0.5d0*rg_nrmlx/((1.0d0-rg_nrmlx)**2))*velocity_c    !normalized v
  real(8),parameter :: nrmlt = nrmlx/nrmlv                                              !normalized t
  real(8),parameter :: nrmlee = nrmlro*nrmlv*nrmlv                                      !normalized ee

  real(8),parameter :: RadCool = nrmlro*sqrt(nrmlte)*nrmlx/(nrmlv**3)*(6.2d0*(10.0d0**20))   !cooling factor
  real(8),parameter :: te_factor = mmw*(nrmlv**2)/Navo/mass_proton/(velocity_c**2)           !factor for te

  real(8),parameter :: rohalo = 1.0d-4
!  real(8),parameter :: rohalo = 1.0d-6

! test1028


!----------------------------------------------------------------------|
!>>> for MPI >>>!
  include 'mpif.h'

  integer :: mpi_comm_cart3d, mpi_err
  integer :: npe, topo_coords(n_mpi_dims), n_tiles(n_mpi_dims)
  type(cans_mpi_rank) :: rank
!<<< for MPI <<<!

  
!======================================================================|
!     prologue
!======================================================================|

  integer :: mcont,ndi

!----------------------------------------------------------------------|
!  initialize counters

  !integer :: nd,ns,merr,ns1,ns2,nscount
  integer :: nd,ns,ns1,ns2,nscount
  real(8) :: time,timep,swtch_t

!----------------------------------------------------------------------|
!   time control parameters
!     nstop : number of total time steps for the run

  real(8) :: tend,dtout,dtout1
  integer :: nstop

  integer :: mwflag,mw,nt1,nt2,mw1
  integer :: nt3,nt4,nd1

  real(8) :: dt,hdt ! time step
  real(8) :: dtg

  real(8) :: safety ! CFL number
  real(8) :: dtmin ! minimum time step

  integer :: idf

  integer :: i,j,k,n

  integer :: mf_qq, err_flg,nanflg

  real(8) :: xin
  real(8) :: te_limit
 
  real(8) :: temp
!======================================================================|
!     prologue
!======================================================================|
  floor = 1.d-6

  mcont=0
  ndi=1000

  nd = 1
  nd1 = 1



!----------------------------------------------------------------------|
!>>> for MPI >>>

  !-- set array indicating MPI decomposition in Cartesian topology --!
  n_tiles(1) = 3
  n_tiles(2) = 3
  n_tiles(3) = 3

  !-- initialize the MPI execution environment --!
  call mpi_init(mpi_err)

  !-- get the size of the group associated with a communictor --!
  call mpi_comm_size(mpi_comm_world, npe, mpi_err)

  !-- check the consistency of npe and n_tiles --!
  if (npe /= product(n_tiles)) then
     call mpi_finalize(mpi_err)
     write(*,*) 'npe (',npe,') does not match n_tiles (',n_tiles,')'
     stop 
  end if

  !-- create a new communicator 
  !                 to which topology information has been attached --!
  call mpi_cart_create &
       &  (mpi_comm_world, n_mpi_dims, n_tiles, periods, reorder &
       & , mpi_comm_cart3d, mpi_err)

  !-- get rank of the calling process in the communicator [nnn] --!
  call mpi_comm_rank(mpi_comm_cart3d, rank%nnn, mpi_err)

  !-- get coords of calling process in cartesian topology --!
  call mpi_cart_coords &
       & ( mpi_comm_cart3d, rank%nnn, n_mpi_dims, topo_coords, mpi_err)

  !-- get ranks of neighbor processes (amount 1)
  !   in x-direction (direction 0: 1st coordinate) [dnn, unn]
  !   in y-direction (direction 1: 2nd coordinate) [ndn, nun]
  !   in z-direction (direction 2: 3rd coordinate) [nnd, nnu] --!
  call mpi_cart_shift &
       & ( mpi_comm_cart3d, 0, 1, rank%dnn, rank%unn, mpi_err )
  call mpi_cart_shift &
       & ( mpi_comm_cart3d, 1, 1, rank%ndn, rank%nun, mpi_err )
  call mpi_cart_shift &
       & ( mpi_comm_cart3d, 2, 1, rank%nnd, rank%nnu, mpi_err )


  !-- get ranks of diagonal processes
  !   in xy-plane [uun, udn, dun, ddn], 
  !   in xz-plane [unu, und, dnu, dnd],
  !   in yz-plane [nuu, nud, ndu, ndd],
  !   and at vertex [uuu, uud, udu, udd, duu, dud, ddu, ddd] --!
  !-- Note: There ranks are necessary for CT scheme but not for GLM-MHD --!
  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles,  1,  1,  0, topo_coords &
       & , rank%uun, mpi_err)
  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles,  1, -1,  0, topo_coords &
       & , rank%udn, mpi_err)
  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles, -1,  1,  0, topo_coords &
       & , rank%dun, mpi_err)
  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles, -1, -1,  0, topo_coords &
       & , rank%ddn, mpi_err)

  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles,  1,  0,  1, topo_coords &
       & , rank%unu, mpi_err)
  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles,  1,  0, -1, topo_coords &
       & , rank%und, mpi_err)
  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles, -1,  0,  1, topo_coords &
       & , rank%dnu, mpi_err)
  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles, -1,  0, -1, topo_coords &
       & , rank%dnd, mpi_err)

  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles,  0,  1,  1, topo_coords &
       & , rank%nuu, mpi_err)
  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles,  0,  1, -1, topo_coords &
       & , rank%nud, mpi_err)
  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles,  0, -1,  1, topo_coords &
       & , rank%ndu, mpi_err)
  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles,  0, -1, -1, topo_coords &
       & , rank%ndd, mpi_err)

  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles,  1,  1,  1, topo_coords &
       & , rank%uuu, mpi_err)
  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles,  1,  1, -1, topo_coords &
       & , rank%uud, mpi_err)
  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles,  1, -1,  1, topo_coords &
       & , rank%udu, mpi_err)
  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles,  1, -1, -1, topo_coords &
       & , rank%udd, mpi_err)
  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles, -1,  1,  1, topo_coords &
       & , rank%duu, mpi_err)
  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles, -1,  1, -1, topo_coords &
       & , rank%dud, mpi_err)
  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles, -1, -1,  1, topo_coords &
       & , rank%ddu, mpi_err)
  call cans_mpi_cart_rank &
       & ( mpi_comm_cart3d, n_mpi_dims, n_tiles, -1, -1, -1, topo_coords &
       & , rank%ddd, mpi_err)

  igx = (ix - 2 * margin) * n_tiles(1) + 2 * margin
  jgx = (kx - 2 * margin) * n_tiles(2) + 2 * margin
  kgx = (kx - 2 * margin) * n_tiles(3) + 2 * margin


  !-- display MPI info --!
  call mpi_barrier (mpi_comm_cart3d, mpi_err)
  
  if (rank%nnn == 0) then
     write(*,'(3(a5,i5),a8,i5)') &
          & 'igx=', igx, 'jgx=', jgx, 'kgx=', kgx, 'margin=', margin
     write(*,'(a10,i3,a11,3i3,a13,3l3)') &
          &   'MPI dims:', n_mpi_dims &
          & , 'MPI tiles:', n_tiles &
          & , 'MPI periods:', periods
  endif

  write(*,'(a30,3i3,a12,i3,/,a44,/,3i3,a3,3i3,a3,3i3,/,3i3,a3,3i3,a3,3i3,/,3i3,a3,3i3,a3,3i3,/)') &
       & 'My coords in MPI topology: [', topo_coords, ' ], My rank:', rank%nnn &
       & ,'Ranks of neighbor and diagonal processes: ' &
       & , rank%ddu, rank%ndu, rank%udu, '|', rank%dnu, rank%nnu, rank%unu, '|', rank%duu, rank%nuu, rank%uuu &
       & , rank%ddn, rank%ndn, rank%udn, '|', rank%dnn, rank%nnn, rank%unn, '|', rank%dun, rank%nun, rank%uun &
       & , rank%ddd, rank%ndd, rank%udd, '|', rank%dnd, rank%nnd, rank%und, '|', rank%dud, rank%nud, rank%uud

  call mpi_barrier (mpi_comm_cart3d, mpi_err)
  
  !<<< for MPI <<<

!----------------------------------------------------------------------|
!  initialize counters

  nd=1
  time  = 0.0d0
  timep = 0.0d0
  swtch_t=62.8d0
  ns    = 0
  mpi_err  = 0
  nscount=0  
      
!----------------------------------------------------------------------|
!   time control parameters
!     nstop : number of total time steps for the run
  
  tend = 6.28d0*10.0d0

  dtout = tend/(100.0d0)
  dtout1 = dtout/10.0d0

  tend = 6.28d0*0.20d0

  dt = tend

  !nstop = 2000000
  nstop = 10

!----------------------------------------------------------------------|
!   setup numerical model (grid, initial conditions, etc.)
!

  call model_mpi_cart3d &
       & (igx, jgx, kgx, ix, jx, kx, margin &
       & , ro, pr, vx, vy, vz, bx, by, bz, phi &
       & , gm &
       & , x, dx, y, dy, z, dz &
       & , gx, gz, cr, eta0, vc, eta, xin &
       & , ccx, ccy, ccz, rg_nrmlx, rohalo, te_limit, te_factor, nrmlv &
       & , nrmlte, boltzmann_const, Navo, mmw &
       & , rank, topo_coords)

  call cans_mpi_exchange_cart3d &
       & (ro, pr, vx, vy, vz, bx, by, bz, phi &
       & , ix, jx, kx, margin, mpi_comm_cart3d, npe, rank, mpi_err)


  call bnd_mpi_cart3d &
       & ( ro, pr, vx, vy, vz, bx, by, bz, phi, eta &
       & , margin, ix, jx, kx, rank)


  call saveInit_3d(ix,jx,kx,ro,roi)
  call saveInit_3d(ix,jx,kx,pr,pri)

  call saveInit_3d(ix,jx,kx,vx,vxi)
  call saveInit_3d(ix,jx,kx,vy,vyi)
  call saveInit_3d(ix,jx,kx,vz,vzi)

  call saveInit_3d(ix,jx,kx,bx,bxi)
  call saveInit_3d(ix,jx,kx,by,byi)

!----------------------------------------------------------------------|
!  file open
  call openfileCor(nd, rank%nnn, ix, jx, kx)
  call openfileAll(nd, rank%nnn, ix, jx, kx)
 
  call dacputparamc(mf_params,'comment','model_mpi_cart3d')
  call dacputparamd(mf_params,'dtout',dtout)
  call dacputparamd(mf_params,'tend',tend)
  call dacputparami(mf_params,'ix',ix)
  call dacputparami(mf_params,'jx',jx)
  call dacputparami(mf_params,'kx',kx)
  call dacputparami(mf_params,'igx',igx)
  call dacputparami(mf_params,'jgx',jgx)
  call dacputparami(mf_params,'kgx',kgx)
  call dacputparami(mf_params,'margin',margin)
  call dacputparami(mf_params,'mpi',1)
  call dacputparami(mf_params,'mpisize',npe)
  call dacputparami(mf_params,'mpirank',rank%nnn)
  call dacputparami(mf_params,'mpitopox',topo_coords(1))
  call dacputparami(mf_params,'mpitopoy',topo_coords(2))
  call dacputparami(mf_params,'mpitopoz',topo_coords(3))
  call dacputparami(mf_params,'mpix',n_tiles(1))
  call dacputparami(mf_params,'mpiy',n_tiles(2))
  call dacputparami(mf_params,'mpiz',n_tiles(3))
  call dacputparami(mf_params,'beta',100)
  call dacputparamd(mf_params,'nrmlro',nrmlro)
  call dacputparamd(mf_params,'nrmlte',nrmlte)
  call dacputparamd(mf_params,'nrmlx',nrmlx)
  call dacputparamd(mf_params,'nrmlv',nrmlv)
  call dacputparamd(mf_params,'nrmlt',nrmlt)
  call dacputparamd(mf_params,'nrmlee',nrmlee)
  call dacputparamd(mf_params,'mass_bh',mass_bh)
  call dacputparamd(mf_params,'rg',rg)
  call dacputparamd(mf_params,'rg_nrmlx',rg_nrmlx)
  call dacputparamd(mf_params,'RadCool',RadCool)
  call dacputparamd(mf_params,'te_factor',te_factor)
  call dacputparamd(mf_params,'rohalo',rohalo)
  call dacputparamd(mf_params,'eta0',eta0)
  call dacputparamd(mf_params,'vc',vc)
  call dacputparamd(mf_params,'x(1)',x(1))
  call dacputparamd(mf_params,'y(1)',y(1))
  call dacputparamd(mf_params,'z(1)',z(1))
  call dacputparamd(mf_params,'dx(1)',dx(1))
  call dacputparamd(mf_params,'dy(1)',dy(1))
  call dacputparamd(mf_params,'dz(1)',dz(1))
!  call dacputparamd(mf_params,'te_limit',te_limit)


!----------------------------------------------------------------------|
!     ready
!----------------------------------------------------------------------|
  
  write(mf_x) x
  write(mf_y) y
  write(mf_z) z

  call dacputparamd(mf_params,'gm',gm)

  if(rank%nnn == 0)then
     write(mf_t) time
     write(6,'(a10,a5,i8,a6,e10.3,a5,i3)') &
          & ' write    ','step=',ns,' time=',time,' nd =',nd
     write(*,*) 'dt :: ',dt
  endif

  write(mf_ro) ro
  write(mf_pr) pr
  write(mf_vx) vx
  write(mf_vy) vy
  write(mf_vz) vz
  write(mf_bx) bx
  write(mf_by) by
  write(mf_bz) bz

  write(mf_phi) phi
  write(mf_eta) eta

  write(mf_gx) gx
  write(mf_gz) gz

  call closefileAll()
  call closefileCor()
  close(mf_params)

!----------------------------------------------------------------------|
!  read-data


  ix0 = ix
  jx0 = jx
  kx0 = kx
  nx0 = nd

!!$  open(100,file='readFileNumber.dat')
!!$  read(100,*) nd
!!$  close(100)
!!$  if(nd .ne. 1)then
!!$     call openReadFileAll(nd,mpirank,ix0,jx0,kx0,nx0)
!!$     read(mfi_ro) ro
!!$     read(mfi_pr) pr
!!$
!!$     read(mfi_vx) vx
!!$     read(mfi_vy) vy
!!$     read(mfi_vz) vz
!!$     
!!$     read(mfi_bx) bx
!!$     read(mfi_by) by
!!$     read(mfi_bz) bz
!!$
!!$     read(mfi_phi) phi
!!$     read(mfi_eta) eta
!!$     time = real(nd-1)*dtout
!!$     if(mpirank .eq. 0)then
!!$        write(*,*) 'time :: ',time
!!$     endif
!!$     
!!$     call closeReadFileAll()
!!$ 	 
!!$  endif
  nd=nd+1



! test1028
!if (FLAG_ModAverage .eq. 1) then
!  call openfileModTCR(mpirank,margin,ix,kx,x,xObs1,xObs2,xObs3,xObs4,xObs5)
!  call closefileModTCR()
!  call open_CalPhysValue(mpirank,'averageVal',mf_calPhysValue)
!end if
! test1028

!======================================================================|
!     time integration 
!======================================================================|
1000  continue  
  ns = ns+1
  nscount = nscount+1
  mwflag=0

!----------------------------------------------------------------------|
!     obtain time spacing

  safety=0.2d0
  dtmin=1.d-10
  
!  call getNewdt_glmcyl_res_mpi(margin,safety,dtmin,ix,jx,kx,gm,ro,pr &
!       ,vx,vy,vz,bx,by,bz,x,dx,y,dy,z,dz,eta &
!       ,dt,merr,ch,mpirank)

  temp=0.0d0
!!$  call getNewdt_glmcyl_res_mpi(margin,safety,dtmin,ix,jx,kx,gm,ro,pr &
!!$       ,vx,vy,vz,bx,by,bz,x,dx,y,dy,z,dz,eta &
!!$       ,dt,merr,ch,mpirank,temp)
!!$
!!$  call mpi_allreduce(dt,dtg,1,mpi_double_precision,mpi_min &
!!$       ,mpi_comm_world,merr)

  call getNewdt_glmcyl_res_mpi &
       & ( margin, safety, dtmin, ix, jx, kx, gm &
       & , ro, pr, vx, vy, vz, bx, by, bz &
       & , x, dx, y, dy, z, dz, eta &
       & , dt, mpi_err, ch, rank%nnn, temp)

  call mpi_allreduce &
       & ( dt, dtg, 1, mpi_double_precision, mpi_min, mpi_comm_cart3d, mpi_err)

  dt = dtg
  err_flg = 1 ! error flag
  if ( dt.lt.dtmin) goto 9999
  err_flg = 2 ! error flag
  if (mpi_err /= 0) goto 9999
  timep = time
  time = time+dt

  mw=0
  mw1=0
!----- check output

! dtout
  nt1=int(timep/dtout)
  nt2=int(time/dtout)
  if (nt1.lt.nt2) mw=1

! hst
  nt3=int(timep/dtout1)
  nt4=int(time/dtout1)
  if (nt3.lt.nt4) mw1=1

!---- integrate

!!$!---1st
!!$  call integrate_cyl_mp5_1(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,dt &
!!$       ,gx,gz,floor,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
!!$       ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
!!$       ,eta0,vc,eta,ccx,ccy,ccz)
!!$
!!$  call exchangeMpixz2(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1 &
!!$       ,phi1,merr)
!!$

!!$  call bnd(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1,eta,x)
!!$
!!$  call bnd_absoub(ix,jx,kx,x,z,xin,roi,pri,vxi,vyi,vzi,bxi,byi,bzi &
!!$       ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1)
!!$
!!$!---2nd
!!$  call integrate_cyl_mp5_2(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,dt &
!!$       ,gx,gz,floor,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
!!$       ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
!!$       ,eta0,vc,eta,ccx,ccy,ccz)
!!$
!!$  call exchangeMpixz2(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1 &
!!$       ,phi1,merr)
!!$
!!$  call bnd(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1,eta,x)
!!$
!!$  call bnd_absoub(ix,jx,kx,x,z,xin,roi,pri,vxi,vyi,vzi,bxi,byi,bzi &
!!$       ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1)
!!$
!!$!----3rd
!!$  call integrate_cyl_mp5_3(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,dt &
!!$       ,gx,gz,floor,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
!!$       ,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
!!$       ,eta0,vc,eta,ccx,ccy,ccz)

!------2stage version

  hdt = 0.5d0*dt

  call integrate_cyl_mp5_1st(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,hdt &
       ,gx,gz,floor,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
       ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
       ,eta0,vc,ccx,ccy,ccz,RadCool,te_factor,time,rohalo,swtch_t,xin)
  
  call cans_mpi_exchange_cart3d &
       & (ro1, pr1, vx1, vy1, vz1, bx1, by1, bz1, phi1 &
       & , ix, jx, kx, margin, mpi_comm_cart3d, npe, rank, mpi_err)

  call bnd_mpi_cart3d &
       & ( ro1, pr1, vx1, vy1, vz1, bx1, by1, bz1, phi1, eta &
       & , margin, ix, jx, kx, rank)

  call bnd_absorb &
       & (ix, jx, kx, x, z, xin &
       & , roi, pri, vxi, vyi, vzi, bxi, byi, bzi &
       & , ro1, pr1, vx1, vy1, vz1, bx1, by1, bz1)

  call integrate_cyl_mp5_2nd(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,dt &
       ,gx,gz,floor,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
       ,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
       ,eta0,vc,eta,ccx,ccy,ccz,RadCool,te_factor,time,rohalo,swtch_t,xin)

  call cans_mpi_exchange_cart3d &
       & (ro, pr, vx, vy, vz, bx, by, bz, phi &
       & , ix, jx, kx, margin, mpi_comm_cart3d, npe, rank, mpi_err)

  call bnd_mpi_cart3d &
       & ( ro, pr, vx, vy, vz, bx, by, bz, phi, eta &
       & , margin, ix, jx, kx, rank)

  call bnd_absorb &
       & (ix, jx, kx, x, z, xin &
       & , roi, pri, vxi, vyi, vzi, bxi, byi, bzi &
       & , ro,  pr,  vx,  vy,  vz,  bx,  by,  bz)

!------2stage version
!----------------------------------------------------------------------|
!     data output 

  
  if(rank%nnn == 0)then
!   if(nscount .ge. 100)then
     print '(" step=",i8," time=",e10.5," dt=",e10.5)', ns, time, dt
!     nscount=0
!   endif
  end if

!-----test1028----------------------
!   ModAverage test

!if (FLAG_ModAverage .eq. 1) then
!  if (mw1.ne.0) then  
!     do k=1,kx
!        do j=1,jx
!           do i=1,ix
!              pb(i,j,k) = 0.5d0*(bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2)
!              rx(i,j,k) = ro(i,j,k)*vx(i,j,k)
!              ry(i,j,k) = ro(i,j,k)*vy(i,j,k)
!              rz(i,j,k) = ro(i,j,k)*vz(i,j,k)
!              am(i,j,k) = vy(i,j,k)*x(i)
!              ee(i,j,k) = pb(i,j,k) &
!                   +0.5d0*ro(i,j,k)*(vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2) &
!                   +pr(i,j,k)/(gm-1.0d0)
!           enddo
!        enddo
!     enddo
!     
!     call calSumBase(margin,ix,jx,kx,x,y,z,xmin_MAV,xmax_MAV,zmin_MAV,zmax_MAV &
!          ,baseAV)
!     
!     call calSum(margin,ix,jx,kx,x,y,z,xmin_MAV,xmax_MAV,zmin_MAV,zmax_MAV &
!          ,pb,pbAV)
!     call calSum(margin,ix,jx,kx,x,y,z,xmin_MAV,xmax_MAV,zmin_MAV,zmax_MAV &
!          ,pr,prAV)
!     call calSum(margin,ix,jx,kx,x,y,z,xmin_MAV,xmax_MAV,zmin_MAV,zmax_MAV &
!          ,ro,roAV)
!     call calSum(margin,ix,jx,kx,x,y,z,xmin_MAV,xmax_MAV,zmin_MAV,zmax_MAV &
!          ,ee,eeAV)
!     call calSum(margin,ix,jx,kx,x,y,z,xmin_MAV,xmax_MAV,zmin_MAV,zmax_MAV &
!          ,am,amAV)
!
!     avInterVal(1) = baseAv
!     avInterVal(2) = pbAv
!     avInterVal(3) = prAv
!     avInterVal(4) = roAv
!     avInterVal(5) = eeAv
!     avInterVal(6) = amAv
!
!     call exchangeSumVal(mpirank,iInterNum,avInterVal)
!     
!     if(mpirank == 0)then
!        pbAv = avInterVal(2)/avInterVal(1)
!        prAv = avInterVal(3)/avInterVal(1)
!        roAv = avInterVal(4)/avInterVal(1)
!        eeAv = avInterVal(5)/avInterVal(1)
!        amAv = avInterVal(6)/avInterVal(1)
!
!        avVal(1) = pbAv
!        avVal(2) = prAv
!        avVal(3) = roAv
!        avVal(4) = eeAv
!        avVal(5) = amAv
!
!        call outputCalPhysValue(mpirank,iCalNum,avVal,ns,time &
!             ,'averageVal',mf_calPhysValue)
!     endif
!!   TimeChnageav@Rad test
!     call openfileModTCRsOld(mpirank,margin,ix,x,xObs1,xObs2,xObs3,xObs4,xObs5)
!
!     call TimeChangeVal(margin,ix,jx,kx,x,y,z,xObs1 &
!          ,ro,pb,pr,rx,ry,rz,bx,by,bz &
!          ,baseTCR,roTCR,pbTCR,prTCR,rxTCR,ryTCR,rzTCR,bxTCR,byTCR,bzTCR)
!     
!     if(isInRadius(margin,ix,x,xObs1) == 1)then
!        write(mtr1_bx) bxTCR
!        write(mtr1_by) byTCR
!        write(mtr1_bz) bzTCR
!        write(mtr1_rx) rxTCR
!        write(mtr1_ry) ryTCR
!        write(mtr1_rz) rzTCR
!        write(mtr1_ro) roTCR
!        write(mtr1_pb) pbTCR
!        write(mtr1_pr) prTCR
!     end if
!
!     call TimeChangeVal(margin,ix,jx,kx,x,y,z,xObs2 &
!          ,ro,pb,pr,rx,ry,rz,bx,by,bz &
!          ,baseTCR,roTCR,pbTCR,prTCR,rxTCR,ryTCR,rzTCR,bxTCR,byTCR,bzTCR)
!     
!     if(isInRadius(margin,ix,x,xObs2) == 1)then
!        write(mtr2_bx) bxTCR
!        write(mtr2_by) byTCR
!        write(mtr2_bz) bzTCR
!        write(mtr2_rx) rxTCR
!        write(mtr2_ry) ryTCR
!        write(mtr2_rz) rzTCR
!        write(mtr2_ro) roTCR
!        write(mtr2_pb) pbTCR
!        write(mtr2_pr) prTCR
!     end if
!
!     call TimeChangeVal(margin,ix,jx,kx,x,y,z,xObs3 &
!          ,ro,pb,pr,rx,ry,rz,bx,by,bz &
!          ,baseTCR,roTCR,pbTCR,prTCR,rxTCR,ryTCR,rzTCR,bxTCR,byTCR,bzTCR)
!     
!     if(isInRadius(margin,ix,x,xObs3) == 1)then
!        write(mtr3_bx) bxTCR
!        write(mtr3_by) byTCR
!        write(mtr3_bz) bzTCR
!        write(mtr3_rx) rxTCR
!        write(mtr3_ry) ryTCR
!        write(mtr3_rz) rzTCR
!        write(mtr3_ro) roTCR
!        write(mtr3_pb) pbTCR
!        write(mtr3_pr) prTCR
!     end if
!
!     call TimeChangeVal(margin,ix,jx,kx,x,y,z,xObs4 &
!          ,ro,pb,pr,rx,ry,rz,bx,by,bz &
!          ,baseTCR,roTCR,pbTCR,prTCR,rxTCR,ryTCR,rzTCR,bxTCR,byTCR,bzTCR)
!     
!     if(isInRadius(margin,ix,x,xObs4) == 1)then
!        write(mtr4_bx) bxTCR
!        write(mtr4_by) byTCR
!        write(mtr4_bz) bzTCR
!        write(mtr4_rx) rxTCR
!        write(mtr4_ry) ryTCR
!        write(mtr4_rz) rzTCR
!        write(mtr4_ro) roTCR
!        write(mtr4_pb) pbTCR
!        write(mtr4_pr) prTCR
!     end if
!
!     call TimeChangeVal(margin,ix,jx,kx,x,y,z,xObs5 &
!          ,ro,pb,pr,rx,ry,rz,bx,by,bz &
!          ,baseTCR,roTCR,pbTCR,prTCR,rxTCR,ryTCR,rzTCR,bxTCR,byTCR,bzTCR)
!     
!     if(isInRadius(margin,ix,x,xObs5) == 1)then
!        write(mtr5_bx) bxTCR
!        write(mtr5_by) byTCR
!        write(mtr5_bz) bzTCR
!        write(mtr5_rx) rxTCR
!        write(mtr5_ry) ryTCR
!        write(mtr5_rz) rzTCR
!        write(mtr5_ro) roTCR
!        write(mtr5_pb) pbTCR
!        write(mtr5_pr) prTCR
!     end if
!
!     call closefileModTCR()
!
!     if(mpirank .eq. 0)then
!        write(6,913) ns,time,nd1
!        write(*,*) '[TCR][AV] dt :: ',dt
!        write(*,*) 'nd1: ',nd1
!     end if
!     nd1 = nd1+1
!  endif
!
!endif
!-----test1028----------------------

  if (mw.ne.0) then
     call openfileAll(nd, rank%nnn, ix, jx, kx)

     write(mf_t) time
     write(mf_ro) ro
     write(mf_pr) pr
     write(mf_vx) vx
     write(mf_vy) vy
     write(mf_vz) vz
     write(mf_bx) bx
     write(mf_by) by
     write(mf_bz) bz

     write(mf_phi) phi
     write(mf_eta) eta

     call closefileAll()

     if(rank%nnn == 0)then
        write(6,'(a10,a5,i8,a6,e10.3,a5,i3)') &
             & ' write    ','step=',ns,' time=',time,' nd =',nd
        write(*,*) '[NORMAL] dt :: ',dt
        write(*,*) 'nd: ',nd
     end if

     open(100,file='readFileNumber.dat')
     write(100,*) nd
     close(100)

     nd=nd+1
     mwflag=1
  endif

!----------------------------------------------------------------------|
!     loop test

  if (ns .lt. nstop .and. time .lt. tend ) goto 1000

!======================================================================|
!     epilogue
!======================================================================|
9999 continue
  if (mpi_err /= 0) then
     write(6,*) '#### abnormal stop ####'
  endif

!----------------------------------------------------------------------|
!  data output
  if (mwflag.eq.0) then
     if(rank%nnn == 0)then
        write(6,'(a10,a5,i8,a6,e10.3,a5,i3)') &
             & ' write    ','step=',ns,' time=',time,' nd =',nd
     end if

     call openfileAll(nd, rank%nnn, ix, jx, kx)

     write(mf_t) time
     write(mf_ro) ro
     write(mf_pr) pr
     write(mf_vx) vx
     write(mf_vy) vy
     write(mf_vz) vz
     write(mf_bx) bx
     write(mf_by) by
     write(mf_bz) bz

     write(mf_phi) phi
     write(mf_eta) eta

     call closefileAll()
  endif

! test1028
!  call closefileModTCR()
! test1028
!======================================================================|
!        enddo
!======================================================================|

  close(mf_t)

  if (rank%nnn == 0) then
     write(6,'(a10,a5,i8,a6,e10.3)') ' stop     ','step=',ns,' time=',time
     if (mpi_err == 0) then
        write(6,*) '  ### normal stop ###'
     else
        write(6,*) '  ### abnormal stop ###'
     endif
  endif

  call mpi_finalize(mpi_err)  
  stop

end program main
