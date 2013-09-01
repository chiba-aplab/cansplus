module init
  use mpi_domain_xz
  use const
  implicit none
!  private

  public :: initial

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
! initial variables
  real(8),dimension(ix,jx,kx) :: roi,pri
  real(8),dimension(ix,jx,kx) :: vxi,vyi,vzi
  real(8),dimension(ix,jx,kx) :: bxi,byi,bzi

  real(8) :: ch

!----------------------------------------------------------------------|
!  initialize counters

  integer :: nd,ns,merr,ns1,ns2,nscount
  real(8) :: time,timep

!----------------------------------------------------------------------|
!   time control parameters
!     nstop : number of total time steps for the run

  integer :: mwflag,mw,nt1,nt2

  real(8) :: dt,hdt
  real(8) :: dtg

  integer :: i,j,k,n

  integer :: mf_qq, err_flg

!-----------------------------------------------------------  
! MPI
  integer,public :: mpisize
  integer,public :: mpirank
! MPI type
  type(mpidomain) :: mpid



contains

subroutine initial()
  use model
  implicit none

!----------------------------------------------------------------------|
!   for MPI
!  
  call  mpi_hoge(mpid,merr)
!----------------------------------------------------------------------|
!   setup numerical model (grid, initial conditions, etc.)
!
 call model_machida(ro,pr,vx,vy,vz,bx,by,bz,phi &
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
               
end subroutine initial


subroutine mpi_hoge(mpid,merr)
  use const
  use mpi_domain_xz
  implicit none
  include 'mpif.h'

  integer :: merr
! MPI type
  type(mpidomain) :: mpid
!----------------------------------------------------------------------|
!   for MPI
!  
  call mpi_init(merr)
  call mpi_comm_size(mpi_comm_world,mpisize,merr)
  call mpi_comm_rank(mpi_comm_world,mpirank,merr)
  
  mpid%mpirank = mpirank
  mpid%mpisize = mpisize
  
  ! Core Numbers, 1: r(i)-direction, 2: z(k)-direction
  mpid%mpisize_2d(1) = mpisize_x !1
  mpid%mpisize_2d(2) = mpisize_z !16
  
!  igx = ix*mpid%mpisize_2d(1)-2*margin*(mpid%mpisize_2d(1)-1)
!  jgx = jx
!  kgx = kx*mpid%mpisize_2d(2)-2*margin*(mpid%mpisize_2d(2)-1)

  ! determine mpid%mpirank_3d
  call setmy2drank(mpid,merr)
  call setmpiboundary(mpid)
  if (mpirank.eq.0) then
  write(6,*) 'igx=',igx,'kgx=',kgx
  write(*,*) 'mpirank',mpirank,'mpirank_2d',mpid%mpirank_2d
  endif

end subroutine mpi_hoge

end module init

