module init
  use const
  use mpi_domain_xz
  implicit none
  private

  public :: mpi_hoge

!  real(8),dimension(ix,jx,kx),public :: ro,pr
!  real(8),dimension(ix,jx,kx),public :: vx,vy,vz
!  real(8),dimension(ix,jx,kx),public :: bx,by,bz,phi
!  real(8),dimension(ix,jx,kx),public :: eta

!  real(8),dimension(ix,jx,kx) :: ro1,pr1
!  real(8),dimension(ix,jx,kx) :: vx1,vy1,vz1
!  real(8),dimension(ix,jx,kx) :: bx1,by1,bz1,phi1

!! gravitation
!  real(8),dimension(ix,jx,kx) :: gx,gy,gz

!  real(8),dimension(ix,jx,kx) :: roi,pri
!  real(8),dimension(ix,jx,kx) :: vxi,vyi,vzi
!  real(8),dimension(ix,jx,kx) :: bxi,byi,bzi

! MPI
  integer,public :: mpisize
  integer,public :: mpirank
! MPI type
  type(mpidomain) :: mpid



contains

!subroutine init~~
  
!  call model()
!  call bnd(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,x,z &
!               ,xin,roi,pri,vxi,vyi,vzi,bxi,byi,bzi)

!  call exchangeMpixz(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz &
!         ,phi,merr)

!  end subroutine init~~


subroutine mpi_hoge(mpid,merr,igx,jgx,kgx)
  use const
  use mpi_domain_xz
  implicit none
  include 'mpif.h'

!  integer :: margin,
  integer :: merr
!  integer :: ix,jx,kx 
  integer :: igx,jgx,kgx 
!  real(8) :: floor 
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
  
  igx = ix*mpid%mpisize_2d(1)-2*margin*(mpid%mpisize_2d(1)-1)
  jgx = jx
  kgx = kx*mpid%mpisize_2d(2)-2*margin*(mpid%mpisize_2d(2)-1)

  ! determine mpid%mpirank_3d
  call setmy2drank(mpid,merr)
  call setmpiboundary(mpid)
  if (mpirank.eq.0) then
  write(6,*) 'igx=',igx,'kgx=',kgx
  write(*,*) 'mpirank',mpirank,'mpirank_2d',mpid%mpirank_2d
  endif

end subroutine mpi_hoge








  end module init

