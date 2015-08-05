module mpi_setup

  implicit none
  private

  include 'mpif.h'

  public :: mpi_setup__init, mpi_setup__init_cyl

  type :: mpivars
     integer :: mpisize, mpirank
     integer, dimension(3) :: mpirank_3d
     integer :: l,r,f,b,t,d,ax
  end type mpivars

  type(mpivars), public :: mpid
  integer, public            :: merr, mcomw = MPI_COMM_WORLD, mstat(MPI_STATUS_SIZE)
  integer, public, parameter :: mnull = MPI_PROC_NULL
  integer, public, parameter :: mdp   = MPI_DOUBLE_PRECISION
  integer, public, parameter :: mmin  = MPI_MIN
  logical, save :: lcheck=.true.
  integer, save , allocatable :: ptable(:,:,:)
 

contains


  subroutine mpi_setup__init(mpisize_x,mpisize_y,mpisize_z,pbcheck)

    logical, intent(in) :: pbcheck(3)
    integer, intent(in) :: mpisize_x, mpisize_y, mpisize_z
    integer :: mpisize, mpirank
    integer :: i, j, k, irank

    allocate(ptable(-1:mpisize_x,-1:mpisize_y,-1:mpisize_z))

    call mpi_init(merr)
    call mpi_comm_size(mcomw,mpisize,merr)
    call mpi_comm_rank(mcomw,mpirank,merr)

    if(mpisize /= (mpisize_x*mpisize_y*mpisize_z)) then
       stop 'error in proc no.'
       call MPI_ABORT(mcomw, 9, merr)
       call MPI_FINALIZE(merr)
    endif

    mpid%mpirank = mpirank
    mpid%mpisize = mpisize
  
    !rank table
    do k=-1,mpisize_z
    do j=-1,mpisize_y
    do i=-1,mpisize_x
       ptable(i,j,k) = MPI_PROC_NULL
    enddo
    enddo
    enddo
    irank = 0
    do k=0,mpisize_z-1
    do j=0,mpisize_y-1
    do i=0,mpisize_x-1
       ptable(i,j,k) = irank
       if(mpirank == irank)then
          mpid%mpirank_3d(1) = i
          mpid%mpirank_3d(2) = j
          mpid%mpirank_3d(3) = k
       endif
       irank = irank+1
    enddo
    enddo
    enddo

    !Periodic boundary condition if pbcheck=true
    if(pbcheck(1))then
       ptable(-1,-1:mpisize_y,-1:mpisize_z) = ptable(mpisize_x-1,-1:mpisize_y,-1:mpisize_z)
       ptable(mpisize_x,-1:mpisize_y,-1:mpisize_z) = ptable(0,-1:mpisize_y,-1:mpisize_z)
    endif
    if(pbcheck(2))then
       ptable(-1:mpisize_x,-1,-1:mpisize_z) = ptable(-1:mpisize_x,mpisize_y-1,-1:mpisize_z)
       ptable(-1:mpisize_x,mpisize_y,-1:mpisize_z) = ptable(-1:mpisize_x,0,-1:mpisize_z)
    endif
    if(pbcheck(3))then
       ptable(-1:mpisize_x,-1:mpisize_y,-1) = ptable(-1:mpisize_x,-1:mpisize_y,mpisize_z-1)
       ptable(-1:mpisize_x,-1:mpisize_y,mpisize_z) = ptable(-1:mpisize_x,-1:mpisize_y,0)
    endif

    !For MPI_SENDRECV
    mpid%r = ptable(mpid%mpirank_3d(1)+1,mpid%mpirank_3d(2),  mpid%mpirank_3d(3)  )
    mpid%l = ptable(mpid%mpirank_3d(1)-1,mpid%mpirank_3d(2),  mpid%mpirank_3d(3)  )
    mpid%f = ptable(mpid%mpirank_3d(1),  mpid%mpirank_3d(2)+1,mpid%mpirank_3d(3)  )
    mpid%b = ptable(mpid%mpirank_3d(1),  mpid%mpirank_3d(2)-1,mpid%mpirank_3d(3)  )
    mpid%t = ptable(mpid%mpirank_3d(1),  mpid%mpirank_3d(2),  mpid%mpirank_3d(3)+1)
    mpid%d = ptable(mpid%mpirank_3d(1),  mpid%mpirank_3d(2),  mpid%mpirank_3d(3)-1)
    mpid%ax = MPI_PROC_NULL

    lcheck = .false.

  end subroutine mpi_setup__init


  subroutine mpi_setup__init_cyl(mpisize_y)

    integer, intent(in) :: mpisize_y

    if(lcheck) stop 'Call mpi_setup__init first.'
    if(mod(mpisize_y,2) /= 0)then
       stop 'Proc. size in the theta direction (mpisize_y) should be an even number.'
       call MPI_ABORT(mcomw, 9, merr)
       call MPI_FINALIZE(merr)
    endif

    mpid%ax = ptable(mpid%mpirank_3d(1),mod(mpid%mpirank_3d(2)+mpisize_y/2,mpisize_y),mpid%mpirank_3d(3))

  end subroutine mpi_setup__init_cyl


end module mpi_setup
