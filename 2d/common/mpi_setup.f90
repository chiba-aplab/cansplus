module mpi_setup

  implicit none
  private

  include 'mpif.h'

  public :: mpi_setup__init

  type :: mpivars
     integer :: mpisize, mpirank
     integer, dimension(2) :: mpirank_2d
     integer :: l,r,f,b
  end type mpivars

  type(mpivars), public :: mpid
  integer, public            :: merr, mcomw = MPI_COMM_WORLD, mstat(MPI_STATUS_SIZE)
  integer, public, parameter :: mnull = MPI_PROC_NULL
  integer, public, parameter :: mdp   = MPI_DOUBLE_PRECISION
  integer, public, parameter :: mmin  = MPI_MIN
  logical, save :: lcheck=.true.
  integer, save , allocatable :: ptable(:,:)
 

contains


  subroutine mpi_setup__init(mpisize_x,mpisize_y,pbcheck)

    logical, intent(in) :: pbcheck(2)
    integer, intent(in) :: mpisize_x, mpisize_y
    integer :: mpisize, mpirank
    integer :: i, j, k, irank

    allocate(ptable(-1:mpisize_x,-1:mpisize_y))

    call mpi_init(merr)
    call mpi_comm_size(mcomw,mpisize,merr)
    call mpi_comm_rank(mcomw,mpirank,merr)

    if(mpisize /= (mpisize_x*mpisize_y)) then
       stop 'error in proc no.'
       call MPI_ABORT(mcomw, 9, merr)
       call MPI_FINALIZE(merr)
    endif

    mpid%mpirank = mpirank
    mpid%mpisize = mpisize
  
    !rank table
    do j=-1,mpisize_y
    do i=-1,mpisize_x
       ptable(i,j) = MPI_PROC_NULL
    enddo
    enddo
    irank = 0
    do j=0,mpisize_y-1
    do i=0,mpisize_x-1
       ptable(i,j) = irank
       if(mpirank == irank)then
          mpid%mpirank_2d(1) = i
          mpid%mpirank_2d(2) = j
       endif
       irank = irank+1
    enddo
    enddo

    !Periodic boundary condition if pbcheck=true
    if(pbcheck(1))then
       ptable(-1,-1:mpisize_y) = ptable(mpisize_x-1,-1:mpisize_y)
       ptable(mpisize_x,-1:mpisize_y) = ptable(0,-1:mpisize_y)
    endif
    if(pbcheck(2))then
       ptable(-1:mpisize_x,-1) = ptable(-1:mpisize_x,mpisize_y-1)
       ptable(-1:mpisize_x,mpisize_y) = ptable(-1:mpisize_x,0)
    endif

    !For MPI_SENDRECV
    mpid%r = ptable(mpid%mpirank_2d(1)+1,mpid%mpirank_2d(2)  )
    mpid%l = ptable(mpid%mpirank_2d(1)-1,mpid%mpirank_2d(2)  )
    mpid%f = ptable(mpid%mpirank_2d(1),  mpid%mpirank_2d(2)+1)
    mpid%b = ptable(mpid%mpirank_2d(1),  mpid%mpirank_2d(2)-1)

    lcheck = .false.

  end subroutine mpi_setup__init


end module mpi_setup
