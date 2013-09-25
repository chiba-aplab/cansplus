module mpi_setup

  implicit none
  private

  include 'mpif.h'

  public :: mpi_setup__init

  type :: mpivars
     integer :: mpisize
     integer :: mpirank
     integer, dimension(3) :: mpirank_3d
     integer :: l,r,f,b,t,d
  end type mpivars

  type(mpivars), public :: mpid
  integer, public            :: merr, mcomw = MPI_COMM_WORLD, mstat(MPI_STATUS_SIZE)
  integer, public, parameter :: mnull = MPI_PROC_NULL
  integer, public, parameter :: mdp   = MPI_DOUBLE_PRECISION
  integer, public, parameter :: mmin  = MPI_MIN
 

contains


  subroutine mpi_setup__init(mpisize_x,mpisize_y,mpisize_z)

    integer, intent(in) :: mpisize_x, mpisize_y, mpisize_z
    integer :: ptable(-1:mpisize_x,-1:mpisize_y,-1:mpisize_z)
    integer :: mpisize, mpirank
    integer :: i, j, k, irank

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
    do i=0,mpisize_x-1
    do j=0,mpisize_y-1
    do k=0,mpisize_z-1
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

    !Periodic boundary condition in y
    ptable(-1:mpisize_x,-1,-1:mpisize_z) = ptable(-1:mpisize_x,mpisize_y-1,-1:mpisize_z)
    ptable(-1:mpisize_x,mpisize_y,-1:mpisize_z) = ptable(-1:mpisize_x,0,-1:mpisize_z)

    !For MPI_SENDRECV
    mpid%r = ptable(mpid%mpirank_3d(1)+1,mpid%mpirank_3d(2),  mpid%mpirank_3d(3)  )
    mpid%l = ptable(mpid%mpirank_3d(1)-1,mpid%mpirank_3d(2),  mpid%mpirank_3d(3)  )
    mpid%f = ptable(mpid%mpirank_3d(1),  mpid%mpirank_3d(2)+1,mpid%mpirank_3d(3)  )
    mpid%b = ptable(mpid%mpirank_3d(1),  mpid%mpirank_3d(2)-1,mpid%mpirank_3d(3)  )
    mpid%t = ptable(mpid%mpirank_3d(1),  mpid%mpirank_3d(2),  mpid%mpirank_3d(3)+1)
    mpid%d = ptable(mpid%mpirank_3d(1),  mpid%mpirank_3d(2),  mpid%mpirank_3d(3)-1)

    write(*,*) 'myrank :: ',mpid%mpirank,mpid%mpirank_3d

  end subroutine mpi_setup__init


end module mpi_setup
