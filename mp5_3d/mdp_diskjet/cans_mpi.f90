module cans_mpi
!======================================================================!
! NAME: cans_mpi
!
! PURPOSE: Set MPI configureation
!
! mpi_comm_cort3d [integer] communicator of 
!                              a virtual Cartesian topology 
! n_mpi_dims      [integer] number of dimension of domain decomposition
! periods         [lgt]     specifier whether the grid is 
!                              periodic (.true.) or not (.false.) 
!                              in each dimension
! reorder         [lgt]     ranking may reordered (.true.) 
!                              or not (.false.)
!
! CONTENTS:
!
! NOTE: in type cans_mpi_rank, 
!       the location of rank in each direction is denoted with
!          n: Nutral (calling process) location
!          u: Up direction 
!          d: Down direction 
!       e.g., nnn --> calling process rank, 
!             unn --> Up direction in x-coords
!             udn --> Up in x-coords and Down in y-coords
!     
!
! UPDATE LOG:
! 2012/12/30 initial coding (by H. ODA)
!======================================================================!

  use cans_type

  implicit none

  type cans_mpi_rank
     !-- local rank (calling process) --!
     integer :: nnn
     !-- ranks of 6 neighbor processes in [x|y|z] direction --!
     integer :: unn, dnn, nun, ndn, nnu, nnd
     !-- ranks of 4 diagonal processes in xy plane --!
     integer :: uun, udn, dun, ddn
     !-- ranks of 4 diagonal processes in xz plane --!
     integer :: unu, und, dnu, dnd
     !-- ranks of 4 diagonal processes in yz plane --!
     integer :: nuu, nud, ndu, ndd
     !-- ranks of 8 diagonal processes at vertex --!
     integer :: uuu, uud, udu, udd, duu, dud, ddu, ddd
  end type cans_mpi_rank

  integer,parameter :: n_mpi_dims=3
  logical(lgt),parameter :: periods(n_mpi_dims)=(/.false., .true., .false. /)
  logical(lgt),parameter :: reorder=.false.

contains
  
  subroutine cans_mpi_cart_rank &
       & (  mpi_comm, n_mpi_dims, n_tiles, displ_x, displ_y, displ_z, topo_coords &
       &  , rank, err)
    
    implicit none
    include 'mpif.h'
    
    integer,intent(in) :: mpi_comm, n_mpi_dims, n_tiles(n_mpi_dims) &
         & , displ_x, displ_y, displ_z, topo_coords(n_mpi_dims)
    integer,intent(out) :: rank
    integer,intent(inout) :: err

    integer :: i, buf_topo_coords(n_mpi_dims)

    
    buf_topo_coords(1)  = topo_coords(1) + displ_x
    buf_topo_coords(2)  = topo_coords(2) + displ_y
    buf_topo_coords(3)  = topo_coords(3) + displ_z

    if ( all(buf_topo_coords < n_tiles) .and. all(buf_topo_coords >= 0) ) then
       call mpi_cart_rank &
            & ( mpi_comm, buf_topo_coords, rank, err)
    else
       rank = mpi_proc_null
    end if

    
  end subroutine cans_mpi_cart_rank
  


end module cans_mpi
