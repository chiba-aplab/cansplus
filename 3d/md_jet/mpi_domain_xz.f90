module mpi_domain_xz
!======================================================================
! example ::
!
!   16,17,18,19 
!   12,13,14,15 
! k  8, 9,10,11 
! ^  4, 5, 6, 7 
! |  0, 1, 2, 3 
!   --> i
!    npe_x = 4, npe_z = 5
! #surface exchange
!     8, 9,10
!     4, 5, 6
!     0, 1, 2
!     5 => l :: 4, r :: 6
!       => t :: 9, d :: 1
!
! #point exchange
!     8, 9,10
!     4, 5, 6
!     0, 1, 2
!     5 => tl :: 8, tr :: 10
!       => dl :: 0, dr :: 2
!======================================================================

  implicit none
  private

  public :: mpi_setup

  type :: mpidomain
     integer :: mpisize
     integer :: mpirank
     integer,dimension(2) :: mpisize_2d
     integer,dimension(2) :: mpirank_2d

     integer :: l,r,t,d

     integer :: tl,tr,dl,dr
  end type mpidomain

  type(mpidomain), public :: mpid

contains

subroutine mpi_setup(mpisize_x,mpisize_z)

  include 'mpif.h'

  integer, intent(in) :: mpisize_x, mpisize_z
  integer :: mpisize,mpirank,merr

  call mpi_init(merr)
  call mpi_comm_size(mpi_comm_world,mpisize,merr)
  call mpi_comm_rank(mpi_comm_world,mpirank,merr)
  
  mpid%mpirank = mpirank
  mpid%mpisize = mpisize
  
  ! Core Numbers, 1: r(i)-direction, 2: z(k)-direction
  mpid%mpisize_2d(1) = mpisize_x !1
  mpid%mpisize_2d(2) = mpisize_z !16
  
  ! determine mpid%mpirank_3d
  call setmy2drank(mpid,merr)
  call setmpiboundary(mpid)

end subroutine mpi_setup

!======================================================================
  subroutine setmy2drank(mpid,merr)
!======================================================================
! Name :: initmy3drank
!         set mpisize3d & mpirank_3d
!======================================================================
    implicit none

    type(mpidomain) :: mpid
    integer :: mpirank_c

    integer :: merr

    integer :: i,j,k

    if(mpid%mpisize_2d(1)*mpid%mpisize_2d(2) .eq. mpid%mpisize)then

       mpid%mpirank_2d(1) = mod(mpid%mpirank,mpid%mpisize_2d(1))
       mpid%mpirank_2d(2) = int(mpid%mpirank/mpid%mpisize_2d(1))
       
       write(*,*) 'myrank :: ',mpid%mpirank,mpid%mpirank_2d
       merr = 0
    else
       merr = 1
    end if
    return
  end subroutine setmy2drank

!======================================================================
  subroutine setmpiboundary(mpid)
!======================================================================
! Name :: setmpiboundary
!         set mpi boundary (xleft, xright, yleft, yright, .. , lu, ld, .. lul
!======================================================================
    implicit none
    include 'mpif.h'

    type(mpidomain) :: mpid,mpid_next
    integer,dimension(2) :: mpir_2d_n

    integer :: neighbor
    
    integer :: i,j,k

!=====================================================================
! surface exchange
!=====================================================================
!---Step 1a.----
! x-direction
!

!-------
! right

    mpir_2d_n(1) = mpid%mpirank_2d(1)
    mpir_2d_n(2) = mpid%mpirank_2d(2)

    mpir_2d_n(1) = mpir_2d_n(1)+1

    call sectorank(mpid,mpir_2d_n,neighbor)
    
    call chkrank(mpid, mpir_2d_n, neighbor)

    mpid%r = neighbor

!-------
! left

    mpir_2d_n(1) = mpid%mpirank_2d(1)
    mpir_2d_n(2) = mpid%mpirank_2d(2)

    mpir_2d_n(1) = mpir_2d_n(1)-1

    call sectorank(mpid,mpir_2d_n,neighbor)
    
    call chkrank(mpid, mpir_2d_n, neighbor)

    mpid%l = neighbor

!---Step 1b.----
! z-direction
!

!-------
! top

    mpir_2d_n(1) = mpid%mpirank_2d(1)
    mpir_2d_n(2) = mpid%mpirank_2d(2)

    mpir_2d_n(2) = mpir_2d_n(2)+1

    call sectorank(mpid,mpir_2d_n,neighbor)
    
    call chkrank(mpid, mpir_2d_n, neighbor)

    mpid%t = neighbor

!-------
! down

    mpir_2d_n(1) = mpid%mpirank_2d(1)
    mpir_2d_n(2) = mpid%mpirank_2d(2)

    mpir_2d_n(2) = mpir_2d_n(2)-1

    call sectorank(mpid,mpir_2d_n,neighbor)
    
    call chkrank(mpid, mpir_2d_n, neighbor)

    mpid%d = neighbor

!---Step 2a.----
! top-right
!
    mpir_2d_n(1) = mpid%mpirank_2d(1)
    mpir_2d_n(2) = mpid%mpirank_2d(2)

    mpir_2d_n(1) = mpir_2d_n(1)+1
    mpir_2d_n(2) = mpir_2d_n(2)+1

    call sectorank(mpid,mpir_2d_n,neighbor)
    
    call chkrank(mpid, mpir_2d_n, neighbor)

    mpid%tr = neighbor

!---Step 2b.----
! top-left
!
    mpir_2d_n(1) = mpid%mpirank_2d(1)
    mpir_2d_n(2) = mpid%mpirank_2d(2)

    mpir_2d_n(1) = mpir_2d_n(1)-1
    mpir_2d_n(2) = mpir_2d_n(2)+1

    call sectorank(mpid,mpir_2d_n,neighbor)
    
    call chkrank(mpid, mpir_2d_n, neighbor)

    mpid%tl = neighbor

!---Step 2a.----
! down-right
!
    mpir_2d_n(1) = mpid%mpirank_2d(1)
    mpir_2d_n(2) = mpid%mpirank_2d(2)

    mpir_2d_n(1) = mpir_2d_n(1)+1
    mpir_2d_n(2) = mpir_2d_n(2)-1

    call sectorank(mpid,mpir_2d_n,neighbor)
    
    call chkrank(mpid, mpir_2d_n, neighbor)

    mpid%dr = neighbor

!---Step 2b.----
! down-left
!
    mpir_2d_n(1) = mpid%mpirank_2d(1)
    mpir_2d_n(2) = mpid%mpirank_2d(2)

    mpir_2d_n(1) = mpir_2d_n(1)-1
    mpir_2d_n(2) = mpir_2d_n(2)-1

    call sectorank(mpid,mpir_2d_n,neighbor)
    
    call chkrank(mpid, mpir_2d_n, neighbor)

    mpid%dl = neighbor
    
  contains
!======================================================================
    subroutine sectorank(mpid,rank_2d,rank)
!======================================================================
! Name :: thdtorank
!         rank3d => rank
!======================================================================
      implicit none
      
      type(mpidomain) :: mpid
    
      integer,dimension(2) :: rank_2d
      integer :: rank

      rank = 0
      
      rank = mpid%mpisize_2d(1)*rank_2d(2)
      rank = rank + rank_2d(1)
    
      return
    end subroutine sectorank

!======================================================================
    subroutine chkrank(mpid,rank_2d,rank)
!======================================================================
! Name :: chkrank
!         check MPI boundary
!======================================================================
      implicit none
      include 'mpif.h'
      
      type(mpidomain) :: mpid
      
      integer,dimension(2) :: rank_2d
      integer :: rank

      if (rank_2d(1) < 0) then
         rank = mpi_proc_null
      endif
      if (rank_2d(1) > (mpid%mpisize_2d(1)-1))then
         rank = mpi_proc_null
      end if

      if (rank_2d(2) < 0) then
         rank = mpi_proc_null
      endif
      if (rank_2d(2) > (mpid%mpisize_2d(2)-1)) then
         rank = mpi_proc_null
      endif

      return
    end subroutine chkrank
  end subroutine setmpiboundary
end module mpi_domain_xz
