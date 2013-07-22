module mpi_domain
!=====================================================================
! example ::
!      
!   16,17,18,19  36,37,38,39  56,57,58,59
!   12,13,14,15  32,33,34,35  52,53,54,55
! j  8, 9,10,11  28,29,30,31  48,49,50,51
! ^  4, 5, 6, 7  24,25,26,27  44,45,46,47
! |  0, 1, 2, 3  20,21,22,23  40,41,42,43
!   --> i
!    npe_x = 4, npe_y = 5, npe_z = 3
! #surface exchange   
!   @ x-exchange
!      mleft = myrank-1
!      mright = myrank+1
!      @ x-boundary
!        left boundary = mod(myrank,npe_x) = 0
!       right boundary = mod(myrank,npe_x) = npe_x-1
!   @ y-exchange
!      mleft = myrank-npe_x
!      mright = myrank+npe_x
!      @ y-boundary
!        left boundary = 0 <= mod(myrank,npe_x*npe_y) <= npe_x-1
!       right boundary = npe_x*(npe_y-1) <= mod(myrank,npe_x*npe_y) <= npe_x*npe_y-1
!   @ z-exhange
!      mleft = myrank-npe_x*npe_y
!      mright = myrank+npe_x*npe_y
!      @ z-boundary
!        left boundary = 0 <= mod(myrank,npe_x*npe_y*npez) <= npe_x*npe_y-1
!       right boundary = npe_x*npe_y*(npe_z-1) <= mod(myrank,npe_x*npe_y*npe_z) <= npe_x*npe_y*npe_z-1
! #side exchange
!     8, 9,10  28,29,30  48,49,50
!     4, 5, 6  24,25,26  44,45,46
!     0, 1, 2  20,21,22  40,41,42
!        25 => lu:9  ld: 1 ll: 4 lr:6
!           => tu:49 td:41 tl:44 tr:46
!           => mul:28 mur:30 mdl:20 mdr:22
!
!               myrank             cell boundary
!                                    send  recv
!        lu = yright(zleft(myrank))   td -> lu
!        ld = yleft (zleft(myrank))   tu -> ld
!        ll = xleft (zleft(myrank))   tr -> ll
!        lr = xright(zleft(myrank))   tl -> lr
!
!       mul = xleft (yright(myrank)) mdr -> mul
!       mur = xright(yright(myrank)) mdl -> mur
!       mdl = xleft (yleft (myrank)) mur -> mdl
!       mdr = xright(yleft (myrank)) mul -> mdr
!
!        tu = yright(zright(myrank))  ld -> tu
!        td = yleft (zright(myrank))  lu -> td
!        tl = xleft (zright(myrank))  lr -> tl
!        tr = xright(zright(myrank))  ll -> tr
!
! #point exchange
!     8, 9,10  28,29,30  48,49,50
!     4, 5, 6  24,25,26  44,45,46
!     0, 1, 2  20,21,22  40,41,42
!        25 => lul:8  lur:10 ldl:0  ldr:2
!              tul:48 tur:50 tdl:40 tdr:42
!
!              myrank                           cell boundary
!                                               send    recv
!        lul = xleft (yright(zleft(myrank)))     tdr -> lul
!        lur = xright(yright(zleft(myrank)))     tdl -> lur
!        ldl = xleft (yleft (zleft(myrank)))     tur -> ldl
!        ldr = xright(yleft (zleft(myrank)))     tul -> ldr
!
!        tul = xleft (yright(zright(myrank)))    ldr -> tul
!        tur = xright(yright(zright(myrank)))    ldl -> tur
!        tdl = xleft (yleft (zright(myrank)))    lur -> tdl
!        tdr = xright(yleft (zright(myrank)))    lul -> tdr
!=====================================================================

  implicit none

  type mpidomain
     integer :: mpisize
     integer :: mpirank
     integer,dimension(3) :: mpisize_3d
     integer,dimension(3) :: mpirank_3d

!---note --------------
! left is lowwer,right is upper
! surface exchange
!--------------
     integer :: xleft,xright
     integer :: yleft,yright
     integer :: zleft,zright

!---note --------------
! side exchange
!
     integer ::  lu,  ld,  ll,  lr
     integer :: mul, mur, mdl, mdr
     integer ::  tu,  td,  tl,  tr

!---note --------------
! point exchange
!
     integer :: lul,lur,ldl,ldr
     integer :: tul,tur,tdl,tdr

  end type mpidomain
contains

!======================================================================
  subroutine setmy3drank(mpid,merr)
!======================================================================
! Name :: initmy3drank
!         set mpisize3d & mpirank_3d
!======================================================================
    implicit none

    type(mpidomain) :: mpid
    integer :: mpirank_c

    integer :: merr
    integer :: i,j,k

    if(mpid%mpisize_3d(1)*mpid%mpisize_3d(2)*mpid%mpisize_3d(3) .eq. mpid%mpisize)then
!------------------------
! set mpirank_z
!
       mpirank_c = mod(mpid%mpirank,(mpid%mpisize_3d(1)*mpid%mpisize_3d(2)))
       mpid%mpirank_3d(3) = int(mpid%mpirank/(mpid%mpisize_3d(1)*mpid%mpisize_3d(2)))

!------------------------
! set mpirank_y
!
       mpid%mpirank_3d(2) = int(mpirank_c/mpid%mpisize_3d(1))

!------------------------
! set mpirank_x
!
       mpid%mpirank_3d(1) = mod(mpirank_c,mpid%mpisize_3d(1))

    else
       merr = 1
    endif

  end subroutine setmy3drank

!======================================================================
  subroutine setmpiboundary(mpid)
!======================================================================
! Name :: setmpiboundary
!         set mpi boundary (xleft, xright, yleft, yright, .. , lu, ld, .. lul
!======================================================================
    implicit none

    include 'mpif.h'
    
    type(mpidomain) :: mpid,mpid_next
    integer,dimension(3) :: mpir_3d_n
    
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

    call copyrank(mpid,mpir_3d_n)
    mpir_3d_n(1) = mpir_3d_n(1) + 1
    
    call thdtorank(mpid,mpir_3d_n,neighbor)
    
    mpid%xright = neighbor
    
!-------
! left

    call copyrank(mpid,mpir_3d_n)
    mpir_3d_n(1) = mpir_3d_n(1) - 1
    
    call thdtorank(mpid,mpir_3d_n,neighbor)
  
    mpid%xleft  = neighbor

    if(mpid%mpirank_3d(1) .eq. 0) then
       mpid%xleft = mpi_proc_null
    endif
    if(mpid%mpirank_3d(1) .eq. mpid%mpisize_3d(1)-1) then
       mpid%xright = mpi_proc_null
    endif

!---Step 1b.----
! y-direction
!

!-------
! right
    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(2) = mpir_3d_n(2) + 1
    
    call thdtorank(mpid,mpir_3d_n,neighbor)
  
    mpid%yright = neighbor

!-------
! left

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(2) = mpir_3d_n(2) - 1

    call thdtorank(mpid,mpir_3d_n,neighbor)

    mpid%yleft  = neighbor

    if(mpid%mpirank_3d(2) .eq. 0)then
       mpid%yleft = mpi_proc_null
    endif
    if(mpid%mpirank_3d(2) .eq. mpid%mpisize_3d(2)-1) then
       mpid%yright = mpi_proc_null
    endif

!---Step 1c.----
! z-direction
!

!-------
! right

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) + 1
    
    call thdtorank(mpid,mpir_3d_n,neighbor)

    mpid%zright = neighbor

!-------
! left

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) - 1

    call thdtorank(mpid,mpir_3d_n,neighbor)

    mpid%zleft  = neighbor

    if(mpid%mpirank_3d(3) .eq. 0) then
       mpid%zleft = mpi_proc_null
    endif
    if(mpid%mpirank_3d(3) .eq. mpid%mpisize_3d(3)-1)then
       mpid%zright = mpi_proc_null
    endif

!=====================================================================
! line exchange
!=====================================================================
!---Step 2a.----
! lower
!

!---------
! lu

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) - 1
    mpir_3d_n(2) = mpir_3d_n(2) + 1
    
    call thdtorank(mpid,mpir_3d_n,neighbor)
    
    call chkrank(mpid, mpir_3d_n, neighbor)

    mpid%lu = neighbor

!---------
! ld

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) - 1
    mpir_3d_n(2) = mpir_3d_n(2) - 1

    call thdtorank(mpid,mpir_3d_n,neighbor)
    
    call chkrank(mpid, mpir_3d_n, neighbor)

    mpid%ld = neighbor

!---------
! ll

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) - 1
    mpir_3d_n(1) = mpir_3d_n(1) - 1

    call thdtorank(mpid,mpir_3d_n,neighbor)
    
    call chkrank(mpid, mpir_3d_n, neighbor)

    mpid%ll = neighbor

!---------
! lr

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) - 1
    mpir_3d_n(1) = mpir_3d_n(1) + 1

    call thdtorank(mpid,mpir_3d_n,neighbor)

    call chkrank(mpid, mpir_3d_n, neighbor)

    mpid%lr = neighbor

!---Step 2b.----
! middle
!

!---------
! mul

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(2) = mpir_3d_n(2) + 1
    mpir_3d_n(1) = mpir_3d_n(1) - 1

    call thdtorank(mpid,mpir_3d_n,neighbor)

    call chkrank(mpid, mpir_3d_n, neighbor)

    mpid%mul = neighbor

!---------
! mur

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(2) = mpir_3d_n(2) + 1
    mpir_3d_n(1) = mpir_3d_n(1) + 1

    call thdtorank(mpid,mpir_3d_n,neighbor)

    call chkrank(mpid, mpir_3d_n, neighbor)

    mpid%mur = neighbor

!---------
! mdl

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(2) = mpir_3d_n(2) - 1
    mpir_3d_n(1) = mpir_3d_n(1) - 1

    call thdtorank(mpid,mpir_3d_n,neighbor)

    call chkrank(mpid, mpir_3d_n, neighbor)

    mpid%mdl = neighbor

!---------
! mdr

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(2) = mpir_3d_n(2) - 1
    mpir_3d_n(1) = mpir_3d_n(1) + 1

    call thdtorank(mpid,mpir_3d_n,neighbor)

    call chkrank(mpid, mpir_3d_n, neighbor)

    mpid%mdr = neighbor

!---Step 2c.----
! top
!

!---------
! tu

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) + 1
    mpir_3d_n(2) = mpir_3d_n(2) + 1

    call thdtorank(mpid,mpir_3d_n,neighbor)
    
    call chkrank(mpid, mpir_3d_n, neighbor)
    
    mpid%tu = neighbor

!---------
! td

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) + 1
    mpir_3d_n(2) = mpir_3d_n(2) - 1

    call thdtorank(mpid,mpir_3d_n,neighbor)

    call chkrank(mpid, mpir_3d_n, neighbor)

    mpid%td = neighbor

!---------
! tl

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) + 1
    mpir_3d_n(1) = mpir_3d_n(1) - 1

    call thdtorank(mpid,mpir_3d_n,neighbor)
    
    call chkrank(mpid, mpir_3d_n, neighbor)
    
    mpid%tl = neighbor

!---------
! tr

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) + 1
    mpir_3d_n(1) = mpir_3d_n(1) + 1

    call thdtorank(mpid,mpir_3d_n,neighbor)
    
    call chkrank(mpid, mpir_3d_n, neighbor)
    
    mpid%tr = neighbor

!=====================================================================
! point exchange
!=====================================================================
!---Step 3a.----
! low
!---------
! lul

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) - 1
    mpir_3d_n(2) = mpir_3d_n(2) + 1
    mpir_3d_n(1) = mpir_3d_n(1) - 1

    call thdtorank(mpid,mpir_3d_n,neighbor)

    call chkrank(mpid, mpir_3d_n, neighbor)

    mpid%lul = neighbor

!---------
! lur

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) - 1
    mpir_3d_n(2) = mpir_3d_n(2) + 1
    mpir_3d_n(1) = mpir_3d_n(1) + 1

    call thdtorank(mpid,mpir_3d_n,neighbor)

    call chkrank(mpid, mpir_3d_n, neighbor)

    mpid%lur = neighbor
  
!---------
! ldl

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) - 1
    mpir_3d_n(2) = mpir_3d_n(2) - 1
    mpir_3d_n(1) = mpir_3d_n(1) - 1

    call thdtorank(mpid,mpir_3d_n,neighbor)

    call chkrank(mpid, mpir_3d_n, neighbor)

    mpid%ldl = neighbor

!---------
! ldr

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) - 1
    mpir_3d_n(2) = mpir_3d_n(2) - 1
    mpir_3d_n(1) = mpir_3d_n(1) + 1

    call thdtorank(mpid,mpir_3d_n,neighbor)

    call chkrank(mpid, mpir_3d_n, neighbor)

    mpid%ldr = neighbor

!---Step 3b.----
! top
!---------
! tul

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) + 1
    mpir_3d_n(2) = mpir_3d_n(2) + 1
    mpir_3d_n(1) = mpir_3d_n(1) - 1
    
    call thdtorank(mpid,mpir_3d_n,neighbor)
  
    call chkrank(mpid, mpir_3d_n, neighbor)

    mpid%tul = neighbor

!---------
! tur

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) + 1
    mpir_3d_n(2) = mpir_3d_n(2) + 1
    mpir_3d_n(1) = mpir_3d_n(1) + 1

    call thdtorank(mpid,mpir_3d_n,neighbor)
    
    call chkrank(mpid, mpir_3d_n, neighbor)
    
    mpid%tur = neighbor
  
!---------
! tdl

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) + 1
    mpir_3d_n(2) = mpir_3d_n(2) - 1
    mpir_3d_n(1) = mpir_3d_n(1) - 1
    
    call thdtorank(mpid,mpir_3d_n,neighbor)

    call chkrank(mpid, mpir_3d_n, neighbor)

    mpid%tdl = neighbor

!---------
! tdr

    mpir_3d_n = mpid%mpirank_3d
    mpir_3d_n(3) = mpir_3d_n(3) + 1
    mpir_3d_n(2) = mpir_3d_n(2) - 1
    mpir_3d_n(1) = mpir_3d_n(1) + 1

    call thdtorank(mpid,mpir_3d_n,neighbor)

    call chkrank(mpid, mpir_3d_n, neighbor)

    mpid%tdr = neighbor

    return
  contains

!======================================================================
    subroutine copyrank(mpid,rank_3d)
!======================================================================
! Name :: copyrank
!         rank copy
!======================================================================

      implicit none

      type(mpidomain) :: mpid

      integer,dimension(3) :: rank_3d

      integer :: i

      do i=1,3
         rank_3d(i) = mpid%mpirank_3d(i)
      enddo
      return
    end subroutine copyrank

!======================================================================
    subroutine thdtorank(mpid,rank_3d,rank)
!======================================================================
! Name :: thdtorank
!         rank3d => rank
!======================================================================
      implicit none
      
      type(mpidomain) :: mpid
    
      integer,dimension(3) :: rank_3d
      integer :: rank

      rank = 0
      
      rank = mpid%mpisize_3d(2)*mpid%mpisize_3d(1)*rank_3d(3)
      rank = rank + mpid%mpisize_3d(1)*rank_3d(2)
      rank = rank + rank_3d(1)
    
      return
    end subroutine thdtorank

!======================================================================
    subroutine chkrank(mpid, rank_3d, rank)
!======================================================================
! Name :: chkrank
!         check MPI boundary
!======================================================================
      implicit none
      
      include 'mpif.h'
      
      type(mpidomain) :: mpid
      
      integer,dimension(3) :: rank_3d
      integer :: rank

      if (rank_3d(1) < 0) then
         rank = mpi_proc_null
      endif
      if (rank_3d(1) > (mpid%mpisize_3d(1) -1)) then
         rank = mpi_proc_null
      endif
      
      if (rank_3d(2) < 0) then
         rank = mpi_proc_null
      endif
      if (rank_3d(2) > (mpid%mpisize_3d(2) -1)) then
         rank = mpi_proc_null
      endif

      if (rank_3d(3) < 0) then
         rank = mpi_proc_null
      endif
      if (rank_3d(3) > (mpid%mpisize_3d(3) -1)) then
         rank = mpi_proc_null
      endif

      return
    end subroutine chkrank
  end subroutine setmpiboundary
end module mpi_domain
