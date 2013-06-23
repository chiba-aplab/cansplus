!======================================================================|
subroutine cans_mpi_exchange_cart3d(ro, pr, vx, vy, vz, bx, by, bz, phi &
     & , ix, jx, kx, margin, mpi_comm, npe, rank, mpi_err)
!======================================================================|
! NAME : cans_mpi_exchange_cart3d
!
! PURPOSE: MPI data exchane
!          use a virtual Cartesian topology for the domain decompositon
!          (communicator: mpi_comm_cart3d)
!
! INPUTS:
!    margin   [integer]             margin outside the boundary
!    ix,jx,kx [integer]             grid number
!    mpi_comm [integer]             MPI communicator
!    npe      [integer]             number of processes
!    rank     [type(cans_mpi_rank)] structure for MPI rank
! 
! INOUTPUTS:
!    ro (ix, jx, kx) [real(8)] : density
!    pr (ix, jx, kx) [real(8)] : pressure
!    vx (ix, jx, kx) [real(8)] : x-velocity 
!    vy (ix, jx, kx) [real(8)] : y-velocity
!    vz (ix, jx, kx) [real(8)] : z-velocity
!    bx (ix, jx, kx) [real(8)] : x-magnetic field
!    by (ix, jx, kx) [real(8)] : y-magnetic field
!    bz (ix, jx, kx) [real(8)] : z-magnetic field
!    phi(ix, jx, kx) [real(8)] : 
!    
! OUTPUTS:
!    mpi_err [integer] error counter
!
! INTERNAL VARIABLES
!    n_element_var [integer]         number of variables for buffer
!    bufsnd_?      [real(real_kind)] buffer for send data in ?-coord
!    bufrcv_?      [real(real_kind)] buffer for receive data  in ?-coord
!    mstatus       [integer]         status object for MPI_sendrecv
!    n_element_buf [integer]         total number of element of buffer
!
! UPDATE LOG:
! 2012/12/28 initial coding  (by H. ODA)
!            implement 3D MPI parallelization in 3D Cartesian topology
! 
!======================================================================!
  use cans_type
  use cans_mpi

  implicit none

  include 'mpif.h'

  integer,intent(in) :: ix,jx,kx,margin

  !type(mpidomain),intent(in) :: mpid
  integer,intent(in) :: mpi_comm, npe
  type(cans_mpi_rank),intent(in) :: rank

  integer,dimension(mpi_status_size) :: mstatus

  real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz,phi
  
  integer,intent(out) :: mpi_err

  !integer,parameter :: mx = 9
  integer,parameter :: n_element_var = 9

!-- surface exchange --!
  real(8),dimension(margin, jx, kx, n_element_var) :: bufsnd_x, bufrcv_x
  real(8),dimension(ix, margin, kx, n_element_var) :: bufsnd_y, bufrcv_y
  real(8),dimension(ix, jx, margin, n_element_var) :: bufsnd_z, bufrcv_z
!----------------------------------------------------------------------|


!-- line exchange --!
!  real(8),dimension(margin, jx, margin, n_element_var) :: bufsnd_xz, bufrcv_xz

  !integer :: mmx, msend, mrecv
  integer :: n_element_buf, rank_send, rank_recv

  !integer :: merr

  integer :: i,j,k

!=== Step 1 ==========================================================
! surface exchange
!=====================================================================

!----------------------------------------------------------------------|
!-- Step 1-1: data exchange in downward surface in x-coordinate 
!----------------------------------------------------------------------|

  !-- set number of elements of buffer and ranks of send/recv --!
  n_element_buf = margin * jx * kx * n_element_var
  rank_send = rank%dnn
  rank_recv = rank%unn

  !-- set send buffer --!
  do k=1,kx
     do j=1,jx
        do i=1,margin
           bufsnd_x(i,j,k,1) = ro(margin+i,j,k)
           bufsnd_x(i,j,k,2) = pr(margin+i,j,k)

           bufsnd_x(i,j,k,3) = vx(margin+i,j,k)
           bufsnd_x(i,j,k,4) = vy(margin+i,j,k)
           bufsnd_x(i,j,k,5) = vz(margin+i,j,k)

           bufsnd_x(i,j,k,6) = bx(margin+i,j,k)
           bufsnd_x(i,j,k,7) = by(margin+i,j,k)
           bufsnd_x(i,j,k,8) = bz(margin+i,j,k)

           bufsnd_x(i,j,k,9) = phi(margin+i,j,k)
        enddo
     end do
  end do

  !-- MPI send receive --!
  call mpi_sendrecv &
       & (bufsnd_x, n_element_buf, mpi_double_precision, rank_send, 0 &
       & ,bufrcv_x, n_element_buf, mpi_double_precision, rank_recv, 0 &
       & ,mpi_comm, mstatus, mpi_err)

  !-- store receive buffer data into arrays for variables --!
  if (rank_recv /= mpi_proc_null) then
     do k = 1, kx
        do j = 1, jx
           do i = 1, margin
              ro(ix-margin+i,j,k) = bufrcv_x(i,j,k,1)
              pr(ix-margin+i,j,k) = bufrcv_x(i,j,k,2)

              vx(ix-margin+i,j,k) = bufrcv_x(i,j,k,3)
              vy(ix-margin+i,j,k) = bufrcv_x(i,j,k,4)
              vz(ix-margin+i,j,k) = bufrcv_x(i,j,k,5)

              bx(ix-margin+i,j,k) = bufrcv_x(i,j,k,6)
              by(ix-margin+i,j,k) = bufrcv_x(i,j,k,7)
              bz(ix-margin+i,j,k) = bufrcv_x(i,j,k,8)

              phi(ix-margin+i,j,k) = bufrcv_x(i,j,k,9)
           enddo
        enddo
     enddo
  endif

!----------------------------------------------------------------------|
!-- Step 1-2: data exchange in upward surface in x-coordinate 
!----------------------------------------------------------------------|

  !-- set number of elements of buffer and ranks of send/recv --!
  n_element_buf = margin * jx * kx * n_element_var
  rank_send = rank%unn
  rank_recv = rank%dnn

  !-- set send buffer --!
  do k=1,kx
     do j=1,jx
        do i=1,margin
           bufsnd_x(i,j,k,1) = ro(ix-2*margin+i,j,k)
           bufsnd_x(i,j,k,2) = pr(ix-2*margin+i,j,k)

           bufsnd_x(i,j,k,3) = vx(ix-2*margin+i,j,k)
           bufsnd_x(i,j,k,4) = vy(ix-2*margin+i,j,k)
           bufsnd_x(i,j,k,5) = vz(ix-2*margin+i,j,k)

           bufsnd_x(i,j,k,6) = bx(ix-2*margin+i,j,k)
           bufsnd_x(i,j,k,7) = by(ix-2*margin+i,j,k)
           bufsnd_x(i,j,k,8) = bz(ix-2*margin+i,j,k)

           bufsnd_x(i,j,k,9) = phi(ix-2*margin+i,j,k)
        end do
     end do
  end do

  !-- MPI send receive --!
  call mpi_sendrecv &
       & (bufsnd_x, n_element_buf, mpi_double_precision, rank_send, 0 &
       & ,bufrcv_x, n_element_buf, mpi_double_precision, rank_recv, 0 &
       & ,mpi_comm, mstatus, mpi_err)

  !-- store receive buffer data into arrays for variables --!
  if (rank_recv /= mpi_proc_null) then
     do k = 1, kx
        do j = 1, jx
           do i = 1, margin
              ro(i,j,k) = bufrcv_x(i,j,k,1)
              pr(i,j,k) = bufrcv_x(i,j,k,2)

              vx(i,j,k) = bufrcv_x(i,j,k,3)
              vy(i,j,k) = bufrcv_x(i,j,k,4)
              vz(i,j,k) = bufrcv_x(i,j,k,5)

              bx(i,j,k) = bufrcv_x(i,j,k,6)
              by(i,j,k) = bufrcv_x(i,j,k,7)
              bz(i,j,k) = bufrcv_x(i,j,k,8)

              phi(i,j,k) = bufrcv_x(i,j,k,9)
           end do
        end do
     end do
  end if

!----------------------------------------------------------------------|
!-- Step 1-3: data exchange in downward surface in y-coordinate 
!----------------------------------------------------------------------|

  !-- set number of elements of buffer and ranks of send/recv --!
  n_element_buf = ix * margin * kx * n_element_var
  rank_send = rank%ndn
  rank_recv = rank%nun

!!$  !-- set send buffer --!
!!$  do k = 1, kx
!!$     do j = 1, margin
!!$        do i = 1, ix
!!$           bufsnd_y(i,j,k,1) = ro(i,margin+j,k)
!!$           bufsnd_y(i,j,k,2) = pr(i,margin+j,k)
!!$
!!$           bufsnd_y(i,j,k,3) = vx(i,margin+j,k)
!!$           bufsnd_y(i,j,k,4) = vy(i,margin+j,k)
!!$           bufsnd_y(i,j,k,5) = vz(i,margin+j,k)
!!$
!!$           bufsnd_y(i,j,k,6) = bx(i,margin+j,k)
!!$           bufsnd_y(i,j,k,7) = by(i,margin+j,k)
!!$           bufsnd_y(i,j,k,8) = bz(i,margin+j,k)
!!$
!!$           bufsnd_y(i,j,k,9) = phi(i,margin+j,k)
!!$        enddo
!!$     end do
!!$  end do
!!$
!!$  !-- MPI send receive --!
!!$  call mpi_sendrecv &
!!$       & (bufsnd_y, n_element_buf, mpi_double_precision, rank_send, 0 &
!!$       & ,bufrcv_y, n_element_buf, mpi_double_precision, rank_recv, 0 &
!!$       & ,mpi_comm, mstatus, mpi_err)
!!$
!!$  !-- store receive buffer data into arrays for variables --!
!!$  if (rank_recv /= mpi_proc_null) then
!!$     do k = 1, kx
!!$        do j = 1, margin
!!$           do i = 1, ix
!!$              ro(i,jx-margin+j,k) = bufrcv_y(i,j,k,1)
!!$              pr(i,jx-margin+j,k) = bufrcv_y(i,j,k,2)
!!$
!!$              vx(i,jx-margin+j,k) = bufrcv_y(i,j,k,3)
!!$              vy(i,jx-margin+j,k) = bufrcv_y(i,j,k,4)
!!$              vz(i,jx-margin+j,k) = bufrcv_y(i,j,k,5)
!!$
!!$              bx(i,jx-margin+j,k) = bufrcv_y(i,j,k,6)
!!$              by(i,jx-margin+j,k) = bufrcv_y(i,j,k,7)
!!$              bz(i,jx-margin+j,k) = bufrcv_y(i,j,k,8)
!!$
!!$              phi(i,jx-margin+j,k) = bufrcv_y(i,j,k,9)
!!$           enddo
!!$        enddo
!!$     enddo
!!$  endif

!----------------------------------------------------------------------|
!-- Step 1-4: data exchange in upward surface in y-coordinate 
!----------------------------------------------------------------------|

  !-- set number of elements of buffer and ranks of send/recv --!
  n_element_buf = ix * margin * kx * n_element_var
  rank_send = rank%nun
  rank_recv = rank%ndn

!!$  !-- set send buffer --!
!!$  do k = 1, kx
!!$     do j = 1, margin
!!$        do i = 1, ix
!!$           bufsnd_y(i,j,k,1) = ro(i,jx-2*margin+j,k)
!!$           bufsnd_y(i,j,k,2) = pr(i,jx-2*margin+j,k)
!!$
!!$           bufsnd_y(i,j,k,3) = vx(i,jx-2*margin+j,k)
!!$           bufsnd_y(i,j,k,4) = vy(i,jx-2*margin+j,k)
!!$           bufsnd_y(i,j,k,5) = vz(i,jx-2*margin+j,k)
!!$
!!$           bufsnd_y(i,j,k,6) = bx(i,jx-2*margin+j,k)
!!$           bufsnd_y(i,j,k,7) = by(i,jx-2*margin+j,k)
!!$           bufsnd_y(i,j,k,8) = bz(i,jx-2*margin+j,k)
!!$
!!$           bufsnd_y(i,j,k,9) = phi(i,jx-2*margin+j,k)
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  !-- MPI send receive --!
!!$  call mpi_sendrecv &
!!$       & (bufsnd_y, n_element_buf, mpi_double_precision, rank_send, 0 &
!!$       & ,bufrcv_y, n_element_buf, mpi_double_precision, rank_recv, 0 &
!!$       & ,mpi_comm, mstatus, mpi_err)
!!$
!!$  !-- store receive buffer data into arrays for variables --!
!!$  if (rank_recv /= mpi_proc_null) then
!!$     do k = 1, kx
!!$        do j = 1, margin
!!$           do i = 1, ix
!!$              ro(i,j,k) = bufrcv_y(i,j,k,1)
!!$              pr(i,j,k) = bufrcv_y(i,j,k,2)
!!$
!!$              vx(i,j,k) = bufrcv_y(i,j,k,3)
!!$              vy(i,j,k) = bufrcv_y(i,j,k,4)
!!$              vz(i,j,k) = bufrcv_y(i,j,k,5)
!!$
!!$              bx(i,j,k) = bufrcv_y(i,j,k,6)
!!$              by(i,j,k) = bufrcv_y(i,j,k,7)
!!$              bz(i,j,k) = bufrcv_y(i,j,k,8)
!!$
!!$              phi(i,j,k) = bufrcv_y(i,j,k,9)
!!$           end do
!!$        end do
!!$     end do
!!$  end if

!----------------------------------------------------------------------|
!-- Step 1-5: data exchange in downward surface in z-coordinate 
!----------------------------------------------------------------------|

  !-- set number of elements of buffer and ranks of send/recv --!
  n_element_buf = ix * jx * margin * n_element_var
  rank_send = rank%nnd
  rank_recv = rank%nnu

  !-- set send buffer --!
  do k = 1, margin
     do j = 1, jx
        do i = 1, ix
           bufsnd_z(i,j,k,1) = ro(i,j,margin+k)
           bufsnd_z(i,j,k,2) = pr(i,j,margin+k)

           bufsnd_z(i,j,k,3) = vx(i,j,margin+k)
           bufsnd_z(i,j,k,4) = vy(i,j,margin+k)
           bufsnd_z(i,j,k,5) = vz(i,j,margin+k)

           bufsnd_z(i,j,k,6) = bx(i,j,margin+k)
           bufsnd_z(i,j,k,7) = by(i,j,margin+k)
           bufsnd_z(i,j,k,8) = bz(i,j,margin+k)

           bufsnd_z(i,j,k,9) = phi(i,j,margin+k)
        end do
     end do
  end do

  !-- MPI send receive --!
  call mpi_sendrecv &
       & (bufsnd_z, n_element_buf, mpi_double_precision, rank_send, 0 &
       & ,bufrcv_z, n_element_buf, mpi_double_precision, rank_recv, 0 &
       & ,mpi_comm, mstatus, mpi_err)

  !-- store receive buffer data into arrays for variables --!
  if(rank_recv /= mpi_proc_null)then
     do k = 1, margin
        do j = 1, jx
           do i = 1, ix
              ro(i,j,kx-margin+k) = bufrcv_z(i,j,k,1)
              pr(i,j,kx-margin+k) = bufrcv_z(i,j,k,2)

              vx(i,j,kx-margin+k) = bufrcv_z(i,j,k,3)
              vy(i,j,kx-margin+k) = bufrcv_z(i,j,k,4)
              vz(i,j,kx-margin+k) = bufrcv_z(i,j,k,5)

              bx(i,j,kx-margin+k) = bufrcv_z(i,j,k,6)
              by(i,j,kx-margin+k) = bufrcv_z(i,j,k,7)
              bz(i,j,kx-margin+k) = bufrcv_z(i,j,k,8)

              phi(i,j,kx-margin+k) = bufrcv_z(i,j,k,9)
           end do
        end do
     end do
  end if

!----------------------------------------------------------------------|
!-- Step 1-6: data exchange in upward surface in z-coordinate 
!----------------------------------------------------------------------|

  !-- set number of elements of buffer and ranks of send/recv --!
  n_element_buf = ix * jx * margin * n_element_var
  rank_send = rank%nnu
  rank_recv = rank%nnd

  !-- set send buffer --!
  do k = 1, margin
     do j = 1, jx
        do i = 1, ix
           bufsnd_z(i,j,k,1) = ro(i,j,kx-2*margin+k)
           bufsnd_z(i,j,k,2) = pr(i,j,kx-2*margin+k)

           bufsnd_z(i,j,k,3) = vx(i,j,kx-2*margin+k)
           bufsnd_z(i,j,k,4) = vy(i,j,kx-2*margin+k)
           bufsnd_z(i,j,k,5) = vz(i,j,kx-2*margin+k)

           bufsnd_z(i,j,k,6) = bx(i,j,kx-2*margin+k)
           bufsnd_z(i,j,k,7) = by(i,j,kx-2*margin+k)
           bufsnd_z(i,j,k,8) = bz(i,j,kx-2*margin+k)

           bufsnd_z(i,j,k,9) = phi(i,j,kx-2*margin+k)
        end do
     end do
  end do

  !-- MPI send receive --!
  call mpi_sendrecv &
       & (bufsnd_z, n_element_buf, mpi_double_precision, rank_send, 0 &
       & ,bufrcv_z, n_element_buf, mpi_double_precision, rank_recv, 0 &
       & ,mpi_comm, mstatus, mpi_err)

  !-- store receive buffer data into arrays for variables --!
  if (rank_recv /= mpi_proc_null) then
     do k = 1, margin
        do j = 1, jx
           do i = 1, ix
              ro(i,j,k) = bufrcv_z(i,j,k,1)
              pr(i,j,k) = bufrcv_z(i,j,k,2)

              vx(i,j,k) = bufrcv_z(i,j,k,3)
              vy(i,j,k) = bufrcv_z(i,j,k,4)
              vz(i,j,k) = bufrcv_z(i,j,k,5)

              bx(i,j,k) = bufrcv_z(i,j,k,6)
              by(i,j,k) = bufrcv_z(i,j,k,7)
              bz(i,j,k) = bufrcv_z(i,j,k,8)

              phi(i,j,k) = bufrcv_z(i,j,k,9)
           end do
        end do
     end do
  end if
  return
end subroutine cans_mpi_exchange_cart3d
