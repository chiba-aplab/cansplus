subroutine exchangeMpixz_res(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz &
     ,phi,eta,merr)
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
  use mpi_domain_xz
  implicit none

  include 'mpif.h'

  integer,intent(in) :: ix,jx,kx,margin

  integer,dimension(mpi_status_size) :: mstatus

  real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz,phi,eta

  integer,parameter :: mx = 10

! surface exchange
  real(8),dimension(margin,jx,kx,mx) :: bufsnd_x,bufrcv_x
  real(8),dimension(ix,jx,margin,mx) :: bufsnd_z,bufrcv_z

! line exchange
  real(8),dimension(margin,jx,margin,mx) :: bufsnd_xz,bufrcv_xz

  integer :: mmx,msend,mrecv

  integer :: merr

  integer :: i,j,k

!=== Step 1.==========================================================
! surface exchange
!=====================================================================

!--- Step 1-1.---
! x-direction
!----------------
! left

  mmx = margin*jx*kx*mx

  msend = mpid%l
  mrecv = mpid%r

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
           bufsnd_x(i,j,k,10) = eta(margin+i,j,k)
        enddo
     end do
  end do

  call mpi_sendrecv &
       (bufsnd_x,mmx,mpi_double_precision,msend, 0 &
       ,bufrcv_x,mmx,mpi_double_precision,mrecv, 0 &
       ,mpi_comm_world,mstatus,merr)

  if (mrecv .ne. mpi_proc_null) then
     do k=1,kx
        do j=1,jx
           do i=1,margin
              ro(ix-margin+i,j,k) = bufrcv_x(i,j,k,1)
              pr(ix-margin+i,j,k) = bufrcv_x(i,j,k,2)

              vx(ix-margin+i,j,k) = bufrcv_x(i,j,k,3)
              vy(ix-margin+i,j,k) = bufrcv_x(i,j,k,4)
              vz(ix-margin+i,j,k) = bufrcv_x(i,j,k,5)

              bx(ix-margin+i,j,k) = bufrcv_x(i,j,k,6)
              by(ix-margin+i,j,k) = bufrcv_x(i,j,k,7)
              bz(ix-margin+i,j,k) = bufrcv_x(i,j,k,8)

              phi(ix-margin+i,j,k) = bufrcv_x(i,j,k,9)
              eta(ix-margin+i,j,k) = bufrcv_x(i,j,k,10)
           enddo
        enddo
     enddo
  endif

!--- Step 1-2.---
! right
!

  mmx = margin*jx*kx*mx

  msend = mpid%r
  mrecv = mpid%l

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
           bufsnd_x(i,j,k,10) = eta(ix-2*margin+i,j,k)
        end do
     end do
  end do

  call mpi_sendrecv &
       (bufsnd_x,mmx,mpi_double_precision,msend, 0 &
       ,bufrcv_x,mmx,mpi_double_precision,mrecv, 0 &
       ,mpi_comm_world,mstatus,merr)

  if (mrecv .ne. mpi_proc_null)then
     do k=1,kx
        do j=1,jx
           do i=1,margin
              ro(i,j,k) = bufrcv_x(i,j,k,1)
              pr(i,j,k) = bufrcv_x(i,j,k,2)

              vx(i,j,k) = bufrcv_x(i,j,k,3)
              vy(i,j,k) = bufrcv_x(i,j,k,4)
              vz(i,j,k) = bufrcv_x(i,j,k,5)

              bx(i,j,k) = bufrcv_x(i,j,k,6)
              by(i,j,k) = bufrcv_x(i,j,k,7)
              bz(i,j,k) = bufrcv_x(i,j,k,8)

              phi(i,j,k) = bufrcv_x(i,j,k,9)
              eta(i,j,k) = bufrcv_x(i,j,k,10)
           end do
        end do
     end do
  end if

!--- Step 1-3.---
! z-direction
!----------------
! down
              
  mmx = ix*jx*margin*mx

  msend = mpid%d
  mrecv = mpid%t

  do k=1,margin
     do j=1,jx
        do i=1,ix
           bufsnd_z(i,j,k,1) = ro(i,j,margin+k)
           bufsnd_z(i,j,k,2) = pr(i,j,margin+k)

           bufsnd_z(i,j,k,3) = vx(i,j,margin+k)
           bufsnd_z(i,j,k,4) = vy(i,j,margin+k)
           bufsnd_z(i,j,k,5) = vz(i,j,margin+k)

           bufsnd_z(i,j,k,6) = bx(i,j,margin+k)
           bufsnd_z(i,j,k,7) = by(i,j,margin+k)
           bufsnd_z(i,j,k,8) = bz(i,j,margin+k)

           bufsnd_z(i,j,k,9) = phi(i,j,margin+k)
           bufsnd_z(i,j,k,10) = eta(i,j,margin+k)
        end do
     end do
  end do

  call mpi_sendrecv &
       (bufsnd_z,mmx,mpi_double_precision,msend, 0 &
       ,bufrcv_z,mmx,mpi_double_precision,mrecv, 0 &
       ,mpi_comm_world,mstatus,merr)

  if(mrecv .ne. mpi_proc_null)then
     do k=1,margin
        do j=1,jx
           do i=1,ix
              ro(i,j,kx-margin+k) = bufrcv_z(i,j,k,1)
              pr(i,j,kx-margin+k) = bufrcv_z(i,j,k,2)

              vx(i,j,kx-margin+k) = bufrcv_z(i,j,k,3)
              vy(i,j,kx-margin+k) = bufrcv_z(i,j,k,4)
              vz(i,j,kx-margin+k) = bufrcv_z(i,j,k,5)

              bx(i,j,kx-margin+k) = bufrcv_z(i,j,k,6)
              by(i,j,kx-margin+k) = bufrcv_z(i,j,k,7)
              bz(i,j,kx-margin+k) = bufrcv_z(i,j,k,8)

              phi(i,j,kx-margin+k) = bufrcv_z(i,j,k,9)
              eta(i,j,kx-margin+k) = bufrcv_z(i,j,k,10)
           end do
        end do
     end do
  end if

!--- Step 1-4.---
! z-direction
!----------------
! top
!

  mmx = ix*jx*margin*mx

  msend = mpid%t
  mrecv = mpid%d

  do k=1,margin
     do j=1,jx
        do i=1,ix
           bufsnd_z(i,j,k,1) = ro(i,j,kx-2*margin+k)
           bufsnd_z(i,j,k,2) = pr(i,j,kx-2*margin+k)

           bufsnd_z(i,j,k,3) = vx(i,j,kx-2*margin+k)
           bufsnd_z(i,j,k,4) = vy(i,j,kx-2*margin+k)
           bufsnd_z(i,j,k,5) = vz(i,j,kx-2*margin+k)

           bufsnd_z(i,j,k,6) = bx(i,j,kx-2*margin+k)
           bufsnd_z(i,j,k,7) = by(i,j,kx-2*margin+k)
           bufsnd_z(i,j,k,8) = bz(i,j,kx-2*margin+k)

           bufsnd_z(i,j,k,9) = phi(i,j,kx-2*margin+k)
           bufsnd_z(i,j,k,10) = eta(i,j,kx-2*margin+k)
        end do
     end do
  end do

  call mpi_sendrecv &
       (bufsnd_z,mmx,mpi_double_precision,msend, 0 &
       ,bufrcv_z,mmx,mpi_double_precision,mrecv, 0 &
       ,mpi_comm_world,mstatus,merr)

  if (mrecv .ne. mpi_proc_null)then
     do k=1,margin
        do j=1,jx
           do i=1,ix
              ro(i,j,k) = bufrcv_z(i,j,k,1)
              pr(i,j,k) = bufrcv_z(i,j,k,2)

              vx(i,j,k) = bufrcv_z(i,j,k,3)
              vy(i,j,k) = bufrcv_z(i,j,k,4)
              vz(i,j,k) = bufrcv_z(i,j,k,5)

              bx(i,j,k) = bufrcv_z(i,j,k,6)
              by(i,j,k) = bufrcv_z(i,j,k,7)
              bz(i,j,k) = bufrcv_z(i,j,k,8)

              phi(i,j,k) = bufrcv_z(i,j,k,9)
              eta(i,j,k) = bufrcv_z(i,j,k,10)
           end do
        end do
     end do
  end if
  return
end subroutine exchangeMpixz_res
