subroutine exchangeMpi(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bxi,byi,bzi &
     ,phi,merr)
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
!     8, 9,10  28,29,30  48,49,50
!     4, 5, 6  24,25,26  44,45,46
!     0, 1, 2  20,21,22  40,41,42
!        25 => xleft:24 xright:26
!           => yleft:21 yright:29
!           => zleft: 5 zright:45  
!
!               myrank             cell boundary
!                                    send  recv
!        xleft = xleft (myrank)    xleft -> xright Step 1-1
!       xright = xright(myrank)   xright -> xleft  Step 1-2
!        yleft = yleft (myrank)    yleft -> yright Step 1-3
!       yright = yright(myrank)   yright -> yleft  Step 1-4
!        zleft = zleft (myrank)    zleft -> zright Step 1-5
!       zright = zright(myrank)   zright -> zleft  Step 1-6
!
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
!        lu = yright(zleft(myrank))   lu -> td  Step 2-8
!        ld = yleft (zleft(myrank))   ld -> tu  Step 2-6
!        ll = xleft (zleft(myrank))   ll -> tr  Step 2-10
!        lr = xright(zleft(myrank))   lr -> tl  Step 2-12
!
!       mul = xleft (yright(myrank)) mul -> mdr Step 2-1
!       mur = xright(yright(myrank)) mur -> mdl Step 2-3
!       mdl = xleft (yleft (myrank)) mdl -> mur Step 2-4
!       mdr = xright(yleft (myrank)) mdr -> mul Step 2-2
!
!        tu = yright(zright(myrank))  tu -> ld  Step 2-5
!        td = yleft (zright(myrank))  td -> lu  Step 2-7
!        tl = xleft (zright(myrank))  tl -> lr  Step 2-11
!        tr = xright(zright(myrank))  tr -> ll  Step 2-9
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
!        lul = xleft (yright(zleft(myrank)))     lul -> tdr Step 3-8
!        lur = xright(yright(zleft(myrank)))     lur -> tdl Step 3-6
!        ldl = xleft (yleft (zleft(myrank)))     ldl -> tur Step 3-2
!        ldr = xright(yleft (zleft(myrank)))     ldr -> tul Step 3-4
!
!        tul = xleft (yright(zright(myrank)))    tul -> ldr Step 3-3
!        tur = xright(yright(zright(myrank)))    tur -> ldl Step 3-1
!        tdl = xleft (yleft (zright(myrank)))    tdl -> lur Step 3-5
!        tdr = xright(yleft (zright(myrank)))    tdr -> lul Step 3-7
!=====================================================================

  use mpi_domain
  implicit none

  include 'mpif.h'

  integer,intent(in) :: ix,jx,kx,margin

  type(mpidomain) :: mpid

  integer,dimension(mpi_status_size) :: mstatus

  real(8),dimension(ix,jx,kx) :: ro,pr
  real(8),dimension(ix,jx,kx) :: vx,vy,vz
  real(8),dimension(ix,jx,kx) :: bxi,byi,bzi,phi

  integer,parameter :: mx = 9

! surface exchange
  real(8),dimension(margin,jx-2*margin,kx-2*margin,mx) :: bufsnd_x,bufrcv_x
  real(8),dimension(ix-2*margin,margin,kx-2*margin,mx) :: bufsnd_y,bufrcv_y
  real(8),dimension(ix-2*margin,jx-2*margin,margin,mx) :: bufsnd_z,bufrcv_z

  integer :: mmx,msend,mrecv

  integer :: merr

  integer :: i,j,k
!=== Step 1.==========================================================
! surface exchange
!=====================================================================

!--- Step 1-1.---
! x-direction
!----------------
! from left

  mmx = margin*(jx-2*margin)*(kx-2*margin)*mx

  msend = mpid%xleft
  mrecv = mpid%xright

  do k=1,kx-2*margin
     do j=1,jx-2*margin
        do i=1,margin
           bufsnd_x(i,j,k,1) = ro(margin+i,j+margin,k+margin)
           bufsnd_x(i,j,k,2) = pr(margin+i,j+margin,k+margin)

           bufsnd_x(i,j,k,3) = vx(margin+i,j+margin,k+margin)
           bufsnd_x(i,j,k,4) = vy(margin+i,j+margin,k+margin)
           bufsnd_x(i,j,k,5) = vz(margin+i,j+margin,k+margin)

           bufsnd_x(i,j,k,6) = bxi(margin+i,j+margin,k+margin)
           bufsnd_x(i,j,k,7) = byi(margin+i,j+margin,k+margin)
           bufsnd_x(i,j,k,8) = bzi(margin+i,j+margin,k+margin)

           bufsnd_x(i,j,k,9) = phi(margin+i,j+margin,k+margin)
        enddo
     enddo
  enddo

  call mpi_sendrecv &
       (bufsnd_x,mmx,mpi_double_precision,msend, 0 &
       ,bufrcv_x,mmx,mpi_double_precision,mrecv,0 &
       ,mpi_comm_world,mstatus,merr)

  if (mrecv .ne. mpi_proc_null) then
     do k=1,kx-2*margin
        do j=1,jx-2*margin
           do i=1,margin
              ro(ix-margin+i,j+margin,k+margin) = bufrcv_x(i,j,k,1)
              pr(ix-margin+i,j+margin,k+margin) = bufrcv_x(i,j,k,2)

              vx(ix-margin+i,j+margin,k+margin) = bufrcv_x(i,j,k,3)
              vy(ix-margin+i,j+margin,k+margin) = bufrcv_x(i,j,k,4)
              vz(ix-margin+i,j+margin,k+margin) = bufrcv_x(i,j,k,5)

              bxi(ix-margin+i,j+margin,k+margin) = bufrcv_x(i,j,k,6)
              byi(ix-margin+i,j+margin,k+margin) = bufrcv_x(i,j,k,7)
              bzi(ix-margin+i,j+margin,k+margin) = bufrcv_x(i,j,k,8)

              phi(ix-margin+i,j+margin,k+margin) = bufrcv_x(i,j,k,9)
           enddo
        enddo
     enddo
  endif

!--- Step 1-2.---
! from right
!

  msend = mpid%xright
  mrecv = mpid%xleft

  do k=1,kx-2*margin
     do j=1,jx-2*margin
        do i=1,margin
           bufsnd_x(i,j,k,1) = ro(ix-2*margin+i,j+margin,k+margin)
           bufsnd_x(i,j,k,2) = pr(ix-2*margin+i,j+margin,k+margin)

           bufsnd_x(i,j,k,3) = vx(ix-2*margin+i,j+margin,k+margin)
           bufsnd_x(i,j,k,4) = vy(ix-2*margin+i,j+margin,k+margin)
           bufsnd_x(i,j,k,5) = vz(ix-2*margin+i,j+margin,k+margin)

           bufsnd_x(i,j,k,6) = bxi(ix-2*margin+i,j+margin,k+margin)
           bufsnd_x(i,j,k,7) = byi(ix-2*margin+i,j+margin,k+margin)
           bufsnd_x(i,j,k,8) = bzi(ix-2*margin+i,j+margin,k+margin)

           bufsnd_x(i,j,k,9) = phi(ix-2*margin+i,j+margin,k+margin)
        enddo
     enddo
  enddo

  call mpi_sendrecv &
       (bufsnd_x,mmx,mpi_double_precision,msend, 1 &
       ,bufrcv_x,mmx,mpi_double_precision,mrecv, 1 &
       ,mpi_comm_world,mstatus,merr)

  if (mrecv .ne. mpi_proc_null)then
     do k=1,kx-2*margin
        do j=1,jx-2*margin
           do i=1,margin
              ro(i,j+margin,k+margin) = bufrcv_x(i,j,k,1)
              pr(i,j+margin,k+margin) = bufrcv_x(i,j,k,2)

              vx(i,j+margin,k+margin) = bufrcv_x(i,j,k,3)
              vy(i,j+margin,k+margin) = bufrcv_x(i,j,k,4)
              vz(i,j+margin,k+margin) = bufrcv_x(i,j,k,5)

              bxi(i,j+margin,k+margin) = bufrcv_x(i,j,k,6)
              byi(i,j+margin,k+margin) = bufrcv_x(i,j,k,7)
              bzi(i,j+margin,k+margin) = bufrcv_x(i,j,k,8)

              phi(i,j+margin,k+margin) = bufrcv_x(i,j,k,9)
           enddo
        enddo
     enddo
  endif

!--- Step 1-3.---
! y-direction
!----------------
! from left

  mmx = (ix-2*margin)*margin*(kx-2*margin)*mx

  msend = mpid%yleft
  mrecv = mpid%yright

  do k=1,kx-2*margin
     do j=1,margin
        do i=1,ix-2*margin
           bufsnd_y(i,j,k,1) = ro(i+margin,margin+j,k+margin)
           bufsnd_y(i,j,k,2) = pr(i+margin,margin+j,k+margin)

           bufsnd_y(i,j,k,3) = vx(i+margin,margin+j,k+margin)
           bufsnd_y(i,j,k,4) = vy(i+margin,margin+j,k+margin)
           bufsnd_y(i,j,k,5) = vz(i+margin,margin+j,k+margin)

           bufsnd_y(i,j,k,6) = bxi(i+margin,margin+j,k+margin)
           bufsnd_y(i,j,k,7) = byi(i+margin,margin+j,k+margin)
           bufsnd_y(i,j,k,8) = bzi(i+margin,margin+j,k+margin)

           bufsnd_y(i,j,k,9) = phi(i+margin,margin+j,k+margin)
        enddo
     enddo
  enddo

  call mpi_sendrecv &
       (bufsnd_y,mmx,mpi_double_precision,msend, 0 &
       ,bufrcv_y,mmx,mpi_double_precision,mrecv, 0 &
       ,mpi_comm_world,mstatus,merr)

  if (mrecv .ne. mpi_proc_null) then
     do k=1,kx-2*margin
        do j=1,margin
           do i=1,ix-2*margin
              ro(i+margin,jx-margin+j,k+margin) = bufrcv_y(i,j,k,1)
              pr(i+margin,jx-margin+j,k+margin) = bufrcv_y(i,j,k,2)

              vx(i+margin,jx-margin+j,k+margin) = bufrcv_y(i,j,k,3)
              vy(i+margin,jx-margin+j,k+margin) = bufrcv_y(i,j,k,4)
              vz(i+margin,jx-margin+j,k+margin) = bufrcv_y(i,j,k,5)

              bxi(i+margin,jx-margin+j,k+margin) = bufrcv_y(i,j,k,6)
              byi(i+margin,jx-margin+j,k+margin) = bufrcv_y(i,j,k,7)
              bzi(i+margin,jx-margin+j,k+margin) = bufrcv_y(i,j,k,8)

              phi(i+margin,jx-margin+j,k+margin) = bufrcv_y(i,j,k,9)
           enddo
        enddo
     enddo
  endif

!--- Step 1-4.---
! from right
!

  msend = mpid%yright
  mrecv = mpid%yleft

  do k=1,kx-2*margin
     do j=1,margin
        do i=1,ix-2*margin
           bufsnd_y(i,j,k,1) = ro(i+margin,jx-2*margin+j,k+margin)
           bufsnd_y(i,j,k,2) = pr(i+margin,jx-2*margin+j,k+margin)

           bufsnd_y(i,j,k,3) = vx(i+margin,jx-2*margin+j,k+margin)
           bufsnd_y(i,j,k,4) = vy(i+margin,jx-2*margin+j,k+margin)
           bufsnd_y(i,j,k,5) = vz(i+margin,jx-2*margin+j,k+margin)

           bufsnd_y(i,j,k,6) = bxi(i+margin,jx-2*margin+j,k+margin)
           bufsnd_y(i,j,k,7) = byi(i+margin,jx-2*margin+j,k+margin)
           bufsnd_y(i,j,k,8) = bzi(i+margin,jx-2*margin+j,k+margin)

           bufsnd_y(i,j,k,9) = phi(i+margin,jx-2*margin+j,k+margin)
        enddo
     enddo
  enddo

  call mpi_sendrecv &
       (bufsnd_y,mmx,mpi_double_precision,msend, 1 &
       ,bufrcv_y,mmx,mpi_double_precision,mrecv, 1 &
       ,mpi_comm_world,mstatus,merr)

  if (mrecv .ne. mpi_proc_null)then
     do k=1,kx-2*margin
        do j=1,margin
           do i=1,ix-2*margin
              ro(i+margin,j,k+margin) = bufrcv_y(i,j,k,1)
              pr(i+margin,j,k+margin) = bufrcv_y(i,j,k,2)

              vx(i+margin,j,k+margin) = bufrcv_y(i,j,k,3)
              vy(i+margin,j,k+margin) = bufrcv_y(i,j,k,4)
              vz(i+margin,j,k+margin) = bufrcv_y(i,j,k,5)

              bxi(i+margin,j,k+margin) = bufrcv_y(i,j,k,6)
              byi(i+margin,j,k+margin) = bufrcv_y(i,j,k,7)
              bzi(i+margin,j,k+margin) = bufrcv_y(i,j,k,8)

              phi(i+margin,j,k+margin) = bufrcv_y(i,j,k,9)
           enddo
        enddo
     enddo
  endif

!--- Step 1-5.---
! z-direction
!----------------
! from left

  mmx = (ix-2*margin)*(jx-2*margin)*margin*mx

  msend = mpid%zleft
  mrecv = mpid%zright

  do k=1,margin
     do j=1,jx-2*margin
        do i=1,ix-2*margin
           bufsnd_z(i,j,k,1) = ro(i+margin,j+margin,margin+k)
           bufsnd_z(i,j,k,2) = pr(i+margin,j+margin,margin+k)

           bufsnd_z(i,j,k,3) = vx(i+margin,j+margin,margin+k)
           bufsnd_z(i,j,k,4) = vy(i+margin,j+margin,margin+k)
           bufsnd_z(i,j,k,5) = vz(i+margin,j+margin,margin+k)

           bufsnd_z(i,j,k,6) = bxi(i+margin,j+margin,margin+k)
           bufsnd_z(i,j,k,7) = byi(i+margin,j+margin,margin+k)
           bufsnd_z(i,j,k,8) = bzi(i+margin,j+margin,margin+k)

           bufsnd_z(i,j,k,9) = phi(i+margin,j+margin,margin+k)
        enddo
     enddo
  enddo

  call mpi_sendrecv &
       (bufsnd_z,mmx,mpi_double_precision,msend, 0 &
       ,bufrcv_z,mmx,mpi_double_precision,mrecv, 0 &
       ,mpi_comm_world,mstatus,merr)

  if (mrecv .ne. mpi_proc_null) then
     do k=1,margin
        do j=1,jx-2*margin
           do i=1,ix-2*margin
              ro(i+margin,j+margin,kx-margin+k) = bufrcv_z(i,j,k,1)
              pr(i+margin,j+margin,kx-margin+k) = bufrcv_z(i,j,k,2)

              vx(i+margin,j+margin,kx-margin+k) = bufrcv_z(i,j,k,3)
              vy(i+margin,j+margin,kx-margin+k) = bufrcv_z(i,j,k,4)
              vz(i+margin,j+margin,kx-margin+k) = bufrcv_z(i,j,k,5)

              bxi(i+margin,j+margin,kx-margin+k) = bufrcv_z(i,j,k,6)
              byi(i+margin,j+margin,kx-margin+k) = bufrcv_z(i,j,k,7)
              bzi(i+margin,j+margin,kx-margin+k) = bufrcv_z(i,j,k,8)

              phi(i+margin,j+margin,kx-margin+k) = bufrcv_z(i,j,k,9)
           enddo
        enddo
     enddo
  endif

!--- Step 1-6.---
! from right
!

  msend = mpid%zright
  mrecv = mpid%zleft

  do k=1,margin
     do j=1,jx-2*margin
        do i=1,ix-2*margin
           bufsnd_z(i,j,k,1) = ro(i+margin,j+margin,kx-2*margin+k)
           bufsnd_z(i,j,k,2) = pr(i+margin,j+margin,kx-2*margin+k)

           bufsnd_z(i,j,k,3) = vx(i+margin,j+margin,kx-2*margin+k)
           bufsnd_z(i,j,k,4) = vy(i+margin,j+margin,kx-2*margin+k)
           bufsnd_z(i,j,k,5) = vz(i+margin,j+margin,kx-2*margin+k)

           bufsnd_z(i,j,k,6) = bxi(i+margin,j+margin,kx-2*margin+k)
           bufsnd_z(i,j,k,7) = byi(i+margin,j+margin,kx-2*margin+k)
           bufsnd_z(i,j,k,8) = bzi(i+margin,j+margin,kx-2*margin+k)

           bufsnd_z(i,j,k,9) = phi(i+margin,j+margin,kx-2*margin+k)
        enddo
     enddo
  enddo

  call mpi_sendrecv &
       (bufsnd_z,mmx,mpi_double_precision,msend, 1 &
       ,bufrcv_z,mmx,mpi_double_precision,mrecv, 1 &
       ,mpi_comm_world,mstatus,merr)

  if (mrecv .ne. mpi_proc_null)then
     do k=1,margin
        do j=1,jx-2*margin
           do i=1,ix-2*margin
              ro(i+margin,j+margin,k) = bufrcv_z(i,j,k,1)
              pr(i+margin,j+margin,k) = bufrcv_z(i,j,k,2)

              vx(i+margin,j+margin,k) = bufrcv_z(i,j,k,3)
              vy(i+margin,j+margin,k) = bufrcv_z(i,j,k,4)
              vz(i+margin,j+margin,k) = bufrcv_z(i,j,k,5)

              bxi(i+margin,j+margin,k) = bufrcv_z(i,j,k,6)
              byi(i+margin,j+margin,k) = bufrcv_z(i,j,k,7)
              bzi(i+margin,j+margin,k) = bufrcv_z(i,j,k,8)

              phi(i+margin,j+margin,k) = bufrcv_z(i,j,k,9)
           enddo
        enddo
     enddo
  endif

  return
end subroutine exchangeMpi
