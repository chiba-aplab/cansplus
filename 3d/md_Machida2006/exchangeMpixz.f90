subroutine exchangeMpixz(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz &
                        ,phi)

  use mpi_setup
  implicit none

  integer,intent(in) :: ix,jx,kx,margin
  real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz,phi

  integer,parameter :: mx = 9
  integer :: i,j,k,mmx,msend,mrecv
  real(8),dimension(margin,jx,kx,mx) :: bufsnd_x,bufrcv_x
  real(8),dimension(ix,margin,kx,mx) :: bufsnd_y,bufrcv_y
  real(8),dimension(ix,jx,margin,mx) :: bufsnd_z,bufrcv_z

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
        enddo
     enddo
  enddo

  call mpi_sendrecv               &
       (bufsnd_x,mmx,mdp,msend, 0 &
       ,bufrcv_x,mmx,mdp,mrecv, 0 &
       ,mcomw,mstat,merr)

  if(mrecv /= mnull)then
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
           enddo
        enddo
     enddo
  endif

!--- Step 1-2.---
! right

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
        enddo
     enddo
  enddo

  call mpi_sendrecv               &
       (bufsnd_x,mmx,mdp,msend, 0 &
       ,bufrcv_x,mmx,mdp,mrecv, 0 &
       ,mcomw,mstat,merr)

  if(mrecv /= mnull)then
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
           enddo
        enddo
     enddo
  end if

!--- Step 1-3.---
! y-direction
!----------------
! back
              
  mmx = ix*kx*margin*mx
  msend = mpid%b
  mrecv = mpid%f

  do k=1,kx
     do j=1,margin
        do i=1,ix
           bufsnd_y(i,j,k,1) = ro(i,margin+j,k)
           bufsnd_y(i,j,k,2) = pr(i,margin+j,k)
           bufsnd_y(i,j,k,3) = vx(i,margin+j,k)
           bufsnd_y(i,j,k,4) = vy(i,margin+j,k)
           bufsnd_y(i,j,k,5) = vz(i,margin+j,k)
           bufsnd_y(i,j,k,6) = bx(i,margin+j,k)
           bufsnd_y(i,j,k,7) = by(i,margin+j,k)
           bufsnd_y(i,j,k,8) = bz(i,margin+j,k)
           bufsnd_y(i,j,k,9) = phi(i,margin+j,k)
        enddo
     enddo
  enddo

  call mpi_sendrecv               &
       (bufsnd_y,mmx,mdp,msend, 0 &
       ,bufrcv_y,mmx,mdp,mrecv, 0 &
       ,mcomw,mstat,merr)

  if(mrecv /= mnull)then
     do k=1,kx
        do j=1,margin
           do i=1,ix
              ro(i,jx-margin+j,k) = bufrcv_y(i,j,k,1)
              pr(i,jx-margin+j,k) = bufrcv_y(i,j,k,2)
              vx(i,jx-margin+j,k) = bufrcv_y(i,j,k,3)
              vy(i,jx-margin+j,k) = bufrcv_y(i,j,k,4)
              vz(i,jx-margin+j,k) = bufrcv_y(i,j,k,5)
              bx(i,jx-margin+j,k) = bufrcv_y(i,j,k,6)
              by(i,jx-margin+j,k) = bufrcv_y(i,j,k,7)
              bz(i,jx-margin+j,k) = bufrcv_y(i,j,k,8)
              phi(i,jx-margin+j,k) = bufrcv_y(i,j,k,9)
           enddo
        enddo
     enddo
  end if

!--- Step 1-4.---
! y-direction
!----------------
! forth

  mmx = ix*kx*margin*mx
  msend = mpid%f
  mrecv = mpid%b

  do k=1,kx
     do j=1,margin
        do i=1,ix
           bufsnd_y(i,j,k,1) = ro(i,jx-2*margin+j,k)
           bufsnd_y(i,j,k,2) = pr(i,jx-2*margin+j,k)
           bufsnd_y(i,j,k,3) = vx(i,jx-2*margin+j,k)
           bufsnd_y(i,j,k,4) = vy(i,jx-2*margin+j,k)
           bufsnd_y(i,j,k,5) = vz(i,jx-2*margin+j,k)
           bufsnd_y(i,j,k,6) = bx(i,jx-2*margin+j,k)
           bufsnd_y(i,j,k,7) = by(i,jx-2*margin+j,k)
           bufsnd_y(i,j,k,8) = bz(i,jx-2*margin+j,k)
           bufsnd_y(i,j,k,9) = phi(i,jx-2*margin+j,k)
        enddo
     enddo
  enddo

  call mpi_sendrecv               &
       (bufsnd_y,mmx,mdp,msend, 0 &
       ,bufrcv_y,mmx,mdp,mrecv, 0 &
       ,mcomw,mstat,merr)

  if(mrecv /= mnull)then
     do k=1,kx
        do j=1,margin
           do i=1,ix
              ro(i,j,k) = bufrcv_y(i,j,k,1)
              pr(i,j,k) = bufrcv_y(i,j,k,2)
              vx(i,j,k) = bufrcv_y(i,j,k,3)
              vy(i,j,k) = bufrcv_y(i,j,k,4)
              vz(i,j,k) = bufrcv_y(i,j,k,5)
              bx(i,j,k) = bufrcv_y(i,j,k,6)
              by(i,j,k) = bufrcv_y(i,j,k,7)
              bz(i,j,k) = bufrcv_y(i,j,k,8)
              phi(i,j,k) = bufrcv_y(i,j,k,9)
           enddo
        enddo
     enddo
  end if

!--- Step 1-5.---
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
        enddo
     enddo
  enddo

  call mpi_sendrecv               &
       (bufsnd_z,mmx,mdp,msend, 0 &
       ,bufrcv_z,mmx,mdp,mrecv, 0 &
       ,mcomw,mstat,merr)

  if(mrecv /= mnull)then
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
           enddo
        enddo
     enddo
  end if

!--- Step 1-6.---
! z-direction
!----------------
! top

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
        enddo
     enddo
  enddo

  call mpi_sendrecv               &
       (bufsnd_z,mmx,mdp,msend, 0 &
       ,bufrcv_z,mmx,mdp,mrecv, 0 &
       ,mcomw,mstat,merr)

  if(mrecv /= mnull)then
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
           enddo
        enddo
     enddo
  end if

end subroutine exchangeMpixz
