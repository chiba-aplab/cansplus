module boundary

  implicit none
  private

  public :: bd_consx, bd_consy, bd_frex, bd_frey, &
            bd_inix, bd_iniy,                     &
            bd_synnx_car, bd_synny_car,           &
            bd_synpx_car, bd_synpy_car,           &
            boundary__mpi


contains


!=====================================================================
! INPUT
!      mbnd :: margin flag
!      ix,jx,kx :: array size
!      margin :: margin
! INPUT&OUTPUT
!      qq
!=====================================================================

subroutine bd_consx(mbnd,margin,cons,qq,ix,jx)
  implicit none

  integer,intent(in) :: ix,jx
  integer,intent(in) :: margin
  integer,intent(in) :: mbnd
  real(8),intent(in) :: cons
  real(8),dimension(ix,jx),intent(inout) :: qq

  integer :: ibnd
  integer :: i,j

  if (mbnd .eq. 0)then
     ibnd = 1+margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,jx
        do i=1,margin
           qq(ibnd-i,j) = cons
        enddo
     enddo
     !$OMP END PARALLEL DO
  else
     ibnd = ix-margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,jx
        do i=1,margin
           qq(ibnd+i,j) = cons
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

  return
end subroutine bd_consx


subroutine bd_consy(mbnd,margin,cons,qq,ix,jx)
  implicit none

  integer,intent(in) :: mbnd ! boundary flag                                       
  integer,intent(in) :: margin ! margin                                            
  integer,intent(in) :: ix,jx! array size                                      
  real(8),intent(in) :: cons
  real(8),dimension(ix,jx),intent(inout) :: qq

  integer :: jbnd
  integer :: i,j

  if( mbnd .eq. 0)then
     jbnd = 1+margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,margin
        do i=1,ix
           qq(i,jbnd-j) = cons
        enddo
     enddo
     !$OMP END PARALLEL DO
  else
     jbnd = jx-margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,margin
        do i=1,ix
           qq(i,jbnd+j) = cons
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

  return
end subroutine bd_consy


subroutine bd_frex(mbnd,margin,qq,ix,jx)
  implicit none

  integer,intent(in) :: mbnd ! boundary flag
  integer,intent(in) :: margin ! margin
  integer,intent(in) :: ix,jx ! array size
  real(8),dimension(ix,jx),intent(inout) :: qq

  integer :: i,j

  if (mbnd .eq. 0)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,jx
        do i=1,margin
           qq(i,j) = qq(margin+1,j)
        enddo
     enddo
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,jx
        do i=1,margin
           qq(ix-margin+i,j) = qq(ix-margin,j)
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

  return
end subroutine bd_frex


subroutine bd_frey(mbnd,margin,qq,ix,jx)
  implicit none

  integer,intent(in) :: mbnd ! boundary flag
  integer,intent(in) :: margin ! margin
  integer,intent(in) :: ix,jx ! array size
  real(8),dimension(ix,jx),intent(inout) :: qq

  integer :: jbnd
  integer :: i,j

  if( mbnd .eq. 0)then
     jbnd = 1+margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,margin
        do i=1,ix
           qq(i,jbnd-j) = qq(i,jbnd)
        enddo
     enddo
     !$OMP END PARALLEL DO
  else
     jbnd = jx-margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,margin
        do i=1,ix
           qq(i,jbnd+j) = qq(i,jbnd)
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

  return
end subroutine bd_frey


subroutine bd_inix(mbnd,margin,qq,qqini,ix,jx)
  implicit none

  integer,intent(in) :: mbnd ! boundary flag
  integer,intent(in) :: margin ! margin
  integer,intent(in) :: ix,jx ! array size
  real(8),dimension(ix,jx),intent(in) :: qqini ! initial quantity
  real(8),dimension(ix,jx),intent(inout) :: qq

  integer :: ibnd
  integer :: i,j

  if (mbnd .eq. 0)then
     ibnd = 1+margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,jx
        do i=1,margin
           qq(ibnd-i,j) = qqini(ibnd-i,j)
        enddo
     enddo
     !$OMP END PARALLEL DO
  else
     ibnd = ix-margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,jx
        do i=1,margin
           qq(ibnd+i,j) = qqini(ibnd+i,j)
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

  return
end subroutine bd_inix


subroutine bd_iniy(mbnd,margin,qq,qqini,ix,jx)
  implicit none

  integer,intent(in) :: mbnd ! boundary flag
  integer,intent(in) :: margin ! margin
  integer,intent(in) :: ix,jx ! array size
  real(8),dimension(ix,jx),intent(in) :: qqini ! initial quantity
  real(8),dimension(ix,jx),intent(inout) :: qq

  integer :: jbnd
  integer :: i,j,k

  if( mbnd .eq. 0)then
     jbnd = 1+margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,margin
        do i=1,ix
           qq(i,jbnd-j) = qqini(i,jbnd-j)
        enddo
     enddo
     !$OMP END PARALLEL DO
  else
     jbnd = jx-margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,margin
        do i=1,ix
           qq(i,jbnd+j) = qqini(i,jbnd+j)
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

  return
end subroutine bd_iniy


subroutine bd_synnx_car(mbnd,margin,qq,ix,jx)
  implicit none

  integer,intent(in) :: margin,ix,jx
  integer,intent(in) :: mbnd
  real(8),dimension(ix,jx),intent(inout) :: qq

  integer :: i,j

  if(mbnd .eq. 0)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,jx
        do i=1,margin
           qq(margin-i+1,j) = -(qq(margin+i,j))
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,jx
        do i=1,margin
           qq(ix-margin+i,j) = -qq(ix-margin-i+1,j)
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  return
end subroutine bd_synnx_car


subroutine bd_synny_car(mbnd,margin,qq,ix,jx)
  implicit none

  integer,intent(in) :: margin,ix,jx
  integer,intent(in) :: mbnd

  real(8),dimension(ix,jx),intent(inout) :: qq

  integer :: i,j

  if(mbnd .eq. 0)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,margin
        do i=1,ix
           qq(i,margin-j+1) = -(qq(i,margin+j))
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,margin
        do i=1,ix
           qq(i,jx-margin+j) = -qq(i,jx-margin-j+1)
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  return
end subroutine bd_synny_car


subroutine bd_synpx_car(mbnd,margin,qq,ix,jx)
  implicit none

  integer,intent(in) :: margin,ix,jx
  integer,intent(in) :: mbnd
  real(8),dimension(ix,jx),intent(inout) :: qq

  integer :: i,j

  if(mbnd .eq. 0)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,jx
        do i=1,margin
           qq(margin-i+1,j) = (qq(margin+i,j))
        enddo
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,jx
        do i=1,margin
           qq(ix-margin+i,j) = qq(ix-margin-i+1,j)
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  return
end subroutine bd_synpx_car


subroutine bd_synpy_car(mbnd,margin,qq,ix,jx)
  implicit none

  integer,intent(in) :: margin,ix,jx
  integer,intent(in) :: mbnd
  real(8),dimension(ix,jx),intent(inout) :: qq

  integer :: i,j

  if(mbnd .eq. 0)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,margin
        do i=1,ix
           qq(i,margin-j+1) = (qq(i,margin+j))
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,margin
        do i=1,ix
           qq(i,jx-margin+j) = qq(i,jx-margin-j+1)
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  return
end subroutine bd_synpy_car


subroutine boundary__mpi(margin,ix,jx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta)

  use mpi_setup

  integer,intent(in) :: ix,jx,margin
  real(8),dimension(ix,jx),intent(inout) :: ro,pr
  real(8),dimension(ix,jx),intent(inout) :: vx,vy,vz
  real(8),dimension(ix,jx),intent(inout) :: bx,by,bz,phi
  real(8),dimension(ix,jx),intent(inout) :: eta

  integer,parameter :: mx = 10
  integer :: i,j,mmx,msend,mrecv
  real(8),dimension(margin,jx,mx) :: bufsnd_x,bufrcv_x
  real(8),dimension(ix,margin,mx) :: bufsnd_y,bufrcv_y

!=== Step 1.==========================================================
! surface exchange
!=====================================================================
!--- Step 1-1.---
! x-direction
!----------------
! left

  mmx = margin*jx*mx
  msend = mpid%l
  mrecv = mpid%r

  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j)
  do j=1,jx
     do i=1,margin
        bufsnd_x(i,j,1) = ro(margin+i,j)
        bufsnd_x(i,j,2) = pr(margin+i,j)
        bufsnd_x(i,j,3) = vx(margin+i,j)
        bufsnd_x(i,j,4) = vy(margin+i,j)
        bufsnd_x(i,j,5) = vz(margin+i,j)
        bufsnd_x(i,j,6) = bx(margin+i,j)
        bufsnd_x(i,j,7) = by(margin+i,j)
        bufsnd_x(i,j,8) = bz(margin+i,j)
        bufsnd_x(i,j,9) = phi(margin+i,j)
        bufsnd_x(i,j,10) = eta(margin+i,j)
     enddo
  enddo
  !$OMP END PARALLEL DO

  call mpi_sendrecv               &
       (bufsnd_x,mmx,mdp,msend, 0 &
       ,bufrcv_x,mmx,mdp,mrecv, 0 &
       ,mcomw,mstat,merr)

  if(mrecv /= mnull)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,jx
        do i=1,margin
           ro(ix-margin+i,j) = bufrcv_x(i,j,1)
           pr(ix-margin+i,j) = bufrcv_x(i,j,2)
           vx(ix-margin+i,j) = bufrcv_x(i,j,3)
           vy(ix-margin+i,j) = bufrcv_x(i,j,4)
           vz(ix-margin+i,j) = bufrcv_x(i,j,5)
           bx(ix-margin+i,j) = bufrcv_x(i,j,6)
           by(ix-margin+i,j) = bufrcv_x(i,j,7)
           bz(ix-margin+i,j) = bufrcv_x(i,j,8)
           phi(ix-margin+i,j) = bufrcv_x(i,j,9)
           eta(ix-margin+i,j) = bufrcv_x(i,j,10)
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

!--- Step 1-2.---
! right

  mmx = margin*jx*mx
  msend = mpid%r
  mrecv = mpid%l

  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j)
  do j=1,jx
     do i=1,margin
        bufsnd_x(i,j,1) = ro(ix-2*margin+i,j)
        bufsnd_x(i,j,2) = pr(ix-2*margin+i,j)
        bufsnd_x(i,j,3) = vx(ix-2*margin+i,j)
        bufsnd_x(i,j,4) = vy(ix-2*margin+i,j)
        bufsnd_x(i,j,5) = vz(ix-2*margin+i,j)
        bufsnd_x(i,j,6) = bx(ix-2*margin+i,j)
        bufsnd_x(i,j,7) = by(ix-2*margin+i,j)
        bufsnd_x(i,j,8) = bz(ix-2*margin+i,j)
        bufsnd_x(i,j,9) = phi(ix-2*margin+i,j)
        bufsnd_x(i,j,10) = eta(ix-2*margin+i,j)
     enddo
  enddo
  !$OMP END PARALLEL DO

  call mpi_sendrecv               &
       (bufsnd_x,mmx,mdp,msend, 0 &
       ,bufrcv_x,mmx,mdp,mrecv, 0 &
       ,mcomw,mstat,merr)

  if(mrecv /= mnull)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,jx
        do i=1,margin
           ro(i,j) = bufrcv_x(i,j,1)
           pr(i,j) = bufrcv_x(i,j,2)
           vx(i,j) = bufrcv_x(i,j,3)
           vy(i,j) = bufrcv_x(i,j,4)
           vz(i,j) = bufrcv_x(i,j,5)
           bx(i,j) = bufrcv_x(i,j,6)
           by(i,j) = bufrcv_x(i,j,7)
           bz(i,j) = bufrcv_x(i,j,8)
           phi(i,j) = bufrcv_x(i,j,9)
           eta(i,j) = bufrcv_x(i,j,10)
        enddo
     enddo
     !$OMP END PARALLEL DO
  end if

!--- Step 1-3.---
! y-direction
!----------------
! back
              
  mmx = ix*margin*mx
  msend = mpid%b
  mrecv = mpid%f

  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j)
  do j=1,margin
     do i=1,ix
        bufsnd_y(i,j,1) = ro(i,margin+j)
        bufsnd_y(i,j,2) = pr(i,margin+j)
        bufsnd_y(i,j,3) = vx(i,margin+j)
        bufsnd_y(i,j,4) = vy(i,margin+j)
        bufsnd_y(i,j,5) = vz(i,margin+j)
        bufsnd_y(i,j,6) = bx(i,margin+j)
        bufsnd_y(i,j,7) = by(i,margin+j)
        bufsnd_y(i,j,8) = bz(i,margin+j)
        bufsnd_y(i,j,9) = phi(i,margin+j)
        bufsnd_y(i,j,10) = eta(i,margin+j)
     enddo
  enddo
  !$OMP END PARALLEL DO

  call mpi_sendrecv               &
       (bufsnd_y,mmx,mdp,msend, 0 &
       ,bufrcv_y,mmx,mdp,mrecv, 0 &
       ,mcomw,mstat,merr)

  if(mrecv /= mnull)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,margin
        do i=1,ix
           ro(i,jx-margin+j) = bufrcv_y(i,j,1)
           pr(i,jx-margin+j) = bufrcv_y(i,j,2)
           vx(i,jx-margin+j) = bufrcv_y(i,j,3)
           vy(i,jx-margin+j) = bufrcv_y(i,j,4)
           vz(i,jx-margin+j) = bufrcv_y(i,j,5)
           bx(i,jx-margin+j) = bufrcv_y(i,j,6)
           by(i,jx-margin+j) = bufrcv_y(i,j,7)
           bz(i,jx-margin+j) = bufrcv_y(i,j,8)
           phi(i,jx-margin+j) = bufrcv_y(i,j,9)
           eta(i,jx-margin+j) = bufrcv_y(i,j,10)
        enddo
     enddo
     !$OMP END PARALLEL DO
  end if

!--- Step 1-4.---
! y-direction
!----------------
! forth

  mmx = ix*margin*mx
  msend = mpid%f
  mrecv = mpid%b

  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j)
  do j=1,margin
     do i=1,ix
        bufsnd_y(i,j,1) = ro(i,jx-2*margin+j)
        bufsnd_y(i,j,2) = pr(i,jx-2*margin+j)
        bufsnd_y(i,j,3) = vx(i,jx-2*margin+j)
        bufsnd_y(i,j,4) = vy(i,jx-2*margin+j)
        bufsnd_y(i,j,5) = vz(i,jx-2*margin+j)
        bufsnd_y(i,j,6) = bx(i,jx-2*margin+j)
        bufsnd_y(i,j,7) = by(i,jx-2*margin+j)
        bufsnd_y(i,j,8) = bz(i,jx-2*margin+j)
        bufsnd_y(i,j,9) = phi(i,jx-2*margin+j)
        bufsnd_y(i,j,10) = eta(i,jx-2*margin+j)
     enddo
  enddo
  !$OMP END PARALLEL DO

  call mpi_sendrecv               &
       (bufsnd_y,mmx,mdp,msend, 0 &
       ,bufrcv_y,mmx,mdp,mrecv, 0 &
       ,mcomw,mstat,merr)

  if(mrecv /= mnull)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do j=1,margin
        do i=1,ix
           ro(i,j) = bufrcv_y(i,j,1)
           pr(i,j) = bufrcv_y(i,j,2)
           vx(i,j) = bufrcv_y(i,j,3)
           vy(i,j) = bufrcv_y(i,j,4)
           vz(i,j) = bufrcv_y(i,j,5)
           bx(i,j) = bufrcv_y(i,j,6)
           by(i,j) = bufrcv_y(i,j,7)
           bz(i,j) = bufrcv_y(i,j,8)
           phi(i,j) = bufrcv_y(i,j,9)
           eta(i,j) = bufrcv_y(i,j,10)
        enddo
     enddo
     !$OMP END PARALLEL DO
  end if

end subroutine boundary__mpi


end module boundary
