module boundary

  implicit none
  private

  public :: bd_consx, bd_consy, bd_consz, bd_frex, bd_frey, bd_frez, &
            bd_inix, bd_iniy, bd_iniz,                               &
            bd_synnx_car, bd_synny_car, bd_synnz_car, bd_synnx,      &
            bd_synpx_car, bd_synpy_car, bd_synpz_car, bd_synpx,      &
            boundary__mpi, boundary__mpi_cyl


contains


!=====================================================================
! INPUT
!      mbnd :: margin flag
!      ix,jx,kx :: array size
!      margin :: margin
! INPUT&OUTPUT
!      qq
!=====================================================================

subroutine bd_consx(mbnd,margin,cons,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: ix,jx,kx
  integer,intent(in) :: margin
  integer,intent(in) :: mbnd
  real(8),intent(in) :: cons

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: ibnd
  integer :: i,j,k

  if (mbnd .eq. 0)then
     ibnd = 1+margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ibnd-i,j,k) = cons
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  else
     ibnd = ix-margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ibnd+i,j,k) = cons
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

  return
end subroutine bd_consx

subroutine bd_consy(mbnd,margin,cons,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: mbnd ! boundary flag                                       
  integer,intent(in) :: margin ! margin                                            
  integer,intent(in) :: ix,jx,kx ! array size                                      
  real(8),intent(in) :: cons

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: jbnd
  integer :: i,j,k

  if( mbnd .eq. 0)then
     jbnd = 1+margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,jbnd-j,k) = cons
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  else
     jbnd = jx-margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,jbnd+j,k) = cons
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

  return
end subroutine bd_consy

subroutine bd_consz(mbnd,margin,cons,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: mbnd ! boundary flag
  integer,intent(in) :: margin ! margin
  integer,intent(in) :: ix,jx,kx ! array size
  real(8),intent(in) :: cons

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: kbnd
  integer :: i,j,k

  if (mbnd .eq. 0)then
     kbnd = 1+margin
     do k=1,margin
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(i)
        do j=1,jx
           do i=1,ix
              qq(i,j,kbnd-k) = cons
           enddo
        enddo
        !$OMP END PARALLEL DO
     enddo
  else
     kbnd = kx-margin
     do k=1,margin
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(i)
        do j=1,jx
           do i=1,ix
              qq(i,j,kbnd+k) = cons
           enddo
        enddo
        !$OMP END PARALLEL DO
     enddo
  endif

  return
end subroutine bd_consz

subroutine bd_frex(mbnd,margin,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: mbnd ! boundary flag
  integer,intent(in) :: margin ! margin
  integer,intent(in) :: ix,jx,kx ! array size

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: i,j,k

  if (mbnd .eq. 0)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(i,j,k) = qq(margin+1,j,k)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ix-margin+i,j,k) = qq(ix-margin,j,k)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

  return
end subroutine bd_frex

subroutine bd_frey(mbnd,margin,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: mbnd ! boundary flag
  integer,intent(in) :: margin ! margin
  integer,intent(in) :: ix,jx,kx ! array size

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: jbnd
  integer :: i,j,k

  if( mbnd .eq. 0)then
     jbnd = 1+margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,jbnd-j,k) = qq(i,jbnd,k)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  else
     jbnd = jx-margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,jbnd+j,k) = qq(i,jbnd,k)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

  return
end subroutine bd_frey

subroutine bd_frez(mbnd,margin,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: mbnd ! boundary flag
  integer,intent(in) :: margin ! margin
  integer,intent(in) :: ix,jx,kx ! array size

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: i,j,k

  if (mbnd .eq. 0)then
     do k=1,margin
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(i)
        do j=1,jx
           do i=1,ix
              qq(i,j,k) = qq(i,j,margin+1)
           enddo
        enddo
        !$OMP END PARALLEL DO
     enddo
  else
     do k=1,margin
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(i)
        do j=1,jx
           do i=1,ix
              qq(i,j,kx-margin+k) = qq(i,j,kx-margin)
           enddo
        enddo
        !$OMP END PARALLEL DO
     enddo
  endif

  return
end subroutine bd_frez

subroutine bd_inix(mbnd,margin,qq,qqini,ix,jx,kx)
  implicit none

  integer,intent(in) :: mbnd ! boundary flag
  integer,intent(in) :: margin ! margin
  integer,intent(in) :: ix,jx,kx ! array size
  real(8),dimension(ix,jx,kx),intent(in) :: qqini ! initial quantity

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: ibnd
  integer :: i,j,k

  if (mbnd .eq. 0)then
     ibnd = 1+margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ibnd-i,j,k) = qqini(ibnd-i,j,k)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  else
     ibnd = ix-margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ibnd+i,j,k) = qqini(ibnd+i,j,k)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

  return
end subroutine bd_inix

subroutine bd_iniy(mbnd,margin,qq,qqini,ix,jx,kx)
  implicit none

  integer,intent(in) :: mbnd ! boundary flag
  integer,intent(in) :: margin ! margin
  integer,intent(in) :: ix,jx,kx ! array size
  real(8),dimension(ix,jx,kx),intent(in) :: qqini ! initial quantity

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: jbnd
  integer :: i,j,k

  if( mbnd .eq. 0)then
     jbnd = 1+margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,jbnd-j,k) = qqini(i,jbnd-j,k)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  else
     jbnd = jx-margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,jbnd+j,k) = qqini(i,jbnd+j,k)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

  return
end subroutine bd_iniy

subroutine bd_iniz(mbnd,margin,qq,qqini,ix,jx,kx)
  implicit none

  integer,intent(in) :: mbnd ! boundary flag
  integer,intent(in) :: margin ! margin
  integer,intent(in) :: ix,jx,kx ! array size
  real(8),dimension(ix,jx,kx),intent(in) :: qqini ! initial quantity

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: kbnd
  integer :: i,j,k

  if (mbnd .eq. 0)then
     kbnd = 1+margin
     do k=1,margin
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(i)
        do j=1,jx
           do i=1,ix
              qq(i,j,kbnd-k) = qqini(i,j,kbnd-k)
           enddo
        enddo
        !$OMP END PARALLEL DO
     enddo
  else
     kbnd = kx-margin
     do k=1,margin
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(i)
        do j=1,jx
           do i=1,ix
              qq(i,j,kbnd+k) = qqini(i,j,kbnd+k)
           enddo
        enddo
        !$OMP END PARALLEL DO
     enddo
  endif

  return
end subroutine bd_iniz

subroutine bd_synnx_car(mbnd,margin,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  integer,intent(in) :: mbnd

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: i,j,k

  if(mbnd .eq. 0)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(margin-i+1,j,k) = -(qq(margin+i,j,k))
           enddo
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ix-margin+i,j,k) = -qq(ix-margin-i+1,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  return
end subroutine bd_synnx_car

subroutine bd_synny_car(mbnd,margin,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  integer,intent(in) :: mbnd

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: i,j,k

  if(mbnd .eq. 0)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,margin-j+1,k) = -(qq(i,margin+j,k))
           enddo
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,jx-margin+j,k) = -qq(i,jx-margin-j+1,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  return
end subroutine bd_synny_car

subroutine bd_synnz_car(mbnd,margin,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  integer,intent(in) :: mbnd

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: i,j,k

  if(mbnd .eq. 0)then
     do k=1,margin
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(i)
        do j=1,jx
           do i=1,ix
              qq(i,j,margin-k+1) = -(qq(i,j,margin+k))
           enddo
        end do
        !$OMP END PARALLEL DO
     end do
  else
     do k=1,margin
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(i)
        do j=1,jx
           do i=1,ix
              qq(i,j,kx-margin+k) = -qq(i,j,kx-margin-k+1)
           end do
        end do
        !$OMP END PARALLEL DO
     end do
  end if

  return
end subroutine bd_synnz_car

subroutine bd_synnx(mbnd,margin,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  integer,intent(in) :: mbnd

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: jpPi,jnPi,hjx,jxmm,jxm
  integer :: i,j,k
  real(8) :: signj

  jxm = jx-margin
  jxmm = jx-2*margin
  hjx = (jx-2*margin)/2

  if(mbnd .eq. 0)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,jpPi,jnPi,signj)
     do k=1,kx
        do j=1,jx
           do i=1,margin
              jpPi = j+hjx
              jnPi = j-hjx

              signj = sign(1.0d0,jxm-jpPi+0.00001d0)

              jpPi = min(jx,jpPi)
              jpPi = max(1,jpPi)
              jnPi = max(1,jnPi)
              jnPi = min(jx,jnPi)

              qq(margin-i+1,j,k) = max(0.0d0,signj)*(-qq(margin+i,jpPi,k)) &
                   +max(0.0d0,-signj)*(-qq(margin+i,jnPi,k))
           enddo
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ix-margin+i,j,k) = -qq(ix-margin-i+1,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  return
end subroutine bd_synnx

subroutine bd_synpx_car(mbnd,margin,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  integer,intent(in) :: mbnd

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: i,j,k

  if(mbnd .eq. 0)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(margin-i+1,j,k) = (qq(margin+i,j,k))
           enddo
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ix-margin+i,j,k) = qq(ix-margin-i+1,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  return
end subroutine bd_synpx_car

subroutine bd_synpy_car(mbnd,margin,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  integer,intent(in) :: mbnd

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: i,j,k

  if(mbnd .eq. 0)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,margin-j+1,k) = (qq(i,margin+j,k))
           enddo
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,jx-margin+j,k) = qq(i,jx-margin-j+1,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  return
end subroutine bd_synpy_car

subroutine bd_synpz_car(mbnd,margin,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  integer,intent(in) :: mbnd

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: i,j,k

  if(mbnd .eq. 0)then
     do k=1,margin
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(i)
        do j=1,jx
           do i=1,ix
              qq(i,j,margin-k+1) = (qq(i,j,margin+k))
           enddo
        end do
        !$OMP END PARALLEL DO
     end do
  else
     do k=1,margin
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(i)
        do j=1,jx
           do i=1,ix
              qq(i,j,kx-margin+k) = qq(i,j,kx-margin-k+1)
           end do
        end do
        !$OMP END PARALLEL DO
     end do
  end if

  return
end subroutine bd_synpz_car

subroutine bd_synpx(mbnd,margin,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  integer,intent(in) :: mbnd

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: jpPi,jnPi,hjx,jxmm,jxm
  integer :: i,j,k
  real(8) :: signj

  jxm = jx-margin
  jxmm = jx-2*margin
  hjx = (jx-2*margin)/2

  if(mbnd .eq. 0)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,jpPi,jnPi,signj)
     do k=1,kx
        do j=1,jx
           do i=1,margin
              jpPi = j+hjx
              jnPi = j-hjx

              signj = sign(1.0d0,jxm-jpPi+0.00001d0)

              jpPi = min(jx,jpPi)
              jpPi = max(1,jpPi)
              jnPi = max(1,jnPi)
              jnPi = min(jx,jnPi)

              qq(margin-i+1,j,k) = max(0.0d0,signj)*(qq(margin+i,jpPi,k)) &
                   +max(0.0d0,-signj)*(qq(margin+i,jnPi,k))
           enddo
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ix-margin+i,j,k) = qq(ix-margin-i+1,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  end subroutine bd_synpx


  subroutine boundary__mpi(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta)

  use mpi_setup

  integer,intent(in) :: ix,jx,kx,margin
  real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz,phi,eta

  integer,parameter :: mx = 10
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

  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j)
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
     !$OMP END PARALLEL DO
  endif

!--- Step 1-2.---
! right

  mmx = margin*jx*kx*mx
  msend = mpid%r
  mrecv = mpid%l

  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j)
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
        enddo
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
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  end if

!--- Step 1-3.---
! y-direction
!----------------
! back
              
  mmx = ix*kx*margin*mx
  msend = mpid%b
  mrecv = mpid%f
  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j)
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
           bufsnd_y(i,j,k,10) = eta(i,margin+j,k)
        enddo
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
              eta(i,jx-margin+j,k) = bufrcv_y(i,j,k,10)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  end if

!--- Step 1-4.---
! y-direction
!----------------
! forth

  mmx = ix*kx*margin*mx
  msend = mpid%f
  mrecv = mpid%b
  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j)
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
           bufsnd_y(i,j,k,10) = eta(i,jx-2*margin+j,k)
        enddo
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
              eta(i,j,k) = bufrcv_y(i,j,k,10)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  end if

!--- Step 1-5.---
! z-direction
!----------------
! down
              
  mmx = ix*jx*margin*mx
  msend = mpid%d
  mrecv = mpid%t

  do k=1,margin
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i)
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
        enddo
     enddo
     !$OMP END PARALLEL DO
  enddo

  call mpi_sendrecv               &
       (bufsnd_z,mmx,mdp,msend, 0 &
       ,bufrcv_z,mmx,mdp,mrecv, 0 &
       ,mcomw,mstat,merr)

  if(mrecv /= mnull)then
     do k=1,margin
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(i)
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
           enddo
        enddo
        !$OMP END PARALLEL DO
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
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i)
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
        enddo
     enddo
     !$OMP END PARALLEL DO
  enddo

  call mpi_sendrecv               &
       (bufsnd_z,mmx,mdp,msend, 0 &
       ,bufrcv_z,mmx,mdp,mrecv, 0 &
       ,mcomw,mstat,merr)

  if(mrecv /= mnull)then
     do k=1,margin
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(i)
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
           enddo
        enddo
        !$OMP END PARALLEL DO
     enddo
  end if

  end subroutine boundary__mpi


  subroutine boundary__mpi_cyl(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta)

  use mpi_setup

  integer,intent(in) :: ix,jx,kx,margin
  real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz,phi,eta

  integer,parameter :: mx = 10
  integer :: i,j,k,mmx,msend,mrecv
  real(8),dimension(margin,jx,kx,mx) :: bufsnd_x,bufrcv_x

!=============================================================
! surface exchange
!=============================================================
! r-direction
!----------------
! across axis

  mmx = margin*jx*kx*mx
  msend = mpid%ax
  mrecv = mpid%ax

  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j)
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
     enddo
  enddo
  !$OMP END PARALLEL DO

  call mpi_sendrecv              &
       (bufsnd_x,mmx,mdp,msend,0 &
       ,bufrcv_x,mmx,mdp,mrecv,0 &
       ,mcomw,mstat,merr)

  if(mpid%l == mnull)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=1,kx
        do j=1,jx
           do i=1,margin
              ro(margin+1-i,j,k)  = +bufrcv_x(i,j,k,1)
              pr(margin+1-i,j,k)  = +bufrcv_x(i,j,k,2)
              vx(margin+1-i,j,k)  = -bufrcv_x(i,j,k,3)
              vy(margin+1-i,j,k)  = -bufrcv_x(i,j,k,4)
              vz(margin+1-i,j,k)  = +bufrcv_x(i,j,k,5)
              bx(margin+1-i,j,k)  = -bufrcv_x(i,j,k,6)
              by(margin+1-i,j,k)  = -bufrcv_x(i,j,k,7)
              bz(margin+1-i,j,k)  = +bufrcv_x(i,j,k,8)
              phi(margin+1-i,j,k) = +bufrcv_x(i,j,k,9)
              eta(margin+1-i,j,k) = +bufrcv_x(i,j,k,10)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

  end subroutine boundary__mpi_cyl


end module boundary
