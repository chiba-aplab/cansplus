module boundary

  implicit none
  private

  public :: bd_consx, bd_consy, bd_consz, bd_frex, bd_frey, bd_frez, &
            bd_inix, bd_iniy, bd_iniz, bd_perx, bd_pery, bd_perz,    &
            bd_synnx_car, bd_synny_car, bd_synnz_car, bd_synnx,      &
            bd_synpx_car, bd_synpy_car, bd_synpz_car, bd_synpx


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
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ibnd-i,j,k) = cons
           enddo
        enddo
     enddo
  else
     ibnd = ix-margin
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ibnd+i,j,k) = cons
           enddo
        enddo
     enddo
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
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,jbnd-j,k) = cons
           enddo
        enddo
     enddo
  else
     jbnd = jx-margin
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,jbnd+j,k) = cons
           enddo
        enddo
     enddo
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
        do j=1,jx
           do i=1,ix
              qq(i,j,kbnd-k) = cons
           enddo
        enddo
     enddo
  else
     kbnd = kx-margin
     do k=1,margin
        do j=1,jx
           do i=1,ix
              qq(i,j,kbnd+k) = cons
           enddo
        enddo
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
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(i,j,k) = qq(margin+1,j,k)
           enddo
        enddo
     enddo
  else
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ix-margin+i,j,k) = qq(ix-margin,j,k)
           enddo
        enddo
     enddo
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
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,jbnd-j,k) = qq(i,jbnd,k)
           enddo
        enddo
     enddo
  else
     jbnd = jx-margin
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,jbnd+j,k) = qq(i,jbnd,k)
           enddo
        enddo
     enddo
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
        do j=1,jx
           do i=1,ix
              qq(i,j,k) = qq(i,j,margin+1)
           enddo
        enddo
     enddo
  else
     do k=1,margin
        do j=1,jx
           do i=1,ix
              qq(i,j,kx-margin+k) = qq(i,j,kx-margin)
           enddo
        enddo
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
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ibnd-i,j,k) = qqini(ibnd-i,j,k)
           enddo
        enddo
     enddo
  else
     ibnd = ix-margin
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ibnd+i,j,k) = qqini(ibnd+i,j,k)
           enddo
        enddo
     enddo
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
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,jbnd-j,k) = qqini(i,jbnd-j,k)
           enddo
        enddo
     enddo
  else
     jbnd = jx-margin
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,jbnd+j,k) = qqini(i,jbnd+j,k)
           enddo
        enddo
     enddo
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
        do j=1,jx
           do i=1,ix
              qq(i,j,kbnd-k) = qqini(i,j,kbnd-k)
           enddo
        enddo
     enddo
  else
     kbnd = kx-margin
     do k=1,margin
        do j=1,jx
           do i=1,ix
              qq(i,j,kbnd+k) = qqini(i,j,kbnd+k)
           enddo
        enddo
     enddo
  endif

  return
end subroutine bd_iniz

subroutine bd_perx(margin1,margin2,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: ix,jx,kx
  integer,intent(in) :: margin1,margin2

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: i,j,k

  do k=1,kx
     do j=1,jx
        do i=1,margin1
           qq(i,j,k) = qq(ix-margin2-margin1+i,j,k)
        enddo
     enddo
  enddo

  do k=1,kx
     do j=1,jx
        do i=1,margin2
           qq(ix-margin2+i,j,k) = qq(margin1+i,j,k)
        enddo
     enddo
  enddo

  return
end subroutine bd_perx

subroutine bd_pery(margin,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: ix,jx,kx
  integer,intent(in) :: margin

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: i,j,k

  do k=1,kx
     do j=1,margin
        do i=1,ix
           qq(i,j,k) = qq(i,jx-2*margin+j,k)
        enddo
     enddo
  enddo
  do k=1,kx
     do j=1,margin
        do i=1,ix
           qq(i,jx-margin+j,k) = qq(i,margin+j,k)
        enddo
     enddo
  enddo

  return
end subroutine bd_pery

subroutine bd_perz(margin1,margin2,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: ix,jx,kx
  integer,intent(in) :: margin1,margin2

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: i,j,k

  do k=1,margin1
     do j=1,jx
        do i=1,ix
           qq(i,j,k) = qq(i,j,kx-margin2-margin1+k)
        enddo
     enddo
  enddo

  do k=1,margin2
     do j=1,jx
        do i=1,ix
           qq(i,j,kx-margin2+k) = qq(i,j,margin1+k)
        enddo
     enddo
  enddo

  return
end subroutine bd_perz

subroutine bd_synnx_car(mbnd,margin,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  integer,intent(in) :: mbnd

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: jpPi,jnPi,hjx,jxmm,jxm
  integer :: i,j,k
  real(8) :: signj

  if(mbnd .eq. 0)then
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(margin-i+1,j,k) = -(qq(margin+i,j,k))
           enddo
        end do
     end do
  else
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ix-margin+i,j,k) = -qq(ix-margin-i+1,j,k)
           end do
        end do
     end do
  end if

  return
end subroutine bd_synnx_car

subroutine bd_synny_car(mbnd,margin,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  integer,intent(in) :: mbnd

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: jpPi,jnPi,hjx,jxmm,jxm
  integer :: i,j,k
  real(8) :: signj

  if(mbnd .eq. 0)then
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,margin-j+1,k) = -(qq(i,margin+j,k))
           enddo
        end do
     end do
  else
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,jx-margin+j,k) = -qq(i,jx-margin-j+1,k)
           end do
        end do
     end do
  end if

  return
end subroutine bd_synny_car

subroutine bd_synnz_car(mbnd,margin,qq,ix,jx,kx)
  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  integer,intent(in) :: mbnd

  real(8),dimension(ix,jx,kx),intent(inout) :: qq

  integer :: jpPi,jnPi,hjx,jxmm,jxm
  integer :: i,j,k
  real(8) :: signj

  if(mbnd .eq. 0)then
     do k=1,margin
        do j=1,jx
           do i=1,ix
              qq(i,j,margin-k+1) = -(qq(i,j,margin+k))
           enddo
        end do
     end do
  else
     do k=1,margin
        do j=1,jx
           do i=1,ix
              qq(i,j,kx-margin+k) = -qq(i,j,kx-margin-k+1)
           end do
        end do
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
  else
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ix-margin+i,j,k) = -qq(ix-margin-i+1,j,k)
           end do
        end do
     end do
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
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(margin-i+1,j,k) = (qq(margin+i,j,k))
           enddo
        end do
     end do
  else
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ix-margin+i,j,k) = qq(ix-margin-i+1,j,k)
           end do
        end do
     end do
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
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,margin-j+1,k) = (qq(i,margin+j,k))
           enddo
        end do
     end do
  else
     do k=1,kx
        do j=1,margin
           do i=1,ix
              qq(i,jx-margin+j,k) = qq(i,jx-margin-j+1,k)
           end do
        end do
     end do
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
        do j=1,jx
           do i=1,ix
              qq(i,j,margin-k+1) = (qq(i,j,margin+k))
           enddo
        end do
     end do
  else
     do k=1,margin
        do j=1,jx
           do i=1,ix
              qq(i,j,kx-margin+k) = qq(i,j,kx-margin-k+1)
           end do
        end do
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
  else
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ix-margin+i,j,k) = qq(ix-margin-i+1,j,k)
           end do
        end do
     end do
  end if

  return
end subroutine bd_synpx

end module boundary
