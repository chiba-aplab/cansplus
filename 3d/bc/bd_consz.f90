!=====================================================================
subroutine bd_consz(mbnd,margin,cons,qq,ix,jx,kx)
!=====================================================================
!
! NAME :: bd_frez
!
! INPUT
!      mbnd :: margin flag
!      margin :: margin
!      ix,jx,kx :: array size
! INPUT&OUTPUT
!      qq :: variables
!
!=====================================================================

  implicit none

!---->INPUT
  integer,intent(in) :: mbnd ! boundary flag
  integer,intent(in) :: margin ! margin
  integer,intent(in) :: ix,jx,kx ! array size
  real(8),intent(in) :: cons
!----<INPUT

!---->INPUT&OUTPUT
  real(8),dimension(ix,jx,kx) :: qq
!----<INPUT&OUTPUT

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
