!=====================================================================
subroutine bd_inix(mbnd,margin,qq,qqini,ix,jx,kx)
!=====================================================================
!
! NAME :: bd_inix
!
! INPUT
!      mbnd :: margin flag
!      margin :: margin
!      ix,jx,kx :: array size
!      qqini :: initial quantity
! INPUT&OUTPUT
!      qq :: variables
!
!=====================================================================

  implicit none

!---->INPUT
  integer,intent(in) :: mbnd ! boundary flag
  integer,intent(in) :: margin ! margin
  integer,intent(in) :: ix,jx,kx ! array size
!----<INPUT

!---->INPUT&OUTPUT
  real(8),dimension(ix,jx,kx) :: qq,qqini
!----<INPUT&OUTPUT

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
