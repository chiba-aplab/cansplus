!=====================================================================
subroutine bd_consy(mbnd,margin,cons,qq,ix,jx,kx)
!=====================================================================
!
! NAME :: bd_consy
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
