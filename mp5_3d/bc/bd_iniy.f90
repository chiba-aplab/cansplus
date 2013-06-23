!=====================================================================
subroutine bd_iniy(mbnd,margin,qq,qqini,ix,jx,kx)
!=====================================================================
!
! NAME :: bd_iniy
!
! INPUT
!      mbnd :: margin flag
!      margin :: margin
!      ix,jx,kx :: array size
!      qqini :: initial
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
