!=====================================================================
subroutine bd_iniz(mbnd,margin,qq,qqini,ix,jx,kx)
!=====================================================================
!
! NAME :: bd_iniz
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
