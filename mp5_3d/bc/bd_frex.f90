!=====================================================================
subroutine bd_frex(mbnd,margin,qq,ix,jx,kx)
!=====================================================================
!
! NAME :: bd_frex
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
!----<INPUT

!---->INPUT&OUTPUT
  real(8),dimension(ix,jx,kx) :: qq
!----<INPUT&OUTPUT

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
