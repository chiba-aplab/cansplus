!=====================================================================
subroutine bd_pery(margin,qq,ix,jx,kx)
!=====================================================================
!
! Name :: bd_pery
!
! INPUT
!      margin1,margin2 :: margin
!      ix,jx,kx :: array size
! INPUT&OUTPUT
!      qq :: variables
!=====================================================================

  implicit none

  integer,intent(in) :: ix,jx,kx
  integer,intent(in) :: margin

  real(8),dimension(ix,jx,kx) :: qq

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
