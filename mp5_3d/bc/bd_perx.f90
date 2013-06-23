!=====================================================================
subroutine bd_perx(margin1,margin2,qq,ix,jx,kx)
!=====================================================================
!
! Name :: bd_perx
!
! INPUT
!      margin1,margin2 :: margin
!      ix,jx,kx :: array size
! INPUT&OUTPUT
!      qq :: variables
!=====================================================================

  implicit none

  integer,intent(in) :: ix,jx,kx
  integer,intent(in) :: margin1,margin2

  real(8),dimension(ix,jx,kx) :: qq

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
