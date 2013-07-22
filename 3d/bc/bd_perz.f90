!=====================================================================
subroutine bd_perz(margin1,margin2,qq,ix,jx,kx)
!=====================================================================
!
! Name :: bd_perz
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
