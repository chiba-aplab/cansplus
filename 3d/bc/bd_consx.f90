!=====================================================================
subroutine bd_consx(mbnd,margin,cons,qq,ix,jx,kx)
!=====================================================================
!
! NAME :: bd_consx
!
! INPUT
!      ix,jx,kx :: array size
!      margin :: margin
!      cons :: constant number
! INPUT&OUTPUT
!      qq
!=====================================================================

  implicit none

  integer,intent(in) :: ix,jx,kx
  integer,intent(in) :: margin
  integer,intent(in) :: mbnd

  real(8),intent(in) :: cons

  real(8),dimension(ix,jx,kx) :: qq

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
