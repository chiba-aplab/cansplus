! Cell surface boundary 
! 1 | 2 | 3 | 4 | 5 | 6 | ..
! q(3) = q(4), q(2) = q(5) , q(1) = q(6)
subroutine bd_synpz_car(mbnd,margin,qq,ix,jx,kx)

  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  integer,intent(in) :: mbnd

  real(8),dimension(ix,jx,kx) :: qq

  integer :: jpPi,jnPi,hjx,jxmm,jxm
  integer :: i,j,k
  real(8) :: signj


  if(mbnd .eq. 0)then
     do k=1,margin
        do j=1,jx
           do i=1,ix
              qq(i,j,margin-k+1) = (qq(i,j,margin+k))
           enddo
        end do
     end do
  else
     do k=1,margin
        do j=1,jx
           do i=1,ix
              qq(i,j,kx-margin+k) = qq(i,j,kx-margin-k+1)
           end do
        end do
     end do
  end if

  return
end subroutine bd_synpz_car
