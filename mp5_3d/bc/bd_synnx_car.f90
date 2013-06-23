! Cell surface boundary 
! 1 | 2 | 3 | 4 | 5 | 6 | ..
! q(3) = q(4), q(2) = q(5) , q(1) = q(6)
subroutine bd_synnx_car(mbnd,margin,qq,ix,jx,kx)

  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  integer,intent(in) :: mbnd

  real(8),dimension(ix,jx,kx) :: qq

  integer :: jpPi,jnPi,hjx,jxmm,jxm
  integer :: i,j,k
  real(8) :: signj


  if(mbnd .eq. 0)then
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(margin-i+1,j,k) = -(qq(margin+i,j,k))
           enddo
        end do
     end do
  else
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ix-margin+i,j,k) = -qq(ix-margin-i+1,j,k)
           end do
        end do
     end do
  end if

  return
end subroutine bd_synnx_car
