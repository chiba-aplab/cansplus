subroutine cal_surface(ix,jx,kx,bx,by,bz &
     ,bxs,bys,bzs)
  implicit none

  integer,intent(in) :: ix,jx,kx

  real(8),dimension(ix,jx,kx) :: bx,by,bz
  real(8),dimension(ix,jx,kx) :: bxs,bys,bzs

  integer :: i,j,k

  do k=1,kx
     do j=1,jx
        do i=1,ix-1
           bxs(i,j,k) = 0.5d0*(bx(i,j,k)+bx(i+1,j,k))
        end do
     end do
  end do

  do k=1,kx
     do j=1,jx-1
        do i=1,ix
           bys(i,j,k) = 0.5d0*(by(i,j,k)+by(i,j+1,k))
        enddo
     enddo
  enddo

  do k=1,kx-1
     do j=1,jx
        do i=1,ix
           bzs(i,j,k) = 0.5d0*(bz(i,j,k)+bz(i,j,k+1))
        end do
     end do
  end do

  return
end subroutine cal_surface

