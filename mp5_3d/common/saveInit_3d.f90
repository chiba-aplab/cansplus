subroutine saveInit_3d(ix,jx,kx,qq,qqi)
  implicit none

  integer :: ix,jx,kx

  real(8),dimension(ix,jx,kx) :: qq, qqi

  integer :: i,j,k

  do k=1,kx
     do j=1,jx
        do i=1,ix
           qqi(i,j,k) = qq(i,j,k)
        enddo
     enddo
  enddo

  return
end subroutine saveInit_3d
