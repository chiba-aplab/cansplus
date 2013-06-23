subroutine lr_state_1D(mdir,ix,jx,kx,qq &
     ,qqw)
  implicit none

  integer,intent(in) :: mdir
  integer,intent(in) :: ix,jx,kx

  real(8),dimension(ix,jx,kx) :: qq

  real(8),dimension(ix,jx,kx,2) :: qqw

  real(8) :: dqqx,dqqy,dqqz

  real(8) :: temp

  integer :: i,j,k

  if(mdir .eq. 1)then
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
              qqw(i,j,k,1) = qq(i,j,k)
              qqw(i-1,j,k,2) = qq(i,j,k)
           enddo
        enddo
     enddo
  else if(mdir .eq. 2)then
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
              qqw(i,j,k,1) = qq(i,j,k)
              qqw(i,j-1,k,2) = qq(i,j,k)
           end do
        end do
     end do
  else
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
              qqw(i,j,k,1) = qq(i,j,k)
              qqw(i,j,k-1,2) = qq(i,j,k)
           end do
        end do
     end do
  endif
  return
end subroutine lr_state_1D
