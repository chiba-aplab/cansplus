! 
! Average eta*current at cell curface
!
subroutine cal_resflux_fee(mdir,ix,jx,kx,fee,curx,cury,curz,bx,by,bz,eta &
     ,fee_res)

  implicit none

  integer,intent(in) :: ix,jx,kx,mdir

! +1 or -1, consistency between Electric fields and numerical flux

  real(8),dimension(ix,jx,kx) :: fee,curx,cury,curz,bx,by,bz,eta
  real(8),dimension(ix,jx,kx) :: fee_res

  real(8) :: fres

  integer :: i,j,k


! Flux at i+1/2
  if(mdir .eq. 1)then
     do k=2,kx-2
        do j=2,jx-2
           do i=2,ix-2
              fres =0.5d0*(eta(i,j,k)*(cury(i,j,k)*bz(i,j,k)-curz(i,j,k)*by(i,j,k)) &
                   +eta(i+1,j,k)*(cury(i+1,j,k)*bz(i+1,j,k)-curz(i+1,j,k)*by(i+1,j,k)))
              fee_res(i,j,k) = fee(i,j,k)+fres
           end do
        end do
     end do
! Flux at j+1/2
  else if(mdir .eq. 2)then
     do k=2,kx-2
        do j=2,jx-2
           do i=2,ix-2
              fres =0.5d0*(eta(i,j,k)*(curz(i,j,k)*bx(i,j,k)-curx(i,j,k)*bz(i,j,k)) &
                   +eta(i,j+1,k)*(curz(i,j+1,k)*bx(i,j+1,k)-curx(i,j+1,k)*bz(i,j+1,k)))
              fee_res(i,j,k) = fee(i,j,k)+fres
           end do
        end do
     end do
! Flux at k+1/2
  else
     do k=2,kx-2
        do j=2,jx-2
           do i=2,ix-2
              fres =0.5d0*(eta(i,j,k)*(curx(i,j,k)*by(i,j,k)-cury(i,j,k)*bx(i,j,k)) &
                   +eta(i,j,k+1)*(curx(i,j,k+1)*by(i,j,k+1)-cury(i,j,k+1)*bx(i,j,k+1)))
              fee_res(i,j,k) = fee(i,j,k)+fres
           end do
        end do
     end do
  end if

  return
end subroutine cal_resflux_fee
