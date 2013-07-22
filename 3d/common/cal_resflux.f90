! 
! Average eta*current at cell curface
!
subroutine cal_resflux(mdir,ix,jx,kx,fbx,curx,eta,pm &
     ,fbx_res)

  implicit none

  integer,intent(in) :: ix,jx,kx,mdir

! +1 or -1, consistency between Electric fields and numerical flux
  real(8),intent(in) :: pm

  real(8),dimension(ix,jx,kx) :: fbx,curx,eta
  real(8),dimension(ix,jx,kx) :: fbx_res

  real(8) :: fres

  integer :: i,j,k


! Flux at i+1/2
  if(mdir .eq. 1)then
     do k=2,kx-2
        do j=2,jx-2
           do i=2,ix-2
              fres = 0.5d0*pm*(eta(i,j,k)*curx(i,j,k) &
                   +eta(i+1,j,k)*curx(i+1,j,k))

              fbx_res(i,j,k) = fbx(i,j,k)+fres
           end do
        end do
     end do
! Flux at j+1/2
  else if(mdir .eq. 2)then
     do k=2,kx-2
        do j=2,jx-2
           do i=2,ix-2
              fres = 0.5d0*pm*(eta(i,j,k)*curx(i,j,k) &
                   +eta(i,j+1,k)*curx(i,j+1,k))
              fbx_res(i,j,k) = fbx(i,j,k)+fres
           end do
        end do
     end do
! Flux at k+1/2
  else
     do k=2,kx-2
        do j=2,jx-2
           do i=2,ix-2
              fres = 0.5d0*pm*(eta(i,j,k)*curx(i,j,k) &
                   +eta(i,j,k+1)*curx(i,j,k+1))
              fbx_res(i,j,k) = fbx(i,j,k)+fres
           end do
        end do
     end do
  end if

  return
end subroutine cal_resflux
