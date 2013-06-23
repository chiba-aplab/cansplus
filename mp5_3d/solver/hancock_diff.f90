subroutine hancock_diff(mdir,ix,jx,kx,qq &
     ,hDiff)
  implicit none

!---Input
  integer,intent(in) :: mdir
  integer,intent(in) :: ix,jx,kx

  real(8),dimension(ix,jx,kx),intent(in) :: qq
  
  real(8),dimension(ix,jx,kx),intent(inout) :: hDiff

  real(8) :: temp

  integer :: i,j,k

  if(mdir .eq. 1)then
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
              hDiff(i,j,k) = ave(qq(i,j,k)-qq(i-1,j,k),qq(i+1,j,k)-qq(i,j,k))
           end do
        end do
     end do
  elseif(mdir .eq. 2)then
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
              hDiff(i,j,k) = ave(qq(i,j,k)-qq(i,j-1,k),qq(i,j+1,k)-qq(i,j,k))
           end do
        end do
     end do
  else
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
              hDiff(i,j,k) = ave(qq(i,j,k)-qq(i,j,k-1),qq(i,j,k+1)-qq(i,j,k))
           end do
        end do
     end do
  end if
  return
contains
  
  function minmod_limiter(qqr,qql)
    implicit none
    real(8) :: qqr, qql
    real(8) :: minmod_limiter

!    minmod_limiter = max(0.0d0, min(dabs(qqr),qql*sign(1.0d0,qqr)))*sign(1.0d0,qqr)
    if(qqr > 0.0d0)then
       minmod_limiter = max(0.0d0,min(qql,qqr))
    else
       minmod_limiter = min(0.0d0,max(qql,qqr))
    endif

    return
  end function minmod_limiter    

  function ave(qqr,qql)
    implicit none
    real(8) :: qqr, qql
    real(8) :: minmod_lr
    real(8) :: ave

!    ave = minmod_limiter(qqr,qql)
    if(qql*qqr > 0.0d0)then
       minmod_lr = minmod_limiter(qql,qqr)
       ave = minmod_limiter(0.5d0*(qql+qqr),minmod_lr)
    else
       ave = 0.0d0
    end if
    return
  end function ave

end subroutine hancock_diff
