subroutine lr_state_minmod(mdir,ix,jx,kx,qq &
     ,qqw)
!======================================================================
! Name :: lr_state_minmod
!         flux limiter :: minmod
! Input :: 
!         mdir :: direction 1:x, 2:y, 3:z
!         ix,jx,kx :: array size
!         qq :: variables
! Output :: 
!         qqw :: qqw(1) :: left state
!                qqw(2) :: right state
!
!======================================================================
  implicit none

!---Input
  integer,intent(in) :: mdir
  integer,intent(in) :: ix,jx,kx

  real(8),dimension(ix,jx,kx),intent(in) :: qq

!---Output
  real(8),dimension(ix,jx,kx,2) :: qqw

!---Temporary
  real(8) :: dqqx,dqqy,dqqz

  real(8) :: temp

  integer :: i,j,k

!-----Step 1a.-------------------------------------------------|
! dqq
!
!-----Step 1b.-------------------------------------------------|
! mdir = 1 :: x-direction
!

  if(mdir .eq. 1)then
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
              dqqx = minmod_limiter2((qq(i+1,j,k)-qq(i,j,k)) &
                   ,(qq(i,j,k)-qq(i-1,j,k)))
              
              qqw(i,j,k,1) = qq(i,j,k) + 0.5d0*dqqx
              
              qqw(i-1,j,k,2) = qq(i,j,k) - 0.5d0*dqqx
              
           enddo
        enddo
     enddo
!-----Step 1b.-------------------------------------------------|
! mdir = 2 :: y-direction
!
  else if(mdir .eq. 2)then
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
              dqqy = minmod_limiter2((qq(i,j+1,k)-qq(i,j,k)) &
                   ,(qq(i,j,k)-qq(i,j-1,k)))
              
              qqw(i,j,k,1) = qq(i,j,k) + 0.5d0*dqqy
              
              qqw(i,j-1,k,2) = qq(i,j,k) - 0.5d0*dqqy
              
           end do
        end do
     end do
  else
!-----Step 3a.-------------------------------------------------|
! dqq
!
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
              dqqz = minmod_limiter2((qq(i,j,k+1)-qq(i,j,k)) &
                   ,(qq(i,j,k)-qq(i,j,k-1)))
           
!-----Step 3b.-------------------------------------------------|
! mdir = 3 :: z-direction
!
              qqw(i,j,k,1) = qq(i,j,k) + 0.5d0*dqqz
              
              qqw(i,j,k-1,2) = qq(i,j,k) - 0.5d0*dqqz
           end do
        end do
     end do
  endif
  return

contains
  function minmod_limiter(qqr,qql)
    implicit none
    
    real(8) :: qqr, qql
    real(8) :: minmod_limiter

    minmod_limiter = max(0.0d0, min(dabs(qqr),qql*sign(1.0d0,qqr)))*sign(1.0d0,qqr)
    return
  end function minmod_limiter

  function minmod_limiter2(qqr,qql)
    implicit none

    real(8) :: qqr, qql
    real(8) :: minmod_limiter2

    if(qqr > 0.0d0)then
       minmod_limiter2 = max(0.0d0,min(qql,qqr))
    else
       minmod_limiter2 = min(0.0d0,max(qql,qqr))
    endif

    return
  end function minmod_limiter2
  
  function MC_limiter(qqc,qql,qqr)
    implicit none

    real(8) :: qqr,qql,qqc
    real(8) :: minmod_lr
    real(8) :: MC_limiter

    minmod_lr = minmod_limiter(2.0d0*qql,2.0d0*qqr)
    MC_limiter = minmod_limiter(0.5d0*qqc,minmod_lr)
    
    return
  end function MC_limiter
end subroutine lr_state_minmod
