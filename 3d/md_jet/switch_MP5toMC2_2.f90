subroutine switchMP5toMC2_2(mdir,ix,jx,kx,ro,qq &
     ,qqw,floor,limit,x,dx,y,dy,z,dz)
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
!密度勾配が急なところで精度切り替え
!
!======================================================================
  implicit none

!---Input
  integer,intent(in) :: mdir
  integer,intent(in) :: ix,jx,kx

  real(8),dimension(ix,jx,kx) :: ro,qq

!---Output
  real(8),dimension(ix,jx,kx,2) :: qqw

  real(8),intent(in) :: floor,limit

  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: y,dy
  real(8),dimension(kx),intent(in) :: z,dz

  real(8) :: floor2

!---Temporary
  real(8) :: dqqx,dqqy,dqqz

  real(8) :: temp,dif1,dif2,dif3,dif4

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
           do i=3,ix-2
				temp=dabs((ro(i+2,j,k)-ro(i+1,j,k))/dx(i+1))
				dif1=dlog10(temp) 
				temp=dabs((ro(i+1,j,k)-ro(i,j,k))/dx(i))
				dif2=dlog10(temp) 
				temp=dabs((ro(i,j,k)-ro(i-1,j,k))/dx(i-1))
				dif3=dlog10(temp) 
				temp=dabs((ro(i-1,j,k)-ro(i-2,j,k))/dx(i-2))
				dif4=dlog10(temp) 
		   temp=max(dif1,dif2,dif3,dif4) 
           if(temp > limit) then
              dqqx = MC_limiter((qq(i+1,j,k)-qq(i-1,j,k)) &
                   ,(qq(i+1,j,k)-qq(i,j,k)) &
                   ,(qq(i,j,k)-qq(i-1,j,k)))
              
              qqw(i,j,k,1) = qq(i,j,k) + 0.5d0*dqqx
              
              qqw(i-1,j,k,2) = qq(i,j,k) - 0.5d0*dqqx
            endif  
           enddo
        enddo
     enddo
!-----Step 1b.-------------------------------------------------|
! mdir = 2 :: y-direction
!
  else if(mdir .eq. 2)then
     do k=2,kx-1
        do j=3,jx-2
           do i=2,ix-1
				temp=dabs((ro(i,j+2,k)-ro(i,j+1,k))/(x(i)*dy(j+1)))
				dif1=dlog10(temp) 
				temp=dabs((ro(i,j+1,k)-ro(i,j,k))/(x(i)*dy(j)))
				dif2=dlog10(temp) 
				temp=dabs((ro(i,j,k)-ro(i,j-1,k))/(x(i)*dy(j-1)))
				dif3=dlog10(temp) 
				temp=dabs((ro(i,j-1,k)-ro(i,j-2,k))/(x(i)*dy(j-2)))
				dif4=dlog10(temp)
		   temp=max(dif1,dif2,dif3,dif4) 
           if(temp > limit) then
              dqqy = MC_limiter((qq(i,j+1,k)-qq(i,j-1,k)) &
                   ,(qq(i,j+1,k)-qq(i,j,k)) &
                   ,(qq(i,j,k)-qq(i,j-1,k)))
              
              qqw(i,j,k,1) = qq(i,j,k) + 0.5d0*dqqy
              
              qqw(i,j-1,k,2) = qq(i,j,k) - 0.5d0*dqqy
            endif  
           end do
        end do
     end do
  else
!-----Step 3a.-------------------------------------------------|
! dqq
!
     do k=3,kx-2
        do j=2,jx-1
           do i=2,ix-1
				temp=dabs((ro(i,j,k+2)-ro(i,j,k+1))/dz(k+1))
				dif1=dlog10(temp) 
				temp=dabs((ro(i,j,k+1)-ro(i,j,k))/dz(k))
				dif2=dlog10(temp) 
				temp=dabs((ro(i,j,k)-ro(i,j,k-1))/dz(k-1))
				dif3=dlog10(temp) 
				temp=dabs((ro(i,j,k-1)-ro(i,j,k-2))/dz(k-2))
				dif4=dlog10(temp) 
		   temp=max(dif1,dif2,dif3,dif4) 
           if(temp > limit) then
              dqqz = MC_limiter((qq(i,j,k+1)-qq(i,j,k-1)) &
                   ,(qq(i,j,k+1)-qq(i,j,k)) &
                   ,(qq(i,j,k)-qq(i,j,k-1)))
           
!-----Step 3b.-------------------------------------------------|
! mdir = 3 :: z-direction
!
              qqw(i,j,k,1) = qq(i,j,k) + 0.5d0*dqqz
              
              qqw(i,j,k-1,2) = qq(i,j,k) - 0.5d0*dqqz
            endif
           end do
        end do
     end do
  endif
  return

contains
  function minmod_limiter(qqr,qql)
    implicit none
    
    real(8) :: qqr,qql
    real(8) :: minmod_limiter

    minmod_limiter = max(0.0d0, min(dabs(qqr),qql*sign(1.0d0,qqr)))*sign(1.0d0,qqr)
    return
  end function minmod_limiter

  function minmod_limiter2(qqr,qql)
    implicit none

    real(8) :: qqr, qql
    real(8) :: minmod_limiter2

    real(8) :: signr,signlr

    signr = sign(1.0d0,qqr)
    signlr = sign(1.0d0,(qqr*qql))

    minmod_limiter2 = max(signlr,0.0d0)*( &
         max(signr,0.0d0)*max(0.0d0,min(qql,qqr)) &
         -min(signr,0.0d0)*min(0.0d0,max(qql,qqr)))
!!$    if(qqr > 0.0d0)then
!!$       minmod_limiter2 = max(0.0d0,min(qql,qqr))
!!$    else
!!$       minmod_limiter2 = min(0.0d0,max(qql,qqr))
!!$    endif

    return
  end function minmod_limiter2
  
  function MC_limiter(qqc,qql,qqr)
    implicit none

    real(8) :: qqr,qql,qqc
    real(8) :: minmod_lr
    real(8) :: MC_limiter

    real(8) :: signlr,signMC
!!$    minmod_lr = minmod_limiter(2.0d0*qql,2.0d0*qqr)
!    MC_limiter = minmod_limiter2(0.5d0*qqc,minmod_lr)

    minmod_lr = minmod_limiter2(qql,qqr)

    signlr = sign(1.0d0,(qqc*minmod_lr))
    signMC = sign(1.0d0,(dabs(0.5d0*qqc)-dabs(2.0d0*minmod_lr)))

    MC_limiter = max(signlr,0.0d0)*( &
         max(signMC,0.0d0)*minmod_lr &
         -min(signMC,0.0d0)*(0.5d0*qqc))
    
!!$    if((qqc*minmod_lr) < 0.0d0)then
!!$       MC_limiter = 0.0d0
!!$    else
!!$       if(dabs(0.5d0*qqc) < dabs(2.0d0*minmod_lr))then
!!$          MC_limiter = 0.5d0*qqc
!!$       else
!!$          MC_limiter = minmod_lr
!!$       endif
!!$    endif
    return
  end function MC_limiter
end subroutine switchMP5toMC2_2