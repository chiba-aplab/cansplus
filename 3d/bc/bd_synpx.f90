! Cell surface boundary 
! 1 | 2 | 3 | 4 | 5 | 6 | ..
! q(3) = q(4), q(2) = q(5) , q(1) = q(6)
subroutine bd_synpx(mbnd,margin,qq,ix,jx,kx)

  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  integer,intent(in) :: mbnd

  real(8),dimension(ix,jx,kx) :: qq

  integer :: jpPi,jnPi,hjx,jxmm,jxm
  integer :: i,j,k
  real(8) :: signj

  jxm = jx-margin
  jxmm = jx-2*margin
  hjx = (jx-2*margin)/2

  if(mbnd .eq. 0)then
     do k=1,kx
        do j=1,jx
           do i=1,margin
              jpPi = j+hjx
              jnPi = j-hjx

              signj = sign(1.0d0,jxm-jpPi+0.00001d0)

              jpPi = min(jx,jpPi)
              jpPi = max(1,jpPi)
              jnPi = max(1,jnPi)
              jnPi = min(jx,jnPi)

              qq(margin-i+1,j,k) = max(0.0d0,signj)*(qq(margin+i,jpPi,k)) &
                   +max(0.0d0,-signj)*(qq(margin+i,jnPi,k))
           enddo
        end do
     end do
  else
     do k=1,kx
        do j=1,jx
           do i=1,margin
              qq(ix-margin+i,j,k) = qq(ix-margin-i+1,j,k)
           end do
        end do
     end do
  end if

  return
end subroutine bd_synpx
