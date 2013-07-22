subroutine TimeChangeVal(margin,ix,jx,kx,x,y,z,xObs &
     ,ro,pb,pr,rx,ry,rz,bx,by,bz &
     ,baseTCR,roTCR,pbTCR,prTCR,rxTCR,ryTCR,rzTCR,bxTCR,byTCR,bzTCR)
  use ModTimeChangeAtRad
  implicit none
  
  integer,intent(in) :: margin,ix,jx,kx
  
  real(8),dimension(ix),intent(in) :: x
  real(8),dimension(jx),intent(in) :: y
  real(8),dimension(kx),intent(in) :: z

  real(8),dimension(ix,jx,kx),intent(in) :: ro,pb,pr
  real(8),dimension(ix,jx,kx),intent(in) :: rx,ry,rz
  real(8),dimension(ix,jx,kx),intent(in) :: bx,by,bz
  
  real(8),intent(in) :: xObs

  real(8),dimension(kx) :: baseTCR

  real(8),dimension(kx) :: bxTCR,byTCR,bzTCR
  real(8),dimension(kx) :: rxTCR,ryTCR,rzTCR
  real(8),dimension(kx) :: roTCR,pbTCR,prTCR

  integer :: k

  call calSumAtRadius(margin,ix,jx,kx,x,y,z,ro,xObs &
       ,roTCR)
  call calSumAtRadius(margin,ix,jx,kx,x,y,z,pb,xObs &
       ,pbTCR)
  call calSumAtRadius(margin,ix,jx,kx,x,y,z,pr,xObs &
       ,prTCR)
  
  call calSumAtRadius(margin,ix,jx,kx,x,y,z,rx,xObs &
       ,rxTCR)
  call calSumAtRadius(margin,ix,jx,kx,x,y,z,ry,xObs &
       ,ryTCR)
  call calSumAtRadius(margin,ix,jx,kx,x,y,z,rz,xObs &
       ,rzTCR)
  
  call calSumAtRadius(margin,ix,jx,kx,x,y,z,bx,xObs &
       ,bxTCR)
  call calSumAtRadius(margin,ix,jx,kx,x,y,z,by,xObs &
       ,byTCR)
  call calSumAtRadius(margin,ix,jx,kx,x,y,z,bz,xObs &
       ,bzTCR)
  
  call calSumAtRadiusOfBase(margin,ix,jx,kx,x,y,z,xObs &
       ,baseTCR)
  
  do k=1,kx
     if(baseTCR(k) .ne. 0.0d0)then
        roTCR(k) = roTCR(k)/baseTCR(k)
        prTCR(k) = prTCR(k)/baseTCR(k)
        pbTCR(k) = pbTCR(k)/baseTCR(k)
        
        bxTCR(k) = bxTCR(k)/baseTCR(k)
        byTCR(k) = byTCR(k)/baseTCR(k)
        bzTCR(k) = bzTCR(k)/baseTCR(k)
        rxTCR(k) = rxTCR(k)/baseTCR(k)
        ryTCR(k) = ryTCR(k)/baseTCR(k)
        rzTCR(k) = rzTCR(k)/baseTCR(k)
     endif
  enddo
  
  return
end subroutine TimeChangeVal
