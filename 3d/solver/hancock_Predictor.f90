subroutine hancock_Predictor(mdir,ix,jx,kx,dt,gm,ch,dx,dy,dz &
     ,roHDiff,prHdiff,vxHDiff,vyHDiff,vzHDiff &
     ,bxHDiff,byHDiff,bzHDiff,phiHDiff,ro,pr,vx,vy,vz,bx,by,bz,phi &
     ,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw)

  implicit none

  integer,intent(in) :: ix,jx,kx
  integer,intent(in) :: mdir

  real(8),intent(in) :: gm,ch,dt

  real(8),dimension(ix) :: dx
  real(8),dimension(jx) :: dy
  real(8),dimension(kx) :: dz

  real(8),dimension(ix,jx,kx),intent(in) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(in) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(in) :: bx,by,bz,phi

  real(8),dimension(ix,jx,kx),intent(in) :: roHDiff,prHDiff
  real(8),dimension(ix,jx,kx),intent(in) :: vxHDiff,vyHDiff,vzHDiff
  real(8),dimension(ix,jx,kx),intent(in) :: bxHDiff,byHDiff,bzHDiff
  real(8),dimension(ix,jx,kx),intent(in) :: phiHDiff

  real(8),dimension(ix,jx,kx,2),intent(inout) :: row,prw,vxw,vyw,vzw
  real(8),dimension(ix,jx,kx,2),intent(inout) :: bxw,byw,bzw,phiw

  integer,parameter :: waveNum = 9

  real(8),dimension(waveNum,waveNum) :: em
  real(8),dimension(waveNum) :: dw,ww,wwl,wwr,wl,wr,wwh
  
  real(8) :: hdtodx,hdtody,hdtodz
  integer :: i,j,k,l,m

  if(mdir .eq. 1)then
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
!-----Step 0.----------------------------------------------------------|
! cal. eigenMatrix
!
              call calGLMEigenMatrix(waveNum,gm,ch,ro(i,j,k) &
                   ,vx(i,j,k),vy(i,j,k),vz(i,j,k) &
                   ,bx(i,j,k),by(i,j,k),bz(i,j,k) &
                   ,pr(i,j,k),phi(i,j,k) &
                   ,em)
! set dw
              dw(1) = roHDiff(i,j,k)
              dw(2) = vxHDiff(i,j,k)
              dw(3) = vyHDiff(i,j,k)
              dw(4) = vzHDiff(i,j,k)

              dw(5) = bxHDiff(i,j,k)
              dw(6) = byHDiff(i,j,k)
              dw(7) = bzHDiff(i,j,k)

              dw(8) = prHDiff(i,j,k)
              dw(9) = phiHDiff(i,j,k)

!-----Step 1.----------------------------------------------------------|
! cal. predictor step
!
              hdtodx = 0.5d0*dt/dx(i)

              ww(1) = ro(i,j,k)
              ww(2) = vx(i,j,k)
              ww(3) = vy(i,j,k)
              ww(4) = vz(i,j,k)
              ww(5) = bx(i,j,k)
              ww(6) = by(i,j,k)
              ww(7) = bz(i,j,k)
              ww(8) = pr(i,j,k)
              ww(9) = phi(i,j,k)

              wwl(1) = ro(i-1,j,k)
              wwl(2) = vx(i-1,j,k)
              wwl(3) = vy(i-1,j,k)
              wwl(4) = vz(i-1,j,k)
              wwl(5) = bx(i-1,j,k)
              wwl(6) = by(i-1,j,k)
              wwl(7) = bz(i-1,j,k)
              wwl(8) = pr(i-1,j,k)
              wwl(9) = phi(i-1,j,k)

              wwr(1) = ro(i+1,j,k)
              wwr(2) = vx(i+1,j,k)
              wwr(3) = vy(i+1,j,k)
              wwr(4) = vz(i+1,j,k)
              wwr(5) = bx(i+1,j,k)
              wwr(6) = by(i+1,j,k)
              wwr(7) = bz(i+1,j,k)
              wwr(8) = pr(i+1,j,k)
              wwr(9) = phi(i+1,j,k)

              do l=1,waveNum
                 do m=1,waveNum
                    wwh(l) = ww(l)-hdtodx*em(l,m)*dw(m)
                 end do
              end do

!-----Step 2.----------------------------------------------------------|
! cal. surface value
!
              do l=1,waveNum
                 wr(l) = wwh(l) + 0.5d0*dw(l)
                 wl(l) = wwh(l) - 0.5d0*dw(l)
              end do

              do l=1,waveNum
                 wl(l) = max(min(ww(l),wwl(l)),wl(l))
                 wl(l) = min(max(ww(l),wwl(l)),wl(l))
                 
                 wr(l) = max(min(ww(l),wwr(l)),wr(l))
                 wr(l) = min(max(ww(l),wwr(l)),wr(l))
              end do

              row(i,j,k,1) = wr(1)
              row(i-1,j,k,2) = wl(1)

              vxw(i,j,k,1) = wr(2)
              vxw(i-1,j,k,2) = wl(2)
              vyw(i,j,k,1) = wr(3)
              vyw(i-1,j,k,2) = wl(3)
              vzw(i,j,k,1) = wr(4)
              vzw(i-1,j,k,2) = wl(4)

              bxw(i,j,k,1) = wr(5)
              bxw(i-1,j,k,2) = wl(5)
              byw(i,j,k,1) = wr(6)
              byw(i-1,j,k,2) = wl(6)
              bzw(i,j,k,1) = wr(7)
              bzw(i-1,j,k,2) = wl(7)

              prw(i,j,k,1) = wr(8)
              prw(i-1,j,k,2) = wl(8) 
              phiw(i,j,k,1) = wr(9)
              phiw(i-1,j,k,2) = wl(9)
           end do
        end do
     end do
  else if(mdir .eq. 2)then
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
!-----Step 0.----------------------------------------------------------|
! cal. eigenMatrix
!
              call calGLMEigenMatrix(waveNum,gm,ch,ro(i,j,k) &
                   ,vy(i,j,k),vz(i,j,k),vx(i,j,k) &
                   ,by(i,j,k),bz(i,j,k),bx(i,j,k) &
                   ,pr(i,j,k),phi(i,j,k) &
                   ,em)
! set dw
              dw(1) = roHDiff(i,j,k)

              dw(2) = vyHDiff(i,j,k)
              dw(3) = vzHDiff(i,j,k)
              dw(4) = vxHDiff(i,j,k)

              dw(5) = byHDiff(i,j,k)
              dw(6) = bzHDiff(i,j,k)
              dw(7) = bxHDiff(i,j,k)

              dw(8) = prHDiff(i,j,k)
              dw(9) = phiHDiff(i,j,k)

!-----Step 1.----------------------------------------------------------|
! cal. predictor step
!
              hdtody = 0.5d0*dt/dy(j)

              ww(1) = ro(i,j,k)
              ww(2) = vy(i,j,k)
              ww(3) = vz(i,j,k)
              ww(4) = vx(i,j,k)
              ww(5) = by(i,j,k)
              ww(6) = bz(i,j,k)
              ww(7) = bx(i,j,k)
              ww(8) = pr(i,j,k)
              ww(9) = phi(i,j,k)

              wwl(1) = ro(i,j-1,k)
              wwl(2) = vy(i,j-1,k)
              wwl(3) = vz(i,j-1,k)
              wwl(4) = vx(i,j-1,k)
              wwl(5) = by(i,j-1,k)
              wwl(6) = bz(i,j-1,k)
              wwl(7) = bx(i,j-1,k)
              wwl(8) = pr(i,j-1,k)
              wwl(9) = phi(i,j-1,k)

              wwr(1) = ro(i,j+1,k)
              wwr(2) = vy(i,j+1,k)
              wwr(3) = vz(i,j+1,k)
              wwr(4) = vx(i,j+1,k)
              wwr(5) = by(i,j+1,k)
              wwr(6) = bz(i,j+1,k)
              wwr(7) = bx(i,j+1,k)
              wwr(8) = pr(i,j+1,k)
              wwr(9) = phi(i,j+1,k)

              do l=1,waveNum
                 do m=1,waveNum
                    wwh(l) = ww(l)-hdtody*em(l,m)*dw(m)
                 end do
              end do

!-----Step 2.----------------------------------------------------------|
! cal. surface value
!
              do l=1,waveNum
                 wr(l) = wwh(l) + 0.5d0*dw(l)
                 wl(l) = wwh(l) - 0.5d0*dw(l)
              end do

              do l=1,waveNum
                 wl(l) = max(min(ww(l),wwl(l)),wl(l))
                 wl(l) = min(max(ww(l),wwl(l)),wl(l))
                 
                 wr(l) = max(min(ww(l),wwr(l)),wr(l))
                 wr(l) = min(max(ww(l),wwr(l)),wr(l))
              end do

              row(i,j,k,1) = wr(1)
              row(i,j-1,k,2) = wl(1)

              vyw(i,j,k,1) = wr(2)
              vyw(i,j-1,k,2) = wl(2)
              vzw(i,j,k,1) = wr(3)
              vzw(i,j-1,k,2) = wl(3)
              vxw(i,j,k,1) = wr(4) 
              vxw(i,j-1,k,2) = wl(4)

              byw(i,j,k,1) = wr(5)
              byw(i,j-1,k,2) = wl(5)
              bzw(i,j,k,1) = wr(6)
              bzw(i,j-1,k,2) = wl(6)
              bxw(i,j,k,1) = wr(7)
              bxw(i,j-1,k,2) = wl(7)

              prw(i,j,k,1) = wr(8)
              prw(i,j-1,k,2) = wl(8)
              phiw(i,j,k,1) = wr(9)
              phiw(i,j-1,k,2) = wl(9)
           end do
        end do
     end do
  else
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
!-----Step 0.----------------------------------------------------------|
! cal. eigenMatrix
!
              call calGLMEigenMatrix(waveNum,gm,ch,ro(i,j,k) &
                   ,vz(i,j,k),vx(i,j,k),vy(i,j,k) &
                   ,bz(i,j,k),bx(i,j,k),by(i,j,k) &
                   ,pr(i,j,k),phi(i,j,k) &
                   ,em)
! set dw
              dw(1) = roHDiff(i,j,k)

              dw(2) = vzHDiff(i,j,k)
              dw(3) = vxHDiff(i,j,k)
              dw(4) = vyHDiff(i,j,k)

              dw(5) = bzHDiff(i,j,k)
              dw(6) = bxHDiff(i,j,k)
              dw(7) = byHDiff(i,j,k)

              dw(8) = prHDiff(i,j,k)
              dw(9) = phiHDiff(i,j,k)

!-----Step 1.----------------------------------------------------------|
! cal. predictor step
!
              hdtodz = 0.5d0*dt/dz(k)

              ww(1) = ro(i,j,k)
              ww(2) = vz(i,j,k)
              ww(3) = vx(i,j,k)
              ww(4) = vy(i,j,k)
              ww(5) = bz(i,j,k)
              ww(6) = bx(i,j,k)
              ww(7) = by(i,j,k)
              ww(8) = pr(i,j,k)
              ww(9) = phi(i,j,k)

              wwl(1) = ro(i,j,k-1)
              wwl(2) = vz(i,j,k-1)
              wwl(3) = vx(i,j,k-1)
              wwl(4) = vy(i,j,k-1)
              wwl(5) = bz(i,j,k-1)
              wwl(6) = bx(i,j,k-1)
              wwl(7) = by(i,j,k-1)
              wwl(8) = pr(i,j,k-1)
              wwl(9) = phi(i,j,k-1)

              wwr(1) = ro(i,j,k+1)
              wwr(2) = vz(i,j,k+1)
              wwr(3) = vx(i,j,k+1)
              wwr(4) = vy(i,j,k+1)
              wwr(5) = bz(i,j,k+1)
              wwr(6) = bx(i,j,k+1)
              wwr(7) = by(i,j,k+1)
              wwr(8) = pr(i,j,k+1)
              wwr(9) = phi(i,j,k+1)

              do l=1,waveNum
                 do m=1,waveNum
                    wwh(l) = ww(l)-hdtodz*em(l,m)*dw(m)
                 end do
              end do

!-----Step 2.----------------------------------------------------------|
! cal. surface value
!
              do l=1,waveNum
                 wr(l) = wwh(l) + 0.5d0*dw(l)
                 wl(l) = wwh(l) - 0.5d0*dw(l)
              end do

              do l=1,waveNum
                 wl(l) = max(min(ww(l),wwl(l)),wl(l))
                 wl(l) = min(max(ww(l),wwl(l)),wl(l))
                 
                 wr(l) = max(min(ww(l),wwr(l)),wr(l))
                 wr(l) = min(max(ww(l),wwr(l)),wr(l))
              end do

              row(i,j,k,1) = wr(1)
              row(i,j,k-1,2) = wl(1)

              vzw(i,j,k,1) = wr(2)
              vzw(i,j,k-1,2) = wl(2)
              vxw(i,j,k,1) = wr(3)
              vxw(i,j,k-1,2) = wl(3)
              vyw(i,j,k,1) = wr(4)
              vyw(i,j,k-1,2) = wl(4)

              bzw(i,j,k,1) = wr(5)
              bzw(i,j,k-1,2) = wl(5)
              bxw(i,j,k,1) = wr(6)
              bxw(i,j,k-1,2) = wl(6)
              byw(i,j,k,1) = wr(7)
              byw(i,j,k-1,2) = wl(7)

              prw(i,j,k,1) = wr(8)
              prw(i,j,k-1,2) = wl(8)
              phiw(i,j,k,1) = wr(9)
              phiw(i,j,k-1,2) = wl(9)
           end do
        end do
     end do
  end if
  return 
contains
  subroutine calGLMEigenMatrix(waveNum,gm,ch,ro,vx,vy,vz,bx,by,bz,pr,phi &
       ,em)
    implicit none

    integer,intent(in) :: waveNum
    real(8),intent(in) :: gm,ch

    real(8),intent(in) :: ro,pr
    real(8),intent(in) :: vx,vy,vz
    real(8),intent(in) :: bx,by,bz,phi

    real(8),dimension(waveNum,waveNum),intent(inout) :: em

    integer :: i,j
    do j=1,waveNum
       do i=1,waveNum
          em(i,j) = 0.0d0
       end do
    end do

    em(1,1) = vx
    em(1,2) = ro
    
    em(2,2) = vx
    em(2,5) = -bx/ro
    em(2,6) = by/ro
    em(2,7) = bz/ro
    em(2,8) = 1.0d0/ro

    em(3,3) = vx
    em(3,5) = -by/ro
    em(3,6) = -bx/ro

    em(4,4) = vx
    em(4,5) = -bz/ro
    em(4,7) = -bx/ro

    em(5,9) = 1.0d0

    em(6,2) = by
    em(6,3) = -bx
    em(6,5) = -vy
    em(6,6) = vx
    
    em(7,2) = bz
    em(7,4) = -bx
    em(7,5) = -vz
    em(7,7) = vx

    em(8,2) = gm*pr
    em(8,5) = (gm-1.0d0)*(vx*bx+vy*by+vz*bz)
    em(8,8) = vx
    em(8,9) = (1.0d0-gm)*bx
    
    em(9,5) = ch**2

    return
  end subroutine calGLMEigenMatrix
end subroutine hancock_Predictor
