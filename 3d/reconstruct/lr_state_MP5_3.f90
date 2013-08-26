subroutine lr_state_MP5_3(mdir,ix,jx,kx,ro,pr &
     ,vx,vy,vz,bx,by,bz,phi &
     ,ch,gm,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw,ccx,ccy,ccz)
  implicit none

  integer,intent(in) :: ix,jx,kx,mdir
  real(8),intent(in) :: gm,ch

  integer,parameter :: nwave = 9

  real(8),dimension(ix,jx,kx),intent(in) :: ro,pr,vx,vy,vz,bx,by,bz,phi
  real(8),dimension(ix,jx,kx,2) :: row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw

  real(8),dimension(5,2,ix),intent(in) :: ccx
  real(8),dimension(5,2,jx),intent(in) :: ccy
  real(8),dimension(5,2,kx),intent(in) :: ccz

  real(8),dimension(nwave) :: wwc_w,qql,qqr
  ! ww(*,1) = qqm2, ..
  real(8),dimension(nwave,5) :: ww,wwc,tmpq
  real(8),dimension(nwave,nwave) :: lem,rem,tmpD1,tmpD2

  real(8) :: wwor
  real(8) :: ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1

!----parameter
  real(8),parameter :: B1 = 0.016666666667
  real(8),parameter :: B2 = 1.333333333333
  real(8),parameter :: Alpha = 4.0d0
  real(8),parameter :: Epsm = 0.0000000001d0

!----function
  integer :: flag

  integer :: i,j,k,l,m,n

!----function
  real(8) :: minmod4,d1,d2,d3,d4
  real(8) :: minmod,x,y
  real(8) :: median
  
  real(8) :: djm1,dj,djp1,dm4jph,dm4jmh,djpp1,dm4jpph
  real(8) :: qqul,qqav,qqmd,qqlc,qqmin,qqmax,qqmin1,qqmax1
  real(8) :: djp2,qqlr

!----esystem
  real(8) :: roi,btsq,vaxsq,asq
  real(8) :: ct2,tsum,tdif,cf2_cs2
  real(8) :: cfsq,cf,cssq,cs
  real(8) :: bt,bet2,bet3
  real(8) :: alpha_f,alpha_s

  real(8) :: na,qf,qs,af_prm,as_prm

  real(8) :: sqrtro,s,a,af,as
  real(8) :: vax

  real(8) :: temp,ich

  minmod4(d1,d2,d3,d4) = 0.125d0*(sign(1.0d0,d1)+sign(1.0d0,d2))* &
       dabs((sign(1.0d0,d1) + sign(1.0d0,d3))* &
       (sign(1.0d0,d1)+sign(1.0d0,d4))) &
       *min(dabs(d1),dabs(d2),dabs(d3),dabs(d4))

  minmod(x,y) = 0.5d0*(sign(1.0d0,x)+sign(1.0d0,y))*min(dabs(x),dabs(y))

  if(mdir .eq. 1)then
     do k=3,kx-2
        do j=3,jx-2
           do i=3,ix-2
              do n=1,5
                 ww(1,n) = ro(i-3+n,j,k)
                 ww(2,n) = vx(i-3+n,j,k)
                 ww(3,n) = vy(i-3+n,j,k)
                 ww(4,n) = vz(i-3+n,j,k)
                 ww(5,n) = pr(i-3+n,j,k)
                 ww(6,n) = bx(i-3+n,j,k)
                 ww(7,n) = by(i-3+n,j,k)
                 ww(8,n) = bz(i-3+n,j,k)
                 ww(9,n) = phi(i-3+n,j,k)
              end do
              ! left state
              
              ! primitive ro characteristic
              
              ro1 = ro(i,j,k)
              pr1 = pr(i,j,k)
              vx1 = vx(i,j,k)
              vy1 = vy(i,j,k)
              vz1 = vz(i,j,k)
              bx1 = bx(i,j,k)
              by1 = by(i,j,k)
              bz1 = bz(i,j,k)
              phi1 = phi(i,j,k)

              !======Esys
              roi = 1.0d0/ro1
              btsq = by1**2+bz1**2
              vaxsq = (bx1**2)*roi
              asq = gm*pr1*roi

              ! fast and slow speed
              ct2 = btsq*roi
              tsum = vaxsq + ct2 + asq
              tdif = vaxsq + ct2 - asq
              cf2_cs2 = sqrt(tdif**2 + 4.0d0*asq*ct2)
              
              cfsq = 0.5d0*(tsum + cf2_cs2)
              cf = sqrt(cfsq)
              
              cssq = asq*vaxsq/cfsq
              cs = sqrt(cssq)
              
              ! compute beta
              
              bt = sqrt(btsq)
              if ( bt .eq. 0.0d0) then
                 bet2 = 1.0d0/sqrt(2.0d0)
                 bet3 = 1.0d0/sqrt(2.0d0)
              else
                 bet2 = by1/bt
                 bet3 = bz1/bt
              endif
              
              ! compute alpha
              
              if(cf2_cs2 .eq. 0.0d0) then
                 alpha_f = 1.0d0
                 alpha_s = 0.0d0
              else if ((asq -cssq) <= 0.0d0)then
                 alpha_f = 0.0d0
                 alpha_s = 1.0d0
              else if ((cfsq - asq) <= 0.0d0)then
                 alpha_f = 1.0d0
                 alpha_s = 0.0d0
              else
                 alpha_f = sqrt((asq -cssq)/cf2_cs2)
                 alpha_s = sqrt((cfsq -asq)/cf2_cs2)
              endif
              
              ! compute Q(s),A(s)
              
              sqrtro = sqrt(ro1)
              s = sign(1.0d0,bx1)
              a = sqrt(asq)
              qf = cf*alpha_f*s
              qs = cs*alpha_s*s
              af = a*alpha_f*sqrtro
              as = a*alpha_s*sqrtro
              
              ! Right-eigenvector
              !  ich = 1.0d0/(ch)
              
              rem(1,1) = ro1*alpha_f
              rem(1,2) = 0.0d0
              rem(1,3) = ro1*alpha_s
              rem(1,4) = 1.0d0
              rem(1,5) = rem(1,3)
              rem(1,6) = 0.0d0
              rem(1,7) = 0.0d0
              rem(1,8) = rem(1,1)
              rem(1,9) = 0.0d0
              
              rem(2,1) = -cf*alpha_f
              rem(2,2) = 0.0d0
              rem(2,3) = -cs*alpha_s
              rem(2,4) = 0.0d0
              rem(2,5) = -rem(2,3)
              rem(2,6) = 0.0d0
              rem(2,7) = 0.0d0
              rem(2,8) = -rem(2,1)
              rem(2,9) = 0.0d0
              
              rem(3,1) = qs*bet2
              rem(3,2) = -bet3
              rem(3,3) = -qf*bet2
              rem(3,4) = 0.0d0
              rem(3,5) = -rem(3,3)
              rem(3,6) = 0.0d0
              rem(3,7) = bet3
              rem(3,8) = -rem(3,1)
              rem(3,9) = 0.0d0
              
              rem(4,1) = qs*bet3
              rem(4,2) = bet2
              rem(4,3) = -qf*bet3
              rem(4,4) = 0.0d0
              rem(4,5) = -rem(4,3)
              rem(4,6) = 0.0d0
              rem(4,7) = -bet2
              rem(4,8) = -rem(4,1)
              rem(4,9) = 0.0d0
              
              rem(5,1) = ro1*asq*alpha_f
              rem(5,2) = 0.0d0
              rem(5,3) = ro1*asq*alpha_s
              rem(5,4) = 0.0d0
              rem(5,5) = rem(5,3) 
              rem(5,6) = 0.0d0
              rem(5,7) = 0.0d0
              rem(5,8) = rem(5,1)
              rem(5,9) = 0.0d0
              
              rem(6,1) = 0.0d0
              rem(6,2) = 0.0d0
              rem(6,3) = 0.0d0
              rem(6,4) = 0.0d0
              rem(6,5) = 0.0d0
              rem(6,6) = 1.0d0
              rem(6,7) = 0.0d0
              rem(6,8) = 0.0d0
              rem(6,9) = 1.0d0
              
              rem(7,1) = as*bet2
              rem(7,2) = -bet3*s*sqrtro
              rem(7,3) = -af*bet2
              rem(7,4) = 0.0d0
              rem(7,5) = rem(7,3)
              rem(7,6) = 0.0d0
              rem(7,7) = rem(7,2) 
              rem(7,8) = rem(7,1)
              rem(7,9) = 0.0d0
              
              rem(8,1) = as*bet3
              rem(8,2) = bet2*s*sqrtro
              rem(8,3) = -af*bet3
              rem(8,4) = 0.0d0
              rem(8,5) = rem(8,3)
              rem(8,6) = 0.0d0
              rem(8,7) = rem(8,2)
              rem(8,8) = rem(8,1)
              rem(8,9) = 0.0d0
              
              rem(9,1) = 0.0d0
              rem(9,2) = 0.0d0
              rem(9,3) = 0.0d0
              rem(9,4) = 0.0d0
              rem(9,5) = 0.0d0
              rem(9,6) = -ch
              rem(9,7) = 0.0d0
              rem(9,8) = 0.0d0
              rem(9,9) = ch
              
              ! Left eigenvector
              
              na = 0.5d0/asq
              qf = na*qf
              qs = na*qs
              af_prm = na*af*roi
              as_prm = na*as*roi
              
              lem(1,1) = 0.0d0
              lem(1,2) = -na*cf*alpha_f
              lem(1,3) = qs*bet2
              lem(1,4) = qs*bet3
              lem(1,5) = na*alpha_f*roi
              lem(1,6) = 0.0d0
              lem(1,7) = as_prm*bet2
              lem(1,8) = as_prm*bet3
              lem(1,9) = 0.0d0
              
              lem(2,1) = 0.0d0
              lem(2,2) = 0.0d0
              lem(2,3) = -0.5d0*bet3
              lem(2,4) = 0.5d0*bet2
              lem(2,5) = 0.0d0
              lem(2,6) = 0.0d0
              lem(2,7) = -0.5d0*bet3*s/sqrtro
              lem(2,8) = 0.5d0*bet2*s/sqrtro
              lem(2,9) = 0.0d0
              
              lem(3,1) = 0.0d0
              lem(3,2) = -na*cs*alpha_s
              lem(3,3) = -qf*bet2
              lem(3,4) = -qf*bet3
              lem(3,5) = na*alpha_s*roi
              lem(3,6) = 0.0d0
              lem(3,7) = -af_prm*bet2
              lem(3,8) = -af_prm*bet3
              lem(3,9) = 0.0d0
              
              lem(4,1) = 1.0d0
              lem(4,2) = 0.0d0
              lem(4,3) = 0.0d0
              lem(4,4) = 0.0d0
              lem(4,5) = -1.0d0/asq
              lem(4,6) = 0.0d0
              lem(4,7) = 0.0d0
              lem(4,8) = 0.0d0
              lem(4,9) = 0.0d0
              
              lem(5,1) = 0.0d0
              lem(5,2) = -lem(3,2)
              lem(5,3) = -lem(3,3)
              lem(5,4) = -lem(3,4)
              lem(5,5) = lem(3,5) 
              lem(5,6) = 0.0d0
              lem(5,7) = lem(3,7)
              lem(5,8) = lem(3,8)
              lem(5,9) = 0.0d0
              
              lem(6,1) = 0.0d0
              lem(6,2) = 0.0d0
              lem(6,3) = 0.0d0
              lem(6,4) = 0.0d0
              lem(6,5) = 0.0d0
              lem(6,6) = 0.5d0
              lem(6,7) = 0.0d0
              lem(6,8) = 0.0d0
              lem(6,9) = -0.5d0/ch
              
              lem(7,1) = 0.0d0
              lem(7,2) = 0.0d0
              lem(7,3) = -lem(2,3)
              lem(7,4) = -lem(2,4)
              lem(7,5) = 0.0d0
              lem(7,6) = 0.0d0
              lem(7,7) = lem(2,7)
              lem(7,8) = lem(2,8)
              lem(7,9) = 0.0d0
              
              lem(8,1) = 0.0d0
              lem(8,2) = -lem(1,2)
              lem(8,3) = -lem(1,3)
              lem(8,4) = -lem(1,4)
              lem(8,5) = lem(1,5)
              lem(8,6) = 0.0d0
              lem(8,7) = lem(1,7)
              lem(8,8) = lem(1,8)
              lem(8,9) = 0.0d0
              
              lem(9,1) = 0.0d0
              lem(9,2) = 0.0d0
              lem(9,3) = 0.0d0
              lem(9,4) = 0.0d0
              lem(9,5) = 0.0d0
              lem(9,6) = 0.5d0
              lem(9,7) = 0.0d0
              lem(9,8) = 0.0d0
              lem(9,9) = 0.5d0/ch
              
              do l=1,5
                 do n=1,nwave
                    wwc(n,l) = lem(n,1)*ww(1,l)
                    do m=2,nwave
                       wwc(n,l) = wwc(n,l)+lem(n,m)*ww(m,l)
                    enddo
                 enddo
              end do
              
              ! mp5
              do n=1,nwave
                 wwor = B1*(ccx(1,2,i)*wwc(n,1)+ccx(2,2,i)*wwc(n,2) &
                      + ccx(3,2,i)*wwc(n,3) + ccx(4,2,i)*wwc(n,4) &
                      + ccx(5,2,i)*wwc(n,5))
                 djm1 = wwc(n,1)-2.0d0*wwc(n,2)+wwc(n,3)
                 dj = wwc(n,2)-2.0d0*wwc(n,3)+wwc(n,4)
                 djp1 = wwc(n,3)-2.0d0*wwc(n,4)+wwc(n,5)
                 
                 dm4jph = minmod4(4.0d0*dj-djp1,4.0d0*djp1-dj,dj,djp1)
                 dm4jmh = minmod4(4.0d0*dj-djm1,4.0d0*djm1-dj,dj,djm1)
                 
                 qqul = wwc(n,3)+Alpha*(wwc(n,3)-wwc(n,2))
                 qqlr = wwc(n,4)+Alpha*(wwc(n,4)-wwc(n,5))
                 
                 qqav = 0.5d0*(wwc(n,3)+wwc(n,4))
                 qqmd = qqav - 0.5d0*dm4jph
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,2)) + B2*dm4jmh
                 
                 qqmin = max(min(wwc(n,3),wwc(n,4),qqmd),min(wwc(n,3),qqul,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,4),qqmd),max(wwc(n,3),qqul,qqlc))

                 wwc_w(n) = wwor + minmod((qqmin-wwor),(qqmax-wwor))
              end do
              
              ! characteristic to primitive
              do n=1,nwave
                 qqr(n) = wwc_w(1)*rem(n,1)
                 do m=2,nwave
                    qqr(n) = qqr(n)+wwc_w(m)*rem(n,m)
                 enddo
              enddo
              
              ! right state
              ! mp5
              do n=1,nwave
                 wwor = B1*(ccx(5,1,i-1)*wwc(n,5)+ccx(4,1,i-1)*wwc(n,4) &
                      + ccx(3,1,i-1)*wwc(n,3) + ccx(2,1,i-1)*wwc(n,2) &
                      + ccx(1,1,i-1)*wwc(n,1))
                 djm1 = wwc(n,1)-2.0d0*wwc(n,2)+wwc(n,3)
                 dj = wwc(n,2)-2.0d0*wwc(n,3)+wwc(n,4)
                 djp1 = wwc(n,3)-2.0d0*wwc(n,4)+wwc(n,5)
                 
                 dm4jph = minmod4(4.0d0*dj-djp1,4.0d0*djp1-dj,dj,djp1)
                 dm4jmh = minmod4(4.0d0*dj-djm1,4.0d0*djm1-dj,dj,djm1)
                 
                 qqul = wwc(n,2)+Alpha*(wwc(n,2)-wwc(n,1))
                 qqlr = wwc(n,3)+Alpha*(wwc(n,3)-wwc(n,4))
                 
                 qqav = 0.5d0*(wwc(n,3)+wwc(n,2))
                 qqmd = qqav - 0.5d0*dm4jmh
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,4)) + B2*dm4jph
                 
                 qqmin = max(min(wwc(n,3),wwc(n,2),qqmd),min(wwc(n,3),qqlr,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,2),qqmd),max(wwc(n,3),qqlr,qqlc))
                 
                 wwc_w(n) = wwor + minmod((qqmin-wwor),(qqmax-wwor))
              end do
              
              ! characteristic to primitive
              do n=1,nwave
                 qql(n) = wwc_w(1)*rem(n,1)
                 do m=2,nwave
                    qql(n) = qql(n)+wwc_w(m)*rem(n,m)
                 enddo
              enddo

              qql(1) = max(min(ro(i,j,k),ro(i-1,j,k)),qql(1))
              qql(1) = min(max(ro(i,j,k),ro(i-1,j,k)),qql(1))
              qqr(1) = max(min(ro(i,j,k),ro(i+1,j,k)),qqr(1))
              qqr(1) = min(max(ro(i,j,k),ro(i+1,j,k)),qqr(1))
              
              qql(2) = max(min(vx(i,j,k),vx(i-1,j,k)),qql(2))
              qql(2) = min(max(vx(i,j,k),vx(i-1,j,k)),qql(2))
              qqr(2) = max(min(vx(i,j,k),vx(i+1,j,k)),qqr(2))
              qqr(2) = min(max(vx(i,j,k),vx(i+1,j,k)),qqr(2))
              
              qql(3) = max(min(vy(i,j,k),vy(i-1,j,k)),qql(3))
              qql(3) = min(max(vy(i,j,k),vy(i-1,j,k)),qql(3))
              qqr(3) = max(min(vy(i,j,k),vy(i+1,j,k)),qqr(3))
              qqr(3) = min(max(vy(i,j,k),vy(i+1,j,k)),qqr(3))
              
              qql(4) = max(min(vz(i,j,k),vz(i-1,j,k)),qql(4))
              qql(4) = min(max(vz(i,j,k),vz(i-1,j,k)),qql(4))
              qqr(4) = max(min(vz(i,j,k),vz(i+1,j,k)),qqr(4))
              qqr(4) = min(max(vz(i,j,k),vz(i+1,j,k)),qqr(4))
              
              qql(5) = max(min(pr(i,j,k),pr(i-1,j,k)),qql(5))
              qql(5) = min(max(pr(i,j,k),pr(i-1,j,k)),qql(5))
              qqr(5) = max(min(pr(i,j,k),pr(i+1,j,k)),qqr(5))
              qqr(5) = min(max(pr(i,j,k),pr(i+1,j,k)),qqr(5))
              
              qql(6) = max(min(bx(i,j,k),bx(i-1,j,k)),qql(6))
              qql(6) = min(max(bx(i,j,k),bx(i-1,j,k)),qql(6))
              qqr(6) = max(min(bx(i,j,k),bx(i+1,j,k)),qqr(6))
              qqr(6) = min(max(bx(i,j,k),bx(i+1,j,k)),qqr(6))
              
              qql(7) = max(min(by(i,j,k),by(i-1,j,k)),qql(7))
              qql(7) = min(max(by(i,j,k),by(i-1,j,k)),qql(7))
              qqr(7) = max(min(by(i,j,k),by(i+1,j,k)),qqr(7))
              qqr(7) = min(max(by(i,j,k),by(i+1,j,k)),qqr(7))
              
              qql(8) = max(min(bz(i,j,k),bz(i-1,j,k)),qql(8))
              qql(8) = min(max(bz(i,j,k),bz(i-1,j,k)),qql(8))
              qqr(8) = max(min(bz(i,j,k),bz(i+1,j,k)),qqr(8))
              qqr(8) = min(max(bz(i,j,k),bz(i+1,j,k)),qqr(8))
              
              qql(9) = max(min(phi(i,j,k),phi(i-1,j,k)),qql(9))
              qql(9) = min(max(phi(i,j,k),phi(i-1,j,k)),qql(9))
              qqr(9) = max(min(phi(i,j,k),phi(i+1,j,k)),qqr(9))
              qqr(9) = min(max(phi(i,j,k),phi(i+1,j,k)),qqr(9))
              
              row(i-1,j,k,2) = qql(1)
              row(i,j,k,1) = qqr(1)
              
              vxw(i-1,j,k,2) = qql(2)
              vxw(i,j,k,1) = qqr(2)
              vyw(i-1,j,k,2) = qql(3)
              vyw(i,j,k,1) = qqr(3)
              vzw(i-1,j,k,2) = qql(4)
              vzw(i,j,k,1) = qqr(4)
              
              prw(i-1,j,k,2) = qql(5)
              prw(i,j,k,1) = qqr(5)
              
              bxw(i-1,j,k,2) = qql(6)
              bxw(i,j,k,1) = qqr(6)
              byw(i-1,j,k,2) = qql(7)
              byw(i,j,k,1) = qqr(7)
              bzw(i-1,j,k,2) = qql(8)
              bzw(i,j,k,1) = qqr(8)
              
              phiw(i-1,j,k,2) = qql(9)
              phiw(i,j,k,1) = qqr(9)
           end do
        end do
     end do
  else if(mdir .eq. 2)then
     do k=3,kx-2
        do j=3,jx-2
           do i=3,ix-2
              do n=1,5
                 ww(1,n) = ro(i,j-3+n,k)
                 ww(2,n) = vx(i,j-3+n,k)
                 ww(3,n) = vy(i,j-3+n,k)
                 ww(4,n) = vz(i,j-3+n,k)
                 ww(5,n) = pr(i,j-3+n,k)
                 ww(6,n) = bx(i,j-3+n,k)
                 ww(7,n) = by(i,j-3+n,k)
                 ww(8,n) = bz(i,j-3+n,k)
                 ww(9,n) = phi(i,j-3+n,k)
              end do
              ! left state
              
              ! primitive ro characteristic
              
              ro1 = ro(i,j,k)
              pr1 = pr(i,j,k)
              vx1 = vx(i,j,k)
              vy1 = vy(i,j,k)
              vz1 = vz(i,j,k)
              bx1 = bx(i,j,k)
              by1 = by(i,j,k)
              bz1 = bz(i,j,k)
              phi1 = phi(i,j,k)
              
              !======Esys
              roi = 1.0d0/ro1
              btsq = by1**2+bz1**2
              vaxsq = (bx1**2)*roi
              asq = gm*pr1*roi

              ! fast and slow speed
              ct2 = btsq*roi
              tsum = vaxsq + ct2 + asq
              tdif = vaxsq + ct2 - asq
              cf2_cs2 = sqrt(tdif**2 + 4.0d0*asq*ct2)
              
              cfsq = 0.5d0*(tsum + cf2_cs2)
              cf = sqrt(cfsq)
              
              cssq = asq*vaxsq/cfsq
              cs = sqrt(cssq)
              
              ! compute beta
              
              bt = sqrt(btsq)
              if ( bt .eq. 0.0d0) then
                 bet2 = 1.0d0/sqrt(2.0d0)
                 bet3 = 1.0d0/sqrt(2.0d0)
              else
                 bet2 = by1/bt
                 bet3 = bz1/bt
              endif
              
              ! compute alpha
              
              if(cf2_cs2 .eq. 0.0d0) then
                 alpha_f = 1.0d0
                 alpha_s = 0.0d0
              else if ((asq -cssq) <= 0.0d0)then
                 alpha_f = 0.0d0
                 alpha_s = 1.0d0
              else if ((cfsq - asq) <= 0.0d0)then
                 alpha_f = 1.0d0
                 alpha_s = 0.0d0
              else
                 alpha_f = sqrt((asq -cssq)/cf2_cs2)
                 alpha_s = sqrt((cfsq -asq)/cf2_cs2)
              endif
              
              ! compute Q(s),A(s)
              
              sqrtro = sqrt(ro1)
              s = sign(1.0d0,bx1)
              a = sqrt(asq)
              qf = cf*alpha_f*s
              qs = cs*alpha_s*s
              af = a*alpha_f*sqrtro
              as = a*alpha_s*sqrtro
              
              ! Right-eigenvector
              !  ich = 1.0d0/(ch)
              
              rem(1,1) = ro1*alpha_f
              rem(1,2) = 0.0d0
              rem(1,3) = ro1*alpha_s
              rem(1,4) = 1.0d0
              rem(1,5) = rem(1,3)
              rem(1,6) = 0.0d0
              rem(1,7) = 0.0d0
              rem(1,8) = rem(1,1)
              rem(1,9) = 0.0d0
              
              rem(2,1) = -cf*alpha_f
              rem(2,2) = 0.0d0
              rem(2,3) = -cs*alpha_s
              rem(2,4) = 0.0d0
              rem(2,5) = -rem(2,3)
              rem(2,6) = 0.0d0
              rem(2,7) = 0.0d0
              rem(2,8) = -rem(2,1)
              rem(2,9) = 0.0d0
              
              rem(3,1) = qs*bet2
              rem(3,2) = -bet3
              rem(3,3) = -qf*bet2
              rem(3,4) = 0.0d0
              rem(3,5) = -rem(3,3)
              rem(3,6) = 0.0d0
              rem(3,7) = bet3
              rem(3,8) = -rem(3,1)
              rem(3,9) = 0.0d0
              
              rem(4,1) = qs*bet3
              rem(4,2) = bet2
              rem(4,3) = -qf*bet3
              rem(4,4) = 0.0d0
              rem(4,5) = -rem(4,3)
              rem(4,6) = 0.0d0
              rem(4,7) = -bet2
              rem(4,8) = -rem(4,1)
              rem(4,9) = 0.0d0
              
              rem(5,1) = ro1*asq*alpha_f
              rem(5,2) = 0.0d0
              rem(5,3) = ro1*asq*alpha_s
              rem(5,4) = 0.0d0
              rem(5,5) = rem(5,3) 
              rem(5,6) = 0.0d0
              rem(5,7) = 0.0d0
              rem(5,8) = rem(5,1)
              rem(5,9) = 0.0d0
              
              rem(6,1) = 0.0d0
              rem(6,2) = 0.0d0
              rem(6,3) = 0.0d0
              rem(6,4) = 0.0d0
              rem(6,5) = 0.0d0
              rem(6,6) = 1.0d0
              rem(6,7) = 0.0d0
              rem(6,8) = 0.0d0
              rem(6,9) = 1.0d0
              
              rem(7,1) = as*bet2
              rem(7,2) = -bet3*s*sqrtro
              rem(7,3) = -af*bet2
              rem(7,4) = 0.0d0
              rem(7,5) = rem(7,3)
              rem(7,6) = 0.0d0
              rem(7,7) = rem(7,2) 
              rem(7,8) = rem(7,1)
              rem(7,9) = 0.0d0
              
              rem(8,1) = as*bet3
              rem(8,2) = bet2*s*sqrtro
              rem(8,3) = -af*bet3
              rem(8,4) = 0.0d0
              rem(8,5) = rem(8,3)
              rem(8,6) = 0.0d0
              rem(8,7) = rem(8,2)
              rem(8,8) = rem(8,1)
              rem(8,9) = 0.0d0
              
              rem(9,1) = 0.0d0
              rem(9,2) = 0.0d0
              rem(9,3) = 0.0d0
              rem(9,4) = 0.0d0
              rem(9,5) = 0.0d0
              rem(9,6) = -ch
              rem(9,7) = 0.0d0
              rem(9,8) = 0.0d0
              rem(9,9) = ch
              
              ! Left eigenvector
              
              na = 0.5d0/asq
              qf = na*qf
              qs = na*qs
              af_prm = na*af*roi
              as_prm = na*as*roi
              
              lem(1,1) = 0.0d0
              lem(1,2) = -na*cf*alpha_f
              lem(1,3) = qs*bet2
              lem(1,4) = qs*bet3
              lem(1,5) = na*alpha_f*roi
              lem(1,6) = 0.0d0
              lem(1,7) = as_prm*bet2
              lem(1,8) = as_prm*bet3
              lem(1,9) = 0.0d0
              
              lem(2,1) = 0.0d0
              lem(2,2) = 0.0d0
              lem(2,3) = -0.5d0*bet3
              lem(2,4) = 0.5d0*bet2
              lem(2,5) = 0.0d0
              lem(2,6) = 0.0d0
              lem(2,7) = -0.5d0*bet3*s/sqrtro
              lem(2,8) = 0.5d0*bet2*s/sqrtro
              lem(2,9) = 0.0d0
              
              lem(3,1) = 0.0d0
              lem(3,2) = -na*cs*alpha_s
              lem(3,3) = -qf*bet2
              lem(3,4) = -qf*bet3
              lem(3,5) = na*alpha_s*roi
              lem(3,6) = 0.0d0
              lem(3,7) = -af_prm*bet2
              lem(3,8) = -af_prm*bet3
              lem(3,9) = 0.0d0
              
              lem(4,1) = 1.0d0
              lem(4,2) = 0.0d0
              lem(4,3) = 0.0d0
              lem(4,4) = 0.0d0
              lem(4,5) = -1.0d0/asq
              lem(4,6) = 0.0d0
              lem(4,7) = 0.0d0
              lem(4,8) = 0.0d0
              lem(4,9) = 0.0d0
              
              lem(5,1) = 0.0d0
              lem(5,2) = -lem(3,2)
              lem(5,3) = -lem(3,3)
              lem(5,4) = -lem(3,4)
              lem(5,5) = lem(3,5) 
              lem(5,6) = 0.0d0
              lem(5,7) = lem(3,7)
              lem(5,8) = lem(3,8)
              lem(5,9) = 0.0d0
              
              lem(6,1) = 0.0d0
              lem(6,2) = 0.0d0
              lem(6,3) = 0.0d0
              lem(6,4) = 0.0d0
              lem(6,5) = 0.0d0
              lem(6,6) = 0.5d0
              lem(6,7) = 0.0d0
              lem(6,8) = 0.0d0
              lem(6,9) = -0.5d0/ch
              
              lem(7,1) = 0.0d0
              lem(7,2) = 0.0d0
              lem(7,3) = -lem(2,3)
              lem(7,4) = -lem(2,4)
              lem(7,5) = 0.0d0
              lem(7,6) = 0.0d0
              lem(7,7) = lem(2,7)
              lem(7,8) = lem(2,8)
              lem(7,9) = 0.0d0
              
              lem(8,1) = 0.0d0
              lem(8,2) = -lem(1,2)
              lem(8,3) = -lem(1,3)
              lem(8,4) = -lem(1,4)
              lem(8,5) = lem(1,5)
              lem(8,6) = 0.0d0
              lem(8,7) = lem(1,7)
              lem(8,8) = lem(1,8)
              lem(8,9) = 0.0d0
              
              lem(9,1) = 0.0d0
              lem(9,2) = 0.0d0
              lem(9,3) = 0.0d0
              lem(9,4) = 0.0d0
              lem(9,5) = 0.0d0
              lem(9,6) = 0.5d0
              lem(9,7) = 0.0d0
              lem(9,8) = 0.0d0
              lem(9,9) = 0.5d0/ch

              
              do l=1,5
                 do n=1,nwave
                    wwc(n,l) = lem(n,1)*ww(1,l)
                    do m=2,nwave
                       wwc(n,l) = wwc(n,l)+lem(n,m)*ww(m,l)
                    enddo
                 enddo
              end do
              
              ! mp5
              do n=1,nwave
                 wwor = B1*(ccy(1,2,j)*wwc(n,1)+ccy(2,2,j)*wwc(n,2) &
                      + ccy(3,2,j)*wwc(n,3) + ccy(4,2,j)*wwc(n,4) &
                      + ccy(5,2,j)*wwc(n,5))
                 djm1 = wwc(n,1)-2.0d0*wwc(n,2)+wwc(n,3)
                 dj = wwc(n,2)-2.0d0*wwc(n,3)+wwc(n,4)
                 djp1 = wwc(n,3)-2.0d0*wwc(n,4)+wwc(n,5)
                 
                 dm4jph = minmod4(4.0d0*dj-djp1,4.0d0*djp1-dj,dj,djp1)
                 dm4jmh = minmod4(4.0d0*dj-djm1,4.0d0*djm1-dj,dj,djm1)
                 
                 qqul = wwc(n,3)+Alpha*(wwc(n,3)-wwc(n,2))
                 qqlr = wwc(n,4)+Alpha*(wwc(n,4)-wwc(n,5))
                 
                 qqav = 0.5d0*(wwc(n,3)+wwc(n,4))
                 qqmd = qqav - 0.5d0*dm4jph
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,2)) + B2*dm4jmh
                 
                 qqmin = max(min(wwc(n,3),wwc(n,4),qqmd),min(wwc(n,3),qqul,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,4),qqmd),max(wwc(n,3),qqul,qqlc))

                 wwc_w(n) = wwor + minmod((qqmin-wwor),(qqmax-wwor))

              end do
              
              ! characteristic to primitive
              do n=1,nwave
                 qqr(n) = wwc_w(1)*rem(n,1)
                 do m=2,nwave
                    qqr(n) = qqr(n)+wwc_w(m)*rem(n,m)
                 enddo
              enddo
              ! right state
              ! mp5
              do n=1,nwave
                 wwor = B1*(ccy(5,1,j-1)*wwc(n,5)+ccy(4,1,j-1)*wwc(n,4) &
                      + ccy(3,1,j-1)*wwc(n,3) + ccy(2,1,j-1)*wwc(n,2) &
                      + ccy(1,1,j-1)*wwc(n,1))
                 djm1 = wwc(n,1)-2.0d0*wwc(n,2)+wwc(n,3)
                 dj = wwc(n,2)-2.0d0*wwc(n,3)+wwc(n,4)
                 djp1 = wwc(n,3)-2.0d0*wwc(n,4)+wwc(n,5)
                 
                 dm4jph = minmod4(4.0d0*dj-djp1,4.0d0*djp1-dj,dj,djp1)
                 dm4jmh = minmod4(4.0d0*dj-djm1,4.0d0*djm1-dj,dj,djm1)
                 
                 qqul = wwc(n,2)+Alpha*(wwc(n,2)-wwc(n,1))
                 qqlr = wwc(n,3)+Alpha*(wwc(n,3)-wwc(n,4))
                 
                 qqav = 0.5d0*(wwc(n,3)+wwc(n,2))
                 qqmd = qqav - 0.5d0*dm4jmh
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,4)) + B2*dm4jph
                 
                 qqmin = max(min(wwc(n,3),wwc(n,2),qqmd),min(wwc(n,3),qqlr,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,2),qqmd),max(wwc(n,3),qqlr,qqlc))
                 
                 wwc_w(n) = wwor + minmod((qqmin-wwor),(qqmax-wwor))
              end do
              
              ! characteristic to primitive
              do n=1,nwave
                 qql(n) = wwc_w(1)*rem(n,1)
                 do m=2,nwave
                    qql(n) = qql(n)+wwc_w(m)*rem(n,m)
                 enddo
              enddo

              qql(1) = max(min(ro(i,j,k),ro(i,j-1,k)),qql(1))
              qql(1) = min(max(ro(i,j,k),ro(i,j-1,k)),qql(1))
              qqr(1) = max(min(ro(i,j,k),ro(i,j+1,k)),qqr(1))
              qqr(1) = min(max(ro(i,j,k),ro(i,j+1,k)),qqr(1))
              
              qql(2) = max(min(vx(i,j,k),vx(i,j-1,k)),qql(2))
              qql(2) = min(max(vx(i,j,k),vx(i,j-1,k)),qql(2))
              qqr(2) = max(min(vx(i,j,k),vx(i,j+1,k)),qqr(2))
              qqr(2) = min(max(vx(i,j,k),vx(i,j+1,k)),qqr(2))
              
              qql(3) = max(min(vy(i,j,k),vy(i,j-1,k)),qql(3))
              qql(3) = min(max(vy(i,j,k),vy(i,j-1,k)),qql(3))
              qqr(3) = max(min(vy(i,j,k),vy(i,j+1,k)),qqr(3))
              qqr(3) = min(max(vy(i,j,k),vy(i,j+1,k)),qqr(3))
              
              qql(4) = max(min(vz(i,j,k),vz(i,j-1,k)),qql(4))
              qql(4) = min(max(vz(i,j,k),vz(i,j-1,k)),qql(4))
              qqr(4) = max(min(vz(i,j,k),vz(i,j+1,k)),qqr(4))
              qqr(4) = min(max(vz(i,j,k),vz(i,j+1,k)),qqr(4))
              
              qql(5) = max(min(pr(i,j,k),pr(i,j-1,k)),qql(5))
              qql(5) = min(max(pr(i,j,k),pr(i,j-1,k)),qql(5))
              qqr(5) = max(min(pr(i,j,k),pr(i,j+1,k)),qqr(5))
              qqr(5) = min(max(pr(i,j,k),pr(i,j+1,k)),qqr(5))
              
              qql(6) = max(min(bx(i,j,k),bx(i,j-1,k)),qql(6))
              qql(6) = min(max(bx(i,j,k),bx(i,j-1,k)),qql(6))
              qqr(6) = max(min(bx(i,j,k),bx(i,j+1,k)),qqr(6))
              qqr(6) = min(max(bx(i,j,k),bx(i,j+1,k)),qqr(6))
              
              qql(7) = max(min(by(i,j,k),by(i,j-1,k)),qql(7))
              qql(7) = min(max(by(i,j,k),by(i,j-1,k)),qql(7))
              qqr(7) = max(min(by(i,j,k),by(i,j+1,k)),qqr(7))
              qqr(7) = min(max(by(i,j,k),by(i,j+1,k)),qqr(7))
              
              qql(8) = max(min(bz(i,j,k),bz(i,j-1,k)),qql(8))
              qql(8) = min(max(bz(i,j,k),bz(i,j-1,k)),qql(8))
              qqr(8) = max(min(bz(i,j,k),bz(i,j+1,k)),qqr(8))
              qqr(8) = min(max(bz(i,j,k),bz(i,j+1,k)),qqr(8))
              
              qql(9) = max(min(phi(i,j,k),phi(i,j-1,k)),qql(9))
              qql(9) = min(max(phi(i,j,k),phi(i,j-1,k)),qql(9))
              qqr(9) = max(min(phi(i,j,k),phi(i,j+1,k)),qqr(9))
              qqr(9) = min(max(phi(i,j,k),phi(i,j+1,k)),qqr(9))
              
              row(i,j-1,k,2) = qql(1)
              row(i,j,k,1) = qqr(1)
              
              vxw(i,j-1,k,2) = qql(2)
              vxw(i,j,k,1) = qqr(2)
              vyw(i,j-1,k,2) = qql(3)
              vyw(i,j,k,1) = qqr(3)
              vzw(i,j-1,k,2) = qql(4)
              vzw(i,j,k,1) = qqr(4)
              
              prw(i,j-1,k,2) = qql(5)
              prw(i,j,k,1) = qqr(5)
              
              bxw(i,j-1,k,2) = qql(6)
              bxw(i,j,k,1) = qqr(6)
              byw(i,j-1,k,2) = qql(7)
              byw(i,j,k,1) = qqr(7)
              bzw(i,j-1,k,2) = qql(8)
              bzw(i,j,k,1) = qqr(8)
              
              phiw(i,j-1,k,2) = qql(9)
              phiw(i,j,k,1) = qqr(9)
           end do
        end do
     end do
  else
     do k=3,kx-2
        do j=3,jx-2
           do i=3,ix-2
              do n=1,5
                 ww(1,n) = ro(i,j,k-3+n)
                 ww(2,n) = vx(i,j,k-3+n)
                 ww(3,n) = vy(i,j,k-3+n)
                 ww(4,n) = vz(i,j,k-3+n)
                 ww(5,n) = pr(i,j,k-3+n)
                 ww(6,n) = bx(i,j,k-3+n)
                 ww(7,n) = by(i,j,k-3+n)
                 ww(8,n) = bz(i,j,k-3+n)
                 ww(9,n) = phi(i,j,k-3+n)
              end do
              ! left state
              
              ! primitive ro characteristic
              
              ro1 = ro(i,j,k)
              pr1 = pr(i,j,k)
              vx1 = vx(i,j,k)
              vy1 = vy(i,j,k)
              vz1 = vz(i,j,k)
              bx1 = bx(i,j,k)
              by1 = by(i,j,k)
              bz1 = bz(i,j,k)
              phi1 = phi(i,j,k)
              
              !======Esys
              roi = 1.0d0/ro1
              btsq = by1**2+bz1**2
              vaxsq = (bx1**2)*roi
              asq = gm*pr1*roi

              ! fast and slow speed
              ct2 = btsq*roi
              tsum = vaxsq + ct2 + asq
              tdif = vaxsq + ct2 - asq
              cf2_cs2 = sqrt(tdif**2 + 4.0d0*asq*ct2)
              
              cfsq = 0.5d0*(tsum + cf2_cs2)
              cf = sqrt(cfsq)
              
              cssq = asq*vaxsq/cfsq
              cs = sqrt(cssq)
              
              ! compute beta
              
              bt = sqrt(btsq)
              if ( bt .eq. 0.0d0) then
                 bet2 = 1.0d0/sqrt(2.0d0)
                 bet3 = 1.0d0/sqrt(2.0d0)
              else
                 bet2 = by1/bt
                 bet3 = bz1/bt
              endif
              
              ! compute alpha
              
              if(cf2_cs2 .eq. 0.0d0) then
                 alpha_f = 1.0d0
                 alpha_s = 0.0d0
              else if ((asq -cssq) <= 0.0d0)then
                 alpha_f = 0.0d0
                 alpha_s = 1.0d0
              else if ((cfsq - asq) <= 0.0d0)then
                 alpha_f = 1.0d0
                 alpha_s = 0.0d0
              else
                 alpha_f = sqrt((asq -cssq)/cf2_cs2)
                 alpha_s = sqrt((cfsq -asq)/cf2_cs2)
              endif
              
              ! compute Q(s),A(s)
              
              sqrtro = sqrt(ro1)
              s = sign(1.0d0,bx1)
              a = sqrt(asq)
              qf = cf*alpha_f*s
              qs = cs*alpha_s*s
              af = a*alpha_f*sqrtro
              as = a*alpha_s*sqrtro
              
              ! Right-eigenvector
              !  ich = 1.0d0/(ch)
              
              rem(1,1) = ro1*alpha_f
              rem(1,2) = 0.0d0
              rem(1,3) = ro1*alpha_s
              rem(1,4) = 1.0d0
              rem(1,5) = rem(1,3)
              rem(1,6) = 0.0d0
              rem(1,7) = 0.0d0
              rem(1,8) = rem(1,1)
              rem(1,9) = 0.0d0
              
              rem(2,1) = -cf*alpha_f
              rem(2,2) = 0.0d0
              rem(2,3) = -cs*alpha_s
              rem(2,4) = 0.0d0
              rem(2,5) = -rem(2,3)
              rem(2,6) = 0.0d0
              rem(2,7) = 0.0d0
              rem(2,8) = -rem(2,1)
              rem(2,9) = 0.0d0
              
              rem(3,1) = qs*bet2
              rem(3,2) = -bet3
              rem(3,3) = -qf*bet2
              rem(3,4) = 0.0d0
              rem(3,5) = -rem(3,3)
              rem(3,6) = 0.0d0
              rem(3,7) = bet3
              rem(3,8) = -rem(3,1)
              rem(3,9) = 0.0d0
              
              rem(4,1) = qs*bet3
              rem(4,2) = bet2
              rem(4,3) = -qf*bet3
              rem(4,4) = 0.0d0
              rem(4,5) = -rem(4,3)
              rem(4,6) = 0.0d0
              rem(4,7) = -bet2
              rem(4,8) = -rem(4,1)
              rem(4,9) = 0.0d0
              
              rem(5,1) = ro1*asq*alpha_f
              rem(5,2) = 0.0d0
              rem(5,3) = ro1*asq*alpha_s
              rem(5,4) = 0.0d0
              rem(5,5) = rem(5,3) 
              rem(5,6) = 0.0d0
              rem(5,7) = 0.0d0
              rem(5,8) = rem(5,1)
              rem(5,9) = 0.0d0
              
              rem(6,1) = 0.0d0
              rem(6,2) = 0.0d0
              rem(6,3) = 0.0d0
              rem(6,4) = 0.0d0
              rem(6,5) = 0.0d0
              rem(6,6) = 1.0d0
              rem(6,7) = 0.0d0
              rem(6,8) = 0.0d0
              rem(6,9) = 1.0d0
              
              rem(7,1) = as*bet2
              rem(7,2) = -bet3*s*sqrtro
              rem(7,3) = -af*bet2
              rem(7,4) = 0.0d0
              rem(7,5) = rem(7,3)
              rem(7,6) = 0.0d0
              rem(7,7) = rem(7,2) 
              rem(7,8) = rem(7,1)
              rem(7,9) = 0.0d0
              
              rem(8,1) = as*bet3
              rem(8,2) = bet2*s*sqrtro
              rem(8,3) = -af*bet3
              rem(8,4) = 0.0d0
              rem(8,5) = rem(8,3)
              rem(8,6) = 0.0d0
              rem(8,7) = rem(8,2)
              rem(8,8) = rem(8,1)
              rem(8,9) = 0.0d0
              
              rem(9,1) = 0.0d0
              rem(9,2) = 0.0d0
              rem(9,3) = 0.0d0
              rem(9,4) = 0.0d0
              rem(9,5) = 0.0d0
              rem(9,6) = -ch
              rem(9,7) = 0.0d0
              rem(9,8) = 0.0d0
              rem(9,9) = ch
              
              ! Left eigenvector
              
              na = 0.5d0/asq
              qf = na*qf
              qs = na*qs
              af_prm = na*af*roi
              as_prm = na*as*roi
              
              lem(1,1) = 0.0d0
              lem(1,2) = -na*cf*alpha_f
              lem(1,3) = qs*bet2
              lem(1,4) = qs*bet3
              lem(1,5) = na*alpha_f*roi
              lem(1,6) = 0.0d0
              lem(1,7) = as_prm*bet2
              lem(1,8) = as_prm*bet3
              lem(1,9) = 0.0d0
              
              lem(2,1) = 0.0d0
              lem(2,2) = 0.0d0
              lem(2,3) = -0.5d0*bet3
              lem(2,4) = 0.5d0*bet2
              lem(2,5) = 0.0d0
              lem(2,6) = 0.0d0
              lem(2,7) = -0.5d0*bet3*s/sqrtro
              lem(2,8) = 0.5d0*bet2*s/sqrtro
              lem(2,9) = 0.0d0
              
              lem(3,1) = 0.0d0
              lem(3,2) = -na*cs*alpha_s
              lem(3,3) = -qf*bet2
              lem(3,4) = -qf*bet3
              lem(3,5) = na*alpha_s*roi
              lem(3,6) = 0.0d0
              lem(3,7) = -af_prm*bet2
              lem(3,8) = -af_prm*bet3
              lem(3,9) = 0.0d0
              
              lem(4,1) = 1.0d0
              lem(4,2) = 0.0d0
              lem(4,3) = 0.0d0
              lem(4,4) = 0.0d0
              lem(4,5) = -1.0d0/asq
              lem(4,6) = 0.0d0
              lem(4,7) = 0.0d0
              lem(4,8) = 0.0d0
              lem(4,9) = 0.0d0
              
              lem(5,1) = 0.0d0
              lem(5,2) = -lem(3,2)
              lem(5,3) = -lem(3,3)
              lem(5,4) = -lem(3,4)
              lem(5,5) = lem(3,5) 
              lem(5,6) = 0.0d0
              lem(5,7) = lem(3,7)
              lem(5,8) = lem(3,8)
              lem(5,9) = 0.0d0
              
              lem(6,1) = 0.0d0
              lem(6,2) = 0.0d0
              lem(6,3) = 0.0d0
              lem(6,4) = 0.0d0
              lem(6,5) = 0.0d0
              lem(6,6) = 0.5d0
              lem(6,7) = 0.0d0
              lem(6,8) = 0.0d0
              lem(6,9) = -0.5d0/ch
              
              lem(7,1) = 0.0d0
              lem(7,2) = 0.0d0
              lem(7,3) = -lem(2,3)
              lem(7,4) = -lem(2,4)
              lem(7,5) = 0.0d0
              lem(7,6) = 0.0d0
              lem(7,7) = lem(2,7)
              lem(7,8) = lem(2,8)
              lem(7,9) = 0.0d0
              
              lem(8,1) = 0.0d0
              lem(8,2) = -lem(1,2)
              lem(8,3) = -lem(1,3)
              lem(8,4) = -lem(1,4)
              lem(8,5) = lem(1,5)
              lem(8,6) = 0.0d0
              lem(8,7) = lem(1,7)
              lem(8,8) = lem(1,8)
              lem(8,9) = 0.0d0
              
              lem(9,1) = 0.0d0
              lem(9,2) = 0.0d0
              lem(9,3) = 0.0d0
              lem(9,4) = 0.0d0
              lem(9,5) = 0.0d0
              lem(9,6) = 0.5d0
              lem(9,7) = 0.0d0
              lem(9,8) = 0.0d0
              lem(9,9) = 0.5d0/ch
              
              do l=1,5
                 do n=1,nwave
                    wwc(n,l) = lem(n,1)*ww(1,l)
                    do m=2,nwave
                       wwc(n,l) = wwc(n,l)+lem(n,m)*ww(m,l)
                    enddo
                 enddo
              end do
              
              do l=1,5
                 do n=1,nwave
                    wwc(n,l) = lem(n,1)*ww(1,l)
                    do m=2,nwave
                       wwc(n,l) = wwc(n,l)+lem(n,m)*ww(m,l)
                    enddo
                 enddo
              end do
              
              ! mp5
              do n=1,nwave
                 wwor = B1*(ccz(1,2,k)*wwc(n,1)+ccz(2,2,k)*wwc(n,2) &
                      + ccz(3,2,k)*wwc(n,3) + ccz(4,2,k)*wwc(n,4) &
                      + ccz(5,2,k)*wwc(n,5))
                 djm1 = wwc(n,1)-2.0d0*wwc(n,2)+wwc(n,3)
                 dj = wwc(n,2)-2.0d0*wwc(n,3)+wwc(n,4)
                 djp1 = wwc(n,3)-2.0d0*wwc(n,4)+wwc(n,5)
                 
                 dm4jph = minmod4(4.0d0*dj-djp1,4.0d0*djp1-dj,dj,djp1)
                 dm4jmh = minmod4(4.0d0*dj-djm1,4.0d0*djm1-dj,dj,djm1)
                 
                 qqul = wwc(n,3)+Alpha*(wwc(n,3)-wwc(n,2))
                 qqlr = wwc(n,4)+Alpha*(wwc(n,4)-wwc(n,5))
                 
                 qqav = 0.5d0*(wwc(n,3)+wwc(n,4))
                 qqmd = qqav - 0.5d0*dm4jph
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,2)) + B2*dm4jmh
                 
                 qqmin = max(min(wwc(n,3),wwc(n,4),qqmd),min(wwc(n,3),qqul,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,4),qqmd),max(wwc(n,3),qqul,qqlc))

                 wwc_w(n) = wwor + minmod((qqmin-wwor),(qqmax-wwor))
              end do
              
              ! characteristic to primitive
              do n=1,nwave
                 qqr(n) = wwc_w(1)*rem(n,1)
                 do m=2,nwave
                    qqr(n) = qqr(n)+wwc_w(m)*rem(n,m)
                 enddo
              enddo
              ! right state
              ! mp5
              do n=1,nwave
                 wwor = B1*(ccz(5,1,k-1)*wwc(n,5)+ccz(4,1,k-1)*wwc(n,4) &
                      + ccz(3,1,k-1)*wwc(n,3) + ccz(2,1,k-1)*wwc(n,2) &
                      + ccz(1,1,k-1)*wwc(n,1))
                 djm1 = wwc(n,1)-2.0d0*wwc(n,2)+wwc(n,3)
                 dj = wwc(n,2)-2.0d0*wwc(n,3)+wwc(n,4)
                 djp1 = wwc(n,3)-2.0d0*wwc(n,4)+wwc(n,5)
                 
                 dm4jph = minmod4(4.0d0*dj-djp1,4.0d0*djp1-dj,dj,djp1)
                 dm4jmh = minmod4(4.0d0*dj-djm1,4.0d0*djm1-dj,dj,djm1)
                 
                 qqul = wwc(n,2)+Alpha*(wwc(n,2)-wwc(n,1))
                 qqlr = wwc(n,3)+Alpha*(wwc(n,3)-wwc(n,4))
                 
                 qqav = 0.5d0*(wwc(n,3)+wwc(n,2))
                 qqmd = qqav - 0.5d0*dm4jmh
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,4)) + B2*dm4jph
                 
                 qqmin = max(min(wwc(n,3),wwc(n,2),qqmd),min(wwc(n,3),qqlr,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,2),qqmd),max(wwc(n,3),qqlr,qqlc))
                 
                 wwc_w(n) = wwor + minmod((qqmin-wwor),(qqmax-wwor))
              end do
              
              ! characteristic to primitive
              do n=1,nwave
                 qql(n) = wwc_w(1)*rem(n,1)
                 do m=2,nwave
                    qql(n) = qql(n)+wwc_w(m)*rem(n,m)
                 enddo
              enddo

              qql(1) = max(min(ro(i,j,k),ro(i,j,k-1)),qql(1))
              qql(1) = min(max(ro(i,j,k),ro(i,j,k-1)),qql(1))
              qqr(1) = max(min(ro(i,j,k),ro(i,j,k+1)),qqr(1))
              qqr(1) = min(max(ro(i,j,k),ro(i,j,k+1)),qqr(1))
              
              qql(2) = max(min(vx(i,j,k),vx(i,j,k-1)),qql(2))
              qql(2) = min(max(vx(i,j,k),vx(i,j,k-1)),qql(2))
              qqr(2) = max(min(vx(i,j,k),vx(i,j,k+1)),qqr(2))
              qqr(2) = min(max(vx(i,j,k),vx(i,j,k+1)),qqr(2))
              
              qql(3) = max(min(vy(i,j,k),vy(i,j,k-1)),qql(3))
              qql(3) = min(max(vy(i,j,k),vy(i,j,k-1)),qql(3))
              qqr(3) = max(min(vy(i,j,k),vy(i,j,k+1)),qqr(3))
              qqr(3) = min(max(vy(i,j,k),vy(i,j,k+1)),qqr(3))
              
              qql(4) = max(min(vz(i,j,k),vz(i,j,k-1)),qql(4))
              qql(4) = min(max(vz(i,j,k),vz(i,j,k-1)),qql(4))
              qqr(4) = max(min(vz(i,j,k),vz(i,j,k+1)),qqr(4))
              qqr(4) = min(max(vz(i,j,k),vz(i,j,k+1)),qqr(4))
              
              qql(5) = max(min(pr(i,j,k),pr(i,j,k-1)),qql(5))
              qql(5) = min(max(pr(i,j,k),pr(i,j,k-1)),qql(5))
              qqr(5) = max(min(pr(i,j,k),pr(i,j,k+1)),qqr(5))
              qqr(5) = min(max(pr(i,j,k),pr(i,j,k+1)),qqr(5))
              
              qql(6) = max(min(bx(i,j,k),bx(i,j,k-1)),qql(6))
              qql(6) = min(max(bx(i,j,k),bx(i,j,k-1)),qql(6))
              qqr(6) = max(min(bx(i,j,k),bx(i,j,k+1)),qqr(6))
              qqr(6) = min(max(bx(i,j,k),bx(i,j,k+1)),qqr(6))
              
              qql(7) = max(min(by(i,j,k),by(i,j,k-1)),qql(7))
              qql(7) = min(max(by(i,j,k),by(i,j,k-1)),qql(7))
              qqr(7) = max(min(by(i,j,k),by(i,j,k+1)),qqr(7))
              qqr(7) = min(max(by(i,j,k),by(i,j,k+1)),qqr(7))
              
              qql(8) = max(min(bz(i,j,k),bz(i,j,k-1)),qql(8))
              qql(8) = min(max(bz(i,j,k),bz(i,j,k-1)),qql(8))
              qqr(8) = max(min(bz(i,j,k),bz(i,j,k+1)),qqr(8))
              qqr(8) = min(max(bz(i,j,k),bz(i,j,k+1)),qqr(8))
              
              qql(9) = max(min(phi(i,j,k),phi(i,j,k-1)),qql(9))
              qql(9) = min(max(phi(i,j,k),phi(i,j,k-1)),qql(9))
              qqr(9) = max(min(phi(i,j,k),phi(i,j,k+1)),qqr(9))
              qqr(9) = min(max(phi(i,j,k),phi(i,j,k+1)),qqr(9))
              
              row(i,j,k-1,2) = qql(1)
              row(i,j,k,1) = qqr(1)
              
              vxw(i,j,k-1,2) = qql(2)
              vxw(i,j,k,1) = qqr(2)
              vyw(i,j,k-1,2) = qql(3)
              vyw(i,j,k,1) = qqr(3)
              vzw(i,j,k-1,2) = qql(4)
              vzw(i,j,k,1) = qqr(4)
              
              prw(i,j,k-1,2) = qql(5)
              prw(i,j,k,1) = qqr(5)
              
              bxw(i,j,k-1,2) = qql(6)
              bxw(i,j,k,1) = qqr(6)
              byw(i,j,k-1,2) = qql(7)
              byw(i,j,k,1) = qqr(7)
              bzw(i,j,k-1,2) = qql(8)
              bzw(i,j,k,1) = qqr(8)
              
              phiw(i,j,k-1,2) = qql(9)
              phiw(i,j,k,1) = qqr(9)
           end do
        end do
     end do
  end if
  return
end subroutine lr_state_MP5_3