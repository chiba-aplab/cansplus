module lr_state

  implicit none
  private

  public :: lr_state__1st, lr_state__2nd, lr_state__MP5, reconstructionConstant


contains


  subroutine lr_state__1st(mdir,ix,jx,kx,qq,qqw)

  integer,intent(in) :: mdir
  integer,intent(in) :: ix,jx,kx
  real(8),dimension(ix,jx,kx), intent(in) :: qq
  real(8),dimension(ix,jx,kx,2), intent(out) :: qqw

  integer :: i,j,k
  real(8) :: dqqx,dqqy,dqqz

  if(mdir == 1)then
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
              qqw(i,j,k,1) = qq(i,j,k)
              qqw(i-1,j,k,2) = qq(i,j,k)
           enddo
        enddo
     enddo
  else if(mdir == 2)then
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

  end subroutine lr_state__1st


  subroutine lr_state__2nd(mdir,ix,jx,kx,dx,dy,dz,qq,qqw)
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

  integer,intent(in) :: mdir
  integer,intent(in) :: ix,jx,kx
  real(8),dimension(ix),intent(in) :: dx
  real(8),dimension(jx),intent(in) :: dy
  real(8),dimension(kx),intent(in) :: dz
  real(8),dimension(ix,jx,kx),intent(in) :: qq
  real(8),dimension(ix,jx,kx,2), intent(out) :: qqw

  integer :: i,j,k
  real(8) :: dqqx,dqqy,dqqz
  real(8) :: dqc,dql,dqr
  real(8) :: dxc,dxl,dxr
  real(8) :: dyc,dyl,dyr
  real(8) :: dzc,dzl,dzr

!-----Step 1a.-------------------------------------------------|
! dqq
!
!-----Step 1b.-------------------------------------------------|
! mdir = 1 :: x-direction
  if(mdir == 1)then
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
              dxl = 0.5d0*(dx(i)+dx(i-1))
              dxr = 0.5d0*(dx(i+1)+dx(i))
              dxc = dxl+dxr

              dql = (qq(i,j,k)-qq(i-1,j,k))/dxl
              dqr = (qq(i+1,j,k)-qq(i,j,k))/dxr
              dqc = (qq(i+1,j,k)-qq(i-1,j,k))/dxc

              dqqx = MC_limiter(dqc,dqr,dql)
              
              qqw(i,j,k,1) = qq(i,j,k) + 0.5d0*dqqx*dx(i)
              qqw(i-1,j,k,2) = qq(i,j,k) - 0.5d0*dqqx*dx(i)
           enddo
        enddo
     enddo
!-----Step 1b.-------------------------------------------------|
! mdir = 2 :: y-direction
  else if(mdir == 2)then
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
              dql = qq(i,j,k)-qq(i,j-1,k)
              dqr = qq(i,j+1,k)-qq(i,j,k)
              dqc = 0.5d0*(qq(i,j+1,k)-qq(i,j-1,k))

              dqqy = MC_limiter(dqc,dqr,dql)
              
              qqw(i,j,k,1) = qq(i,j,k) + 0.5d0*dqqy
              qqw(i,j-1,k,2) = qq(i,j,k) - 0.5d0*dqqy
           end do
        end do
     end do
  else
!-----Step 3a.-------------------------------------------------|
! dqq
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
              dzl = 0.5d0*(dz(k)+dz(k-1))
              dzr = 0.5d0*(dz(k+1)+dz(k))
              dzc = dzl+dzr

              dql = (qq(i,j,k)-qq(i,j,k-1))/dzl
              dqr = (qq(i,j,k+1)-qq(i,j,k))/dzr
              dqc = (qq(i,j,k+1)-qq(i,j,k-1))/dzc

              dqqz = MC_limiter(dqc,dqr,dql)
           
              qqw(i,j,k,1) = qq(i,j,k) + 0.5d0*dqqz*dz(k)
              qqw(i,j,k-1,2) = qq(i,j,k) - 0.5d0*dqqz*dz(k)
           end do
        end do
     end do
  endif

  end subroutine lr_state__2nd


  subroutine lr_state__MP5(mdir,ix,jx,kx,ro,pr &
                          ,vx,vy,vz,bx,by,bz,phi &
                          ,ch,gm,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw,ccx,ccy,ccz)

  integer,intent(in) :: ix,jx,kx,mdir
  real(8),intent(in) :: gm,ch
  real(8),dimension(5,2,ix),intent(in) :: ccx
  real(8),dimension(5,2,jx),intent(in) :: ccy
  real(8),dimension(5,2,kx),intent(in) :: ccz
  real(8),dimension(ix,jx,kx),intent(in) :: ro,pr,vx,vy,vz,bx,by,bz,phi
  real(8),dimension(ix,jx,kx,2),intent(out) :: row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw

  integer,parameter :: nwave = 9
  real(8),dimension(nwave) :: qql,qqr
  real(8),dimension(2,nwave) :: wwc_w
  ! ww(*,1) = qqm2, ..
  real(8),dimension(nwave,5) :: ww,wwc
  real(8),dimension(nwave,nwave) :: lem,rem
  real(8) :: wwor,minvalue,smv,psmv,msmv
  real(8) :: ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1
  real(8) :: temp,temp1,temp2

!----parameter
  real(8),parameter :: B1 = 0.016666666667
  real(8),parameter :: B2 = 1.333333333333
  real(8),parameter :: Alpha = 4.0d0
  real(8),parameter :: Epsm = 0.0000000001d0

!----function
  integer :: i,j,k,l,m,n
  real(8) :: minmod4,d1,d2,d3,d4
  real(8) :: minmod,x,y
  real(8) :: djm1,dj,djp1,dm4jph,dm4jmh,djpp1,dm4jpph
  real(8) :: qqul,qqav,qqmd,qqlc,qqmin,qqmax,qqmin1,qqmax1
  real(8) :: djp2,qqlr

!----Recipe for rarefaction ---
  real(8), parameter :: floor=1d-2
  real(8) :: romaxvalue,prmaxvalue

  minmod4(d1,d2,d3,d4) = 0.125d0*(sign(1.0d0,d1)+sign(1.0d0,d2))* &
       dabs((sign(1.0d0,d1) + sign(1.0d0,d3))* &
       (sign(1.0d0,d1)+sign(1.0d0,d4))) &
       *min(dabs(d1),dabs(d2),dabs(d3),dabs(d4))

  minmod(x,y) = 0.5d0*(sign(1.0d0,x)+sign(1.0d0,y))*min(dabs(x),dabs(y))
  
  if(mdir == 1)then

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
              ro1 = ww(1,3)
              vx1 = ww(2,3)
              vy1 = ww(3,3)
              vz1 = ww(4,3)
              pr1 = ww(5,3)
              bx1 = ww(6,3)
              by1 = ww(7,3)
              bz1 = ww(8,3)
              phi1 = ww(9,3)

              ! primitive to characteristic
              call esystem_glmmhd(lem,rem,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1,ch,gm)
              wwc = 0d0
              do l=1,5
                 do m=1,nwave
                    temp = ww(m,l)
                    do n=1,nwave
                       wwc(n,l) = wwc(n,l)+lem(n,m)*temp
                    enddo
                 enddo
              end do
              
              ! mp5
              do n=1,nwave
                 !left state
                 wwor = B1*(ccx(1,2,i)*wwc(n,1)+ccx(2,2,i)*wwc(n,2) &
                      + ccx(3,2,i)*wwc(n,3) + ccx(4,2,i)*wwc(n,4) &
                      + ccx(5,2,i)*wwc(n,5))

                 djm1 = wwc(n,1)-2.0d0*wwc(n,2)+wwc(n,3)
                 dj = wwc(n,2)-2.0d0*wwc(n,3)+wwc(n,4)
                 djp1 = wwc(n,3)-2.0d0*wwc(n,4)+wwc(n,5)
                 dm4jph = minmod4(4.0d0*dj-djp1,4.0d0*djp1-dj,dj,djp1)
                 dm4jmh = minmod4(4.0d0*dj-djm1,4.0d0*djm1-dj,dj,djm1)

                 qqul = wwc(n,3)+Alpha*(wwc(n,3)-wwc(n,2))
                 qqav = 0.5d0*(wwc(n,3)+wwc(n,4))
                 qqmd = qqav - 0.5d0*dm4jph
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,2)) + B2*dm4jmh
                 
                 qqmin = max(min(wwc(n,3),wwc(n,4),qqmd),min(wwc(n,3),qqul,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,4),qqmd),max(wwc(n,3),qqul,qqlc))
                 wwc_w(1,n) = wwor + minmod((qqmin-wwor),(qqmax-wwor))

                 !right state
                 wwor = B1*(ccx(5,1,i-1)*wwc(n,5)+ccx(4,1,i-1)*wwc(n,4) &
                      + ccx(3,1,i-1)*wwc(n,3) + ccx(2,1,i-1)*wwc(n,2) &
                      + ccx(1,1,i-1)*wwc(n,1))

                 qqlr = wwc(n,3)+Alpha*(wwc(n,3)-wwc(n,4))
                 qqav = 0.5d0*(wwc(n,3)+wwc(n,2))
                 qqmd = qqav - 0.5d0*dm4jmh
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,4)) + B2*dm4jph
                 
                 qqmin = max(min(wwc(n,3),wwc(n,2),qqmd),min(wwc(n,3),qqlr,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,2),qqmd),max(wwc(n,3),qqlr,qqlc))
                 wwc_w(2,n) = wwor + minmod((qqmin-wwor),(qqmax-wwor))
              end do
              
              ! characteristic to primitive
              qql = 0d0
              qqr = 0d0
              do m=1,nwave
                 temp1 = wwc_w(1,m)
                 temp2 = wwc_w(2,m)
                 do n=1,nwave
                    qqr(n) = qqr(n)+temp1*rem(n,m)
                    qql(n) = qql(n)+temp2*rem(n,m)
                 enddo
              enddo

              minvalue = min(qql(1),qqr(1),qql(5),qqr(5))
              romaxvalue = max(ww(1,1),ww(1,2),ww(1,3),ww(1,4),ww(1,5))
              prmaxvalue = max(ww(5,1),ww(5,2),ww(5,3),ww(5,4),ww(5,5))
              minvalue = min(minvalue,ro1-romaxvalue*floor,pr1-prmaxvalue*floor)
              smv = sign(1d0,minvalue)
              psmv = max(0d0,smv)
              msmv = max(0d0,-smv)
           
              row(i-1,j,k,2) = qql(1)*psmv + ro1*msmv
              row(i  ,j,k,1) = qqr(1)*psmv + ro1*msmv
              vxw(i-1,j,k,2) = qql(2)*psmv + vx1*msmv
              vxw(i  ,j,k,1) = qqr(2)*psmv + vx1*msmv
              vyw(i-1,j,k,2) = qql(3)*psmv + vy1*msmv
              vyw(i  ,j,k,1) = qqr(3)*psmv + vy1*msmv
              vzw(i-1,j,k,2) = qql(4)*psmv + vz1*msmv
              vzw(i  ,j,k,1) = qqr(4)*psmv + vz1*msmv
              prw(i-1,j,k,2) = qql(5)*psmv + pr1*msmv
              prw(i  ,j,k,1) = qqr(5)*psmv + pr1*msmv
              bxw(i-1,j,k,2) = qql(6)*psmv + bx1*msmv
              bxw(i  ,j,k,1) = qqr(6)*psmv + bx1*msmv
              byw(i-1,j,k,2) = qql(7)*psmv + by1*msmv
              byw(i  ,j,k,1) = qqr(7)*psmv + by1*msmv
              bzw(i-1,j,k,2) = qql(8)*psmv + bz1*msmv
              bzw(i  ,j,k,1) = qqr(8)*psmv + bz1*msmv
              phiw(i-1,j,k,2) = qql(9)*psmv + phi1*msmv
              phiw(i  ,j,k,1) = qqr(9)*psmv + phi1*msmv
           end do
        end do
     end do
  else if(mdir == 2)then
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
              ro1 = ww(1,3)
              vx1 = ww(2,3)
              vy1 = ww(3,3)
              vz1 = ww(4,3)
              pr1 = ww(5,3)
              bx1 = ww(6,3)
              by1 = ww(7,3)
              bz1 = ww(8,3)
              phi1 = ww(9,3)
              
              ! primitive to characteristic
              call esystem_glmmhd(lem,rem,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1,ch,gm)
              wwc = 0d0
              do l=1,5
                 do m=1,nwave
                    temp = ww(m,l)
                    do n=1,nwave
                       wwc(n,l) = wwc(n,l)+lem(n,m)*temp
                    enddo
                 enddo
              end do
              

              ! mp5
              do n=1,nwave
                 ! left state
                 wwor = B1*(ccy(1,2,j)*wwc(n,1)+ccy(2,2,j)*wwc(n,2) &
                      + ccy(3,2,j)*wwc(n,3) + ccy(4,2,j)*wwc(n,4) &
                      + ccy(5,2,j)*wwc(n,5))

                 djm1 = wwc(n,1)-2.0d0*wwc(n,2)+wwc(n,3)
                 dj = wwc(n,2)-2.0d0*wwc(n,3)+wwc(n,4)
                 djp1 = wwc(n,3)-2.0d0*wwc(n,4)+wwc(n,5)
                 
                 dm4jph = minmod4(4.0d0*dj-djp1,4.0d0*djp1-dj,dj,djp1)
                 dm4jmh = minmod4(4.0d0*dj-djm1,4.0d0*djm1-dj,dj,djm1)
                 
                 qqul = wwc(n,3)+Alpha*(wwc(n,3)-wwc(n,2))
                 qqav = 0.5d0*(wwc(n,3)+wwc(n,4))
                 qqmd = qqav - 0.5d0*dm4jph
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,2)) + B2*dm4jmh
                 
                 qqmin = max(min(wwc(n,3),wwc(n,4),qqmd),min(wwc(n,3),qqul,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,4),qqmd),max(wwc(n,3),qqul,qqlc))

                 wwc_w(1,n) = wwor + minmod((qqmin-wwor),(qqmax-wwor))
                 
                 ! right state
                 wwor = B1*(ccy(5,1,j-1)*wwc(n,5)+ccy(4,1,j-1)*wwc(n,4) &
                      + ccy(3,1,j-1)*wwc(n,3) + ccy(2,1,j-1)*wwc(n,2) &
                      + ccy(1,1,j-1)*wwc(n,1))

                 qqlr = wwc(n,3)+Alpha*(wwc(n,3)-wwc(n,4))
                 qqav = 0.5d0*(wwc(n,3)+wwc(n,2))
                 qqmd = qqav - 0.5d0*dm4jmh
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,4)) + B2*dm4jph
                 
                 qqmin = max(min(wwc(n,3),wwc(n,2),qqmd),min(wwc(n,3),qqlr,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,2),qqmd),max(wwc(n,3),qqlr,qqlc))
                 
                 wwc_w(2,n) = wwor + minmod((qqmin-wwor),(qqmax-wwor))
              end do
              
              ! characteristic to primitive
              qqr = 0d0
              qql = 0d0
              do m=1,nwave
                 temp1 = wwc_w(1,m)
                 temp2 = wwc_w(2,m)
                 do n=1,nwave
                    qqr(n) = qqr(n)+temp1*rem(n,m)
                    qql(n) = qql(n)+temp2*rem(n,m)
                 enddo
              enddo

              minvalue = min(qql(1),qqr(1),qql(5),qqr(5))
              romaxvalue = max(ww(1,1),ww(1,2),ww(1,3),ww(1,4),ww(1,5))
              prmaxvalue = max(ww(5,1),ww(5,2),ww(5,3),ww(5,4),ww(5,5))
              minvalue = min(minvalue,ro1-romaxvalue*floor,pr1-prmaxvalue*floor)
              smv = sign(1d0,minvalue)
              psmv = max(0d0,smv)
              msmv = max(0d0,-smv)

              row(i,j-1,k,2) = qql(1)*psmv + ro1*msmv
              row(i,j  ,k,1) = qqr(1)*psmv + ro1*msmv
              vxw(i,j-1,k,2) = qql(2)*psmv + vx1*msmv
              vxw(i,j  ,k,1) = qqr(2)*psmv + vx1*msmv
              vyw(i,j-1,k,2) = qql(3)*psmv + vy1*msmv
              vyw(i,j  ,k,1) = qqr(3)*psmv + vy1*msmv
              vzw(i,j-1,k,2) = qql(4)*psmv + vz1*msmv
              vzw(i,j  ,k,1) = qqr(4)*psmv + vz1*msmv
              prw(i,j-1,k,2) = qql(5)*psmv + pr1*msmv
              prw(i,j  ,k,1) = qqr(5)*psmv + pr1*msmv
              bxw(i,j-1,k,2) = qql(6)*psmv + bx1*msmv
              bxw(i,j  ,k,1) = qqr(6)*psmv + bx1*msmv
              byw(i,j-1,k,2) = qql(7)*psmv + by1*msmv
              byw(i,j  ,k,1) = qqr(7)*psmv + by1*msmv
              bzw(i,j-1,k,2) = qql(8)*psmv + bz1*msmv
              bzw(i,j  ,k,1) = qqr(8)*psmv + bz1*msmv
              phiw(i,j-1,k,2) = qql(9)*psmv + phi1*msmv
              phiw(i,j  ,k,1) = qqr(9)*psmv + phi1*msmv
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
              ro1 = ww(1,3)
              vx1 = ww(2,3)
              vy1 = ww(3,3)
              vz1 = ww(4,3)
              pr1 = ww(5,3)
              bx1 = ww(6,3)
              by1 = ww(7,3)
              bz1 = ww(8,3)
              phi1 = ww(9,3)
              
              ! primitive to characteristic
              call esystem_glmmhd(lem,rem,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1,ch,gm)
              wwc = 0d0
              do l=1,5
                 do m=1,nwave
                    temp = ww(m,l)
                    do n=1,nwave
                       wwc(n,l) = wwc(n,l)+lem(n,m)*temp
                    enddo
                 enddo
              end do
              
              ! mp5
              do n=1,nwave
                 ! left state
                 wwor = B1*(ccz(1,2,k)*wwc(n,1)+ccz(2,2,k)*wwc(n,2) &
                      + ccz(3,2,k)*wwc(n,3) + ccz(4,2,k)*wwc(n,4) &
                      + ccz(5,2,k)*wwc(n,5))

                 djm1 = wwc(n,1)-2.0d0*wwc(n,2)+wwc(n,3)
                 dj = wwc(n,2)-2.0d0*wwc(n,3)+wwc(n,4)
                 djp1 = wwc(n,3)-2.0d0*wwc(n,4)+wwc(n,5)
                 
                 dm4jph = minmod4(4.0d0*dj-djp1,4.0d0*djp1-dj,dj,djp1)
                 dm4jmh = minmod4(4.0d0*dj-djm1,4.0d0*djm1-dj,dj,djm1)
                 
                 qqul = wwc(n,3)+Alpha*(wwc(n,3)-wwc(n,2))
                 qqav = 0.5d0*(wwc(n,3)+wwc(n,4))
                 qqmd = qqav - 0.5d0*dm4jph
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,2)) + B2*dm4jmh
                 
                 qqmin = max(min(wwc(n,3),wwc(n,4),qqmd),min(wwc(n,3),qqul,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,4),qqmd),max(wwc(n,3),qqul,qqlc))

                 wwc_w(1,n) = wwor + minmod((qqmin-wwor),(qqmax-wwor))

              ! right state
                 wwor = B1*(ccz(5,1,k-1)*wwc(n,5)+ccz(4,1,k-1)*wwc(n,4) &
                      + ccz(3,1,k-1)*wwc(n,3) + ccz(2,1,k-1)*wwc(n,2) &
                      + ccz(1,1,k-1)*wwc(n,1))

                 qqlr = wwc(n,3)+Alpha*(wwc(n,3)-wwc(n,4))
                 qqav = 0.5d0*(wwc(n,3)+wwc(n,2))
                 qqmd = qqav - 0.5d0*dm4jmh
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,4)) + B2*dm4jph
                 
                 qqmin = max(min(wwc(n,3),wwc(n,2),qqmd),min(wwc(n,3),qqlr,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,2),qqmd),max(wwc(n,3),qqlr,qqlc))
                 
                 wwc_w(2,n) = wwor + minmod((qqmin-wwor),(qqmax-wwor))
              end do
              
              ! characteristic to primitive
              qqr = 0d0
              qql = 0d0
              do m=1,nwave
                 temp1 = wwc_w(1,m)
                 temp2 = wwc_w(2,m)
                 do n=1,nwave
                    qqr(n) = qqr(n)+temp1*rem(n,m)
                    qql(n) = qql(n)+temp2*rem(n,m)
                 enddo
              enddo

              minvalue = min(qql(1),qqr(1),qql(5),qqr(5))
              romaxvalue = max(ww(1,1),ww(1,2),ww(1,3),ww(1,4),ww(1,5))
              prmaxvalue = max(ww(5,1),ww(5,2),ww(5,3),ww(5,4),ww(5,5))
              minvalue = min(minvalue,ro1-romaxvalue*floor,pr1-prmaxvalue*floor)
              smv = sign(1d0,minvalue)
              psmv = max(0d0,smv)
              msmv = max(0d0,-smv)

              row(i,j,k-1,2) = qql(1)*psmv + ro1*msmv
              row(i,j,k  ,1) = qqr(1)*psmv + ro1*msmv
              vxw(i,j,k-1,2) = qql(2)*psmv + vx1*msmv
              vxw(i,j,k  ,1) = qqr(2)*psmv + vx1*msmv
              vyw(i,j,k-1,2) = qql(3)*psmv + vy1*msmv
              vyw(i,j,k  ,1) = qqr(3)*psmv + vy1*msmv
              vzw(i,j,k-1,2) = qql(4)*psmv + vz1*msmv
              vzw(i,j,k  ,1) = qqr(4)*psmv + vz1*msmv
              prw(i,j,k-1,2) = qql(5)*psmv + pr1*msmv
              prw(i,j,k  ,1) = qqr(5)*psmv + pr1*msmv
              bxw(i,j,k-1,2) = qql(6)*psmv + bx1*msmv
              bxw(i,j,k  ,1) = qqr(6)*psmv + bx1*msmv
              byw(i,j,k-1,2) = qql(7)*psmv + by1*msmv
              byw(i,j,k  ,1) = qqr(7)*psmv + by1*msmv
              bzw(i,j,k-1,2) = qql(8)*psmv + bz1*msmv
              bzw(i,j,k  ,1) = qqr(8)*psmv + bz1*msmv
              phiw(i,j,k-1,2) = qql(9)*psmv + phi1*msmv
              phiw(i,j,k  ,1) = qqr(9)*psmv + phi1*msmv
           end do
        end do
     end do
  endif

  end subroutine lr_state__MP5


  subroutine reconstructionConstant(margin,ix,x,xm,dx,cc)

  integer,intent(in) :: ix,margin
  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(0:ix),intent(in) :: xm

  integer :: i,m,l,q
  integer,parameter :: k = 5
  integer :: j,r
  real(8),dimension(5,2,ix) :: cc
  real(8) :: tmp1,tmp2,tmp3,tmp4

  do i=3,ix-3
     do r=1,2
        do j=0,4
           tmp3 = 0.0d0
           do m=(j+1),k
              tmp2 = 0.0d0
              do l=0,k
                 if(l .ne. m)then
                    tmp1 = 1.0d0
                    do q=0,k
                       if(q .ne. m)then
                          if(q .ne. l)then
                             tmp1 = tmp1*(xm(i)-xm(i-r+q-1))
                          endif
                       endif
                    enddo
                    tmp2 = tmp2 + tmp1
                 endif
              enddo

              tmp4 = 1.0d0
              do l=0,k
                 if(l .ne. m)then
                    tmp4 = tmp4*(xm(i-r+m-1)-xm(i-r+l-1))
                 end if
              enddo
              tmp3 = tmp3 + tmp2/tmp4
              tmp2 = 0.d0
              tmp4 = 1.0d0
           end do
           cc(j+1,r,i) = 60.0d0*tmp3*dx(i-r+j)
           tmp3 = 0.0d0
        end do
     end do
  end do

  end subroutine reconstructionConstant


  subroutine esystem_glmmhd(lem,rem,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,gm)

  real(8),intent(in) :: gm,ch
  real(8),intent(in) :: ro,pr,vx,vy,vz
  real(8),intent(in) :: bx,by,bz,phi
  real(8),dimension(9,9),intent(out) :: lem,rem

  real(8) :: roi,btsq,vaxsq,asq
  real(8) :: ct2,tsum,tdif,cf2_cs2
  real(8) :: cfsq,cf,cssq,cs
  real(8) :: bt,ibt,sbt,psbt,msbt,bet2,bet3
  real(8) :: sfs,psfs,msfs,sas,saf,eps,ifs
  real(8) :: alpha_f,alpha_s
  real(8) :: na,qf,qs,af_prm,as_prm
  real(8) :: sqrtro,s,a,af,as,sqrt2i
  real(8) :: ich

  sqrt2i = 1.D0/sqrt(2.D0)
  roi = 1.0d0/ro
  btsq = by**2+bz**2
  vaxsq = (bx**2)*roi
  asq = gm*pr*roi

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
  eps = 1d-40
  sbt = sign(1d0,bt-eps)
  psbt = max(0d0,sbt)
  msbt = max(0d0,-sbt)
  ibt = psbt/(bt+msbt)

  bet2 = msbt*sqrt2i + by*ibt
  bet3 = msbt*sqrt2i + bz*ibt
  
! compute alpha
  sfs = sign(1d0,cf2_cs2-eps)
  psfs = max(0d0,sfs)
  msfs = max(0d0,-sfs)
  sas = max(0d0,sign(1d0,asq-cssq))
  saf = max(0d0,sign(1d0,asq-cfsq))
  ifs = 1.0d0/(cf2_cs2+msfs)

  alpha_f = msfs + psfs*((1.0d0-saf)*sqrt(max(asq-cssq,0d0)*ifs) + saf)
  alpha_s =        psfs*(sas*sqrt(max(cfsq-asq,0d0)*ifs) + (1.0d0-sas))

! compute Q(s),A(s)
  sqrtro = sqrt(ro)
  s = sign(1.0d0,bx)
  a = sqrt(asq)
  qf = cf*alpha_f*s
  qs = cs*alpha_s*s
  af = a*alpha_f*sqrtro
  as = a*alpha_s*sqrtro

! Right-eigenvector
!  ich = 1.0d0/(ch)
  rem(1,1) = ro*alpha_f
  rem(1,2) = 0.0d0
  rem(1,3) = ro*alpha_s
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

  rem(5,1) = ro*asq*alpha_f
  rem(5,2) = 0.0d0
  rem(5,3) = ro*asq*alpha_s
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

  end subroutine esystem_glmmhd


  function minmod_limiter2(qqr,qql)

    real(8), intent(in) :: qqr, qql
    real(8) :: minmod_limiter2
    real(8) :: signr,signlr

    signr = sign(1.0d0,qqr)
    signlr = sign(1.0d0,(qqr*qql))

    minmod_limiter2 = max(signlr,0.0d0)*( &
         max(signr,0.0d0)*max(0.0d0,min(qql,qqr)) &
         -min(signr,0.0d0)*min(0.0d0,max(qql,qqr)))

  end function minmod_limiter2
  

  function MC_limiter(qqc,qql,qqr)

    real(8), intent(in) :: qqr,qql,qqc
    real(8) :: minmod_lr
    real(8) :: MC_limiter
    real(8) :: signlr,signMC

    minmod_lr = minmod_limiter2(qql,qqr)

    signlr = sign(1.0d0,(qqc*minmod_lr))
    signMC = sign(1.0d0,(dabs(qqc)-dabs(2.0d0*minmod_lr)))

    MC_limiter = max(signlr,0.0d0)*( &
         max(signMC,0.0d0)*minmod_lr &
         -min(signMC,0.0d0)*(qqc))

  end function MC_limiter


end module lr_state
