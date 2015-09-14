module lr_state

  implicit none
  private

  public :: lr_state__1st, lr_state__MSCL2, lr_state__MP5, reconstructionConstant


contains


  subroutine lr_state__1st(mdir,ix,jx,kx,&
                           ro,pr,vx,vy,vz,bx,by,bz,phi,&
                           row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw)
  integer,intent(in) :: mdir
  integer,intent(in) :: ix,jx,kx
  real(8),dimension(ix,jx,kx),intent(in) :: ro,pr,vx,vy,vz,bx,by,bz,phi
  real(8),dimension(ix,jx,kx,2),intent(out) :: row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw

  integer :: i,j,k

  if(mdir == 1)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
              row(i-1,j,k,2) = ro(i,j,k)
              row(i  ,j,k,1) = ro(i,j,k)
              vxw(i-1,j,k,2) = vx(i,j,k)
              vxw(i  ,j,k,1) = vx(i,j,k)
              vyw(i-1,j,k,2) = vy(i,j,k)
              vyw(i  ,j,k,1) = vy(i,j,k)
              vzw(i-1,j,k,2) = vz(i,j,k)
              vzw(i  ,j,k,1) = vz(i,j,k)
              prw(i-1,j,k,2) = pr(i,j,k)
              prw(i  ,j,k,1) = pr(i,j,k)
              bxw(i-1,j,k,2) = bx(i,j,k)
              bxw(i  ,j,k,1) = bx(i,j,k)
              byw(i-1,j,k,2) = by(i,j,k)
              byw(i  ,j,k,1) = by(i,j,k)
              bzw(i-1,j,k,2) = bz(i,j,k)
              bzw(i  ,j,k,1) = bz(i,j,k)
              phiw(i-1,j,k,2) = phi(i,j,k)
              phiw(i  ,j,k,1) = phi(i,j,k)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  else if(mdir == 2)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
              row(i,j-1,k,2) = ro(i,j,k)
              row(i,j  ,k,1) = ro(i,j,k)
              vxw(i,j-1,k,2) = vx(i,j,k)
              vxw(i,j  ,k,1) = vx(i,j,k)
              vyw(i,j-1,k,2) = vy(i,j,k)
              vyw(i,j  ,k,1) = vy(i,j,k)
              vzw(i,j-1,k,2) = vz(i,j,k)
              vzw(i,j  ,k,1) = vz(i,j,k)
              prw(i,j-1,k,2) = pr(i,j,k)
              prw(i,j  ,k,1) = pr(i,j,k)
              bxw(i,j-1,k,2) = bx(i,j,k)
              bxw(i,j  ,k,1) = bx(i,j,k)
              byw(i,j-1,k,2) = by(i,j,k)
              byw(i,j  ,k,1) = by(i,j,k)
              bzw(i,j-1,k,2) = bz(i,j,k)
              bzw(i,j  ,k,1) = bz(i,j,k)
              phiw(i,j-1,k,2) = phi(i,j,k)
              phiw(i,j  ,k,1) = phi(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j)
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1
              row(i,j,k-1,2) = ro(i,j,k)
              row(i,j,k  ,1) = ro(i,j,k)
              vxw(i,j,k-1,2) = vx(i,j,k)
              vxw(i,j,k  ,1) = vx(i,j,k)
              vyw(i,j,k-1,2) = vy(i,j,k)
              vyw(i,j,k  ,1) = vy(i,j,k)
              vzw(i,j,k-1,2) = vz(i,j,k)
              vzw(i,j,k  ,1) = vz(i,j,k)
              prw(i,j,k-1,2) = pr(i,j,k)
              prw(i,j,k  ,1) = pr(i,j,k)
              bxw(i,j,k-1,2) = bx(i,j,k)
              bxw(i,j,k  ,1) = bx(i,j,k)
              byw(i,j,k-1,2) = by(i,j,k)
              byw(i,j,k  ,1) = by(i,j,k)
              bzw(i,j,k-1,2) = bz(i,j,k)
              bzw(i,j,k  ,1) = bz(i,j,k)
              phiw(i,j,k-1,2) = phi(i,j,k)
              phiw(i,j,k  ,1) = phi(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  endif

  end subroutine lr_state__1st



  subroutine lr_state__MSCL2(mdir,ix,jx,kx,&
                             ro,pr,vx,vy,vz,bx,by,bz,phi,ch,gm,&
                             row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw,dx,dy,dz)

  integer,intent(in) :: ix,jx,kx,mdir
  real(8),intent(in) :: gm,ch
  real(8),dimension(ix),intent(in) :: dx
  real(8),dimension(jx),intent(in) :: dy
  real(8),dimension(kx),intent(in) :: dz
  real(8),dimension(ix,jx,kx),intent(in) :: ro,pr,vx,vy,vz,bx,by,bz,phi
  real(8),dimension(ix,jx,kx,2),intent(out) :: row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw

  integer,parameter :: nwave = 9
  integer :: i,j,k,l,n
  real(8),dimension(nwave) :: qql,qqr
  real(8),dimension(nwave,2) :: wwc_w
  ! ww(*,1) = qqm2, ..
  real(8),dimension(nwave,3) :: ww,wwc
  real(8),dimension(nwave,nwave) :: lem,rem
  real(8) :: minvalue,smv,psmv,msmv
  real(8) :: ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1
  real(8) :: temp1,temp2,temp3,ich
  real(8) :: dqqx,dqqy,dqqz
  real(8) :: dqc,dql,dqr
  real(8) :: dxc,dxl,dxr
  real(8) :: dyc,dyl,dyr
  real(8) :: dzc,dzl,dzr

!----Recipe for rarefaction ---
  real(8), parameter :: floor=1d-2
  real(8) :: romaxvalue,prmaxvalue

  ich = 0.5d0/ch

  if(mdir == 1)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,l,n,ww,ro1,vx1,vy1,vz1,pr1,bx1,by1,bz1,phi1,lem,rem,wwc,&
     !$OMP         temp1,temp2,temp3,dxl,dxr,dxc,dql,dqr,dqc,dqqx,wwc_w,       &
     !$OMP         qql,qqr,minvalue,romaxvalue,prmaxvalue,smv,psmv,msmv)
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1

              do l=1,3
                 ww(1,l) = ro(i-2+l,j,k)
                 ww(2,l) = vx(i-2+l,j,k)
                 ww(3,l) = vy(i-2+l,j,k)
                 ww(4,l) = vz(i-2+l,j,k)
                 ww(5,l) = pr(i-2+l,j,k)
                 ww(6,l) = bx(i-2+l,j,k)
                 ww(7,l) = by(i-2+l,j,k)
                 ww(8,l) = bz(i-2+l,j,k)
                 ww(9,l) = phi(i-2+l,j,k)
              end do
              ro1 = ww(1,2)
              vx1 = ww(2,2)
              vy1 = ww(3,2)
              vz1 = ww(4,2)
              pr1 = ww(5,2)
              bx1 = ww(6,2)
              by1 = ww(7,2)
              bz1 = ww(8,2)
              phi1 = ww(9,2)

              ! primitive to characteristic
              call esystem_glmmhd(lem,rem,ro1,pr1,bx1,by1,bz1,gm)

              do l=1,3
                 wwc(4,l) = ww(1,l)+lem(4,5)*ww(5,l)

                 temp1 = lem(1,2)*ww(2,l) + lem(1,3)*ww(3,l) + lem(1,4)*ww(4,l)
                 temp2 = lem(1,5)*ww(5,l) + lem(1,7)*ww(7,l) + lem(1,8)*ww(8,l)
                 wwc(1,l) = temp1 + temp2
                 wwc(8,l) = -temp1 + temp2

                 temp1 = lem(2,3)*ww(3,l) + lem(2,4)*ww(4,l)
                 temp2 = lem(2,7)*ww(7,l) + lem(2,8)*ww(8,l)
                 wwc(2,l) = temp1 + temp2
                 wwc(7,l) = -temp1 + temp2
                      
                 temp1 = lem(3,2)*ww(2,l) + lem(3,3)*ww(3,l) + lem(3,4)*ww(4,l)
                 temp2 = lem(3,5)*ww(5,l) + lem(3,7)*ww(7,l) + lem(3,8)*ww(8,l)
                 wwc(3,l) = temp1 +temp2
                 wwc(5,l) = -temp1 + temp2
                 
                 temp1 = 0.5d0*ww(6,l)
                 temp2 = ich*ww(9,l) 
                 wwc(6,l) = temp1 - temp2
                 wwc(9,l) = temp1 + temp2
              enddo
              
              !MUSCL
              do n=1,nwave
                 dxl = 0.5d0*(dx(i)+dx(i-1))
                 dxr = 0.5d0*(dx(i+1)+dx(i))
                 dxc = dxl+dxr

                 dql = (wwc(n,2)-wwc(n,1))/dxl
                 dqr = (wwc(n,3)-wwc(n,2))/dxr
                 dqc = (wwc(n,3)-wwc(n,1))/dxc

                 dqqx = MC2(2.*dqr,2.*dql,dqc)

                 !right-hand left-state
                 wwc_w(n,1) = wwc(n,2) + 0.5d0*dqqx*dx(i)
                 !left-hand right-state
                 wwc_w(n,2) = wwc(n,2) - 0.5d0*dqqx*dx(i)
              end do

              ! characteristic to primitive
              temp1 = wwc_w(1,1)-wwc_w(8,1)
              temp2 = wwc_w(2,1)-wwc_w(7,1)
              temp3 = wwc_w(3,1)-wwc_w(5,1)
              qqr(2) = rem(2,1)*temp1 + rem(2,3)*temp3
              qqr(3) = rem(3,1)*temp1 + rem(3,2)*temp2 + rem(3,3)*temp3
              qqr(4) = rem(4,1)*temp1 + rem(4,2)*temp2 + rem(4,3)*temp3

              temp1 = wwc_w(1,1)+wwc_w(8,1)
              temp2 = wwc_w(2,1)+wwc_w(7,1)
              temp3 = wwc_w(3,1)+wwc_w(5,1)
              qqr(1) = rem(1,1)*temp1 + rem(1,3)*temp3 + wwc_w(4,1)
              qqr(5) = rem(5,1)*temp1 + rem(5,3)*temp3
              qqr(6) = wwc_w(6,1) + wwc_w(9,1)
              qqr(7) = rem(7,1)*temp1 + rem(7,2)*temp2 + rem(7,3)*temp3
              qqr(8) = rem(8,1)*temp1 + rem(8,2)*temp2 + rem(8,3)*temp3
              qqr(9) = ch*(-wwc_w(6,1) + wwc_w(9,1)) 

              temp1 = wwc_w(1,2)-wwc_w(8,2)
              temp2 = wwc_w(2,2)-wwc_w(7,2)
              temp3 = wwc_w(3,2)-wwc_w(5,2)
              qql(2) = rem(2,1)*temp1 + rem(2,3)*temp3
              qql(3) = rem(3,1)*temp1 + rem(3,2)*temp2 + rem(3,3)*temp3
              qql(4) = rem(4,1)*temp1 + rem(4,2)*temp2 + rem(4,3)*temp3

              temp1 = wwc_w(1,2)+wwc_w(8,2)
              temp2 = wwc_w(2,2)+wwc_w(7,2)
              temp3 = wwc_w(3,2)+wwc_w(5,2)
              qql(1) = rem(1,1)*temp1 + rem(1,3)*temp3 + wwc_w(4,2)
              qql(5) = rem(5,1)*temp1 + rem(5,3)*temp3
              qql(6) = wwc_w(6,2) + wwc_w(9,2)
              qql(7) = rem(7,1)*temp1 + rem(7,2)*temp2 + rem(7,3)*temp3
              qql(8) = rem(8,1)*temp1 + rem(8,2)*temp2 + rem(8,3)*temp3
              qql(9) = ch*(-wwc_w(6,2) + wwc_w(9,2)) 

              minvalue = min(qql(1),qqr(1),qql(5),qqr(5))
              romaxvalue = max(ww(1,1),ww(1,2),ww(1,3))
              prmaxvalue = max(ww(5,1),ww(5,2),ww(5,3))
              minvalue = min(minvalue,ro1-romaxvalue*floor,pr1-prmaxvalue*floor)
              smv = sign(1d0,minvalue)
              psmv = max(0d0,smv)
              msmv = max(0d0,-smv)

              row(i-1,j,k,2) = qql(1)*psmv + ro1*msmv
              row(i  ,j,k,1) = qqr(1)*psmv + ro1*msmv
              vxw(i-1,j,k,2) = qql(2)*psmv + vx1*msmv
              vxw(i  ,j,k,1) = qqr(2)*psmv + vx1*msmv
              vyw(i-1,j,k,2) = qql(3)
              vyw(i  ,j,k,1) = qqr(3)
              vzw(i-1,j,k,2) = qql(4)
              vzw(i  ,j,k,1) = qqr(4)
              prw(i-1,j,k,2) = qql(5)*psmv + pr1*msmv
              prw(i  ,j,k,1) = qqr(5)*psmv + pr1*msmv
              bxw(i-1,j,k,2) = qql(6)*psmv + bx1*msmv
              bxw(i  ,j,k,1) = qqr(6)*psmv + bx1*msmv
              byw(i-1,j,k,2) = qql(7)
              byw(i  ,j,k,1) = qqr(7)
              bzw(i-1,j,k,2) = qql(8)
              bzw(i  ,j,k,1) = qqr(8)
              phiw(i-1,j,k,2) = qql(9)*psmv + phi1*msmv
              phiw(i  ,j,k,1) = qqr(9)*psmv + phi1*msmv
           end do
        end do
     end do
!$OMP END PARALLEL DO
  else if(mdir == 2)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,n,l,ww,ro1,vx1,vy1,vz1,pr1,bx1,by1,bz1,phi1,lem,rem,wwc,&
     !$OMP         temp1,temp2,temp3,dxl,dxr,dxc,dql,dqr,dqc,dqqx,wwc_w,       &
     !$OMP         qql,qqr,minvalue,romaxvalue,prmaxvalue,smv,psmv,msmv)
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1

              do l=1,3
                 ww(1,l) = ro(i,j-2+l,k)
                 ww(2,l) = vx(i,j-2+l,k)
                 ww(3,l) = vy(i,j-2+l,k)
                 ww(4,l) = vz(i,j-2+l,k)
                 ww(5,l) = pr(i,j-2+l,k)
                 ww(6,l) = bx(i,j-2+l,k)
                 ww(7,l) = by(i,j-2+l,k)
                 ww(8,l) = bz(i,j-2+l,k)
                 ww(9,l) = phi(i,j-2+l,k)
              end do
              ro1 = ww(1,2)
              vx1 = ww(2,2)
              vy1 = ww(3,2)
              vz1 = ww(4,2)
              pr1 = ww(5,2)
              bx1 = ww(6,2)
              by1 = ww(7,2)
              bz1 = ww(8,2)
              phi1 = ww(9,2)
              
              ! primitive to characteristic
              call esystem_glmmhd(lem,rem,ro1,pr1,bx1,by1,bz1,gm)

              do l=1,3
                 wwc(4,l) = ww(1,l)+lem(4,5)*ww(5,l)

                 temp1 = lem(1,2)*ww(2,l) + lem(1,3)*ww(3,l) + lem(1,4)*ww(4,l)
                 temp2 = lem(1,5)*ww(5,l) + lem(1,7)*ww(7,l) + lem(1,8)*ww(8,l)
                 wwc(1,l) = temp1 + temp2
                 wwc(8,l) = -temp1 + temp2

                 temp1 = lem(2,3)*ww(3,l) + lem(2,4)*ww(4,l)
                 temp2 = lem(2,7)*ww(7,l) + lem(2,8)*ww(8,l)
                 wwc(2,l) = temp1 + temp2
                 wwc(7,l) = -temp1 + temp2
                      
                 temp1 = lem(3,2)*ww(2,l) + lem(3,3)*ww(3,l) + lem(3,4)*ww(4,l)
                 temp2 = lem(3,5)*ww(5,l) + lem(3,7)*ww(7,l) + lem(3,8)*ww(8,l)
                 wwc(3,l) = temp1 +temp2
                 wwc(5,l) = -temp1 + temp2
                 
                 temp1 = 0.5d0*ww(6,l)
                 temp2 = ich*ww(9,l) 
                 wwc(6,l) = temp1 - temp2
                 wwc(9,l) = temp1 + temp2
              enddo

              !MUSCL
              do n=1,nwave
                 dyl = 0.5d0*(dy(j)+dy(j-1))
                 dyr = 0.5d0*(dy(j+1)+dy(j))
                 dyc = dyl+dyr

                 dql = (wwc(n,2)-wwc(n,1))/dyl
                 dqr = (wwc(n,3)-wwc(n,2))/dyr
                 dqc = (wwc(n,3)-wwc(n,1))/dyc

                 dqqy = MC2(2.*dqr,2.*dql,dqc)

                 !right-hand left-state
                 wwc_w(n,1) = wwc(n,2) + 0.5d0*dqqy*dy(j)
                 !left-hand right-state
                 wwc_w(n,2) = wwc(n,2) - 0.5d0*dqqy*dy(j)
              end do
              
              ! characteristic to primitive
              temp1 = wwc_w(1,1)-wwc_w(8,1)
              temp2 = wwc_w(2,1)-wwc_w(7,1)
              temp3 = wwc_w(3,1)-wwc_w(5,1)
              qqr(2) = rem(2,1)*temp1 + rem(2,3)*temp3
              qqr(3) = rem(3,1)*temp1 + rem(3,2)*temp2 + rem(3,3)*temp3
              qqr(4) = rem(4,1)*temp1 + rem(4,2)*temp2 + rem(4,3)*temp3

              temp1 = wwc_w(1,1)+wwc_w(8,1)
              temp2 = wwc_w(2,1)+wwc_w(7,1)
              temp3 = wwc_w(3,1)+wwc_w(5,1)
              qqr(1) = rem(1,1)*temp1 + rem(1,3)*temp3 + wwc_w(4,1)
              qqr(5) = rem(5,1)*temp1 + rem(5,3)*temp3
              qqr(6) = wwc_w(6,1) + wwc_w(9,1)
              qqr(7) = rem(7,1)*temp1 + rem(7,2)*temp2 + rem(7,3)*temp3
              qqr(8) = rem(8,1)*temp1 + rem(8,2)*temp2 + rem(8,3)*temp3
              qqr(9) = ch*(-wwc_w(6,1) + wwc_w(9,1)) 

              temp1 = wwc_w(1,2)-wwc_w(8,2)
              temp2 = wwc_w(2,2)-wwc_w(7,2)
              temp3 = wwc_w(3,2)-wwc_w(5,2)
              qql(2) = rem(2,1)*temp1 + rem(2,3)*temp3
              qql(3) = rem(3,1)*temp1 + rem(3,2)*temp2 + rem(3,3)*temp3
              qql(4) = rem(4,1)*temp1 + rem(4,2)*temp2 + rem(4,3)*temp3

              temp1 = wwc_w(1,2)+wwc_w(8,2)
              temp2 = wwc_w(2,2)+wwc_w(7,2)
              temp3 = wwc_w(3,2)+wwc_w(5,2)
              qql(1) = rem(1,1)*temp1 + rem(1,3)*temp3 + wwc_w(4,2)
              qql(5) = rem(5,1)*temp1 + rem(5,3)*temp3
              qql(6) = wwc_w(6,2) + wwc_w(9,2)
              qql(7) = rem(7,1)*temp1 + rem(7,2)*temp2 + rem(7,3)*temp3
              qql(8) = rem(8,1)*temp1 + rem(8,2)*temp2 + rem(8,3)*temp3
              qql(9) = ch*(-wwc_w(6,2) + wwc_w(9,2)) 

              minvalue = min(qql(1),qqr(1),qql(5),qqr(5))
              romaxvalue = max(ww(1,1),ww(1,2),ww(1,3))
              prmaxvalue = max(ww(5,1),ww(5,2),ww(5,3))
              minvalue = min(minvalue,ro1-romaxvalue*floor,pr1-prmaxvalue*floor)
              smv = sign(1d0,minvalue)
              psmv = max(0d0,smv)
              msmv = max(0d0,-smv)

              row(i,j-1,k,2) = qql(1)*psmv + ro1*msmv
              row(i,j  ,k,1) = qqr(1)*psmv + ro1*msmv
              vxw(i,j-1,k,2) = qql(2)*psmv + vx1*msmv
              vxw(i,j  ,k,1) = qqr(2)*psmv + vx1*msmv
              vyw(i,j-1,k,2) = qql(3)
              vyw(i,j  ,k,1) = qqr(3)
              vzw(i,j-1,k,2) = qql(4)
              vzw(i,j  ,k,1) = qqr(4)
              prw(i,j-1,k,2) = qql(5)*psmv + pr1*msmv
              prw(i,j  ,k,1) = qqr(5)*psmv + pr1*msmv
              bxw(i,j-1,k,2) = qql(6)*psmv + bx1*msmv
              bxw(i,j  ,k,1) = qqr(6)*psmv + bx1*msmv
              byw(i,j-1,k,2) = qql(7)
              byw(i,j  ,k,1) = qqr(7)
              bzw(i,j-1,k,2) = qql(8)
              bzw(i,j  ,k,1) = qqr(8)
              phiw(i,j-1,k,2) = qql(9)*psmv + phi1*msmv
              phiw(i,j  ,k,1) = qqr(9)*psmv + phi1*msmv
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,n,l,ww,ro1,vx1,vy1,vz1,pr1,bx1,by1,bz1,phi1,lem,rem,wwc,&
     !$OMP         temp1,temp2,temp3,dxl,dxr,dxc,dql,dqr,dqc,dqqx,wwc_w,       &
     !$OMP         qql,qqr,minvalue,romaxvalue,prmaxvalue,smv,psmv,msmv)
     do k=2,kx-1
        do j=2,jx-1
           do i=2,ix-1

              do l=1,3
                 ww(1,l) = ro(i,j,k-2+l)
                 ww(2,l) = vx(i,j,k-2+l)
                 ww(3,l) = vy(i,j,k-2+l)
                 ww(4,l) = vz(i,j,k-2+l)
                 ww(5,l) = pr(i,j,k-2+l)
                 ww(6,l) = bx(i,j,k-2+l)
                 ww(7,l) = by(i,j,k-2+l)
                 ww(8,l) = bz(i,j,k-2+l)
                 ww(9,l) = phi(i,j,k-2+l)
              end do
              ro1 = ww(1,2)
              vx1 = ww(2,2)
              vy1 = ww(3,2)
              vz1 = ww(4,2)
              pr1 = ww(5,2)
              bx1 = ww(6,2)
              by1 = ww(7,2)
              bz1 = ww(8,2)
              phi1 = ww(9,2)
              
              ! primitive to characteristic
              call esystem_glmmhd(lem,rem,ro1,pr1,bx1,by1,bz1,gm)

              do l=1,3
                 wwc(4,l) = ww(1,l)+lem(4,5)*ww(5,l)

                 temp1 = lem(1,2)*ww(2,l) + lem(1,3)*ww(3,l) + lem(1,4)*ww(4,l)
                 temp2 = lem(1,5)*ww(5,l) + lem(1,7)*ww(7,l) + lem(1,8)*ww(8,l)
                 wwc(1,l) = temp1 + temp2
                 wwc(8,l) = -temp1 + temp2

                 temp1 = lem(2,3)*ww(3,l) + lem(2,4)*ww(4,l)
                 temp2 = lem(2,7)*ww(7,l) + lem(2,8)*ww(8,l)
                 wwc(2,l) = temp1 + temp2
                 wwc(7,l) = -temp1 + temp2
                      
                 temp1 = lem(3,2)*ww(2,l) + lem(3,3)*ww(3,l) + lem(3,4)*ww(4,l)
                 temp2 = lem(3,5)*ww(5,l) + lem(3,7)*ww(7,l) + lem(3,8)*ww(8,l)
                 wwc(3,l) = temp1 +temp2
                 wwc(5,l) = -temp1 + temp2
                 
                 temp1 = 0.5d0*ww(6,l)
                 temp2 = ich*ww(9,l) 
                 wwc(6,l) = temp1 - temp2
                 wwc(9,l) = temp1 + temp2
              enddo

              !MUSCL
              do n=1,nwave
                 dzl = 0.5d0*(dz(k)+dz(k-1))
                 dzr = 0.5d0*(dz(k+1)+dz(k))
                 dzc = dzl+dzr

                 dql = (wwc(n,2)-wwc(n,1))/dzl
                 dqr = (wwc(n,3)-wwc(n,2))/dzr
                 dqc = (wwc(n,3)-wwc(n,1))/dzc

                 dqqz = MC2(2.*dqr,2.*dql,dqc)
              
                 !right-hand left-state
                 wwc_w(n,1) = wwc(n,2) + 0.5d0*dqqz*dz(k)
                 !left-hand right-state
                 wwc_w(n,2) = wwc(n,2) - 0.5d0*dqqz*dz(k)
              end do
              
              ! characteristic to primitive
              temp1 = wwc_w(1,1)-wwc_w(8,1)
              temp2 = wwc_w(2,1)-wwc_w(7,1)
              temp3 = wwc_w(3,1)-wwc_w(5,1)
              qqr(2) = rem(2,1)*temp1 + rem(2,3)*temp3
              qqr(3) = rem(3,1)*temp1 + rem(3,2)*temp2 + rem(3,3)*temp3
              qqr(4) = rem(4,1)*temp1 + rem(4,2)*temp2 + rem(4,3)*temp3

              temp1 = wwc_w(1,1)+wwc_w(8,1)
              temp2 = wwc_w(2,1)+wwc_w(7,1)
              temp3 = wwc_w(3,1)+wwc_w(5,1)
              qqr(1) = rem(1,1)*temp1 + rem(1,3)*temp3 + wwc_w(4,1)
              qqr(5) = rem(5,1)*temp1 + rem(5,3)*temp3
              qqr(6) = wwc_w(6,1) + wwc_w(9,1)
              qqr(7) = rem(7,1)*temp1 + rem(7,2)*temp2 + rem(7,3)*temp3
              qqr(8) = rem(8,1)*temp1 + rem(8,2)*temp2 + rem(8,3)*temp3
              qqr(9) = ch*(-wwc_w(6,1) + wwc_w(9,1)) 

              temp1 = wwc_w(1,2)-wwc_w(8,2)
              temp2 = wwc_w(2,2)-wwc_w(7,2)
              temp3 = wwc_w(3,2)-wwc_w(5,2)
              qql(2) = rem(2,1)*temp1 + rem(2,3)*temp3
              qql(3) = rem(3,1)*temp1 + rem(3,2)*temp2 + rem(3,3)*temp3
              qql(4) = rem(4,1)*temp1 + rem(4,2)*temp2 + rem(4,3)*temp3

              temp1 = wwc_w(1,2)+wwc_w(8,2)
              temp2 = wwc_w(2,2)+wwc_w(7,2)
              temp3 = wwc_w(3,2)+wwc_w(5,2)
              qql(1) = rem(1,1)*temp1 + rem(1,3)*temp3 + wwc_w(4,2)
              qql(5) = rem(5,1)*temp1 + rem(5,3)*temp3
              qql(6) = wwc_w(6,2) + wwc_w(9,2)
              qql(7) = rem(7,1)*temp1 + rem(7,2)*temp2 + rem(7,3)*temp3
              qql(8) = rem(8,1)*temp1 + rem(8,2)*temp2 + rem(8,3)*temp3
              qql(9) = ch*(-wwc_w(6,2) + wwc_w(9,2)) 

              minvalue = min(qql(1),qqr(1),qql(5),qqr(5))
              romaxvalue = max(ww(1,1),ww(1,2),ww(1,3))
              prmaxvalue = max(ww(5,1),ww(5,2),ww(5,3))
              minvalue = min(minvalue,ro1-romaxvalue*floor,pr1-prmaxvalue*floor)
              smv = sign(1d0,minvalue)
              psmv = max(0d0,smv)
              msmv = max(0d0,-smv)

              row(i,j,k-1,2) = qql(1)*psmv + ro1*msmv
              row(i,j,k  ,1) = qqr(1)*psmv + ro1*msmv
              vxw(i,j,k-1,2) = qql(2)*psmv + vx1*msmv
              vxw(i,j,k  ,1) = qqr(2)*psmv + vx1*msmv
              vyw(i,j,k-1,2) = qql(3)
              vyw(i,j,k  ,1) = qqr(3)
              vzw(i,j,k-1,2) = qql(4)
              vzw(i,j,k  ,1) = qqr(4)
              prw(i,j,k-1,2) = qql(5)*psmv + pr1*msmv
              prw(i,j,k  ,1) = qqr(5)*psmv + pr1*msmv
              bxw(i,j,k-1,2) = qql(6)*psmv + bx1*msmv
              bxw(i,j,k  ,1) = qqr(6)*psmv + bx1*msmv
              byw(i,j,k-1,2) = qql(7)
              byw(i,j,k  ,1) = qqr(7)
              bzw(i,j,k-1,2) = qql(8)
              bzw(i,j,k  ,1) = qqr(8)
              phiw(i,j,k-1,2) = qql(9)*psmv + phi1*msmv
              phiw(i,j,k  ,1) = qqr(9)*psmv + phi1*msmv
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  endif

  end subroutine lr_state__MSCL2


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
  real(8),dimension(nwave,2) :: wwc_w
  ! ww(*,1) = qqm2, ..
  real(8),dimension(nwave,5) :: ww,wwc
  real(8),dimension(nwave,nwave) :: lem,rem
  real(8) :: wwor,minvalue,smv,psmv,msmv
  real(8) :: ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1
  real(8) :: temp1,temp2,temp3,ich

!----parameter
  real(8),parameter :: B1 = 0.016666666667
  real(8),parameter :: B2 = 1.333333333333
  real(8),parameter :: Alpha = 4.0d0
  real(8),parameter :: Epsm = 0.0000000001d0

!----function
  integer :: i,j,k,l,n
  real(8) :: djm1,dj,djp1,dm4jph,dm4jmh
  real(8) :: qqul,qqmd,qqlc,qqmin,qqmax
  real(8) :: qqlr

!----Recipe for rarefaction ---
  real(8), parameter :: floor=1d-2
  real(8) :: romaxvalue,prmaxvalue

  ich = 0.5d0/ch

  if(mdir == 1)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,l,n,ww,ro1,vx1,vy1,vz1,pr1,bx1,by1,bz1,phi1,lem,rem,wwc,&
     !$OMP         temp1,temp2,temp3,wwor,djm1,dj,djp1,dm4jph,dm4jmh,qqul,qqmd,qqlc,qqmin,qqmax,wwc_w,qqlr,&
     !$OMP         qqr,qql,minvalue,romaxvalue,prmaxvalue,smv,psmv,msmv)
     do k=3,kx-2
        do j=3,jx-2
           ww(1,2:5) = ro(1:4,j,k)
           ww(2,2:5) = vx(1:4,j,k)
           ww(3,2:5) = vy(1:4,j,k)
           ww(4,2:5) = vz(1:4,j,k)
           ww(5,2:5) = pr(1:4,j,k)
           ww(6,2:5) = bx(1:4,j,k)
           ww(7,2:5) = by(1:4,j,k)
           ww(8,2:5) = bz(1:4,j,k)
           ww(9,2:5) = phi(1:4,j,k)
           do i=3,ix-2
              ww(1:9,1) = ww(1:9,2)
              ww(1:9,2) = ww(1:9,3)
              ww(1:9,3) = ww(1:9,4)
              ww(1:9,4) = ww(1:9,5)

              ww(1,5) = ro(i+2,j,k)
              ww(2,5) = vx(i+2,j,k)
              ww(3,5) = vy(i+2,j,k)
              ww(4,5) = vz(i+2,j,k)
              ww(5,5) = pr(i+2,j,k)
              ww(6,5) = bx(i+2,j,k)
              ww(7,5) = by(i+2,j,k)
              ww(8,5) = bz(i+2,j,k)
              ww(9,5) = phi(i+2,j,k)

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
              call esystem_glmmhd(lem,rem,ro1,pr1,bx1,by1,bz1,gm)

              do l=1,5
                 wwc(4,l) = ww(1,l)+lem(4,5)*ww(5,l)

                 temp1 = lem(1,2)*ww(2,l) + lem(1,3)*ww(3,l) + lem(1,4)*ww(4,l)
                 temp2 = lem(1,5)*ww(5,l) + lem(1,7)*ww(7,l) + lem(1,8)*ww(8,l)
                 wwc(1,l) = temp1 + temp2
                 wwc(8,l) = -temp1 + temp2

                 temp1 = lem(2,3)*ww(3,l) + lem(2,4)*ww(4,l)
                 temp2 = lem(2,7)*ww(7,l) + lem(2,8)*ww(8,l)
                 wwc(2,l) = temp1 + temp2
                 wwc(7,l) = -temp1 + temp2
                      
                 temp1 = lem(3,2)*ww(2,l) + lem(3,3)*ww(3,l) + lem(3,4)*ww(4,l)
                 temp2 = lem(3,5)*ww(5,l) + lem(3,7)*ww(7,l) + lem(3,8)*ww(8,l)
                 wwc(3,l) = temp1 +temp2
                 wwc(5,l) = -temp1 + temp2
                 
                 temp1 = 0.5d0*ww(6,l)
                 temp2 = ich*ww(9,l) 
                 wwc(6,l) = temp1 - temp2
                 wwc(9,l) = temp1 + temp2
              enddo
              
              ! mp5
              do n=1,nwave
                 !right-hand left state
                 wwor = B1*(ccx(1,2,i)*wwc(n,1)+ccx(2,2,i)*wwc(n,2) &
                      + ccx(3,2,i)*wwc(n,3) + ccx(4,2,i)*wwc(n,4) &
                      + ccx(5,2,i)*wwc(n,5))

                 djm1 = wwc(n,1)-2.0d0*wwc(n,2)+wwc(n,3)
                 dj = wwc(n,2)-2.0d0*wwc(n,3)+wwc(n,4)
                 djp1 = wwc(n,3)-2.0d0*wwc(n,4)+wwc(n,5)
                 
                 dm4jph = minmod4(4.0d0*dj-djp1,4.0d0*djp1-dj,dj,djp1)
                 dm4jmh = minmod4(4.0d0*dj-djm1,4.0d0*djm1-dj,dj,djm1)

                 qqul = wwc(n,3)+Alpha*(wwc(n,3)-wwc(n,2))
                 qqmd = 0.5d0*(wwc(n,3)+wwc(n,4) - dm4jph)
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,2)) + B2*dm4jmh
                 
                 qqmin = max(min(wwc(n,3),wwc(n,4),qqmd),min(wwc(n,3),qqul,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,4),qqmd),max(wwc(n,3),qqul,qqlc))
                 wwc_w(n,1) = wwor + minmod((qqmin-wwor),(qqmax-wwor))

                 !left-hand right state
                 wwor = B1*(ccx(5,1,i-1)*wwc(n,5)+ccx(4,1,i-1)*wwc(n,4) &
                      + ccx(3,1,i-1)*wwc(n,3) + ccx(2,1,i-1)*wwc(n,2) &
                      + ccx(1,1,i-1)*wwc(n,1))

                 qqlr = wwc(n,3)+Alpha*(wwc(n,3)-wwc(n,4))
                 qqmd = 0.5d0*(wwc(n,3)+wwc(n,2) - dm4jmh)
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,4)) + B2*dm4jph
                 
                 qqmin = max(min(wwc(n,3),wwc(n,2),qqmd),min(wwc(n,3),qqlr,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,2),qqmd),max(wwc(n,3),qqlr,qqlc))
                 wwc_w(n,2) = wwor + minmod((qqmin-wwor),(qqmax-wwor))
              end do

              ! characteristic to primitive
              temp1 = wwc_w(1,1)-wwc_w(8,1)
              temp2 = wwc_w(2,1)-wwc_w(7,1)
              temp3 = wwc_w(3,1)-wwc_w(5,1)
              qqr(2) = rem(2,1)*temp1 + rem(2,3)*temp3
              qqr(3) = rem(3,1)*temp1 + rem(3,2)*temp2 + rem(3,3)*temp3
              qqr(4) = rem(4,1)*temp1 + rem(4,2)*temp2 + rem(4,3)*temp3

              temp1 = wwc_w(1,1)+wwc_w(8,1)
              temp2 = wwc_w(2,1)+wwc_w(7,1)
              temp3 = wwc_w(3,1)+wwc_w(5,1)
              qqr(1) = rem(1,1)*temp1 + rem(1,3)*temp3 + wwc_w(4,1)
              qqr(5) = rem(5,1)*temp1 + rem(5,3)*temp3
              qqr(6) = wwc_w(6,1) + wwc_w(9,1)
              qqr(7) = rem(7,1)*temp1 + rem(7,2)*temp2 + rem(7,3)*temp3
              qqr(8) = rem(8,1)*temp1 + rem(8,2)*temp2 + rem(8,3)*temp3
              qqr(9) = ch*(-wwc_w(6,1) + wwc_w(9,1)) 

              temp1 = wwc_w(1,2)-wwc_w(8,2)
              temp2 = wwc_w(2,2)-wwc_w(7,2)
              temp3 = wwc_w(3,2)-wwc_w(5,2)
              qql(2) = rem(2,1)*temp1 + rem(2,3)*temp3
              qql(3) = rem(3,1)*temp1 + rem(3,2)*temp2 + rem(3,3)*temp3
              qql(4) = rem(4,1)*temp1 + rem(4,2)*temp2 + rem(4,3)*temp3

              temp1 = wwc_w(1,2)+wwc_w(8,2)
              temp2 = wwc_w(2,2)+wwc_w(7,2)
              temp3 = wwc_w(3,2)+wwc_w(5,2)
              qql(1) = rem(1,1)*temp1 + rem(1,3)*temp3 + wwc_w(4,2)
              qql(5) = rem(5,1)*temp1 + rem(5,3)*temp3
              qql(6) = wwc_w(6,2) + wwc_w(9,2)
              qql(7) = rem(7,1)*temp1 + rem(7,2)*temp2 + rem(7,3)*temp3
              qql(8) = rem(8,1)*temp1 + rem(8,2)*temp2 + rem(8,3)*temp3
              qql(9) = ch*(-wwc_w(6,2) + wwc_w(9,2)) 

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
              vyw(i-1,j,k,2) = qql(3)
              vyw(i  ,j,k,1) = qqr(3)
              vzw(i-1,j,k,2) = qql(4)
              vzw(i  ,j,k,1) = qqr(4)
              prw(i-1,j,k,2) = qql(5)*psmv + pr1*msmv
              prw(i  ,j,k,1) = qqr(5)*psmv + pr1*msmv
              bxw(i-1,j,k,2) = qql(6)*psmv + bx1*msmv
              bxw(i  ,j,k,1) = qqr(6)*psmv + bx1*msmv
              byw(i-1,j,k,2) = qql(7)
              byw(i  ,j,k,1) = qqr(7)
              bzw(i-1,j,k,2) = qql(8)
              bzw(i  ,j,k,1) = qqr(8)
              phiw(i-1,j,k,2) = qql(9)*psmv + phi1*msmv
              phiw(i  ,j,k,1) = qqr(9)*psmv + phi1*msmv
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  else if(mdir == 2)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,l,ww,ro1,vx1,vy1,vz1,pr1,bx1,by1,bz1,phi1,lem,rem,wwc,&
     !$OMP         temp1,temp2,n,wwor,djm1,dj,djp1,dm4jph,dm4jmh,qqul,qqmd,qqlc,qqmin,qqmax,&
     !$OMP         wwc_w,qqlr,temp3,qqr,qql,minvalue,romaxvalue,prmaxvalue,smv,psmv,msmv)
     do k=3,kx-2
        do i=3,ix-2
           ww(1,2:5) = ro(i,1:4,k)
           ww(2,2:5) = vx(i,1:4,k)
           ww(3,2:5) = vy(i,1:4,k)
           ww(4,2:5) = vz(i,1:4,k)
           ww(5,2:5) = pr(i,1:4,k)
           ww(6,2:5) = bx(i,1:4,k)
           ww(7,2:5) = by(i,1:4,k)
           ww(8,2:5) = bz(i,1:4,k)
           ww(9,2:5) = phi(i,1:4,k)
           do j=3,jx-2
              ww(1:9,1) = ww(1:9,2)
              ww(1:9,2) = ww(1:9,3)
              ww(1:9,3) = ww(1:9,4)
              ww(1:9,4) = ww(1:9,5)

              ww(1,5) = ro(i,j+2,k)
              ww(2,5) = vx(i,j+2,k)
              ww(3,5) = vy(i,j+2,k)
              ww(4,5) = vz(i,j+2,k)
              ww(5,5) = pr(i,j+2,k)
              ww(6,5) = bx(i,j+2,k)
              ww(7,5) = by(i,j+2,k)
              ww(8,5) = bz(i,j+2,k)
              ww(9,5) = phi(i,j+2,k)

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
              call esystem_glmmhd(lem,rem,ro1,pr1,bx1,by1,bz1,gm)

              do l=1,5
                 wwc(4,l) = ww(1,l)+lem(4,5)*ww(5,l)

                 temp1 = lem(1,2)*ww(2,l) + lem(1,3)*ww(3,l) + lem(1,4)*ww(4,l)
                 temp2 = lem(1,5)*ww(5,l) + lem(1,7)*ww(7,l) + lem(1,8)*ww(8,l)
                 wwc(1,l) = temp1 + temp2
                 wwc(8,l) = -temp1 + temp2

                 temp1 = lem(2,3)*ww(3,l) + lem(2,4)*ww(4,l)
                 temp2 = lem(2,7)*ww(7,l) + lem(2,8)*ww(8,l)
                 wwc(2,l) = temp1 + temp2
                 wwc(7,l) = -temp1 + temp2
                      
                 temp1 = lem(3,2)*ww(2,l) + lem(3,3)*ww(3,l) + lem(3,4)*ww(4,l)
                 temp2 = lem(3,5)*ww(5,l) + lem(3,7)*ww(7,l) + lem(3,8)*ww(8,l)
                 wwc(3,l) = temp1 +temp2
                 wwc(5,l) = -temp1 + temp2
                 
                 temp1 = 0.5d0*ww(6,l)
                 temp2 = ich*ww(9,l) 
                 wwc(6,l) = temp1 - temp2
                 wwc(9,l) = temp1 + temp2
              enddo

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
                 qqmd = 0.5d0*(wwc(n,3)+wwc(n,4) - dm4jph)
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,2)) + B2*dm4jmh
                 
                 qqmin = max(min(wwc(n,3),wwc(n,4),qqmd),min(wwc(n,3),qqul,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,4),qqmd),max(wwc(n,3),qqul,qqlc))

                 wwc_w(n,1) = wwor + minmod((qqmin-wwor),(qqmax-wwor))
                 
                 ! right state
                 wwor = B1*(ccy(5,1,j-1)*wwc(n,5)+ccy(4,1,j-1)*wwc(n,4) &
                      + ccy(3,1,j-1)*wwc(n,3) + ccy(2,1,j-1)*wwc(n,2) &
                      + ccy(1,1,j-1)*wwc(n,1))

                 qqlr = wwc(n,3)+Alpha*(wwc(n,3)-wwc(n,4))
                 qqmd = 0.5d0*(wwc(n,3)+wwc(n,2) - dm4jmh)
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,4)) + B2*dm4jph
                 
                 qqmin = max(min(wwc(n,3),wwc(n,2),qqmd),min(wwc(n,3),qqlr,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,2),qqmd),max(wwc(n,3),qqlr,qqlc))
                 
                 wwc_w(n,2) = wwor + minmod((qqmin-wwor),(qqmax-wwor))
              end do
              
              ! characteristic to primitive
              temp1 = wwc_w(1,1)-wwc_w(8,1)
              temp2 = wwc_w(2,1)-wwc_w(7,1)
              temp3 = wwc_w(3,1)-wwc_w(5,1)
              qqr(2) = rem(2,1)*temp1 + rem(2,3)*temp3
              qqr(3) = rem(3,1)*temp1 + rem(3,2)*temp2 + rem(3,3)*temp3
              qqr(4) = rem(4,1)*temp1 + rem(4,2)*temp2 + rem(4,3)*temp3

              temp1 = wwc_w(1,1)+wwc_w(8,1)
              temp2 = wwc_w(2,1)+wwc_w(7,1)
              temp3 = wwc_w(3,1)+wwc_w(5,1)
              qqr(1) = rem(1,1)*temp1 + rem(1,3)*temp3 + wwc_w(4,1)
              qqr(5) = rem(5,1)*temp1 + rem(5,3)*temp3
              qqr(6) = wwc_w(6,1) + wwc_w(9,1)
              qqr(7) = rem(7,1)*temp1 + rem(7,2)*temp2 + rem(7,3)*temp3
              qqr(8) = rem(8,1)*temp1 + rem(8,2)*temp2 + rem(8,3)*temp3
              qqr(9) = ch*(-wwc_w(6,1) + wwc_w(9,1)) 

              temp1 = wwc_w(1,2)-wwc_w(8,2)
              temp2 = wwc_w(2,2)-wwc_w(7,2)
              temp3 = wwc_w(3,2)-wwc_w(5,2)
              qql(2) = rem(2,1)*temp1 + rem(2,3)*temp3
              qql(3) = rem(3,1)*temp1 + rem(3,2)*temp2 + rem(3,3)*temp3
              qql(4) = rem(4,1)*temp1 + rem(4,2)*temp2 + rem(4,3)*temp3

              temp1 = wwc_w(1,2)+wwc_w(8,2)
              temp2 = wwc_w(2,2)+wwc_w(7,2)
              temp3 = wwc_w(3,2)+wwc_w(5,2)
              qql(1) = rem(1,1)*temp1 + rem(1,3)*temp3 + wwc_w(4,2)
              qql(5) = rem(5,1)*temp1 + rem(5,3)*temp3
              qql(6) = wwc_w(6,2) + wwc_w(9,2)
              qql(7) = rem(7,1)*temp1 + rem(7,2)*temp2 + rem(7,3)*temp3
              qql(8) = rem(8,1)*temp1 + rem(8,2)*temp2 + rem(8,3)*temp3
              qql(9) = ch*(-wwc_w(6,2) + wwc_w(9,2)) 

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
              vyw(i,j-1,k,2) = qql(3)
              vyw(i,j  ,k,1) = qqr(3)
              vzw(i,j-1,k,2) = qql(4)
              vzw(i,j  ,k,1) = qqr(4)
              prw(i,j-1,k,2) = qql(5)*psmv + pr1*msmv
              prw(i,j  ,k,1) = qqr(5)*psmv + pr1*msmv
              bxw(i,j-1,k,2) = qql(6)*psmv + bx1*msmv
              bxw(i,j  ,k,1) = qqr(6)*psmv + bx1*msmv
              byw(i,j-1,k,2) = qql(7)
              byw(i,j  ,k,1) = qqr(7)
              bzw(i,j-1,k,2) = qql(8)
              bzw(i,j  ,k,1) = qqr(8)
              phiw(i,j-1,k,2) = qql(9)*psmv + phi1*msmv
              phiw(i,j  ,k,1) = qqr(9)*psmv + phi1*msmv
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,l,ww,ro1,vx1,vy1,vz1,pr1,bx1,by1,bz1,phi1,lem,rem,wwc,&
     !$OMP         temp1,temp2,n,wwor,djm1,dj,djp1,dm4jph,dm4jmh,qqul,qqmd,qqlc,qqmin,qqmax,&
     !$OMP         wwc_w,qqlr,temp3,qqr,qql,minvalue,romaxvalue,prmaxvalue,smv,psmv,msmv)
     do j=3,jx-2
        do i=3,ix-2
           ww(1,2:5) = ro(i,j,1:4)
           ww(2,2:5) = vx(i,j,1:4)
           ww(3,2:5) = vy(i,j,1:4)
           ww(4,2:5) = vz(i,j,1:4)
           ww(5,2:5) = pr(i,j,1:4)
           ww(6,2:5) = bx(i,j,1:4)
           ww(7,2:5) = by(i,j,1:4)
           ww(8,2:5) = bz(i,j,1:4)
           ww(9,2:5) = phi(i,j,1:4)
           do k=3,kx-2
              ww(1:9,1) = ww(1:9,2)
              ww(1:9,2) = ww(1:9,3)
              ww(1:9,3) = ww(1:9,4)
              ww(1:9,4) = ww(1:9,5)

              ww(1,5) = ro(i,j,k+2)
              ww(2,5) = vx(i,j,k+2)
              ww(3,5) = vy(i,j,k+2)
              ww(4,5) = vz(i,j,k+2)
              ww(5,5) = pr(i,j,k+2)
              ww(6,5) = bx(i,j,k+2)
              ww(7,5) = by(i,j,k+2)
              ww(8,5) = bz(i,j,k+2)
              ww(9,5) = phi(i,j,k+2)

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
              call esystem_glmmhd(lem,rem,ro1,pr1,bx1,by1,bz1,gm)

              do l=1,5
                 wwc(4,l) = ww(1,l)+lem(4,5)*ww(5,l)

                 temp1 = lem(1,2)*ww(2,l) + lem(1,3)*ww(3,l) + lem(1,4)*ww(4,l)
                 temp2 = lem(1,5)*ww(5,l) + lem(1,7)*ww(7,l) + lem(1,8)*ww(8,l)
                 wwc(1,l) = temp1 + temp2
                 wwc(8,l) = -temp1 + temp2

                 temp1 = lem(2,3)*ww(3,l) + lem(2,4)*ww(4,l)
                 temp2 = lem(2,7)*ww(7,l) + lem(2,8)*ww(8,l)
                 wwc(2,l) = temp1 + temp2
                 wwc(7,l) = -temp1 + temp2
                      
                 temp1 = lem(3,2)*ww(2,l) + lem(3,3)*ww(3,l) + lem(3,4)*ww(4,l)
                 temp2 = lem(3,5)*ww(5,l) + lem(3,7)*ww(7,l) + lem(3,8)*ww(8,l)
                 wwc(3,l) = temp1 +temp2
                 wwc(5,l) = -temp1 + temp2
                 
                 temp1 = 0.5d0*ww(6,l)
                 temp2 = ich*ww(9,l) 
                 wwc(6,l) = temp1 - temp2
                 wwc(9,l) = temp1 + temp2
              enddo

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
                 qqmd = 0.5d0*(wwc(n,3)+wwc(n,4) - dm4jph)
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,2)) + B2*dm4jmh
                 
                 qqmin = max(min(wwc(n,3),wwc(n,4),qqmd),min(wwc(n,3),qqul,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,4),qqmd),max(wwc(n,3),qqul,qqlc))

                 wwc_w(n,1) = wwor + minmod((qqmin-wwor),(qqmax-wwor))

              ! right state
                 wwor = B1*(ccz(5,1,k-1)*wwc(n,5)+ccz(4,1,k-1)*wwc(n,4) &
                      + ccz(3,1,k-1)*wwc(n,3) + ccz(2,1,k-1)*wwc(n,2) &
                      + ccz(1,1,k-1)*wwc(n,1))

                 qqlr = wwc(n,3)+Alpha*(wwc(n,3)-wwc(n,4))
                 qqmd = 0.5d0*(wwc(n,3)+wwc(n,2) - dm4jmh)
                 qqlc = wwc(n,3) + 0.5d0*(wwc(n,3)-wwc(n,4)) + B2*dm4jph
                 
                 qqmin = max(min(wwc(n,3),wwc(n,2),qqmd),min(wwc(n,3),qqlr,qqlc))
                 qqmax = min(max(wwc(n,3),wwc(n,2),qqmd),max(wwc(n,3),qqlr,qqlc))
                 
                 wwc_w(n,2) = wwor + minmod((qqmin-wwor),(qqmax-wwor))
              end do
              
              ! characteristic to primitive
              temp1 = wwc_w(1,1)-wwc_w(8,1)
              temp2 = wwc_w(2,1)-wwc_w(7,1)
              temp3 = wwc_w(3,1)-wwc_w(5,1)
              qqr(2) = rem(2,1)*temp1 + rem(2,3)*temp3
              qqr(3) = rem(3,1)*temp1 + rem(3,2)*temp2 + rem(3,3)*temp3
              qqr(4) = rem(4,1)*temp1 + rem(4,2)*temp2 + rem(4,3)*temp3

              temp1 = wwc_w(1,1)+wwc_w(8,1)
              temp2 = wwc_w(2,1)+wwc_w(7,1)
              temp3 = wwc_w(3,1)+wwc_w(5,1)
              qqr(1) = rem(1,1)*temp1 + rem(1,3)*temp3 + wwc_w(4,1)
              qqr(5) = rem(5,1)*temp1 + rem(5,3)*temp3
              qqr(6) = wwc_w(6,1) + wwc_w(9,1)
              qqr(7) = rem(7,1)*temp1 + rem(7,2)*temp2 + rem(7,3)*temp3
              qqr(8) = rem(8,1)*temp1 + rem(8,2)*temp2 + rem(8,3)*temp3
              qqr(9) = ch*(-wwc_w(6,1) + wwc_w(9,1)) 

              temp1 = wwc_w(1,2)-wwc_w(8,2)
              temp2 = wwc_w(2,2)-wwc_w(7,2)
              temp3 = wwc_w(3,2)-wwc_w(5,2)
              qql(2) = rem(2,1)*temp1 + rem(2,3)*temp3
              qql(3) = rem(3,1)*temp1 + rem(3,2)*temp2 + rem(3,3)*temp3
              qql(4) = rem(4,1)*temp1 + rem(4,2)*temp2 + rem(4,3)*temp3

              temp1 = wwc_w(1,2)+wwc_w(8,2)
              temp2 = wwc_w(2,2)+wwc_w(7,2)
              temp3 = wwc_w(3,2)+wwc_w(5,2)
              qql(1) = rem(1,1)*temp1 + rem(1,3)*temp3 + wwc_w(4,2)
              qql(5) = rem(5,1)*temp1 + rem(5,3)*temp3
              qql(6) = wwc_w(6,2) + wwc_w(9,2)
              qql(7) = rem(7,1)*temp1 + rem(7,2)*temp2 + rem(7,3)*temp3
              qql(8) = rem(8,1)*temp1 + rem(8,2)*temp2 + rem(8,3)*temp3
              qql(9) = ch*(-wwc_w(6,2) + wwc_w(9,2)) 

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
              vyw(i,j,k-1,2) = qql(3)
              vyw(i,j,k  ,1) = qqr(3)
              vzw(i,j,k-1,2) = qql(4)
              vzw(i,j,k  ,1) = qqr(4)
              prw(i,j,k-1,2) = qql(5)*psmv + pr1*msmv
              prw(i,j,k  ,1) = qqr(5)*psmv + pr1*msmv
              bxw(i,j,k-1,2) = qql(6)*psmv + bx1*msmv
              bxw(i,j,k  ,1) = qqr(6)*psmv + bx1*msmv
              byw(i,j,k-1,2) = qql(7)
              byw(i,j,k  ,1) = qqr(7)
              bzw(i,j,k-1,2) = qql(8)
              bzw(i,j,k  ,1) = qqr(8)
              phiw(i,j,k-1,2) = qql(9)*psmv + phi1*msmv
              phiw(i,j,k  ,1) = qqr(9)*psmv + phi1*msmv
           end do
        end do
     end do
     !$OMP END PARALLEL DO
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

  !$OMP PARALLEL DO &
  !$OMP PRIVATE(r,j,tmp3,m,tmp2,l,tmp1,q,tmp4)
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
  !$OMP END PARALLEL DO

  end subroutine reconstructionConstant


  subroutine esystem_glmmhd(lem,rem,ro,pr,bx,by,bz,gm)

  real(8),intent(in) :: gm
  real(8),intent(in) :: ro,pr
  real(8),intent(in) :: bx,by,bz
  real(8),dimension(9,9),intent(out) :: lem,rem

  real(8) :: roi,btsq,vaxsq,asq
  real(8) :: ct2,tsum,tdif,cf2_cs2
  real(8) :: cfsq,cf,cssq,cs
  real(8) :: bt,ibt,sbt,psbt,msbt,bet2,bet3
  real(8) :: sfs,psfs,msfs,sas,saf,eps,ifs
  real(8) :: alpha_f,alpha_s
  real(8) :: na,qf,qs,af_prm,as_prm
  real(8) :: sqrtro,sqrtroi,s,a,af,as,sqrt2i
  real(8) :: temp1,temp2

  sqrt2i = 0.5D0*sqrt(2.D0)
  roi = 1.0d0/ro
  sqrtroi = sqrt(roi)
  btsq = by*by+bz*bz
  vaxsq = bx*bx*roi
  asq = gm*pr*roi

! fast and slow speed
  ct2 = btsq*roi
  tsum = vaxsq + ct2 + asq
  tdif = vaxsq + ct2 - asq
  cf2_cs2 = sqrt(tdif*tdif + 4.0d0*asq*ct2)
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

  temp1 = ro*alpha_f
  temp2 = ro*alpha_s

  rem(1,1) = temp1
  rem(1,3) = temp2

  rem(2,1) = -cf*alpha_f
  rem(2,3) = -cs*alpha_s

  rem(3,1) = qs*bet2
  rem(3,2) = -bet3
  rem(3,3) = -qf*bet2

  rem(4,1) = qs*bet3
  rem(4,2) = bet2
  rem(4,3) = -qf*bet3

  rem(5,1) = asq*temp1
  rem(5,3) = asq*temp2

  temp1 = s*sqrtro
  rem(7,1) = as*bet2
  rem(7,2) = -bet3*temp1
  rem(7,3) = -af*bet2

  rem(8,1) = as*bet3
  rem(8,2) = bet2*temp1
  rem(8,3) = -af*bet3

! Left eigenvector
  na = 0.5d0/asq
  qf = na*qf
  qs = na*qs
  af_prm = na*af*roi
  as_prm = na*as*roi

  temp1 = na*roi
  lem(1,2) = -na*cf*alpha_f
  lem(1,3) = qs*bet2
  lem(1,4) = qs*bet3
  lem(1,5) = alpha_f*temp1
  lem(1,7) = as_prm*bet2
  lem(1,8) = as_prm*bet3

  temp2 = 0.5d0*s*sqrtroi
  lem(2,3) = -0.5d0*bet3
  lem(2,4) = 0.5d0*bet2
  lem(2,7) = -bet3*temp2
  lem(2,8) = bet2*temp2

  lem(3,2) = -na*cs*alpha_s
  lem(3,3) = -qf*bet2
  lem(3,4) = -qf*bet3
  lem(3,5) = alpha_s*temp1
  lem(3,7) = -af_prm*bet2
  lem(3,8) = -af_prm*bet3

  lem(4,5) = -2d0*na

  end subroutine esystem_glmmhd


  real(8) function minmod(x,y)
    real(8), intent(in) :: x,y
   
    minmod = 0.5d0*(sign(1.0d0,x)+sign(1.0d0,y))*min(dabs(x),dabs(y))

  end function minmod


  real(8) function minmod4(d1,d2,d3,d4)
    real(8), intent(in) :: d1,d2,d3,d4
    real(8) :: sign1

    sign1 = sign(1d0,d1)
    minmod4 = 0.125d0*(sign1+sign(1.0d0,d2)) &
              *dabs( (sign1+sign(1.0d0,d3)) &
                    *(sign1+sign(1.0d0,d4))) &
              *min(dabs(d1),dabs(d2),dabs(d3),dabs(d4))

  end function minmod4


  real(8) function MC2(qqr,qql,qqc)

    real(8), intent(in) :: qqr,qql,qqc

    MC2 = minmod(qqc,minmod(qql,qqr))

  end function MC2


end module lr_state
