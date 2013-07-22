function getMP5value(lr,qqm2,qqm1,qq0,qqp1,qqp2,qqor,Alpha)
  implicit none

  integer,intent(in) :: lr

  real(8) :: qqm2,qqm1,qq0,qqp1,qqp2,qqor
  real(8),intent(in) :: Alpha
  
!----function
  real(8) :: minmod4,d1,d2,d3,d4
  real(8) :: minmod,x,y
  real(8) :: median
  
  real(8) :: djm1,dj,djp1,dm4jph,dm4jmh,djpp1,dm4jpph
  real(8) :: qqul,qqav,qqmd,qqlc,qqmin,qqmax,qqmin1,qqmax1
  real(8) :: djp2,qqlr

  real(8) :: getMP5value

  real(8),parameter :: B1 = 0.016666666667
  real(8),parameter :: B2 = 1.333333333333

  minmod4(d1,d2,d3,d4) = 0.125d0*(sign(1.0d0,d1)+sign(1.0d0,d2))* &
       dabs((sign(1.0d0,d1) + sign(1.0d0,d3))* &
       (sign(1.0d0,d1)+sign(1.0d0,d4))) &
       *min(dabs(d1),dabs(d2),dabs(d3),dabs(d4))

  minmod(x,y) = 0.5d0*(sign(1.0d0,x)+sign(1.0d0,y))*min(dabs(x),dabs(y))

  if(lr .eq. 1)then
     djm1 = qqm2-2.0d0*qqm1+qq0
     dj = qqm1-2.0d0*qq0+qqp1
     djp1 = qq0-2.0d0*qqp1+qqp2

     dm4jph = minmod4(4.0d0*dj-djp1,4.0d0*djp1-dj,dj,djp1)
     dm4jmh = minmod4(4.0d0*dj-djm1,4.0d0*djm1-dj,dj,djm1)
     
     qqul = qq0+Alpha*(qq0-qqm1)
     qqlr = qqp1+Alpha*(qqp1-qqp2)
     
     qqav = 0.5d0*(qq0+qqp1)
     qqmd = qqav - 0.5d0*dm4jph
     qqlc = qq0 + 0.5d0*(qq0-qqm1) + B2*dm4jmh
          
     qqmin = max(min(qq0,qqp1,qqmd),min(qq0,qqul,qqlc))
     qqmax = min(max(qq0,qqp1,qqmd),max(qq0,qqul,qqlc))

!     getMP5value = median(qqor,qqmin,qqmax)
     getMP5value = qqor + minmod((qqmin-qqor),(qqmax-qqor))
  else
     djm1 = qqm2-2.0d0*qqm1+qq0
     dj = qqm1-2.0d0*qq0+qqp1
     djp1 = qq0-2.0d0*qqp1+qqp2
     
     dm4jph = minmod4(4.0d0*dj-djp1,4.0d0*djp1-dj,dj,djp1)
     dm4jmh = minmod4(4.0d0*dj-djm1,4.0d0*djm1-dj,dj,djm1)

     qqul = qqm1+Alpha*(qqm1-qqm2)
     qqlr = qq0+Alpha*(qq0-qqp1)

     qqav = 0.5d0*(qq0+qqm1)
     qqmd = qqav - 0.5d0*dm4jmh
     qqlc = qq0 + 0.5d0*(qq0-qqp1) + B2*dm4jph
     
     qqmin = max(min(qq0,qqm1,qqmd),min(qq0,qqlr,qqlc))
     qqmax = min(max(qq0,qqm1,qqmd),max(qq0,qqlr,qqlc))

!     getMP5value = median(qqor,qqmin,qqmax)
     getMP5value = qqor + minmod((qqmin-qqor),(qqmax-qqor))
  end if
end function getMP5value
