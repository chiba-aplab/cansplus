function minmod4(d1,d2,d3,d4)

  implicit none

  real(8) ,intent(in) :: d1,d2,d3,d4

  real(8) :: ss

  real(8) :: min1,min2
  real(8) :: minmod4

!!$  ss = 0.5d0*(sign(1.0d0,d1)+sign(1.0d0,d2)) &
!!$       *dabs(0.5d0*(sign(1.0d0,d1)+sign(1.0d0,d3)) &
!!$       *0.5d0*(sign(1.0d0,d1)+sign(1.0d0,d4)))
!!$
!!$  min1 = min(dabs(d1),dabs(d2))
!!$  min2 = min(min1,dabs(d3))

!  minmod4 = ss*min(min2, dabs(d4))
  minmod4 = 0.125d0*(sign(1.0d0,d1)+sign(1.0d0,d2))* &
       dabs((sign(1.0d0,d1) + sign(1.0d0,d3))* &
       (sign(1.0d0,d1)+sign(1.0d0,d4))) &
       *min(dabs(d1),dabs(d2),dabs(d3),dabs(d4))
end function minmod4

!!$program main
!!$
!!$  implicit none
!!$
!!$  real(8) :: d1,d2,d3,d4
!!$
!!$  real(8) :: minmod4
!!$
!!$  real(8) :: answer
!!$
!!$  write(*,*) "all positive"
!!$
!!$  d1 = 2.0d0
!!$  d2 = 2.0d0
!!$  d3 = 3.0d0
!!$  d4 = 4.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = 2.0d0
!!$  d2 = 3.0d0
!!$  d3 = 4.0d0
!!$  d4 = 2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = 3.0d0
!!$  d2 = 4.0d0
!!$  d3 = 2.0d0
!!$  d4 = 2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = 4.0d0
!!$  d2 = 2.0d0
!!$  d3 = 2.0d0
!!$  d4 = 3.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = 2.0d0
!!$  d2 = 2.0d0
!!$  d3 = 3.0d0
!!$  d4 = 4.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = 2.0d0
!!$  d2 = 2.0d0
!!$  d3 = 3.0d0
!!$  d4 = 2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = 2.0d0
!!$  d2 = 2.0d0
!!$  d3 = 2.0d0
!!$  d4 = 4.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = 2.0d0
!!$  d2 = 2.0d0
!!$  d3 = 2.0d0
!!$  d4 = 4.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = 2.0d0
!!$  d2 = 2.0d0
!!$  d3 = 2.0d0
!!$  d4 = 2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = 2.0d0
!!$  d2 = 2.0d0
!!$  d3 = 3.0d0
!!$  d4 = 2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = 2.0d0
!!$  d2 = 2.0d0
!!$  d3 = 2.0d0
!!$  d4 = 2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  write(*,*) "######all negative"
!!$
!!$  d1 = -2.0d0
!!$  d2 = -2.0d0
!!$  d3 = -3.0d0
!!$  d4 = -4.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = -2.0d0
!!$  d2 = -3.0d0
!!$  d3 = -4.0d0
!!$  d4 = -2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = -3.0d0
!!$  d2 = -4.0d0
!!$  d3 = -2.0d0
!!$  d4 = -2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = -4.0d0
!!$  d2 = -2.0d0
!!$  d3 = -2.0d0
!!$  d4 = -3.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = -2.0d0
!!$  d2 = -2.0d0
!!$  d3 = -3.0d0
!!$  d4 = -4.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = -2.0d0
!!$  d2 = -2.0d0
!!$  d3 = -3.0d0
!!$  d4 = -2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = -2.0d0
!!$  d2 = -2.0d0
!!$  d3 = -2.0d0
!!$  d4 = -4.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = -2.0d0
!!$  d2 = -2.0d0
!!$  d3 = -2.0d0
!!$  d4 = -4.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = -2.0d0
!!$  d2 = -2.0d0
!!$  d3 = -2.0d0
!!$  d4 = -2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = -2.0d0
!!$  d2 = -2.0d0
!!$  d3 = -3.0d0
!!$  d4 = -2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = -2.0d0
!!$  d2 = -2.0d0
!!$  d3 = -2.0d0
!!$  d4 = -2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  write(*,*) "########Mix"
!!$
!!$  d1 = -2.0d0
!!$  d2 = 2.0d0
!!$  d3 = 3.0d0
!!$  d4 = 4.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = 2.0d0
!!$  d2 = 3.0d0
!!$  d3 = 4.0d0
!!$  d4 = -2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = 3.0d0
!!$  d2 = 4.0d0
!!$  d3 = -2.0d0
!!$  d4 = 2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = 4.0d0
!!$  d2 = -2.0d0
!!$  d3 = 2.0d0
!!$  d4 = 3.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = -2.0d0
!!$  d2 = -2.0d0
!!$  d3 = 3.0d0
!!$  d4 = 4.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = -2.0d0
!!$  d2 = 2.0d0
!!$  d3 = 3.0d0
!!$  d4 = -2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = -2.0d0
!!$  d2 = 2.0d0
!!$  d3 = -2.0d0
!!$  d4 = 4.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = 2.0d0
!!$  d2 = 2.0d0
!!$  d3 = 2.0d0
!!$  d4 = -4.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = 2.0d0
!!$  d2 = -2.0d0
!!$  d3 = 2.0d0
!!$  d4 = 2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = 2.0d0
!!$  d2 = 2.0d0
!!$  d3 = -3.0d0
!!$  d4 = 2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$  d1 = -2.0d0
!!$  d2 = 2.0d0
!!$  d3 = 2.0d0
!!$  d4 = 2.0d0
!!$  write(*,*) "d1 = ", d1, " d2 = ", d2, " d3 = ", d3, " d4 = ", d4
!!$  answer = minmod4(d1,d2,d3,d4)
!!$  write(*,*) "answer = " , answer
!!$
!!$end program main
