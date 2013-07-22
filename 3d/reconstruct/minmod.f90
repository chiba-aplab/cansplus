function minmod(x,y)
  implicit none

  real(8),intent(in) :: x,y

  real(8) :: minmod

  minmod = 0.5d0*(sign(1.0d0,x)+sign(1.0d0,y))*min(dabs(x),dabs(y))

end function minmod

!!$program main
!!$
!!$  implicit none
!!$
!!$  real(8) :: x,y
!!$  
!!$  real(8) :: minmod
!!$
!!$  real(8) :: answer
!!$
!!$  x = 1.0d0
!!$  y = 2.0d0
!!$
!!$  write(*,*) "x = ", x , " y = ", y
!!$
!!$  answer = minmod(x,y)
!!$
!!$  write(*,*) "answer = " , answer
!!$
!!$  x = -1.0d0
!!$  y = 2.0d0
!!$
!!$  write(*,*) "x = ", x , " y = ", y
!!$
!!$  answer = minmod(x,y)
!!$
!!$  write(*,*) "answer = " , answer
!!$
!!$  
!!$  x = 1.0d0
!!$  y = -2.0d0
!!$
!!$  write(*,*) "x = ", x , " y = ", y
!!$
!!$  answer = minmod(x,y)
!!$
!!$  write(*,*) "answer = " , answer
!!$
!!$  x = -1.0d0
!!$  y = -2.0d0
!!$
!!$  write(*,*) "x = ", x , " y = ", y
!!$
!!$  answer = minmod(x,y)
!!$
!!$  write(*,*) "answer = " , answer
!!$
!!$end program main
