function median(x,y,z)

  implicit none

  real(8),intent(in) :: x,y,z

  real(8) :: minmod
  real(8) :: median

  median = x + minmod((y-x),(z-x))

end function median
!!$
!!$program main 
!!$
!!$  implicit none
!!$
!!$  real(8) :: x,y,z
!!$
!!$  real(8) :: median
!!$  real(8) :: answer
!!$
!!$  write(*,*) " all positive"
!!$  x = 1.0d0
!!$  y = 2.0d0
!!$  z = 3.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$
!!$  x = 2.0d0
!!$  y = 3.0d0
!!$  z = 1.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$  
!!$  x = 3.0d0
!!$  y = 1.0d0
!!$  z = 2.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$
!!$  x = 1.0d0
!!$  y = 1.0d0
!!$  z = 2.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$
!!$  x = 1.0d0
!!$  y = 1.0d0
!!$  z = 1.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$
!!$  write(*,*) " all negative"
!!$
!!$  x = -1.0d0
!!$  y = -2.0d0
!!$  z = -3.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$
!!$  x = -2.0d0
!!$  y = -3.0d0
!!$  z = -1.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$  
!!$  x = -3.0d0
!!$  y = -1.0d0
!!$  z = -2.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$
!!$  x = -1.0d0
!!$  y = -1.0d0
!!$  z = -2.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$
!!$  x = -1.0d0
!!$  y = -1.0d0
!!$  z = -1.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$
!!$  write(*,*) " mix"
!!$  x = 1.0d0
!!$  y = -2.0d0
!!$  z = -3.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$
!!$  x = -2.0d0
!!$  y = -3.0d0
!!$  z = 1.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$  
!!$  x = -3.0d0
!!$  y = 1.0d0
!!$  z = -2.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$
!!$  x = 1.0d0
!!$  y = 2.0d0
!!$  z = -3.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$
!!$  x = 2.0d0
!!$  y = -3.0d0
!!$  z = 1.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$  
!!$  x = -3.0d0
!!$  y = 1.0d0
!!$  z = 2.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$
!!$  x = 1.0d0
!!$  y = -2.0d0
!!$  z = 3.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$
!!$  x = 2.0d0
!!$  y = -3.0d0
!!$  z = 1.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$  
!!$  x = 3.0d0
!!$  y = 1.0d0
!!$  z = -2.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$
!!$  x = -1.0d0
!!$  y = 2.0d0
!!$  z = 3.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$
!!$  x = -2.0d0
!!$  y = 3.0d0
!!$  z = 1.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$  
!!$  x = 3.0d0
!!$  y = -1.0d0
!!$  z = 2.0d0
!!$  write(*,*) "x = ", x , " y = ", y, " z = ", z
!!$  answer = median(x,y,z)
!!$  write(*,*) "answer = " , answer
!!$
!!$end program main
