subroutine reconstructionConstant(margin,ix,x,xm,dx,cc)
  implicit none

  integer,intent(in) :: ix,margin

  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(0:ix),intent(in) :: xm

  real(8),dimension(5,2,ix) :: cc

  real(8) :: tmp1,tmp2,tmp3,tmp4
  integer :: i,m,l,q
  integer,parameter :: k = 5
  integer :: j,r

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

  return 
end subroutine reconstructionConstant

!!$program main
!!$  implicit none
!!$
!!$  integer,parameter :: ix = 50,margin = 5
!!$  
!!$  real(8),dimension(ix) :: x,dx,xm
!!$  real(8),dimension(ix,5,2) :: cc
!!$
!!$  real(8) :: xmin,xmax,dx0
!!$
!!$  integer :: izero
!!$
!!$  integer :: i,j
!!$
!!$  xmax = 0.5d0
!!$  xmin = -0.5d0
!!$
!!$  !
!!$  ! Set X-coordinate
!!$  !
!!$  dx0 = (xmax-xmin)/real(ix-(2*margin+1))
!!$
!!$  do i=1,ix
!!$     dx(i) = dx0
!!$  enddo
!!$
!!$  izero = margin+1
!!$
!!$  xm(izero) = xmin + 0.5d0*dx(izero)
!!$  do i=izero,ix-1
!!$     xm(i+1) = xm(i) + dx(i+1)
!!$  enddo
!!$
!!$  do i=izero-1,1,-1
!!$     xm(i) = xm(i+1)-dx(i+1)
!!$  enddo
!!$
!!$  do i=2,ix
!!$     x(i) = 0.5d0*(xm(i)+xm(i-1))
!!$  enddo
!!$
!!$  x(1) = 0.5d0*(xm(1) + (xm(1)-dx0))
!!$
!!$  call reconstructionConstant(margin,ix,x,xm,dx,cc)
!!$
!!$  write(*,*) '##############left'
!!$  do i=1,ix
!!$     write(*,*) cc(i,1,1), cc(i,2,1), cc(i,3,1), cc(i,4,1), cc(i,5,1)
!!$  end do
!!$
!!$  write(*,*) '##############right'
!!$  do i=1,ix
!!$     write(*,*) cc(i,1,2), cc(i,2,2), cc(i,3,2), cc(i,4,2), cc(i,5,2)
!!$  end do
!!$
!!$end program main
