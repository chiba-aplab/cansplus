!======================================================================
! NOTE
!  please change 917 format (1x,i8,ES14.6,1(" ",ES14.6)), 
!                                         ^ this number is iCalNum
!  if you change iCalNum
!======================================================================
module ModAverageValue
  implicit none

contains
  subroutine calSum(margin,ix,jx,kx,x,y,z,xmin,xmax,zmin,zmax &
       ,qq,sumVal)
    implicit none

    integer,intent(in) :: ix,jx,kx,margin
    
    real(8),dimension(ix),intent(in) :: x
    real(8),dimension(jx),intent(in) :: y
    real(8),dimension(kx),intent(in) :: z

    real(8),intent(in) :: xmin,xmax
    real(8),intent(in) :: zmin,zmax
    
    real(8),dimension(ix,jx,kx),intent(in) :: qq

    real(8) :: sumVal

    integer :: i,j,k

    sumVal = 0.0d0
       
    do k=margin+1,kx-margin
       do j=margin+1,jx-margin
          do i=margin+1,ix-margin
             if(x(i) > xmax)then
                cycle
             endif
             if(x(i) <= xmin)then
                cycle
             endif
             if(z(k) > zmax)then
                cycle
             endif
             if(z(k) <= zmin)then
                cycle
             endif
             sumVal = sumVal + qq(i,j,k)
          end do
       end do
    end do

    return
  end subroutine calSum

  subroutine calSumBase(margin,ix,jx,kx,x,y,z,xmin,xmax,zmin,zmax &
       ,sumVal)
    implicit none

    integer,intent(in) :: ix,jx,kx,margin
    
    real(8),dimension(ix),intent(in) :: x
    real(8),dimension(jx),intent(in) :: y
    real(8),dimension(kx),intent(in) :: z

    real(8),intent(in) :: xmin,xmax
    real(8),intent(in) :: zmin,zmax
    
    real(8) :: sumVal

    integer :: i,j,k

    sumVal = 0.0d0
       
    do k=margin+1,kx-margin
       do j=margin+1,jx-margin
          do i=margin+1,ix-margin
             if(x(i) > xmax)then
                cycle
             endif
             if(x(i) <= xmin)then
                cycle
             endif
             if(z(k) > zmax)then
                cycle
             endif
             if(z(k) <= zmin)then
                cycle
             endif

             sumVal = sumVal + 1.0d0
          end do
       end do
    end do

    return
  end subroutine calSumBase
  
  subroutine calAverageVal(sumVal, baseVal, averageVal)
    implicit none

    real(8),intent(in) :: sumVal, baseVal
    real(8) :: averageVal 

    if(baseVal == 0.0d0)then
       averageVal = 0.0d0
    else
       averageVal = sumVal/averageVal
    endif
    return
  end subroutine calAverageVal

  subroutine exchangeSumVal(mpirank,iCalNum,calAverageValues)
    implicit none

    include 'mpif.h'

    integer,intent(in) :: iCalNum,mpirank
    real(8),dimension(iCalNum) :: calAverageValues,send,recv

    integer :: merr

    integer :: i

    do i=1,iCalNum
       recv(i) = calAverageValues(i)
       send(i) = calAverageValues(i)
    end do

    call mpi_reduce &
         (send,recv,iCalNum,mpi_double_precision &
         ,mpi_sum,0,mpi_comm_world,merr)

    do i=1,iCalNum
       calAverageValues(i) = recv(i)
    end do

    return
  end subroutine exchangeSumVal

  subroutine outputCalPhysValue(mpirank,iCalNum,calPhysValues,ns,time &
       ,fileName,mf_calPhysValueTest)
    implicit none

    integer,intent(in) :: iCalNum,ns,mpirank
    real(8),dimension(iCalNum) :: calPhysValues

    real(8),intent(in) :: time
    character(10),intent(in) :: fileName
    integer,intent(in) :: mf_calPhysValueTest

    integer :: merr

    if(mpirank .eq. 0)then
       open(mf_calPhysValueTest,file=fileName,status='old',position='append')
       write(mf_calPhysValueTest,917) ns,time,calPhysValues
       close(mf_calPhysValueTest)
    endif

917 format (1x,i8,ES14.6,5(" ",ES14.6))
  end subroutine outputCalPhysValue

  subroutine open_CalPhysValue(mpirank,fileName,mf_calPhysValueTest)
    implicit none

    integer,intent(in) :: mpirank
    integer,intent(in) :: mf_calPhysValueTest
    
    character(10),intent(in) :: fileName

    if(mpirank .eq. 0)then
       open(mf_calPhysValueTest,file=fileName)
       write(mf_calPhysValueTest,*) '1:step,2:time,3:pb,4:pr,5:ro,6:ee,7:am'
       close(mf_calPhysValueTest)
    endif
    return
  end subroutine open_CalPhysValue
end module ModAverageValue
