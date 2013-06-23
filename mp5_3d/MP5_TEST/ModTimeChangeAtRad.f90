module ModTimeChangeAtRad
  implicit none
  
  integer,public :: mtr1_bx,mtr1_by,mtr1_bz
  integer,public :: mtr1_rx,mtr1_ry,mtr1_rz
  integer,public :: mtr1_ro,mtr1_pb,mtr1_pr

  integer,public :: mtr2_bx,mtr2_by,mtr2_bz
  integer,public :: mtr2_rx,mtr2_ry,mtr2_rz
  integer,public :: mtr2_ro,mtr2_pb,mtr2_pr

  integer,public :: mtr3_bx,mtr3_by,mtr3_bz
  integer,public :: mtr3_rx,mtr3_ry,mtr3_rz
  integer,public :: mtr3_ro,mtr3_pb,mtr3_pr

  integer,public :: mtr4_bx,mtr4_by,mtr4_bz
  integer,public :: mtr4_rx,mtr4_ry,mtr4_rz
  integer,public :: mtr4_ro,mtr4_pb,mtr4_pr

  integer,public :: mtr5_bx,mtr5_by,mtr5_bz
  integer,public :: mtr5_rx,mtr5_ry,mtr5_rz
  integer,public :: mtr5_ro,mtr5_pb,mtr5_pr

  character,private :: cno*4

contains
  subroutine calSumAtRadius(margin,ix,jx,kx,x,y,z,qq,xObs &
       , sumVal)
    implicit none

    integer,intent(in) :: ix,jx,kx,margin

    integer :: isInRadius

    real(8),dimension(ix),intent(in) :: x
    real(8),dimension(jx),intent(in) :: y
    real(8),dimension(kx),intent(in) :: z

    real(8),intent(in) :: xObs

    real(8),dimension(ix,jx,kx),intent(in) :: qq

    real(8),dimension(kx) :: sumVal
    integer :: iObs
    integer :: existFlag, overFlag
    
    integer :: i,j,k

    do k=1,kx
       sumVal = 0.0d0
    end do

    existFlag = 0
    overFlag = 0

    existFlag = isInRadius(margin,ix,x,xObs)

    if(existFlag == 1)then
       do i=1,ix
          if(overFlag == 0)then
             if(x(i) >= xObs)then
                iObs = i
                overFlag = 1
             end if
          end if
       end do

       do k=1,kx
          do j=1,jx
             sumVal(k) = sumVal(k)+qq(iObs,j,k)
          end do
       end do
    endif
  end subroutine calSumAtRadius

  subroutine calSumAtRadiusOfBase(margin,ix,jx,kx,x,y,z,xObs &
       , sumVal)
    implicit none

    integer,intent(in) :: ix,jx,kx,margin

    integer :: isInRadius

    real(8),dimension(ix),intent(in) :: x
    real(8),dimension(jx),intent(in) :: y
    real(8),dimension(kx),intent(in) :: z

    real(8),intent(in) :: xObs

    real(8),dimension(kx) :: sumVal
    integer :: iObs
    integer :: existFlag, overFlag
    
    integer :: i,j,k

    do k=1,kx
       sumVal = 0.0d0
    end do

    existFlag = 0
    overFlag = 0

    existFlag = isInRadius(margin,ix,x,xObs)

    if(existFlag == 1)then
       do i=1,ix
          if(overFlag == 0)then
             if(x(i) >= xObs)then
                iObs = i
                overFlag = 1
             end if
          end if
       end do

       do k=1,kx
          do j=1,jx
             sumVal(k) = sumVal(k)+1.0d0
          end do
       end do
    endif
  end subroutine calSumAtRadiusOfBase
  
  subroutine openfileModTCRsOld(mpirank,margin,ix,x,xObs1,xObs2,xObs3,xObs4,xObs5)
    implicit none

    integer,intent(in) :: mpirank
    integer,intent(in) :: ix,margin

    real(8),intent(in) :: xObs1,xObs2,xObs3,xObs4,xObs5
    real(8),dimension(ix),intent(in) :: x

    integer :: flag
integer :: isInRadius

    flag = 0
    flag = isInRadius(margin,ix,x,xObs1)

    if(flag == 1)then
       write(cno,'(i4.4)') mpirank
       open(mtr1_bx,file='mtr1_bx.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr1_by,file='mtr1_by.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr1_bz,file='mtr1_bz.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr1_rx,file='mtr1_rx.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr1_ry,file='mtr1_ry.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr1_rz,file='mtr1_rz.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr1_ro,file='mtr1_ro.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr1_pb,file='mtr1_pb.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr1_pr,file='mtr1_pr.dac.'//cno,form='unformatted',status='old',position='append')
    endif
    flag = 0
    flag = isInRadius(margin,ix,x,xObs2)
    if(flag == 1)then
       write(cno,'(i4.4)') mpirank
       open(mtr2_bx,file='mtr2_bx.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr2_by,file='mtr2_by.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr2_bz,file='mtr2_bz.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr2_rx,file='mtr2_rx.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr2_ry,file='mtr2_ry.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr2_rz,file='mtr2_rz.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr2_ro,file='mtr2_ro.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr2_pb,file='mtr2_pb.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr2_pr,file='mtr2_pr.dac.'//cno,form='unformatted',status='old',position='append')
    endif

    flag = 0
    flag = isInRadius(margin,ix,x,xObs3)
    if(flag == 1)then
       write(cno,'(i4.4)') mpirank
       open(mtr3_bx,file='mtr3_bx.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr3_by,file='mtr3_by.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr3_bz,file='mtr3_bz.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr3_rx,file='mtr3_rx.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr3_ry,file='mtr3_ry.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr3_rz,file='mtr3_rz.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr3_ro,file='mtr3_ro.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr3_pb,file='mtr3_pb.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr3_pr,file='mtr3_pr.dac.'//cno,form='unformatted',status='old',position='append')
    endif

    flag = 0
    flag = isInRadius(margin,ix,x,xObs4)
    if(flag == 1)then
       write(cno,'(i4.4)') mpirank
       open(mtr4_bx,file='mtr4_bx.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr4_by,file='mtr4_by.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr4_bz,file='mtr4_bz.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr4_rx,file='mtr4_rx.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr4_ry,file='mtr4_ry.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr4_rz,file='mtr4_rz.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr4_ro,file='mtr4_ro.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr4_pb,file='mtr4_pb.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr4_pr,file='mtr4_pr.dac.'//cno,form='unformatted',status='old',position='append')
    endif

    flag = 0
    flag = isInRadius(margin,ix,x,xObs5)
    if(flag == 1)then
       write(cno,'(i4.4)') mpirank
       open(mtr5_bx,file='mtr5_bx.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr5_by,file='mtr5_by.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr5_bz,file='mtr5_bz.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr5_rx,file='mtr5_rx.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr5_ry,file='mtr5_ry.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr5_rz,file='mtr5_rz.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr5_ro,file='mtr5_ro.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr5_pb,file='mtr5_pb.dac.'//cno,form='unformatted',status='old',position='append')
       open(mtr5_pr,file='mtr5_pr.dac.'//cno,form='unformatted',status='old',position='append')
    endif

    return
  end subroutine openfileModTCRsOld

  subroutine openfileModTCR(mpirank,margin,ix,kx,x,xObs1,xObs2,xObs3,xObs4,xObs5)
    implicit none

    integer :: isInRadius

    integer,intent(in) :: mpirank

    integer,intent(in) :: ix,kx,margin
    real(8),intent(in) :: xObs1,xObs2,xObs3,xObs4,xObs5
    real(8),dimension(ix),intent(in) :: x

    integer :: flag

    flag = 0
    flag = isInRadius(margin,ix,x,xObs1)
    if(flag == 1)then
       write(cno,'(i4.4)') mpirank
       mtr1_bx=501
       call dacdef1s(mtr1_bx,'mtr1_bx.dac.'//cno,6,kx)
       mtr1_by=502
       call dacdef1s(mtr1_by,'mtr1_by.dac.'//cno,6,kx)
       mtr1_bz=503
       call dacdef1s(mtr1_bz,'mtr1_bz.dac.'//cno,6,kx)

       mtr1_rx=504
       call dacdef1s(mtr1_rx,'mtr1_rx.dac.'//cno,6,kx)
       mtr1_ry=505
       call dacdef1s(mtr1_ry,'mtr1_ry.dac.'//cno,6,kx)
       mtr1_rz=506
       call dacdef1s(mtr1_rz,'mtr1_rz.dac.'//cno,6,kx)

       mtr1_ro=507
       call dacdef1s(mtr1_ro,'mtr1_ro.dac.'//cno,6,kx)
       mtr1_pb=508
       call dacdef1s(mtr1_pb,'mtr1_pb.dac.'//cno,6,kx)
       mtr1_pr=509
       call dacdef1s(mtr1_pr,'mtr1_pr.dac.'//cno,6,kx)
    endif

    flag = 0
    flag = isInRadius(margin,ix,x,xObs2)
    if(flag == 1)then
       write(cno,'(i4.4)') mpirank
       mtr2_bx=511
       call dacdef1s(mtr2_bx,'mtr2_bx.dac.'//cno,6,kx)
       mtr2_by=512
       call dacdef1s(mtr2_by,'mtr2_by.dac.'//cno,6,kx)
       mtr2_bz=513
       call dacdef1s(mtr2_bz,'mtr2_bz.dac.'//cno,6,kx)

       mtr2_rx=514
       call dacdef1s(mtr2_rx,'mtr2_rx.dac.'//cno,6,kx)
       mtr2_ry=515
       call dacdef1s(mtr2_ry,'mtr2_ry.dac.'//cno,6,kx)
       mtr2_rz=516
       call dacdef1s(mtr2_rz,'mtr2_rz.dac.'//cno,6,kx)

       mtr2_ro=517
       call dacdef1s(mtr2_ro,'mtr2_ro.dac.'//cno,6,kx)
       mtr2_pb=518
       call dacdef1s(mtr2_pb,'mtr2_pb.dac.'//cno,6,kx)
       mtr2_pr=519
       call dacdef1s(mtr2_pr,'mtr2_pr.dac.'//cno,6,kx)
    endif

    flag = 0
    flag = isInRadius(margin,ix,x,xObs3)
    if(flag == 1)then
       write(cno,'(i4.4)') mpirank
       mtr3_bx=521
       call dacdef1s(mtr3_bx,'mtr3_bx.dac.'//cno,6,kx)
       mtr3_by=522
       call dacdef1s(mtr3_by,'mtr3_by.dac.'//cno,6,kx)
       mtr3_bz=523
       call dacdef1s(mtr3_bz,'mtr3_bz.dac.'//cno,6,kx)

       mtr3_rx=524
       call dacdef1s(mtr3_rx,'mtr3_rx.dac.'//cno,6,kx)
       mtr3_ry=525
       call dacdef1s(mtr3_ry,'mtr3_ry.dac.'//cno,6,kx)
       mtr3_rz=526
       call dacdef1s(mtr3_rz,'mtr3_rz.dac.'//cno,6,kx)

       mtr3_ro=527
       call dacdef1s(mtr3_ro,'mtr3_ro.dac.'//cno,6,kx)
       mtr3_pb=528
       call dacdef1s(mtr3_pb,'mtr3_pb.dac.'//cno,6,kx)
       mtr3_pr=529
       call dacdef1s(mtr3_pr,'mtr3_pr.dac.'//cno,6,kx)
    endif

    flag = 0
    flag = isInRadius(margin,ix,x,xObs4)
    if(flag == 1)then
       write(cno,'(i4.4)') mpirank
       mtr4_bx=531
       call dacdef1s(mtr4_bx,'mtr4_bx.dac.'//cno,6,kx)
       mtr4_by=532
       call dacdef1s(mtr4_by,'mtr4_by.dac.'//cno,6,kx)
       mtr4_bz=533
       call dacdef1s(mtr4_bz,'mtr4_bz.dac.'//cno,6,kx)

       mtr4_rx=534
       call dacdef1s(mtr4_rx,'mtr4_rx.dac.'//cno,6,kx)
       mtr4_ry=535
       call dacdef1s(mtr4_ry,'mtr4_ry.dac.'//cno,6,kx)
       mtr4_rz=536
       call dacdef1s(mtr4_rz,'mtr4_rz.dac.'//cno,6,kx)

       mtr4_ro=537
       call dacdef1s(mtr4_ro,'mtr4_ro.dac.'//cno,6,kx)
       mtr4_pb=538
       call dacdef1s(mtr4_pb,'mtr4_pb.dac.'//cno,6,kx)
       mtr4_pr=539
       call dacdef1s(mtr4_pr,'mtr4_pr.dac.'//cno,6,kx)
    endif

    flag = 0
    flag = isInRadius(margin,ix,x,xObs5)
    if(flag == 1)then
       write(cno,'(i4.4)') mpirank
       mtr5_bx=541
       call dacdef1s(mtr5_bx,'mtr5_bx.dac.'//cno,6,kx)
       mtr5_by=542
       call dacdef1s(mtr5_by,'mtr5_by.dac.'//cno,6,kx)
       mtr5_bz=543
       call dacdef1s(mtr5_bz,'mtr5_bz.dac.'//cno,6,kx)

       mtr5_rx=544
       call dacdef1s(mtr5_rx,'mtr5_rx.dac.'//cno,6,kx)
       mtr5_ry=545
       call dacdef1s(mtr5_ry,'mtr5_ry.dac.'//cno,6,kx)
       mtr5_rz=546
       call dacdef1s(mtr5_rz,'mtr5_rz.dac.'//cno,6,kx)

       mtr5_ro=547
       call dacdef1s(mtr5_ro,'mtr5_ro.dac.'//cno,6,kx)
       mtr5_pb=548
       call dacdef1s(mtr5_pb,'mtr5_pb.dac.'//cno,6,kx)
       mtr5_pr=549
       call dacdef1s(mtr5_pr,'mtr5_pr.dac.'//cno,6,kx)
    endif

    return
  end subroutine openfileModTCR

  subroutine closefileModTCR()
    implicit none

    close(mtr1_bx)
    close(mtr1_by)
    close(mtr1_bz)
    close(mtr1_rx)
    close(mtr1_ry)
    close(mtr1_rz)
    close(mtr1_ro)
    close(mtr1_pb)
    close(mtr1_pr)

    close(mtr2_bx)
    close(mtr2_by)
    close(mtr2_bz)
    close(mtr2_rx)
    close(mtr2_ry)
    close(mtr2_rz)
    close(mtr2_ro)
    close(mtr2_pb)
    close(mtr2_pr)

    close(mtr3_bx)
    close(mtr3_by)
    close(mtr3_bz)
    close(mtr3_rx)
    close(mtr3_ry)
    close(mtr3_rz)
    close(mtr3_ro)
    close(mtr3_pb)
    close(mtr3_pr)

    close(mtr4_bx)
    close(mtr4_by)
    close(mtr4_bz)
    close(mtr4_rx)
    close(mtr4_ry)
    close(mtr4_rz)
    close(mtr4_ro)
    close(mtr4_pb)
    close(mtr4_pr)

    close(mtr5_bx)
    close(mtr5_by)
    close(mtr5_bz)
    close(mtr5_rx)
    close(mtr5_ry)
    close(mtr5_rz)
    close(mtr5_ro)
    close(mtr5_pb)
    close(mtr5_pr)

    return 
  end subroutine closefileModTCR

!======================================================================
! if x(margin) < xObs < x(ix-margin) => return 1
! else return 0
!======================================================================
  function isInRadius(margin,ix,x,xObs)
    implicit none
    
    integer,intent(in) :: margin,ix
    real(8),dimension(ix),intent(in) :: x

    real(8),intent(in) :: xObs

    integer :: isInRadius

    integer :: existFlag
    integer :: i

    existFlag = 0
    if(x(margin+1) <= xObs)then
       if(x(ix-margin) >= xObs)then
          existFlag = 1
       end if
    end if

    isInRadius = existFlag
  end function isInRadius

end module ModTimeChangeAtRad
