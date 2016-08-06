program main

  use mpi_setup
  use openfile
  use const
  use init
  use getNewdt
  use integrate, only : integrate__TVDRK3

  implicit none

!----------------------------------------------------------------------|
!  initialize
  call initialize
!----------------------------------------------------------------------|

!----------------------------------------------------------------------|
!  file open
   call file_output_param(nd,dtout,tend,ix,jx,igx,jgx,margin &
         ,mpid%mpisize,mpid%mpirank,mpisize_x,mpisize_y &
         ,gm,x,y,dx,dy)
   call file_output(nd,mpid%mpirank,ro,pr,vx,vy,vz,bx,by,bz,phi,eta &
         ,ix,jx)

  if(mpid%mpirank == 0)then
    write(6,913) ns,time,nd
    write(*,*) 'dt :: ',dt
  endif
!----------------------------------------------------------------------|

!----------------------------------------------------------------------|
!  read-data
  if(restart)then
     open(91,file=input_dir//'readFileNumber.dat')
     read(91,*) nd
     close(91)
     call file_input(nd,mpid%mpirank,ro,pr,vx,vy,vz,bx,by,bz,phi,eta &
                    ,ix,jx)
     time = real(nd-1)*dtout
     if(mpid%mpirank == 0)then
        write(*,*) 'restart -> time :: ',time,' nd ::',nd
     endif
  endif
  nd=nd+1
!----------------------------------------------------------------------|

!======================================================================|
!     time integration 
!======================================================================|
  loop: do ns=1,nstop

  mwflag=0
  mw=0

!----------------------------------------------------------------------|
!     obtain time spacing
!----------------------------------------------------------------------|
  call getNewdt__glm(margin,safety,dtmin,ix,jx,gm,ro,pr &
       ,vx,vy,vz,bx,by,bz,x,dx,y,dy,eta &
       ,dt,ch,cp,min_dx)

  if(merr /= 0) exit loop

  timep = time
  time = time+dt

!---- integrate--------------------------------------------------------|
  call integrate__TVDRK3(margin,ix,jx,gm,x,dx,y,dy,dt &
                        ,gy,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cp &
                        ,roi,pri,vxi,vyi,vzi,bxi,byi,bzi,phii &
                        ,eta0,vc,eta,ccx,ccy,time )
!----------------------------------------------------------------------|

!     data output 
  if(mpid%mpirank == 0)then
   if(mod(ns,1000) == 0)then
     print '(" step=",i8," time=",e10.5," dt=",e10.5)', ns, time, dt
   endif
  endif


!----- check output
! dtout
  nt1=int(timep/dtout)
  nt2=int(time/dtout)
  if(nt1 < nt2) mw=1
  if(mw /= 0)then
     call file_output(nd,mpid%mpirank,ro,pr,vx,vy,vz,bx,by,bz,phi,eta &
                     ,ix,jx)
     if(mpid%mpirank == 0)then
        write(6,913) ns,time,nd
        write(*,*) '[NORMAL] dt :: ',dt
        write(*,*) 'nd: ',nd
     endif

    open(91,file=output_dir//'readFileNumber.dat')
    write(91,*) nd
    close(91)

     nd=nd+1
     mwflag=1
  endif

! loop test
  if(time > tend ) exit loop

  enddo loop
!======================================================================|
!  end of time integration 
!======================================================================|

!  data output
  if(mwflag == 0)then
     if(mpid%mpirank == 0)then
        write(6,913) ns,time,nd
     endif
     call file_output(nd,mpid%mpirank,ro,pr,vx,vy,vz,bx,by,bz,phi,eta &
                     ,ix,jx)
  endif

  call mpi_finalize(merr)

  write(6,915) ns,time
  if(merr == 0)then
     write(6,*) '  ### normal stop ###'
  else
     write(6,*) '  ### abnormal stop ###'
  endif
  
913 format (1x,' write    ','step=',i8,' time=',e10.3,' nd =',i3)
915 format (1x,' stop     ','step=',i8,' time=',e10.3)


end program main
