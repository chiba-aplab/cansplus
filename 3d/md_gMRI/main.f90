program main

  use mpi_setup
  use openfile
  use const
  use init
  use getNewdt
  use integrate, only : integrate__TVDRK3
  use omp_lib

  implicit none

  real(8) :: etime, etlim, etime0

  etlim = 23.*60.*60.+30.*60.
!  etlim = 2.*60.*60.+30.*60.
!  etlim = 5.*60.

!----------------------------------------------------------------------|
!  initialize
  call initialize
!----------------------------------------------------------------------|

  etime0 = omp_get_wtime()
  call MPI_BCAST(etime0,1,mdp,0,mcomw,merr)

!----------------------------------------------------------------------|
!  file open
  call file_output(nd,mpid%mpirank,ro,pr,vx,vy,vz,bx,by,bz,phi,eta &
                  ,ix,jx,kx)
  call file_output_param(nd,dtout,tend,ix,jx,kx,igx,jgx,kgx,margin &
       ,mpid%mpisize,mpid%mpirank,mpisize_x,mpisize_y,mpisize_z & 
       ,nrmlro,nrmlte,nrmlx,nrmlv,nrmlt,nrmlee,mass_bh,rg,rg_nrmlx &
       ,RadCool,te_factor,rohalo,eta0,vc,gm,x,y,z,dx,dy,dz,gx,gz)
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
                    ,ix,jx,kx)
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

  if(mpid%mpirank == 0) etime = omp_get_wtime()
  call MPI_BCAST(etime,1,mdp,0,mcomw,merr)
  if(etime-etime0 >= etlim) then
     if(mpid%mpirank == 0) write(*,*) '*** elapse time over ***',ns,etime-etime0
     exit loop
  endif

  mwflag=0
  mw=0

!----------------------------------------------------------------------|
!     obtain time spacing
!----------------------------------------------------------------------|
  call getNewdt__glmcyl(margin,safety,dtmin,ix,jx,kx,gm,ro,pr &
       ,vx,vy,vz,bx,by,bz,x,dx,y,dy,z,dz,eta &
       ,dt,ch,cp,min_dx)

  if(merr /= 0) exit loop

  timep = time
  time = time+dt

!---- integrate--------------------------------------------------------|
  call integrate__TVDRK3(margin,ix,jx,kx,gm,x,xm,dx,y,dy,z,dz,dt &
                        ,gx,gz,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cp &
                        ,roi,pri,vxi,vyi,vzi,bxi,byi,bzi,phii &
                        ,eta0,vc,eta,ccx,ccy,ccz,RadCool,te_factor,time &
                        ,rohalo,swtch_t,xin)
!----------------------------------------------------------------------|

!     data output 
  if(mpid%mpirank == 0)then
   if(mod(ns,100) == 0)then
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
                     ,ix,jx,kx)
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
