!
! $Log: 
!
program main
  use mpi_domain_xz
  use openfile
  use init
  use const
  use model

  implicit none
  include 'mpif.h'

!======================================================================|
!     prologue
!======================================================================|

!----------------------------------------------------------------------|
!  initialize counters

  nd=1
  time  = 0.0d0
  timep = 0.0d0

  ns    = 0
  merr  = 0
  nscount=0  

  dt = tend
  call initial()

!----------------------------------------------------------------------|
!  file open

  call openfileCor(nd,mpirank,ix,jx,kx)
  call openfileAll(nd,mpirank,ix,jx,kx)

!----------------------------------------------------------------------|
!     ready
!----------------------------------------------------------------------|

    call file_output(ro,pr,vx,vy,vz,bx,by,bz,phi,eta,time,ix,jx,kx)
    call file_output_param(dtout,tend,ix,jx,kx,igx,jgx,kgx,margin,mpisize &
       ,mpirank,mpid,nrmlro,nrmlte,nrmlx,nrmlv,nrmlt,nrmlee,mass_bh,rg,rg_nrmlx &
       ,RadCool,te_factor,rohalo,eta0,vc,gm,x,y,z,dx,dy,dz,gx,gz)


    if(mpirank .eq. 0)then
         write(6,913) ns,time,nd
         write(*,*) 'dt :: ',dt
    endif
   call closefileAll()
   call closefileCor()

!----------------------------------------------------------------------|
!  read-data


  ix0 = ix
  jx0 = jx
  kx0 = kx
  nx0 = nd

  if(nd .ne. 1)then
  open(91,file='readFileNumber.dat')
  read(91,*) nd
  close(91)
     call openReadFileAll(nd,mpirank,ix0,jx0,kx0,nx0)
     call file_input(ro,pr,vx,vy,vz,bx,by,bz,phi,eta,ix,jx,kx)

     time = real(nd-1)*dtout
     if(mpirank .eq. 0)then
        write(*,*) 'time :: ',time
     endif
     
     call closeReadFileAll()
 
  endif
  nd=nd+1



!======================================================================|
!     time integration 
!======================================================================|
1000  continue  
  ns = ns+1
  nscount = nscount+1
  mwflag=0

!----------------------------------------------------------------------|
!     obtain time spacing
  call getNewdt_glmcyl(margin,safety,dtmin,ix,jx,kx,gm,ro,pr &
       ,vx,vy,vz,bx,by,bz,x,dx,y,dy,z,dz,eta &
       ,dt,merr,ch)

  call mpi_allreduce(dt,dtg,1,mpi_double_precision,mpi_min &
       ,mpi_comm_world,merr)

  dt = dtg
  err_flg = 1 ! error flag
  if ( dt.lt.dtmin) goto 9999
  err_flg = 2 ! error flag
  if (merr.ne.0) goto 9999
  timep = time
  time = time+dt


!---- integrate

  call integrate_cyl(mpid,margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,dt &
       ,gx,gz,floor,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
       ,roi,pri,vxi,vyi,vzi,bxi,byi,bzi &
       ,eta0,vc,eta,ccx,ccy,ccz,RadCool,te_factor,time,rohalo,swtch_t,xin)
  
!----------------------------------------------------------------------|
!     data output 

  
  if(mpirank .eq. 0)then
!   if(nscount .ge. 100)then
     print '(" step=",i8," time=",e10.5," dt=",e10.5)', ns, time, dt
!     nscount=0
!   endif
  end if

  mw=0
!----- check output

! dtout
  nt1=int(timep/dtout)
  nt2=int(time/dtout)
  if (nt1.lt.nt2) mw=1
  if (mw.ne.0) then
     call openfileAll(nd,mpirank,ix,jx,kx)
     call file_output(ro,pr,vx,vy,vz,bx,by,bz,phi,eta,time,ix,jx,kx)
     call closefileAll()

     if(mpirank .eq. 0)then
        write(6,913) ns,time,nd
        write(*,*) '[NORMAL] dt :: ',dt
        write(*,*) 'nd: ',nd
     end if

    open(91,file='readFileNumber.dat')
    write(91,*) nd
    close(91)

     nd=nd+1
     mwflag=1
  endif
!----------------------------------------------------------------------|
!     loop test

  if (ns .lt. nstop .and. time .lt. tend ) goto 1000

!======================================================================|
!     epilogue
!======================================================================|
9999 continue
  if (merr .ne. 0) then
     write(6,*) '#### abnormal stop ####'
  endif
!----------------------------------------------------------------------|
!  data output
  if (mwflag.eq.0) then
     if(mpirank .eq. 0)then
        write(6,913) ns,time,nd
     end if

     call openfileAll(nd,mpirank,ix,jx,kx)
     call file_output(ro,pr,vx,vy,vz,bx,by,bz,phi,eta,time,ix,jx,kx)
     call closefileAll()
  endif

!======================================================================|
!        enddo
!======================================================================|

  close(mf_t)
  call mpi_finalize(merr)

  write(6,915) ns,time
  if (merr.eq.0) then
     write(6,*) '  ### normal stop ###'
  else
     write(6,*) '  ### abnormal stop ###'
  endif
  
  stop
  
913 format (1x,' write    ','step=',i8,' time=',e10.3,' nd =',i3)
915 format (1x,' stop     ','step=',i8,' time=',e10.3)


end program main
