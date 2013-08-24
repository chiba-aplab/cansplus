!
! $Log: 
!
program main
  use mpi_domain_xz
  use openfile

  implicit none
  include 'mpif.h'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer,parameter :: FLAG_Accuracy=0 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
   integer,parameter :: margin = 3
!   integer,parameter :: margin = 6

  integer,parameter :: ix = 8+2*margin,jx=64+2*margin,kx=8+2*margin
!  integer,parameter :: ix = 256+2*margin,jx=64+2*margin,kx=4+2*margin

  integer :: igx,jgx,kgx
  integer :: ix0,jx0,kx0,nx0

  real(8),dimension(ix) :: x,dx
  real(8),dimension(jx) :: y,dy
  real(8),dimension(kx) :: z,dz
  real(8),dimension(5,2,ix) :: ccx
  real(8),dimension(5,2,jx) :: ccy
  real(8),dimension(5,2,kx) :: ccz

! cell center variables
  real(8),dimension(ix,jx,kx) :: ro,pr
  real(8),dimension(ix,jx,kx) :: vx,vy,vz
  real(8),dimension(ix,jx,kx) :: bx,by,bz,phi
  real(8),dimension(ix,jx,kx) :: eta

  real(8),dimension(ix,jx,kx) :: ro1,pr1
  real(8),dimension(ix,jx,kx) :: vx1,vy1,vz1
  real(8),dimension(ix,jx,kx) :: bx1,by1,bz1,phi1
  
! gravitation
  real(8),dimension(ix,jx,kx) :: gx,gy,gz

  real(8),dimension(ix,jx,kx) :: roi,pri
  real(8),dimension(ix,jx,kx) :: vxi,vyi,vzi
  real(8),dimension(ix,jx,kx) :: bxi,byi,bzi

  real(8) :: ch,cr

  real(8) :: gm ! specific heat retio
  real(8) :: eta0,vc

  real(8) :: floor ! minimum value

! MPI
  integer :: mpisize
  integer :: mpirank
  
! MPI type
  type(mpidomain) :: mpid

!--cooling
  real(8),parameter :: velocity_c = 2.9979d0*(10.0d0**10)                                       !cgs
  real(8),parameter :: mass_bh = 10.0d0                                                         !per mass_solar
  real(8),parameter :: mass_proton = 1.6726d0*(10.0d0**(-24))                                   !cgs
  real(8),parameter :: boltzmann_const = 1.3807d0*(10.0d0**(-16))                               !cgs
  real(8),parameter :: Navo = 6.0221d0*(10.0d0**23)                                             !Avogadro's constant
  real(8),parameter :: mmw = 1.0d0                                                              !mean molecular weight
  real(8),parameter :: rg = 3.0d0*(10.0d0**6)*(mass_bh/10.0d0)                                  !cgs Schwarzschild radius
  real(8),parameter :: rg_nrmlx = 1.0d0/10.0d0                                                  !per nrmlx
!  real(8),parameter :: rg_nrmlx = 1.0d0/1000.0d0                                                  !per nrmlx


  real(8),parameter :: nrmlro = 0.29d0*8.3d0*(10.0d0**(-7))/(mass_bh/10.0d0)            !normalized ro
  real(8),parameter :: nrmlte = mass_proton*(velocity_c**2)/boltzmann_const             !normalized te
  real(8),parameter :: nrmlx = rg/rg_nrmlx                                              !normalized x
  real(8),parameter :: nrmlv = sqrt(0.5d0*rg_nrmlx/((1.0d0-rg_nrmlx)**2))*velocity_c    !normalized v
  real(8),parameter :: nrmlt = nrmlx/nrmlv                                              !normalized t
  real(8),parameter :: nrmlee = nrmlro*nrmlv*nrmlv                                      !normalized ee

  real(8),parameter :: RadCool = nrmlro*sqrt(nrmlte)*nrmlx/(nrmlv**3)*(6.2d0*(10.0d0**20))   !cooling factor
  real(8),parameter :: te_factor = mmw*(nrmlv**2)/Navo/mass_proton/(velocity_c**2)           !factor for te

  real(8),parameter :: rohalo = 1.0d-4
!  real(8),parameter :: rohalo = 1.0d-6

!======================================================================|
!     prologue
!======================================================================|

!----------------------------------------------------------------------|
!  initialize counters

  integer :: nd,ns,merr,ns1,ns2,nscount
  real(8) :: time,timep,swtch_t

!----------------------------------------------------------------------|
!   time control parameters
!     nstop : number of total time steps for the run

  real(8) :: tend,dtout 
  integer :: nstop

  integer :: mwflag,mw,nt1,nt2 

  real(8) :: dt,hdt ! time step
  real(8) :: dtg

  real(8) :: safety ! CFL number
  real(8) :: dtmin ! minimum time step

  integer :: i,j,k,n

  integer :: mf_qq, err_flg   

  real(8) :: xin
  real(8) :: te_limit
 
  real(8) :: temp
!======================================================================|
!     prologue
!======================================================================|
  floor = 1.d-6

  nd = 1
!----------------------------------------------------------------------|
!   for MPI
!  
  call mpi_init(merr)
  call mpi_comm_size(mpi_comm_world,mpisize,merr)
  call mpi_comm_rank(mpi_comm_world,mpirank,merr)

  mpid%mpirank = mpirank
  mpid%mpisize = mpisize

! Core Numbers, 1: r(i)-direction, 2: z(k)-direction
  mpid%mpisize_2d(1) = 32
  mpid%mpisize_2d(2) = 32

  igx = ix*mpid%mpisize_2d(1)-2*margin*(mpid%mpisize_2d(1)-1)
  jgx = jx
  kgx = kx*mpid%mpisize_2d(2)-2*margin*(mpid%mpisize_2d(2)-1)

! determine mpid%mpirank_3d
  call setmy2drank(mpid,merr)
  call setmpiboundary(mpid)

  if (mpirank.eq.0) then
     write(6,*) 'igx=',igx,'kgx=',kgx
     write(*,*) 'mpirank',mpirank,'mpirank_2d',mpid%mpirank_2d
  endif

!----------------------------------------------------------------------|
!  initialize counters

  nd=1
  time  = 0.0d0
  timep = 0.0d0
  swtch_t=62.8d0
  ns    = 0
  merr  = 0
  nscount=0  
      
!----------------------------------------------------------------------|
!   time control parameters
!     nstop : number of total time steps for the run

  tend = 6.28d0*10.0d0
  dtout = tend/(100.0d0)
  tend = 6.28d0*40.0d0

  dt = tend

  nstop = 2000000

!----------------------------------------------------------------------|
!   setup numerical model (grid, initial conditions, etc.)
!

  call model_machida(mpid,igx,jgx,kgx,ix,jx,kx,margin &
      ,ro,pr,vx,vy,vz,bx,by,bz,phi,gm &
      ,roi,pri,vxi,vyi,vzi,bxi,byi,bzi &
      ,x,dx,y,dy,z,dz &
      ,gx,gz,mf_params,dtout,tend,cr,eta0,vc,eta,xin &
      ,ccx,ccy,ccz,rg_nrmlx,rohalo,te_factor,nrmlv &
      ,nrmlte,boltzmann_const,Navo,mmw) 

  call exchangeMpixz(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz &
       ,phi,merr)

  call bnd(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,x,z &
             ,xin,roi,pri,vxi,vyi,vzi,bxi,byi,bzi)

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
   close(mf_params)

!----------------------------------------------------------------------|
!  read-data


  ix0 = ix
  jx0 = jx
  kx0 = kx
  nx0 = nd

  if(nd .ne. 1)then
  open(100,file='readFileNumber.dat')
  read(100,*) nd
  close(100)
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

  safety=0.2d0
  dtmin=1.d-10
  
  temp=0.0d0
  call getNewdt_glmcyl_res_mpi(margin,safety,dtmin,ix,jx,kx,gm,ro,pr &
       ,vx,vy,vz,bx,by,bz,x,dx,y,dy,z,dz,eta &
       ,dt,merr,ch,mpirank,temp)

  call mpi_allreduce(dt,dtg,1,mpi_double_precision,mpi_min &
       ,mpi_comm_world,merr)

  dt = dtg
  err_flg = 1 ! error flag
  if ( dt.lt.dtmin) goto 9999
  err_flg = 2 ! error flag
  if (merr.ne.0) goto 9999
  timep = time
  time = time+dt
  hdt = 0.5d0*dt

  mw=0
!----- check output

! dtout
  nt1=int(timep/dtout)
  nt2=int(time/dtout)
  if (nt1.lt.nt2) mw=1

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

  if (mw.ne.0) then
     call openfileAll(nd,mpirank,ix,jx,kx)
     call file_output(ro,pr,vx,vy,vz,bx,by,bz,phi,eta,time,ix,jx,kx)
     call closefileAll()

     if(mpirank .eq. 0)then
        write(6,913) ns,time,nd
        write(*,*) '[NORMAL] dt :: ',dt
        write(*,*) 'nd: ',nd
     end if

    open(100,file='readFileNumber.dat')
    write(100,*) nd
    close(100)

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
