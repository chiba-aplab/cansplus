!
! $Log: 
!
program main
  use mpi_domain
  use openfile

  implicit none
  include 'mpif.h'

  integer,parameter :: margin = 3

  integer,parameter :: ix = 200+2*margin,jx=200+2*margin,kx=2+2*margin

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

! test1028
  integer,parameter :: iCalNum = 5
  integer,parameter :: iInterNum = 6
  integer,parameter :: mf_calPhysValue = 401

  real(8),dimension(ix,jx,kx) :: pb,ee,am,rx,ry,rz
  real(8),dimension(iCalNum) :: avVal
  real(8),dimensioN(iInterNum) :: avInterVal

  
!======================================================================|
!     prologue
!======================================================================|

  integer :: mcont,ndi

!----------------------------------------------------------------------|
!  initialize counters

  integer :: nd,ns,merr,ns1,ns2
  real(8) :: time,timep

!----------------------------------------------------------------------|
!   time control parameters
!     nstop : number of total time steps for the run

  real(8) :: tend,dtout,dtout1
  integer :: nstop

  integer :: mwflag,mw,nt1,nt2,mw1
  integer :: nt3,nt4,nd1

  real(8) :: dt,hdt ! time step
  real(8) :: dtg

  real(8) :: safety ! CFL number
  real(8) :: dtmin ! minimum time step

  integer :: idf

  integer :: i,j,k,n

  integer :: mf_qq, err_flg,nanflg

  real(8) :: xin

!======================================================================|
!     prologue
!======================================================================|
  floor = 1.d-10

  mcont=0
  ndi=1000

  nd = 1
  nd1 = 1

!----------------------------------------------------------------------|
!   for MPI
!  
  call mpi_init(merr)
  call mpi_comm_size(mpi_comm_world,mpisize,merr)
  call mpi_comm_rank(mpi_comm_world,mpirank,merr)

  mpid%mpirank = mpirank
  mpid%mpisize = mpisize

  mpid%mpisize_3d(1) = 1
  mpid%mpisize_3d(2) = 1
  mpid%mpisize_3d(3) = 128

  igx = ix*mpid%mpisize_3d(1)-2*margin*(mpid%mpisize_3d(1)-1)
  jgx = jx*mpid%mpisize_3d(2)-2*margin*(mpid%mpisize_3d(2)-1)
  kgx = kx*mpid%mpisize_3d(3)-2*margin*(mpid%mpisize_3d(3)-1)

!----------------------
! determine mpid%mpirank_3d

  call setmy3drank(mpid,merr)
  call setmpiboundary(mpid)

  if (mpirank.eq.0) then
     write(6,*) 'igx=',igx,'jgx=',jgx,'kgx=',kgx
     write(*,*) 'mpirank',mpirank,'mpirank_3d',mpid%mpirank_3d
  endif

!----------------------------------------------------------------------|
!  file open

  call openfileCor(nd,mpirank,ix,jx,kx)
  call openfileAll(nd,mpirank,ix,jx,kx)

  call dacputparamc(mf_params,'comment','model_asahina')
  call dacputparami(mf_params,'ix',ix)
  call dacputparami(mf_params,'jx',jx)
  call dacputparami(mf_params,'kx',kx)
  call dacputparami(mf_params,'igx',igx)
  call dacputparami(mf_params,'jgx',jgx)
  call dacputparami(mf_params,'kgx',kgx)
  call dacputparami(mf_params,'margin',margin)
  call dacputparami(mf_params,'mpi',1)
  call dacputparami(mf_params,'mpisize',mpisize)
  call dacputparami(mf_params,'mpirank',mpirank)
  call dacputparami(mf_params,'beta',100)

!----------------------------------------------------------------------|
!  initialize counters

  nd=1
  time  = 0.0d0
  timep = 0.0d0
  ns    = 0
  merr  = 0
        
!----------------------------------------------------------------------|
!   time control parameters
!     nstop : number of total time steps for the run

  
!  tend = 1.0d0
  tend = 0.2d0

!  dtout = 0.5d0

!  tend = 0.5d0
!  dtout = 251.0d0/10.0d0
  dtout = tend/(10.0d0)
  dtout1 = tend*10.0d0

  dt = tend

  nstop = 1000000
!  nstop = 20
!----------------------------------------------------------------------|
!   setup numerical model (grid, initial conditions, etc.)
!

  call model_asahina(idf,mpid,igx,jgx,kgx,ix,jx,kx,margin &
       ,ro,pr,vx,vy,vz,bx,by,bz,phi,gm &
       ,x,dx,y,dy,z,dz &
       ,gx,gy,gz,mf_params,dtout,tend,cr,eta0,vc,eta,xin &
       ,ccx,ccy,ccz)

  call exchangeMpi(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz &
       ,phi,merr)

  call bnd_asahina(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta &
       ,x,y,gm)

!----------------------------------------------------------------------|
!  read-data


  ix0 = ix
  jx0 = jx
  kx0 = kx
  nx0 = nd

  mcont = 0
  if(mcont .eq. 1)then

     nd = 16  ! read number
     call openReadFileAll(nd,mpirank,ix0,jx0,kx0,nx0)

     read(mfi_ro) ro
     read(mfi_pr) pr

     read(mfi_vx) vx
     read(mfi_vy) vy
     read(mfi_vz) vz
     
     read(mfi_bx) bx
     read(mfi_by) by
     read(mfi_bz) bz

     read(mfi_phi) phi
     read(mfi_eta) eta

     time = real(nd-1)*dtout
     if(mpirank .eq. 0)then
        write(*,*) 'time :: ',time
     endif
     
     call closeReadFileAll()
  endif

!----------------------------------------------------------------------|
!     ready
!----------------------------------------------------------------------|

  write(mf_x) x
  write(mf_y) y
  write(mf_z) z
  
  call dacputparamd(mf_params,'gm',gm)
  
  write(mf_t) time
  write(mf_ro) ro
  write(mf_pr) pr
  write(mf_vx) vx
  write(mf_vy) vy
  write(mf_vz) vz
  write(mf_bx) bx
  write(mf_by) by
  write(mf_bz) bz

  write(mf_phi) phi
  write(mf_eta) eta

  write(mf_gx) gx
  write(mf_gz) gz

  if(mpirank .eq. 0)then
     write(6,913) ns,time,nd
     write(*,*) 'dt :: ',dt
  endif

  nd=nd+1

  call closefileAll()
  call closefileCor()
  close(mf_params)

!======================================================================|
!     time integration 
!======================================================================|
1000  continue  
  ns = ns+1
  mwflag=0

!----------------------------------------------------------------------|
!     obtain time spacing

  safety=0.2d0
  dtmin=1.d-10
  
  call getNewdt_glm_res_mpi(margin,safety,dtmin,ix,jx,kx,gm,ro,pr &
       ,vx,vy,vz,bx,by,bz,x,dx,y,dy,z,dz,eta &
       ,dt,merr,ch,mpirank)

  call mpi_allreduce(dt,dtg,1,mpi_double_precision,mpi_min &
       ,mpi_comm_world,merr)

  dt = dtg
  err_flg = 1 ! error flag
  if ( dt.lt.dtmin) goto 9999
  err_flg = 2 ! error flag
  if (merr.ne.0) goto 9999
  timep = time
  time = time+dt

  mw=0
  mw1=0
!----- check output

! dtout
  nt1=int(timep/dtout)
  nt2=int(time/dtout)
  if (nt1.lt.nt2) mw=1

! hst
  nt3=int(timep/dtout1)
  nt4=int(time/dtout1)
  if (nt3.lt.nt4) mw1=1

!---- integrate

!!$!---1st
!!$  call integrate_cyl_mp5_1(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,dt &
!!$       ,gx,gz,floor,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
!!$       ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
!!$       ,eta0,vc,eta,ccx,ccy,ccz)
!!$
!!$  call exchangeMpixz2(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1 &
!!$       ,phi1,merr)
!!$
!!$  call bnd(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1,eta,x)
!!$
!!$  call bnd_absoub(ix,jx,kx,x,z,xin,roi,pri,vxi,vyi,vzi,bxi,byi,bzi &
!!$       ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1)
!!$
!!$!---2nd
!!$  call integrate_cyl_mp5_2(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,dt &
!!$       ,gx,gz,floor,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
!!$       ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
!!$       ,eta0,vc,eta,ccx,ccy,ccz)
!!$
!!$  call exchangeMpixz2(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1 &
!!$       ,phi1,merr)
!!$
!!$  call bnd(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1,eta,x)
!!$
!!$  call bnd_absoub(ix,jx,kx,x,z,xin,roi,pri,vxi,vyi,vzi,bxi,byi,bzi &
!!$       ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1)
!!$
!!$!----3rd
!!$  call integrate_cyl_mp5_3(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,dt &
!!$       ,gx,gz,floor,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
!!$       ,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
!!$       ,eta0,vc,eta,ccx,ccy,ccz)

!------2stage version

  hdt = 0.5d0*dt

  call integrate_mp5_1st(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,hdt &
       ,gx,gz,floor,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
       ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
       ,eta0,vc,ccx,ccy,ccz)
  
  call exchangeMpi(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1 &
       ,phi1,merr)

  call bnd_asahina(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1,eta &
       ,x,y,gm)

  call integrate_mp5_2nd(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,dt &
       ,gx,gz,floor,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
       ,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
       ,eta0,vc,eta,ccx,ccy,ccz)

  call exchangeMpi(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz &
       ,phi,merr)

  call bnd_asahina(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta &
       ,x,y,gm)

!------2stage version
!----------------------------------------------------------------------|
!     data output 

  
  if(mpirank .eq. 0)then
     print '(" step=",i8," time=",e10.5," dt=",e10.5)', ns, time, dt
  end if

  if (mw.ne.0) then
     call openfileAll(nd,mpirank,ix,jx,kx)

     write(mf_t) time
     write(mf_ro) ro
     write(mf_pr) pr
     write(mf_vx) vx
     write(mf_vy) vy
     write(mf_vz) vz
     write(mf_bx) bx
     write(mf_by) by
     write(mf_bz) bz

     write(mf_phi) phi
     write(mf_eta) eta

     call closefileAll()

     if(mpirank .eq. 0)then
        write(6,913) ns,time,nd
        write(*,*) '[NORMAL] dt :: ',dt
        write(*,*) 'nd: ',nd
     end if

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

     write(mf_t) time
     write(mf_ro) ro
     write(mf_pr) pr
     write(mf_vx) vx
     write(mf_vy) vy
     write(mf_vz) vz
     write(mf_bx) bx
     write(mf_by) by
     write(mf_bz) bz

     write(mf_phi) phi
     write(mf_eta) eta

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
