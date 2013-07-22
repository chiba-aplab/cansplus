!
! $Log: 
!
program main
  use mpi_domain_xz
  use openfile

  implicit none
  include 'mpif.h'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer,parameter :: FLAG_ModAverage=0 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
   integer,parameter :: margin = 3
!   integer,parameter :: margin = 6

!  integer,parameter :: ix = 25+2*margin,jx=32+2*margin,kx=30+2*margin
!  integer,parameter :: ix = 200+2*margin,jx=32+2*margin,kx=8+2*margin

!  integer,parameter :: ix = 256+2*margin,jx=64+2*margin,kx=4+2*margin




!  integer,parameter :: ix = 32+2*margin,jx=64+2*margin,kx=16+2*margin
!  integer,parameter :: ix = 16+2*margin,jx=64+2*margin,kx=8+2*margin
  integer,parameter :: ix = 1+2*margin,jx=1+2*margin,kx=5+2*margin




!  integer,parameter :: ix = 256+2*margin,jx=64+2*margin,kx=4+2*margin
!  integer,parameter :: ix = 100+2*margin,jx=32+2*margin,kx=4+2*margin

!  integer,parameter :: ix = 100+2*margin,jx=32+2*margin,kx=100+2*margin
!  integer,parameter :: ix = 50+2*margin,jx=32+2*margin,kx=1+2*margin

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

!--av
  real(8) :: baseAV
  real(8) :: pbAV,prAV
  real(8) :: ebAV,roAv,eeAv,amAv

  real(8),parameter :: xmin_MAV=0.7d0
  real(8),parameter :: xmax_MAV=1.3d0
  real(8),parameter :: zmin_MAV=0.0d0
  real(8),parameter :: zmax_MAV=0.3d0
!--tcr
  real(8),dimension(kx) :: baseTCR

  real(8),dimension(kx) :: bxTCR,byTCR,bzTCR
  real(8),dimension(kx) :: rxTCR,ryTCR,rzTCR
  real(8),dimension(kx) :: roTCR,pbTCR,prTCR

  real(8),parameter :: xObs1 = 0.3d0
  real(8),parameter :: xObs2 = 0.5d0
  real(8),parameter :: xObs3 = 0.75d0
  real(8),parameter :: xObs4 = 1.0d0
  real(8),parameter :: xObs5 = 1.5d0

!--cooling
  real(8),parameter :: velocity_c = 2.9979d0*(10.0d0**10)                                       !cgs
  real(8),parameter :: mass_bh = 10.0d0                                                         !per mass_solar
  real(8),parameter :: mass_proton = 1.6726d0*(10.0d0**(-24))                                   !cgs
  real(8),parameter :: boltzmann_const = 1.3807d0*(10.0d0**(-16))                               !cgs
  real(8),parameter :: Navo = 6.0221d0*(10.0d0**23)                                             !Avogadro's constant
!  real(8),parameter :: mmw = 0.62d0									!mean molecular weight
  real(8),parameter :: mmw = 1.0d0                                                              !mean molecular weight
  real(8),parameter :: rg = 3.0d0*(10.0d0**6)*(mass_bh/10.0d0)                                  !cgs Schwarzschild radius
!  real(8),parameter :: rg_nrmlx = 1.0d0/50.0d0					         	!per nrmlx
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

! test1028
  
!======================================================================|
!     prologue
!======================================================================|

  integer :: mcont,ndi

!----------------------------------------------------------------------|
!  initialize counters

  integer :: nd,ns,merr,ns1,ns2,nscount
  real(8) :: time,timep,swtch_t

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
  real(8) :: te_limit
 
  real(8) :: temp
!======================================================================|
!     prologue
!======================================================================|
  floor = 1.d-6

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

!  mpid%mpisize_2d(1) = 8
!  mpid%mpisize_2d(2) = 16
  mpid%mpisize_2d(1) = 1
  mpid%mpisize_2d(2) = 64
!  mpid%mpisize_2d(1) = 2
!  mpid%mpisize_2d(2) = 8

  igx = ix*mpid%mpisize_2d(1)-2*margin*(mpid%mpisize_2d(1)-1)
  jgx = jx
  kgx = kx*mpid%mpisize_2d(2)-2*margin*(mpid%mpisize_2d(2)-1)

!----------------------
! determine mpid%mpirank_3d

  call setmy2drank(mpid,merr)
  call setmpiboundary(mpid)

  if (mpirank.eq.0) then
     write(6,*) 'igx=',igx,'kgx=',kgx
     write(*,*) 'mpirank',mpirank,'mpirank_2d',mpid%mpirank_2d
  endif

!----------------------------------------------------------------------|
!  file open
write(*,*) nd,mpirank,ix,jx,kx
  call openfileCor(nd,mpirank,ix,jx,kx)
  call openfileAll(nd,mpirank,ix,jx,kx)
 
  call dacputparamc(mf_params,'comment','model_machida,int_2')
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
  call dacputparami(mf_params,'mpix',mpid%mpisize_2d(1))
  call dacputparami(mf_params,'mpiz',mpid%mpisize_2d(2))
  call dacputparami(mf_params,'beta',100)
  call dacputparamd(mf_params,'dtout',dtout)
  call dacputparamd(mf_params,'nrmlro',nrmlro)
  call dacputparamd(mf_params,'nrmlte',nrmlte)
  call dacputparamd(mf_params,'nrmlx',nrmlx)
  call dacputparamd(mf_params,'nrmlv',nrmlv)
  call dacputparamd(mf_params,'nrmlt',nrmlt)
  call dacputparamd(mf_params,'nrmlee',nrmlee)
  call dacputparamd(mf_params,'mass_bh',mass_bh)
  call dacputparamd(mf_params,'rg',rg)
  call dacputparamd(mf_params,'rg_nrmlx',rg_nrmlx)
  call dacputparamd(mf_params,'RadCool',RadCool)
  call dacputparamd(mf_params,'te_factor',te_factor)
!  call dacputparamd(mf_params,'te_limit',te_limit)

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

  
!  tend = 1.0d0
!  tend = 0.4d0
  tend = 0.36d0
!  tend = 6.28d0

!  dtout = 0.5d0

!  tend = 0.5d0
!  dtout = 251.0d0/10.0d0
  dtout = tend/(10.0d0)
  dtout1 = dtout/10.0d0

!  tend = 6.28d0*100.0d0
!  tend = 6.28d0*60.0d0
!  tend = 6.28d0*3.0d0

  dt = tend

  nstop = 2000000
!nstop=10
!  nstop = 100000
!----------------------------------------------------------------------|
!   setup numerical model (grid, initial conditions, etc.)
!

  call model_machida(idf,mpid,igx,jgx,kgx,ix,jx,kx,margin &
       ,ro,pr,vx,vy,vz,bx,by,bz,phi,gm &
       ,x,dx,y,dy,z,dz &
       ,gx,gz,mf_params,dtout,tend,cr,eta0,vc,eta,xin &
       ,ccx,ccy,ccz,rg_nrmlx,rohalo,te_limit,te_factor,nrmlv &
       ,nrmlte,boltzmann_const,Navo,mmw)

  call exchangeMpixz2(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz &
       ,phi,merr)

  call bnd(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,x)

  call saveInit_3d(ix,jx,kx,ro,roi)
  call saveInit_3d(ix,jx,kx,pr,pri)

  call saveInit_3d(ix,jx,kx,vx,vxi)
  call saveInit_3d(ix,jx,kx,vy,vyi)
  call saveInit_3d(ix,jx,kx,vz,vzi)

  call saveInit_3d(ix,jx,kx,bx,bxi)
  call saveInit_3d(ix,jx,kx,by,byi)
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
   call closefileAll()
   call closefileCor()
   close(mf_params)

!----------------------------------------------------------------------|
!  read-data


  ix0 = ix
  jx0 = jx
  kx0 = kx
  nx0 = nd

  open(100,file='readFileNumber.dat')
  read(100,*) nd
  close(100)
  if(nd .ne. 1)then
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
  nd=nd+1



! test1028
!if (FLAG_ModAverage .eq. 1) then
!  call openfileModTCR(mpirank,margin,ix,kx,x,xObs1,xObs2,xObs3,xObs4,xObs5)
!  call closefileModTCR()
!  call open_CalPhysValue(mpirank,'averageVal',mf_calPhysValue)
!end if
! test1028

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
  
!  call getNewdt_glmcyl_res_mpi(margin,safety,dtmin,ix,jx,kx,gm,ro,pr &
!       ,vx,vy,vz,bx,by,bz,x,dx,y,dy,z,dz,eta &
!       ,dt,merr,ch,mpirank)

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

  call integrate_cyl_mp5_1st(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,hdt &
       ,gx,gz,floor,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
       ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
       ,eta0,vc,ccx,ccy,ccz,RadCool,te_factor,time,rohalo,swtch_t,xin)
  
  call exchangeMpixz2(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1 &
       ,phi1,merr)

  call bnd(mpid,margin,ix,jx,kx,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1,eta,x)

  call bnd_absoub(ix,jx,kx,x,z,xin,roi,pri,vxi,vyi,vzi,bxi,byi,bzi &
       ,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1)

  call integrate_cyl_mp5_2nd(margin,ix,jx,kx,gm,x,dx,y,dy,z,dz,dt &
       ,gx,gz,floor,ro1,pr1,vx1,vy1,vz1,bx1,by1,bz1,phi1 &
       ,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cr &
       ,eta0,vc,eta,ccx,ccy,ccz,RadCool,te_factor,time,rohalo,swtch_t,xin)

  call exchangeMpixz2(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz &
       ,phi,merr)

  call bnd(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,x)

  call bnd_absoub(ix,jx,kx,x,z,xin,roi,pri,vxi,vyi,vzi,bxi,byi,bzi &
       ,ro,pr,vx,vy,vz,bx,by,bz)

!------2stage version
!----------------------------------------------------------------------|
!     data output 

  
  if(mpirank .eq. 0)then
!   if(nscount .ge. 100)then
     print '(" step=",i8," time=",e10.5," dt=",e10.5)', ns, time, dt
!     nscount=0
!   endif
  end if

!-----test1028----------------------
!   ModAverage test

!if (FLAG_ModAverage .eq. 1) then
!  if (mw1.ne.0) then  
!     do k=1,kx
!        do j=1,jx
!           do i=1,ix
!              pb(i,j,k) = 0.5d0*(bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2)
!              rx(i,j,k) = ro(i,j,k)*vx(i,j,k)
!              ry(i,j,k) = ro(i,j,k)*vy(i,j,k)
!              rz(i,j,k) = ro(i,j,k)*vz(i,j,k)
!              am(i,j,k) = vy(i,j,k)*x(i)
!              ee(i,j,k) = pb(i,j,k) &
!                   +0.5d0*ro(i,j,k)*(vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2) &
!                   +pr(i,j,k)/(gm-1.0d0)
!           enddo
!        enddo
!     enddo
!     
!     call calSumBase(margin,ix,jx,kx,x,y,z,xmin_MAV,xmax_MAV,zmin_MAV,zmax_MAV &
!          ,baseAV)
!     
!     call calSum(margin,ix,jx,kx,x,y,z,xmin_MAV,xmax_MAV,zmin_MAV,zmax_MAV &
!          ,pb,pbAV)
!     call calSum(margin,ix,jx,kx,x,y,z,xmin_MAV,xmax_MAV,zmin_MAV,zmax_MAV &
!          ,pr,prAV)
!     call calSum(margin,ix,jx,kx,x,y,z,xmin_MAV,xmax_MAV,zmin_MAV,zmax_MAV &
!          ,ro,roAV)
!     call calSum(margin,ix,jx,kx,x,y,z,xmin_MAV,xmax_MAV,zmin_MAV,zmax_MAV &
!          ,ee,eeAV)
!     call calSum(margin,ix,jx,kx,x,y,z,xmin_MAV,xmax_MAV,zmin_MAV,zmax_MAV &
!          ,am,amAV)
!
!     avInterVal(1) = baseAv
!     avInterVal(2) = pbAv
!     avInterVal(3) = prAv
!     avInterVal(4) = roAv
!     avInterVal(5) = eeAv
!     avInterVal(6) = amAv
!
!     call exchangeSumVal(mpirank,iInterNum,avInterVal)
!     
!     if(mpirank == 0)then
!        pbAv = avInterVal(2)/avInterVal(1)
!        prAv = avInterVal(3)/avInterVal(1)
!        roAv = avInterVal(4)/avInterVal(1)
!        eeAv = avInterVal(5)/avInterVal(1)
!        amAv = avInterVal(6)/avInterVal(1)
!
!        avVal(1) = pbAv
!        avVal(2) = prAv
!        avVal(3) = roAv
!        avVal(4) = eeAv
!        avVal(5) = amAv
!
!        call outputCalPhysValue(mpirank,iCalNum,avVal,ns,time &
!             ,'averageVal',mf_calPhysValue)
!     endif
!!   TimeChnageav@Rad test
!     call openfileModTCRsOld(mpirank,margin,ix,x,xObs1,xObs2,xObs3,xObs4,xObs5)
!
!     call TimeChangeVal(margin,ix,jx,kx,x,y,z,xObs1 &
!          ,ro,pb,pr,rx,ry,rz,bx,by,bz &
!          ,baseTCR,roTCR,pbTCR,prTCR,rxTCR,ryTCR,rzTCR,bxTCR,byTCR,bzTCR)
!     
!     if(isInRadius(margin,ix,x,xObs1) == 1)then
!        write(mtr1_bx) bxTCR
!        write(mtr1_by) byTCR
!        write(mtr1_bz) bzTCR
!        write(mtr1_rx) rxTCR
!        write(mtr1_ry) ryTCR
!        write(mtr1_rz) rzTCR
!        write(mtr1_ro) roTCR
!        write(mtr1_pb) pbTCR
!        write(mtr1_pr) prTCR
!     end if
!
!     call TimeChangeVal(margin,ix,jx,kx,x,y,z,xObs2 &
!          ,ro,pb,pr,rx,ry,rz,bx,by,bz &
!          ,baseTCR,roTCR,pbTCR,prTCR,rxTCR,ryTCR,rzTCR,bxTCR,byTCR,bzTCR)
!     
!     if(isInRadius(margin,ix,x,xObs2) == 1)then
!        write(mtr2_bx) bxTCR
!        write(mtr2_by) byTCR
!        write(mtr2_bz) bzTCR
!        write(mtr2_rx) rxTCR
!        write(mtr2_ry) ryTCR
!        write(mtr2_rz) rzTCR
!        write(mtr2_ro) roTCR
!        write(mtr2_pb) pbTCR
!        write(mtr2_pr) prTCR
!     end if
!
!     call TimeChangeVal(margin,ix,jx,kx,x,y,z,xObs3 &
!          ,ro,pb,pr,rx,ry,rz,bx,by,bz &
!          ,baseTCR,roTCR,pbTCR,prTCR,rxTCR,ryTCR,rzTCR,bxTCR,byTCR,bzTCR)
!     
!     if(isInRadius(margin,ix,x,xObs3) == 1)then
!        write(mtr3_bx) bxTCR
!        write(mtr3_by) byTCR
!        write(mtr3_bz) bzTCR
!        write(mtr3_rx) rxTCR
!        write(mtr3_ry) ryTCR
!        write(mtr3_rz) rzTCR
!        write(mtr3_ro) roTCR
!        write(mtr3_pb) pbTCR
!        write(mtr3_pr) prTCR
!     end if
!
!     call TimeChangeVal(margin,ix,jx,kx,x,y,z,xObs4 &
!          ,ro,pb,pr,rx,ry,rz,bx,by,bz &
!          ,baseTCR,roTCR,pbTCR,prTCR,rxTCR,ryTCR,rzTCR,bxTCR,byTCR,bzTCR)
!     
!     if(isInRadius(margin,ix,x,xObs4) == 1)then
!        write(mtr4_bx) bxTCR
!        write(mtr4_by) byTCR
!        write(mtr4_bz) bzTCR
!        write(mtr4_rx) rxTCR
!        write(mtr4_ry) ryTCR
!        write(mtr4_rz) rzTCR
!        write(mtr4_ro) roTCR
!        write(mtr4_pb) pbTCR
!        write(mtr4_pr) prTCR
!     end if
!
!     call TimeChangeVal(margin,ix,jx,kx,x,y,z,xObs5 &
!          ,ro,pb,pr,rx,ry,rz,bx,by,bz &
!          ,baseTCR,roTCR,pbTCR,prTCR,rxTCR,ryTCR,rzTCR,bxTCR,byTCR,bzTCR)
!     
!     if(isInRadius(margin,ix,x,xObs5) == 1)then
!        write(mtr5_bx) bxTCR
!        write(mtr5_by) byTCR
!        write(mtr5_bz) bzTCR
!        write(mtr5_rx) rxTCR
!        write(mtr5_ry) ryTCR
!        write(mtr5_rz) rzTCR
!        write(mtr5_ro) roTCR
!        write(mtr5_pb) pbTCR
!        write(mtr5_pr) prTCR
!     end if
!
!     call closefileModTCR()
!
!     if(mpirank .eq. 0)then
!        write(6,913) ns,time,nd1
!        write(*,*) '[TCR][AV] dt :: ',dt
!        write(*,*) 'nd1: ',nd1
!     end if
!     nd1 = nd1+1
!  endif
!
!endif
!-----test1028----------------------

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

! test1028
!  call closefileModTCR()
! test1028
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
