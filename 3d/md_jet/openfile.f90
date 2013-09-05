!======================================================================
module openfile
!======================================================================
! Name :: openfile
!         file open
!
!======================================================================
  implicit none

  integer,public :: mf_params
  integer,public :: mf_t,mf_x,mf_y,mf_z

  integer,public :: mf_ro,mf_pr       ! primitive variables
  integer,public :: mf_vx,mf_vy,mf_vz !
  integer,public :: mf_bx,mf_by,mf_bz !
  integer,public :: mf_phi,mf_eta

  integer,public :: mfi_x,mfi_y,mfi_z ! coordianates
  integer,public :: mfi_t ! time

  integer,public :: mfi_ro,mfi_pr       ! primitive variables
  integer,public :: mfi_vx,mfi_vy,mfi_vz !
  integer,public :: mfi_bx,mfi_by,mfi_bz !
  integer,public :: mfi_phi,mfi_eta
  
  integer,public :: mf_gx,mf_gy,mf_gz ! graviation

  integer,public :: mf_lambda1,mf_lambda2

  character,private :: cno*4
  character,private :: cnond*4

contains

  subroutine openReadFileAll(nd,mpirank,ix0,jx0,kx0,nx0)
    implicit none

    integer,intent(in) :: nd,mpirank
    integer :: ix0,jx0,kx0,nx0
    integer :: mt
    mt=6
    
    write(cnond,'(i4.4)') nd
    write(cno,'(i4.4)') mpirank
    mfi_ro=70
    call dacopnr3s(mfi_ro,'ro.dac.'//cnond//'.'//cno,mt,ix0,jx0,kx0,nx0)
    mfi_pr=71
    call dacopnr3s(mfi_pr,'pr.dac.'//cnond//'.'//cno,mt,ix0,jx0,kx0,nx0)
    mfi_vx=72
    call dacopnr3s(mfi_vx,'vx.dac.'//cnond//'.'//cno,mt,ix0,jx0,kx0,nx0)
    mfi_vy=73
    call dacopnr3s(mfi_vy,'vy.dac.'//cnond//'.'//cno,mt,ix0,jx0,kx0,nx0)
    mfi_vz=74
    call dacopnr3s(mfi_vz,'vz.dac.'//cnond//'.'//cno,mt,ix0,jx0,kx0,nx0)
    mfi_bx=75
    call dacopnr3s(mfi_bx,'bx.dac.'//cnond//'.'//cno,mt,ix0,jx0,kx0,nx0)
    mfi_by=76
    call dacopnr3s(mfi_by,'by.dac.'//cnond//'.'//cno,mt,ix0,jx0,kx0,nx0)
    mfi_bz=77
    call dacopnr3s(mfi_bz,'bz.dac.'//cnond//'.'//cno,mt,ix0,jx0,kx0,nx0)

    mfi_phi=78
    call dacopnr3s(mfi_phi,'phi.dac.'//cnond//'.'//cno,mt,ix0,jx0,kx0,nx0)
    mfi_eta=79
    call dacopnr3s(mfi_eta,'eta.dac.'//cnond//'.'//cno,mt,ix0,jx0,kx0,nx0)
    
  end subroutine openReadFileAll

  subroutine closeReadFileAll()
    implicit none

    close(mfi_ro)
    close(mfi_pr)

    close(mfi_vx)
    close(mfi_vy)
    close(mfi_vz)

    close(mfi_bx)
    close(mfi_by)
    close(mfi_bz)

    close(mfi_eta)
    close(mfi_phi)

  end subroutine closeReadFileAll

  subroutine file_input(ro,pr,vx,vy,vz,bx,by,bz,phi,eta,ix,jx,kx)
    implicit none

    integer,intent(in) :: ix,jx,kx
    real(8),intent(out),dimension(ix,jx,kx) :: ro,pr,vx,vy,vz,bx,by,bz,phi,eta

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
     
  end subroutine file_input

  subroutine file_output(nd,mpirank,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,time,ix,jx,kx)
    implicit none

    integer,intent(in) :: ix,jx,kx,nd,mpirank
    real(8),intent(in),dimension(ix,jx,kx) :: ro,pr,vx,vy,vz,bx,by,bz,phi,eta
    real(8),intent(in) :: time

    write(cnond,'(i4.4)') nd
    write(cno,'(i4.4)') mpirank

    write(mf_t) time

    mf_ro=20
    call dacdef3s(mf_ro,'ro.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    write(mf_ro) ro
    close(mf_ro) 

    mf_pr=21
    call dacdef3s(mf_pr,'pr.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    write(mf_pr) pr
    close(mf_pr) 

    mf_vx=22
    call dacdef3s(mf_vx,'vx.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    write(mf_vx) vx
    close(mf_vx)

    mf_vy=23
    call dacdef3s(mf_vy,'vy.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    write(mf_vy) vy
    close(mf_vy)

    mf_vz=24
    call dacdef3s(mf_vz,'vz.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    write(mf_vz) vz
    close(mf_vz)

    mf_bx=25
    call dacdef3s(mf_bx,'bx.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    write(mf_bx) bx
    close(mf_bx)

    mf_by=26
    call dacdef3s(mf_by,'by.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    write(mf_by) by
    close(mf_by)

    mf_bz=27
    call dacdef3s(mf_bz,'bz.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    write(mf_bz) bz
    close(mf_bz)
    
    mf_phi=28
    call dacdef3s(mf_phi,'phi.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    write(mf_phi) phi
    close(mf_phi)

    mf_eta=29
    call dacdef3s(mf_eta,'eta.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    write(mf_eta) eta
    close(mf_eta)

  end subroutine file_output

  subroutine file_output_param(nd,dtout,tend,ix,jx,kx,igx,jgx,kgx,margin,mpisize &
       ,mpirank,mpisize_x,mpisize_z,eta0,vc,gm,x,y,z,dx,dy,dz,gx,gz)
    
    implicit none

    integer,intent(in) :: nd,ix,jx,kx,igx,jgx,kgx,margin,mpisize,mpirank
    integer,intent(in) :: mpisize_x,mpisize_z
    real(8),intent(in),dimension(ix) :: x,dx
    real(8),intent(in),dimension(jx) :: y,dy
    real(8),intent(in),dimension(kx) :: z,dz
    real(8),intent(in),dimension(ix,jx,kx) :: gx,gz
    real(8),intent(in) :: dtout,tend
    real(8),intent(in) :: eta0,vc,gm

    write(cnond,'(i4.4)') nd
    write(cno,'(i4.4)') mpirank

    mf_params=9
    call dacdefparam(mf_params,'params.txt.'//cno)
    
    call dacputparamc(mf_params,'comment','model_machida,int_2')
    call dacputparami(mf_params,'ix',ix)
    call dacputparami(mf_params,'jx',jx)
    call dacputparami(mf_params,'kx',kx)
    call dacputparami(mf_params,'igx',igx)
    call dacputparami(mf_params,'jgx',jgx)
    call dacputparami(mf_params,'kgx',kgx)
    call dacputparami(mf_params,'margin',margin)
    call dacputparamd(mf_params,'tend',tend)
    call dacputparami(mf_params,'mpi',1)
    call dacputparamd(mf_params,'dtout',dtout)
    call dacputparami(mf_params,'mpisize',mpisize)
    call dacputparami(mf_params,'mpirank',mpirank)
    call dacputparami(mf_params,'mpix',mpisize_x)
    call dacputparami(mf_params,'mpiz',mpisize_z)
    call dacputparamd(mf_params,'eta0',eta0)
    call dacputparamd(mf_params,'vc',vc)
    call dacputparamd(mf_params,'x(1)',x(1))
    call dacputparamd(mf_params,'y(1)',y(1))
    call dacputparamd(mf_params,'z(1)',z(1))
    call dacputparamd(mf_params,'dx(1)',dx(1))
    call dacputparamd(mf_params,'dy(1)',dy(1))
    call dacputparamd(mf_params,'dz(1)',dz(1))
    call dacputparamd(mf_params,'gm',gm)
    close(mf_params)


    mf_x=11
    call dacdef1d(mf_x,'x.dac.'//cno,6,ix)
    write(mf_x) x
    close(mf_x)

    mf_y=12
    call dacdef1d(mf_y,'y.dac.'//cno,6,jx)
    write(mf_y) y
    close(mf_y)

    mf_z=13
    call dacdef1d(mf_z,'z.dac.'//cno,6,kx)
    write(mf_z) z
    close(mf_z)

    mf_gx=40
    call dacdef3s(mf_gx,'gx.dac.'//cno,6,ix,jx,kx)
    write(mf_gx) gx
    close(mf_gx)

    close(mf_gy)

    mf_gz=42
    call dacdef3s(mf_gz,'gz.dac.'//cno,6,ix,jx,kx)
    write(mf_gz) gz
    close(mf_gz)

    mf_t =10
    call dacdef0s(mf_t,'t.dac.'//cno,6)

  end subroutine file_output_param

end module openfile
    
