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

  subroutine openfileCor(nd,mpirank,ix,jx,kx)
    implicit none

    integer,intent(in) :: ix,jx,kx
    integer,intent(in) :: nd
    integer,intent(in) :: mpirank

    write(cnond,'(i4.4)') nd
    write(cno,'(i4.4)') mpirank

    mf_params=9
    call dacdefparam(mf_params,'params.txt.'//cno)
    mf_t =10
    call dacdef0s(mf_t,'t.dac.'//cno,6)
    mf_x=11
    call dacdef1d(mf_x,'x.dac.'//cno,6,ix)
    mf_y=12
    call dacdef1d(mf_y,'y.dac.'//cno,6,jx)
    mf_z=13
    call dacdef1d(mf_z,'z.dac.'//cno,6,kx)

    mf_gx=40
    call dacdef3s(mf_gx,'gx.dac.'//cno,6,ix,jx,kx)
    mf_gz=42
    call dacdef3s(mf_gz,'gz.dac.'//cno,6,ix,jx,kx)

    return
  end subroutine openfileCor

  subroutine openfileAll(nd,mpirank,ix,jx,kx)
    implicit none

    integer,intent(in) :: ix,jx,kx
    integer,intent(in) :: nd
    integer,intent(in) :: mpirank

    write(cnond,'(i4.4)') nd
    write(cno,'(i4.4)') mpirank

    mf_ro=20
    call dacdef3s(mf_ro,'ro.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    mf_pr=21
    call dacdef3s(mf_pr,'pr.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    mf_vx=22
    call dacdef3s(mf_vx,'vx.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    mf_vy=23
    call dacdef3s(mf_vy,'vy.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    mf_vz=24
    call dacdef3s(mf_vz,'vz.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    mf_bx=25
    call dacdef3s(mf_bx,'bx.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    mf_by=26
    call dacdef3s(mf_by,'by.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    mf_bz=27
    call dacdef3s(mf_bz,'bz.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    
    mf_phi=28
    call dacdef3s(mf_phi,'phi.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    mf_eta=29
    call dacdef3s(mf_eta,'eta.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    
    return
  end subroutine openfileAll

  subroutine openfilelambda(nd,mpirank,ix,jx,kx)
    implicit none

    integer,intent(in) :: ix,jx,kx
    integer,intent(in) :: nd
    integer,intent(in) :: mpirank

    write(cnond,'(i4.4)') nd
    write(cno,'(i4.4)') mpirank

    mf_lambda1=320
    call dacdef3s(mf_lambda1,'lambda1.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    mf_lambda2=321
    call dacdef3s(mf_lambda2,'lambda2.dac.'//cnond//'.'//cno,6,ix,jx,kx)
    return
  end subroutine openfilelambda

  subroutine closefilelambda()
    implicit none

      close(mf_lambda1) 
      close(mf_lambda2) 

    end subroutine closefilelambda

  subroutine closefileCor()
    implicit none

    close(mf_x)
    close(mf_y)
    close(mf_z)

    close(mf_gx)
    close(mf_gy)
    close(mf_gz)
  end subroutine closefileCor

  subroutine closefileAll()
    implicit none

      close(mf_ro) 
      close(mf_pr) 
      close(mf_vx)
      close(mf_vy)
      close(mf_vz)
      close(mf_bx)
      close(mf_by)
      close(mf_bz)
      
      close(mf_phi)
      close(mf_eta)

    end subroutine closefileAll

  subroutine openReadFileAll(nd,mpirank,ix0,jx0,kx0,nx0)
    implicit none

    integer,intent(in) :: nd,mpirank
    integer :: ix0,jx0,kx0,nx0

    write(cnond,'(i4.4)') nd
    write(cno,'(i4.4)') mpirank

    mfi_ro=70
    call dacopnr3s(mfi_ro,'ro.dac.'//cnond//'.'//cno,6,ix0,jx0,kx0,nx0)
    mfi_pr=71
    call dacopnr3s(mfi_pr,'pr.dac.'//cnond//'.'//cno,6,ix0,jx0,kx0,nx0)

    mfi_vx=72
    call dacopnr3s(mfi_vx,'vx.dac.'//cnond//'.'//cno,6,ix0,jx0,kx0,nx0)
    mfi_vy=73
    call dacopnr3s(mfi_vy,'vy.dac.'//cnond//'.'//cno,6,ix0,jx0,kx0,nx0)
    mfi_vz=74
    call dacopnr3s(mfi_vz,'vz.dac.'//cnond//'.'//cno,6,ix0,jx0,kx0,nx0)
    
    mfi_bx=75
    call dacopnr3s(mfi_bx,'bx.dac.'//cnond//'.'//cno,6,ix0,jx0,kx0,nx0)
    mfi_by=76
    call dacopnr3s(mfi_by,'by.dac.'//cnond//'.'//cno,6,ix0,jx0,kx0,nx0)
    mfi_bz=77
    call dacopnr3s(mfi_bz,'bz.dac.'//cnond//'.'//cno,6,ix0,jx0,kx0,nx0)

    mfi_phi=78
    call dacopnr3s(mfi_phi,'phi.dac.'//cnond//'.'//cno,6,ix0,jx0,kx0,nx0)
    mfi_eta=79
    call dacopnr3s(mfi_eta,'eta.dac.'//cnond//'.'//cno,6,ix0,jx0,kx0,nx0)
    
!!$     do n=1,nd
!!$        read(mfi_t) time
!!$     enddo
    
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

end module openfile
    
