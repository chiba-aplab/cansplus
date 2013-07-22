!======================================================================
! 1 : jouleheat, 2 : MassAcc, 3 : Massout, 4 : ee
! 5 : pbx, 6 : pby, 7 : pbz, 8 : kex, 9 : key
! 10 : kez, 11 : pr, 12 : density, 13 : angular_moment
! 14 : alpha, 15 : beta, 16 : vvt
!======================================================================
module Hst2
  implicit none

  integer,public,parameter :: iHst2 = 16
  real(8),public,dimension(iHst2) :: cHst2

  integer,private :: is, ie
  integer,private :: js, je
  integer,private :: ks, ke
  real(8),private,parameter :: x_is=4.0d0, x_ie = 10.0d0
  real(8),private,parameter :: z_is=0.0d0, z_ke = 1.0d0

  integer,private :: chk_ma
  integer,private :: is_ma, ie_ma
  integer,private :: js_ma, je_ma
  integer,private :: ks_ma, ke_ma
  real(8),private,parameter :: x_ima=2.5d0, z_sma=-20.0d0, z_ema=20.0d0

  integer,private :: chk_mo
  integer,private :: is_mo, ie_mo
  integer,private :: js_mo, je_mo
  integer,private :: ks_mo, ke_mo
  real(8),private,parameter :: z_imo=20.0d0, x_smo=2.0d0, x_emo=43.0d0

  integer,parameter :: mf_hst2 = 301
contains
!===================================================================
! Check
!===================================================================

!-------------------------------------------------------------------
! check Output range
!-------------------------------------------------------------------
  subroutine chkHst2_range(mpid,ix,jx,kx,x,y,z,margin,pmargin)
    use mpi_domain_xz

    implicit none

    type(mpidomain) :: mpid

    integer,intent(in) :: ix,jx,kx,margin,pmargin
    real(8),dimension(ix) :: x
    real(8),dimension(jx) :: y
    real(8),dimension(kx) :: z

    integer :: flg1,flg2
    integer :: i,j,k

    flg1 = 0
    flg2 = 0
    if(mpid%mpirank_2d(1) .eq.0)then
       do i=pmargin+1,ix-margin
          if(x(i) < x_is)then
             is = i
          endif
       enddo
    else
       do i=margin+1,ix-margin
          if(x(i) < x_is)then
             is = i
          endif
       enddo
    endif

    if(mpid%mpirank_2d(1) .eq.0)then
       do i=pmargin+1,ix-margin
          if(x(i) > x_ie)then
             if(flg1 .eq. 0)then
                ie = i
                flg1 =1
             endif
          endif
       end do
    else
       do i=margin+1,ix-margin
          if(x(i) > x_ie)then
             if(flg1 .eq. 0)then
                ie = i
                flg1 =1
             endif
          endif
       end do
    endif

    do k=1,margin
       if(z(k) .le. z_is)then
          ks = k
       endif
    enddo
    do k=1,margin
       if(z(k) .ge. z_ke)then
          if(flg2 .eq. 0)then
             ke = k
             flg2 = 1
          endif
       endif
    enddo

    js = margin+1
    je = jx-margin

    return
  end subroutine chkHst2_range
!-------------------------------------------------------------------
! check Mass accretion 
!-------------------------------------------------------------------
  subroutine chkHst2_MassAcc(mpid,ix,jx,kx,margin,pmargin,x,dx,z,dz)
    use mpi_domain_xz
    implicit none

    type(mpidomain) :: mpid

    integer,intent(in) :: ix,jx,kx,margin,pmargin

    real(8),dimension(ix) :: x,dx
    real(8),dimension(kx) :: z,dz

    integer :: i,j,k
    integer :: flg1,flg2,flg3

    chk_ma = 0

    flg1 = 0
    flg2 = 0
    flg3 = 0

    if(mpid%mpirank_2d(1) .eq.0)then
       do i=pmargin+1,ix-margin
          if((x(i)-0.5d0*dx(i)) < x_ima)then
             if((x(i)+0.5d0*dx(i)) > x_ima)then
                is_ma = i
                flg1 = 1
             endif
          endif
       enddo
    else
       do i=margin+1,ix-margin
          if((x(i)-0.5d0*dx(i)) < x_ima)then
             if((x(i)+0.5d0*dx(i)) > x_ima)then
                is_ma = i
                flg1 = 1
             endif
          endif
       enddo
    endif
    do k=margin+1,kx-margin
       if(flg2 .eq. 0)then
          if(z(k) > z_sma)then
             ks_ma = k
             flg2 = 1
          endif
       endif
    enddo

    do k=margin+1,kx-margin
       if(z(k) < z_ema)then
          ke_ma = k
          flg3 = 1
       endif
    enddo

    chk_ma = flg1*flg2*flg3

    return
  end subroutine chkHst2_MassAcc
!-------------------------------------------------------------------
! check Mass out
!-------------------------------------------------------------------
  subroutine chkHst2_MassOut(mpid,ix,jx,kx,margin,pmargin,x,dx,z,dz)
    use mpi_domain_xz
    implicit none

    type(mpidomain) :: mpid

    integer,intent(in) :: ix,jx,kx,margin,pmargin

    real(8),dimension(ix) :: x,dx
    real(8),dimension(kx) :: z,dz

    integer :: i,j,k
    integer :: flg1,flg2,flg3

    chk_mo = 0

    flg1 = 0
    flg2 = 0
    flg3 = 0

    do k=margin+1,kx-margin
       if((z(k)-0.5d0*dz(k)) < z_imo)then
          if((z(k)+0.5d0*dz(k)) > z_imo)then
             ks_mo = k
             flg1 = 1
          endif
       endif
    enddo

    if(mpid%mpirank_2d(1) .eq.0)then
       do i=pmargin+1,ix-margin
          if(flg2 .eq. 0)then
             if(x(i) > x_smo)then
                is_mo = i
                flg2 = 1
             endif
          endif
       enddo
    else
       do i=margin+1,ix-margin
          if(flg2 .eq. 0)then
             if(x(i) > x_smo)then
                is_mo = i
                flg2 = 1
             endif
          endif
       enddo
    endif

    if(mpid%mpirank_2d(1) .eq.0)then
       do i=pmargin+1,ix-margin
          if(x(i) < x_emo)then
             ie_mo = i
             flg3 = 1
          endif
       enddo
    else
       do i=margin+1,ix-margin
          if(x(i) < x_emo)then
             ie_mo = i
             flg3 = 1
          endif
       enddo
    endif
    chk_mo = flg1*flg2*flg3
    return
  end subroutine chkHst2_MassOut
!===================================================================
! Calculate variables
!===================================================================
!-------------------------------------------------------------------
! calculate Mass Accretion
!-------------------------------------------------------------------
  subroutine calHst2_Massacc(ix,jx,kx,x,dx,dy,z,dz,ro,vx)
    implicit none
    
    integer,intent(in) :: ix,jx,kx
    
    real(8),dimension(ix) :: x,dx
    real(8),dimension(jx) :: dy
    real(8),dimension(kx) :: z,dz
    
    real(8),dimension(ix,jx,kx) :: ro,vx
    
    real(8) :: Macc,dMacc
    
    integer :: i,j,k
    
    real(8) :: dv
    
    Macc = 0.0d0
    if(chk_ma .eq. 1)then
       do k=ks_ma,ke_ma
          do j=js,je
             dv = dabs(x(is_ma))*dy(j)*dz(k)
             dMacc = ro(is_ma,j,k)*vx(is_ma,j,k)
             Macc = Macc + dMacc*dv
          enddo
       enddo
    endif

    cHst2(2) = Macc
    return
  end subroutine calHst2_Massacc
!-------------------------------------------------------------------
! calculate Mass Out
!-------------------------------------------------------------------
  subroutine calHst2_Massout(ix,jx,kx,x,dx,dy,z,dz,ro,vz)

    implicit none
    
    integer,intent(in) :: ix,jx,kx
    
    real(8),dimension(ix) :: x,dx
    real(8),dimension(jx) :: dy
    real(8),dimension(kx) :: z,dz
    
    real(8),dimension(ix,jx,kx) :: ro,vz
    
    real(8) :: Mout,dMout
    
    integer :: i,j,k
    
    real(8) :: dv
    
    Mout = 0.0d0
    if(chk_mo .eq. 1)then
       do j=js,je
          do i=is_mo,ie_mo
             dv = dabs(x(i))*dy(j)*dabs(dx(i))
             dMout = ro(i,j,ks_mo)*vz(i,j,ks_mo)
             Mout = Mout + dMout*dv
          enddo
       enddo
    endif

    cHst2(3) = Mout
    return
  end subroutine calHst2_Massout
!-------------------------------------------------------------------
! calculate all variables
!-------------------------------------------------------------------
  subroutine calHst2_all(ix,jx,kx,margin,x,dx,dy,dz,gm &
       ,ro,pr,vx,vy,vz,bxc,byc,bzc)

    implicit none

    integer,intent(in) :: ix,jx,kx,margin
    
    real(8),intent(in) :: gm

    real(8),dimension(ix) :: x,dx
    real(8),dimension(jx) :: dy
    real(8),dimension(kx) :: dz
    
    real(8),dimension(ix,jx,kx) :: ro,pr
    real(8),dimension(ix,jx,kx) :: vx,vy,vz
    real(8),dimension(ix,jx,kx) :: bxc,byc,bzc
    
    real(8) :: eet,prt,rot,amt
    real(8) :: pbxt,pbyt,pbzt
    real(8) :: kext,keyt,kezt
    real(8) :: alphat,betat
    real(8) :: bxbyt

    real(8) :: deet,dprt,drot,damt
    real(8) :: dpbxt,dpbyt,dpbzt
    real(8) :: dkext,dkeyt,dkezt
    real(8) :: dalphat,dbetat
    real(8) :: dbxbyt

    real(8) :: vvt
    real(8) :: dv
    integer :: i,j,k

    vvt = 0.0d0
    rot = 0.0d0
    eet = 0.0d0
    prt = 0.0d0
    amt = 0.0d0
    kext = 0.0d0
    keyt = 0.0d0
    kezt = 0.0d0
    pbxt = 0.0d0
    pbyt = 0.0d0
    pbzt = 0.0d0
    alphat = 0.0d0

    do k=ks,ke
       do j=js,je
          do i=is,ie
             dv = dabs(x(i))*dabs(dx(i))*dy(j)*dz(k)
             vvt = vvt + dv

             drot = ro(i,j,k)*dv
             rot = rot + drot

             dkext = (ro(i,j,k)*(vx(i,j,k)**2))*dv
             kext = kext + dkext
             dkeyt = (ro(i,j,k)*(vy(i,j,k)**2))*dv
             keyt = keyt + dkeyt
             dkezt = (ro(i,j,k)*(vz(i,j,k)**2))*dv
             kezt = kezt + dkezt

             dpbxt = (0.5d0*(bxc(i,j,k)**2))*dv
             pbxt = pbxt + dpbxt
             dpbyt = (0.5d0*(byc(i,j,k)**2))*dv
             pbyt = pbyt + dpbyt
             dpbzt = (0.5d0*(bzc(i,j,k)**2))*dv
             pbzt = pbzt + dpbzt
             
             dprt = (pr(i,j,k)/(gm-1.0d0))*dv
             prt = prt + dprt
             damt = (ro(i,j,k)*vy(i,j,k)*dabs(x(i)))*dv
             amt = amt + damt

             deet = (dpbxt+dpbyt+dpbzt+dkext+dkeyt+dkezt+dprt)*dv
             eet = eet + deet

             dbetat = 2.0d0*pr(i,j,k)/(bxc(i,j,k)**2+byc(i,j,k)**2+bzc(i,j,k)**2)
             betat = betat + dbetat
             dbxbyt = bxc(i,j,k)*byc(i,j,k)
             bxbyt = bxbyt + dbxbyt
          enddo
       enddo
    enddo
    alphat = (gm-1.0d0)*bxbyt/prt
    cHst2(4) = eet
    cHst2(5) = pbxt
    cHst2(6) = pbyt
    cHst2(7) = pbzt
    cHst2(8) = kext
    cHst2(9) = keyt
    cHst2(10) = kezt
    cHst2(11) = prt
    cHst2(12) = rot
    cHst2(13) = amt
    cHst2(14) = alphat/real((ie-is)*(je-js)*(ke-ks))
    cHst2(15) = betat/real((ie-is)*(je-js)*(ke-ks))
    cHst2(16) = vvt
    return
  end subroutine calHst2_all
!===================================================================
! Output
!===================================================================
  subroutine output_Hst2(mpirank,ns,time)
    implicit none

    include 'mpif.h'

    integer,intent(in) :: ns,mpirank
    real(8),dimension(iHst2) :: send,recv

    real(8),intent(in) :: time
    integer :: merr

    integer :: i

    do i=1,iHst2
       recv(i) = cHst2(i)
       send(i) = cHst2(i)
    enddo

    call mpi_reduce &
         (send,recv,iHst2,mpi_double_precision &
         ,mpi_sum,0,mpi_comm_world,merr)

    if(mpirank .eq. 0)then
       open(mf_hst2,file='output.hst2',status='old',position='append')
       write(mf_hst2,917) ns,time,recv
       close(mf_hst2)
    endif

917 format (1x,i8,ES14.6,16(" ",ES14.6))

    return
  end subroutine output_Hst2

  subroutine open_Hst2(mpirank)
    implicit none

    integer,intent(in) :: mpirank

    if(mpirank .eq. 0)then
       open(mf_hst2,file='output.hst2')
       write(mf_hst2,*) '1 : jouleheat, 2 : MassAcc, 3 : Massout, 4 : ee, 5 : pbx, 6 : pby, 7 : pbz, 8 : kex, 9 : key, 10 : kez, 11 : pr, 12 : density, 13 : angular_moment,14 : alpha, 15 : beta'
       close(mf_hst2)
    endif
    return
  end subroutine open_Hst2
end module Hst2
