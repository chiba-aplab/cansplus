!======================================================================
! 1 : jouleheat, 2 : MassAcc, 3 : Massout, 4 : ee
! 5 : pbx, 6 : pby, 7 : pbz, 8 : kex, 9 : key
! 10 : kez, 11 : pr, 12 : density, 13 : angular_moment
! 14 : alpha, 15 : beta, 16 : vvt 17 : beta2 18: vvtt 19:pr2 20:pb2
!======================================================================
module Hst1
  implicit none

  integer,public,parameter :: iHst1 = 23
  real(8),public,dimension(iHst1) :: cHst1

  integer,private :: is, ie
  integer,private :: js, je
  integer,private :: ks, ke

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

  integer,private :: chk_bt
  integer,private :: is_bt, ie_bt
  integer,private :: js_bt, je_bt
  integer,private :: ks_bt, ke_bt
  real(8),private,parameter :: z_sbt=0.0d0, z_ebt=1.0d0, x_sbt=7.0d0, x_ebt=80.0d0

  integer,private :: chk_bt2
  integer,private :: is_bt2, ie_bt2
  integer,private :: js_bt2, je_bt2
  integer,private :: ks_bt2, ke_bt2
  real(8),private,parameter :: z_sbt2=1.0d0, z_ebt2=2.0d0, x_sbt2=7.0d0, x_ebt2=80.0d0

  integer,parameter :: mf_hst1 = 300
contains
!===================================================================
! Check
!===================================================================

!-------------------------------------------------------------------
! check Output range
!-------------------------------------------------------------------
  subroutine chkHst1_range(mpid,ix,jx,kx,margin,pmargin)
    use mpi_domain_xz
    implicit none
    
    type(mpidomain) :: mpid

    integer,intent(in) :: ix,jx,kx,margin,pmargin

    integer :: i,j,k

    if(mpid%mpirank_2d(1) .eq.0)then
       is = pmargin+1
    else
       is = margin+1
    endif

    js = margin+1
    ks = margin+1

    ie = ix-margin
    je = jx-margin
    ke = kx-margin
    return
  end subroutine chkHst1_range
!-------------------------------------------------------------------
! check Mass accretion 
!-------------------------------------------------------------------
  subroutine chkHst1_MassAcc(mpid,ix,jx,kx,margin,pmargin,x,dx,z,dz)
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
  end subroutine chkHst1_MassAcc
!-------------------------------------------------------------------
! check Mass out
!-------------------------------------------------------------------
  subroutine chkHst1_MassOut(mpid,ix,jx,kx,margin,pmargin,x,dx,z,dz)
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
  end subroutine chkHst1_MassOut
!-------------------------------------------------------------------
! check Output beta
!-------------------------------------------------------------------
  subroutine chkHst1_bt(mpid,ix,jx,kx,margin,x,dx,z,dz)
    use mpi_domain_xz
    implicit none

    type(mpidomain) :: mpid
    
    real(8),dimension(ix) :: x,dx
    real(8),dimension(kx) :: z,dz

    integer,intent(in) :: ix,jx,kx,margin
    integer :: i,j,k

    integer :: flg1,flg2

    flg1 = 0
    flg2 = 0

    do i=1,ix
       if(x(i) .ge. x_sbt)then
          if(x(i) .lt. x_ebt)then
             flg1 = 1
             ie_bt = i
          endif
       endif
    enddo

    if(flg1 .eq. 1)then
       do i=ix,1,-1
          if(x(i) .lt. x_ebt)then
             if(x(i) .ge. x_sbt)then
                is_bt = i
             endif
          endif
       enddo
    endif

    do k=1,kx
       if(z(k) .ge. z_sbt)then
          if(z(k) .lt. z_ebt)then
             flg2 = 1
             ke_bt = k
          endif
       endif
    enddo

    if(flg2 .eq. 1)then
       do k=kx,1,-1
          if(z(k) .lt. z_ebt)then
             if(z(k) .ge. z_sbt)then
                ks_bt = k
             endif
          endif
       enddo
    endif

    chk_bt = flg1*flg2
    return
  end subroutine chkHst1_bt

  subroutine chkHst1_bt2(mpid,ix,jx,kx,margin,x,dx,z,dz)
    use mpi_domain_xz
    implicit none

    type(mpidomain) :: mpid
    
    real(8),dimension(ix) :: x,dx
    real(8),dimension(kx) :: z,dz

    integer,intent(in) :: ix,jx,kx,margin
    integer :: i,j,k

    integer :: flg1,flg2

    flg1 = 0
    flg2 = 0

    do i=1,ix
       if(x(i) .ge. x_sbt2)then
          if(x(i) .lt. x_ebt2)then
             flg1 = 1
             ie_bt2 = i
          endif
       endif
    enddo

    if(flg1 .eq. 1)then
       do i=ix,1,-1
          if(x(i) .lt. x_ebt2)then
             if(x(i) .ge. x_sbt2)then
                is_bt2 = i
             endif
          endif
       enddo
    endif

    do k=1,kx
       if(z(k) .ge. z_sbt2)then
          if(z(k) .lt. z_ebt2)then
             flg2 = 1
             ke_bt2 = k
          endif
       endif
    enddo

    if(flg2 .eq. 1)then
       do k=kx,1,-1
          if(z(k) .lt. z_ebt2)then
             if(z(k) .ge. z_sbt2)then
                ks_bt2 = k
             endif
          endif
       enddo
    endif

    chk_bt2 = flg1*flg2
    return
  end subroutine chkHst1_bt2

!===================================================================
! Calculate variables
!===================================================================
!-------------------------------------------------------------------
! calculate Mass Accretion
!-------------------------------------------------------------------
  subroutine calHst1_Massacc(ix,jx,kx,x,dx,dy,z,dz,ro,vx)
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

    cHst1(2) = Macc
    return
  end subroutine calHst1_Massacc
!-------------------------------------------------------------------
! calculate Mass Out
!-------------------------------------------------------------------
  subroutine calHst1_Massout(ix,jx,kx,x,dx,dy,z,dz,ro,vz)

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

    cHst1(3) = Mout
    return
  end subroutine calHst1_Massout
!-------------------------------------------------------------------
! calculate all variables
!-------------------------------------------------------------------
  subroutine calHst1_all(ix,jx,kx,margin,x,dx,dy,dz,gm &
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
    real(8) :: pbxt,pbyt,pbzt,pbt
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

    real(8) :: dv_beta,vvt_beta
    real(8) :: inbeta

    integer :: flg_beta
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
    betat = 0.0d0

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

    cHst1(4) = eet
    cHst1(5) = pbxt
    cHst1(6) = pbyt
    cHst1(7) = pbzt
    cHst1(8) = kext
    cHst1(9) = keyt
    cHst1(10) = kezt
    cHst1(11) = prt
    cHst1(12) = rot
    cHst1(13) = amt
    cHst1(14) = alphat/real((ie-is)*(je-js)*(ke-ks))
    cHst1(15) = betat/real((ie-is)*(je-js)*(ke-ks))
    cHst1(16) = vvt

    return
  end subroutine calHst1_all

  subroutine calbeta2(ix,jx,kx,bxc,byc,bzc,pr,x,dx,dy,z,dz)
    implicit none

    integer :: ix,jx,kx

    real(8),dimension(ix) :: x,dx
    real(8),dimension(jx) :: dy
    real(8),dimension(kx) :: z,dz

    real(8),dimension(ix,jx,kx) :: bxc,byc,bzc
    real(8),dimension(ix,jx,kx) :: pr

    real(8) :: betat2,dbetat2
    real(8) :: vvtt 

    real(8) :: dpb2,pb2,dpr2,pr2
    real(8) :: dv,pbt

    integer :: i,j,k

    betat2 = 0.0d0
    vvtt = 0.0d0
    pr2 = 0.0d0
    pb2 = 0.0d0

    if(chk_bt .eq. 1)then
       do k=ks_bt,ke_bt
          do j=js_bt,je_bt
             do i=is_bt,ie_bt
                dv = dabs(x(i))*dabs(dx(i))*dy(j)*dz(k)
                vvtt = vvtt + dv

                dbetat2 = 2.0d0*pr(i,j,k)/(bxc(i,j,k)**2+byc(i,j,k)**2+bzc(i,j,k)**2)
                betat2 = betat2 + dbetat2

                dpr2 = pr(i,j,k)*dv
                pr2 = pr2 + dpr2

                dpb2 = 0.5d0*(bxc(i,j,k)**2+byc(i,j,k)**2+bzc(i,j,k)**2)
                pb2 = pb2 + dpb2
             enddo
          enddo
       enddo
    endif

    cHst1(17) = betat2
    cHst1(18) = vvtt
    cHst1(19) = pr2
    cHst1(20) = pb2
    return
  end subroutine calbeta2

  subroutine calbeta3(ix,jx,kx,bxc,byc,bzc,pr,x,dx,dy,z,dz)
    implicit none

    integer :: ix,jx,kx

    real(8),dimension(ix) :: x,dx
    real(8),dimension(jx) :: dy
    real(8),dimension(kx) :: z,dz

    real(8),dimension(ix,jx,kx) :: bxc,byc,bzc
    real(8),dimension(ix,jx,kx) :: pr

    real(8) :: betat2,dbetat2
    real(8) :: vvtt 

    real(8) :: dpb2,pb2,dpr2,pr2
    real(8) :: dv,pbt

    integer :: i,j,k

    betat2 = 0.0d0
    vvtt = 0.0d0
    pr2 = 0.0d0
    pb2 = 0.0d0

    if(chk_bt2 .eq. 1)then
       do k=ks_bt2,ke_bt2
          do j=js_bt2,je_bt2
             do i=is_bt2,ie_bt2
                dv = dabs(x(i))*dabs(dx(i))*dy(j)*dz(k)
                vvtt = vvtt + dv

                dbetat2 = 2.0d0*pr(i,j,k)/(bxc(i,j,k)**2+byc(i,j,k)**2+bzc(i,j,k)**2)
                betat2 = betat2 + dbetat2

                dpr2 = pr(i,j,k)*dv
                pr2 = pr2 + dpr2

                dpb2 = 0.5d0*(bxc(i,j,k)**2+byc(i,j,k)**2+bzc(i,j,k)**2)
                pb2 = pb2 + dpb2
             enddo
          enddo
       enddo
    endif

    cHst1(21) = vvtt
    cHst1(22) = pr2
    cHst1(23) = pb2
    return
  end subroutine calbeta3

!===================================================================
! Output
!===================================================================
  subroutine output_Hst1(mpirank,ns,time)
    implicit none

    include 'mpif.h'

    integer,intent(in) :: ns,mpirank
    real(8),dimension(iHst1) :: send,recv

    real(8),intent(in) :: time
    integer :: merr

    integer :: i

    do i=1,iHst1
       recv(i) = cHst1(i)
       send(i) = cHst1(i)
    enddo

    call mpi_reduce &
         (send,recv,iHst1,mpi_double_precision &
         ,mpi_sum,0,mpi_comm_world,merr)

    if(mpirank .eq. 0)then
       open(mf_hst1,file='output.hst1',status='old',position='append')
       write(mf_hst1,917) ns,time,recv
       close(mf_hst1)
    endif

917 format (1x,i8,ES14.6,23(" ",ES14.6))

    return
  end subroutine output_Hst1

  subroutine open_Hst1(mpirank)
    implicit none

    integer,intent(in) :: mpirank

    if(mpirank .eq. 0)then
       open(mf_hst1,file='output.hst1')
       write(mf_hst1,*) '1: step, 2: time, 3: jouleheat, 4: MassAcc, 5: Massout, 6: ee, 7: pbx, 8: pby, 9: pbz, 10: kex, 11: key, 12: kez'
       write(mf_hst1,*) '13: pr, 14: density, 15: angular_moment, 16: alpha, 17: beta, 18: vvt, 19: beta2, 20:vvtt, 21: pr2, 22: pb2, 23:vvtt2, 24: pr2, 25: pb2'
       close(mf_hst1)
    endif
    return
  end subroutine open_Hst1
end module Hst1
