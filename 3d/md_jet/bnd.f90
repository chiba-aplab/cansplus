subroutine bnd(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,x,z &
           ,xin,roi,pri,vxi,vyi,vzi,bxi,byi,bzi)

  use mpi_domain_xz
  use boundary
  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  type(mpidomain) :: mpid

  real(8),dimension(ix),intent(in) :: x
  real(8),dimension(kx),intent(in) :: z

  real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(inout) :: phi,eta

  real(8),dimension(ix,jx,kx),intent(in) :: roi,pri
  real(8),dimension(ix,jx,kx),intent(in) :: vxi,vyi,vzi
  real(8),dimension(ix,jx,kx),intent(in) :: bxi,byi,bzi

  real(8),intent(in) :: xin
  
  integer :: i,j,k

  call bd_pery(margin,ro,ix,jx,kx)
  call bd_pery(margin,pr,ix,jx,kx)
  
  call bd_pery(margin,vx,ix,jx,kx)
  call bd_pery(margin,vy,ix,jx,kx)
  call bd_pery(margin,vz,ix,jx,kx)
  
  call bd_pery(margin,phi,ix,jx,kx)
  call bd_pery(margin,eta,ix,jx,kx)
  
  call bd_pery(margin,bx,ix,jx,kx)
  call bd_pery(margin,by,ix,jx,kx)
  call bd_pery(margin,bz,ix,jx,kx)

!======================================================================
! inner x-boundary

  if(mpid%mpirank_2d(1) .eq.0)then

     call bd_synpx(0,margin,ro,ix,jx,kx)
     call bd_synpx(0,margin,pr,ix,jx,kx)
     call bd_synnx(0,margin,vx,ix,jx,kx)
     call bd_synnx(0,margin,vy,ix,jx,kx)
     call bd_synpx(0,margin,vz,ix,jx,kx)
     
     call bd_synnx(0,margin,bx,ix,jx,kx)
     call bd_synnx(0,margin,by,ix,jx,kx)
     call bd_synpx(0,margin,bz,ix,jx,kx)
     
     call bd_synpx(0,margin,phi,ix,jx,kx)
     call bd_synpx(0,margin,eta,ix,jx,kx)

  end if
!======================================================================
! outer x-boudary

  if(mpid%mpirank_2d(1) .eq. (mpid%mpisize_2d(1)-1))then

     call bd_frex(1,margin,ro,ix,jx,kx)
     call bd_frex(1,margin,pr,ix,jx,kx)
     
     call bd_frex(1,margin,vx,ix,jx,kx)
     call bd_frex(1,margin,vy,ix,jx,kx)
     call bd_frex(1,margin,vz,ix,jx,kx)

     call bd_frex(1,margin,bx,ix,jx,kx)
     call bd_frex(1,margin,by,ix,jx,kx)
     call bd_frex(1,margin,bz,ix,jx,kx)

     call bd_frex(1,margin,phi,ix,jx,kx)
     call bd_frex(1,margin,eta,ix,jx,kx)
  end if
!----------------------------------------------------------------------|
! outer z-boundary
!


  if(mpid%mpirank_2d(2) .eq. (mpid%mpisize_2d(2)-1)) then  
     call bd_frez(1,margin,ro,ix,jx,kx)
     call bd_frez(1,margin,pr,ix,jx,kx)
     
     call bd_frez(1,margin,vx,ix,jx,kx)
     call bd_frez(1,margin,vy,ix,jx,kx)
     call bd_frez(1,margin,vz,ix,jx,kx)
     
     call bd_frez(1,margin,bx,ix,jx,kx)
     call bd_frez(1,margin,by,ix,jx,kx)
     call bd_frez(1,margin,bz,ix,jx,kx)

     call bd_frez(1,margin,phi,ix,jx,kx)
     call bd_frez(1,margin,eta,ix,jx,kx)
  end if

!----------------------------------------------------------------------|
! inner z-boundary
!

  if(mpid%mpirank_2d(2) .eq. 0)then

     call bd_frez(0,margin,ro,ix,jx,kx)
     call bd_frez(0,margin,pr,ix,jx,kx)
     
     call bd_frez(0,margin,vx,ix,jx,kx)
     call bd_frez(0,margin,vy,ix,jx,kx)
     call bd_frez(0,margin,vz,ix,jx,kx)
     
     call bd_frez(0,margin,bx,ix,jx,kx)
     call bd_frez(0,margin,by,ix,jx,kx)
     call bd_frez(0,margin,bz,ix,jx,kx)

     call bd_frez(0,margin,phi,ix,jx,kx)
     call bd_frez(0,margin,eta,ix,jx,kx)
  end if

  do k=1,margin
     do j=1,jx
        do i=1,ix
           if (x(i).le.1d0)then
              ro(i,j,k) = ro_jet
              pr(i,j,k) = pr_jet
              vx(i,j,k) = vx_jet
              vy(i,j,k) = vy_jet
              vz(i,j,k) = vz_jet
              bx(i,j,k) = bx_jet
              by(i,j,k) = by_jet
              bz(i,j,k) = bz_jet
           endif
        enddo
     enddo
  enddo

return
end subroutine bnd

subroutine absorb(ix,jx,kx,x,z,xin,qq,qqi)
  implicit none

  integer,intent(in) :: ix,jx,kx

  real(8),intent(in) :: xin

  real(8),dimension(ix) :: x
  real(8),dimension(kx) :: z

  real(8),dimension(ix,jx,kx) :: qq,qqi

  real(8) :: ai,halfx,ss,dx0

  integer :: i,j,k

  dx0 = 0.01d0
  do k=1,kx
    do j=1,jx
      do i=1,ix
         ss = sqrt(x(i)**2+z(k)**2)
         if(ss .le. xin)then
         ai = 0.1d0*(1.0d0-tanh((ss-xin+5.0d0*dx0)/(2.0d0*dx0)))
         qq(i,j,k) = (1.0d0-ai)*qq(i,j,k)+ai*qqi(i,j,k)
         endif
      enddo
    enddo
  enddo

return
end subroutine absorb

