subroutine bnd_asahina(mpid,margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta &
     ,x,y,gm)
  use mpi_domain
  implicit none

  integer,intent(in) :: margin,ix,jx,kx
  real(8),intent(in) :: gm
  type(mpidomain) :: mpid

  real(8),dimension(ix) :: x,y

  real(8),dimension(ix,jx,kx) :: ro,pr
  real(8),dimension(ix,jx,kx) :: vx,vy,vz
  real(8),dimension(ix,jx,kx) :: bx,by,bz
  real(8),dimension(ix,jx,kx) :: phi,eta

  integer :: i,j,k
  real :: r,pi,si,co,br,bp,vr,vp,f

  pi=acos(-1.0d0)
  if(mpid%mpirank_3d(2) .eq.0)then  
     call bd_frey(0,margin,ro,ix,jx,kx)
     call bd_frey(0,margin,pr,ix,jx,kx)
     
     call bd_frey(0,margin,vx,ix,jx,kx)
     call bd_frey(0,margin,vy,ix,jx,kx)
     call bd_frey(0,margin,vz,ix,jx,kx)
     
     call bd_frey(0,margin,phi,ix,jx,kx)
     call bd_frey(0,margin,eta,ix,jx,kx)
     
     call bd_frey(0,margin,bx,ix,jx,kx)
     call bd_frey(0,margin,by,ix,jx,kx)
     call bd_frey(0,margin,bz,ix,jx,kx)
     
  end if
  if(mpid%mpirank_3d(2) .eq. (mpid%mpisize_3d(2)-1))then
     call bd_frey(1,margin,ro,ix,jx,kx)
     call bd_frey(1,margin,pr,ix,jx,kx)
     
     call bd_frey(1,margin,vx,ix,jx,kx)
     call bd_frey(1,margin,vy,ix,jx,kx)
     call bd_frey(1,margin,vz,ix,jx,kx)
     
     call bd_frey(1,margin,phi,ix,jx,kx)
     call bd_frey(1,margin,eta,ix,jx,kx)
     
     call bd_frey(1,margin,bx,ix,jx,kx)
     call bd_frey(1,margin,by,ix,jx,kx)
     call bd_frey(1,margin,bz,ix,jx,kx)
  end if
!======================================================================
! inner x-boundary

  if(mpid%mpirank_3d(1) .eq.0)then

     call bd_frex(0,margin,ro,ix,jx,kx)
     call bd_frex(0,margin,pr,ix,jx,kx)
     call bd_frex(0,margin,vx,ix,jx,kx)
     call bd_frex(0,margin,vy,ix,jx,kx)
     call bd_frex(0,margin,vz,ix,jx,kx)
     
     call bd_frex(0,margin,bx,ix,jx,kx)
     call bd_frex(0,margin,by,ix,jx,kx)
     call bd_frex(0,margin,bz,ix,jx,kx)
     
     call bd_frex(0,margin,phi,ix,jx,kx)
     call bd_frex(0,margin,eta,ix,jx,kx)
  end if
!======================================================================
! outer x-boudary

  if(mpid%mpirank_3d(1) .eq. (mpid%mpisize_3d(1)-1))then

     call bd_frex(1,margin,ro,ix,jx,kx)
     call bd_frex(1,margin,pr,ix,jx,kx)
     
     call bd_frex(1,margin,vx,ix,jx,kx)
     call bd_frex(1,margin,vy,ix,jx,kx)
     call bd_frex(1,margin,vz,ix,jx,kx)
     
     call bd_frex(1,margin,bx,ix,jx,kx)
     call bd_frex(1,margin,by,ix,jx,kx)
     call bd_frex(1,margin,bz,ix,jx,kx)

!!$  call bd_openbx_cyl(1,margin,bxs,bzs,bx,by,bz &
!!$       ,ix,jx,kx,x,dx,dz)

     call bd_frex(1,margin,phi,ix,jx,kx)
     call bd_frex(1,margin,eta,ix,jx,kx)
  end if
!----------------------------------------------------------------------|
! outer z-boundary
!
  if(mpid%mpirank_3d(3) .eq. (mpid%mpisize_3d(3)-1)) then
  
     call bd_frez(1,margin,ro,ix,jx,kx)
     call bd_frez(1,margin,pr,ix,jx,kx)
     
     call bd_frez(1,margin,vx,ix,jx,kx)
     call bd_frez(1,margin,vy,ix,jx,kx)
     call bd_frez(1,margin,vz,ix,jx,kx)
     
     call bd_frez(1,margin,bx,ix,jx,kx)
     call bd_frez(1,margin,by,ix,jx,kx)
     call bd_frez(1,margin,bz,ix,jx,kx)

!!$  call bd_openbz_cyl(1,margin,bxs,bzs,bx,by,bz &
!!$       ,ix,jx,kx,x,dx,dz)

     call bd_frez(1,margin,phi,ix,jx,kx)
     call bd_frez(1,margin,eta,ix,jx,kx)
  end if

!----------------------------------------------------------------------|
! inner z-boundary
!

  if(mpid%mpirank_3d(3) .eq. 0)then

!     call bd_frez(0,margin,ro,ix,jx,kx)
!     call bd_frez(0,margin,pr,ix,jx,kx)
    
!     call bd_frez(0,margin,vx,ix,jx,kx)
!     call bd_frez(0,margin,vy,ix,jx,kx)
!     call bd_frez(0,margin,vz,ix,jx,kx)
     
!     call bd_frez(0,margin,bx,ix,jx,kx)
!     call bd_frez(0,margin,by,ix,jx,kx)
!     call bd_frez(0,margin,bz,ix,jx,kx)

!!$  call bd_openbz_cyl(0,margin,bxs,bzs,bx,by,bz &
!!$       ,ix,jx,kx,x,dx,dz)

!     call bd_frez(0,margin,phi,ix,jx,kx)
!     call bd_frez(0,margin,eta,ix,jx,kx)

     call bd_synpz_car(0,margin,ro,ix,jx,kx)
     call bd_synpz_car(0,margin,pr,ix,jx,kx)
!     call bd_synpz(0,margin,vx,ix,jx,kx)
!     call bd_synpz(0,margin,vy,ix,jx,kx)
     call bd_synnz_car(0,margin,vz,ix,jx,kx)

!     call bd_synnz(0,margin,bx,ix,jx,kx)
!     call bd_synnz(0,margin,by,ix,jx,kx)
     call bd_synnz_car(0,margin,bz,ix,jx,kx)
     call bd_synpz_car(0,margin,phi,ix,jx,kx)
     call bd_synpz_car(0,margin,eta,ix,jx,kx)

     do k=1,margin
        do j=1,jx
           do i=1,ix
              si=y(j)/max(sqrt(x(i)**2+y(j)**2),1d-10)
              co=x(i)/max(sqrt(x(i)**2+y(j)**2),1d-10)
              br=bx(i,j,2*margin-k+1)*co+by(i,j,2*margin-k+1)*si
              bp=-bx(i,j,2*margin-k+1)*si+by(i,j,2*margin-k+1)*co
              vr=vx(i,j,2*margin-k+1)*co+vy(i,j,2*margin-k+1)*si
              vp=-vx(i,j,2*margin-k+1)*si+vy(i,j,2*margin-k+1)*co
              bx(i,j,k)=br*co+bp*si
              by(i,j,k)=br*si-bp*co
              vx(i,j,k)=vr*co+vp*si
              vy(i,j,k)=vr*si-vp*co
              
!           if (vz(i,k) > 0.0d0)then
!              vz(i,k) = 0.0d0
!           endif
!  vz(i,j,k)=min(vz(i,j,k),0d0)
              r=max(sqrt(x(i)**2+y(j)**2),1d-10)
              if (r.le.3d18) then
!              f = 0.5d0*(1.0d0+tanh((1.0d0-r*r/9d36)/(0.21d0*0.21d0+0.42d0)))
                 f = 1.0d0
                 ro(i,j,k)=0.1d0*1.6d-24*0.1478709138d0*(10d0-9d0*f)
                 pr(i,j,k)=6.881726147d0*1.38d-16*2d2
                 vx(i,j,k)=0d0
                 vy(i,j,k)=0d0
!              if(time.le.10.0d0)then
                 vz(i,j,k)=6.0d0*sqrt(gm*6.881726147d0*1.38d-16*2d2/(0.1d0*1.6d-24*0.1478709138d0))*f
!              else
!                 vz(i,k)=jt_bnd%vz*exp(-time+10.0d0)
!              endif
!           r=sqrt(x(i)*x(i) + z(k)*z(k)+ep**2)
                 bx(i,j,k)=-sqrt(2d0*pr(i,j,k)/1000d0)*y(j)/r*sin(r*pi/3d18)**4
                 by(i,j,k)=sqrt(2d0*pr(i,j,k)/1000d0)*x(i)/r*sin(r*pi/3d18)**4
                 bz(i,j,k)=1d-30
!---
! bphi
!!$              by(i,k) = jt_bnd%b0*sqrt(0.5d0*jt_bnd%aa*jt_bnd%dd*x(i)**2 &
!!$                   /(x(i)+0.5d0*jt_bnd%dd)**3)
!!$              bz(i,k) = jt_bnd%b0*sqrt(1.0d0-jt_bnd%aa*(x(i)**2) &
!!$                   *(x(i)+jt_bnd%dd)/(x(i)+0.5d0*jt_bnd%dd)**3)
              endif
           enddo
        enddo
     enddo
  end if
  return
end subroutine bnd_asahina
