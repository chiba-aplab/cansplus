subroutine bnd( )

  use mpi_domain_xz
  use const
  use init
  use boundary
  implicit none

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

     call bd_synpz_car(0,margin,ro,ix,jx,kx)
     call bd_synpz_car(0,margin,pr,ix,jx,kx)
     
     call bd_synpz_car(0,margin,vx,ix,jx,kx)
     call bd_synpz_car(0,margin,vy,ix,jx,kx)
     call bd_synnz_car(0,margin,vz,ix,jx,kx)
     
     call bd_synpz_car(0,margin,bx,ix,jx,kx)
     call bd_synpz_car(0,margin,by,ix,jx,kx)
     call bd_synnz_car(0,margin,bz,ix,jx,kx)

     call bd_synpz_car(0,margin,phi,ix,jx,kx)
     call bd_synpz_car(0,margin,eta,ix,jx,kx)


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
end if
return
end subroutine bnd
