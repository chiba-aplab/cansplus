module getNewdt

  use mpi_setup
  implicit none
  private

  public :: getNewdt__glm, getNewdt__glmcyl


contains


  subroutine getNewdt__glm(margin,safety,dtmin,ix,jx,kx,gm,ro,pr &
                          ,vx,vy,vz,bx,by,bz,x,dx,y,dy,z,dz,eta &
                          ,dt,ch,min_dx)

  integer,intent(in) :: ix,jx,kx,margin
  real(8),intent(in) :: safety ! CFL number
  real(8),intent(in) :: dtmin ! minimum dt
  real(8),intent(in) :: min_dx
  real(8),intent(in) :: gm ! specific heat retio
  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: y,dy
  real(8),dimension(kx),intent(in) :: z,dz  
  real(8),dimension(ix,jx,kx),intent(in) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(in) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(in) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(in) :: eta

  real(8), intent(inout) :: dt
! GLM-MHD's wave velocity
! this wave velocity shuld be maximum velocity in system.
  real(8), intent(out) :: ch 
  
  integer :: imin,jmin,kmin
  integer :: i,j,k
  real(8) :: roinverse,vsq,cssq,bsq
  real(8) :: cfxsq,cfysq,cfzsq
  real(8) :: b1,b2,b3,dtmaxi
  real(8) :: v1,v2,v3
  real(8) :: temp, temp2
  real(8) :: temp_sum,temp_dif
  real(8) :: beforedt,dtg

! calculate time step dt in resistive glmmhd @ cartesian grid
  dtmaxi = 0.0d0
  beforedt = dt

  do k=1,kx
     do j=1,jx
        do i=1,ix
           roinverse = 1.0d0/ro(i,j,k)
           vsq = vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2
           
           v1 = vx(i,j,k)
           v2 = vy(i,j,k)
           v3 = vz(i,j,k)
           
           b1 = dabs(bx(i,j,k))
           b2 = dabs(by(i,j,k))
           b3 = dabs(bz(i,j,k))
           bsq = b1**2 + b2**2 + b3**2
           
           cssq = gm*pr(i,j,k)*roinverse
           
           temp_sum = bsq*roinverse + cssq
           temp_dif = bsq*roinverse - cssq

! calculate wave velocity in each direction
           cfxsq = 0.5d0*(temp_sum + sqrt(temp_dif*temp_dif &
                + 4.0d0*cssq*(b2*b2 + b3*b3)*roinverse))
           cfysq = 0.5d0*(temp_sum + sqrt(temp_dif*temp_dif &
                + 4.0d0*cssq*(b3*b3 + b1*b1)*roinverse))
           cfzsq = 0.5d0*(temp_sum + sqrt(temp_dif*temp_dif &
                + 4.0d0*cssq*(b1*b1 + b2*b2)*roinverse))
           
! absolute wave velocity
           temp = max((dabs(v1)+sqrt(cfxsq))/dx(i) &
                , (dabs(v2)+sqrt(cfysq))/dy(j) &
                , (dabs(v3)+sqrt(cfzsq))/dz(k))
           
! diffusion velocity
           temp2 = max(eta(i,j,k)/dx(i)**2 &
                , eta(i,j,k)/dy(j)**2 &
                , eta(i,j,k)/dz(k)**2)

           temp = max(temp,temp2)
           if ( temp > dtmaxi) then
              dtmaxi = temp
              imin = i
              jmin = j
              kmin = k
           endif
        end do
     end do
  end do

  dtg = min(2.0d0*beforedt,safety/dtmaxi)
  call mpi_allreduce(dtg,dt,1,mdp,mmin,mcomw,merr)
  ch = safety*min_dx/dt

! Exception print

  if (dt < dtmin) then
     write(6,*) '  ### stop due to small dt, less than dtmin ###'
     write(6,*) dt,dtmin,imin,jmin,kmin
620  format('   dt = ',1pe10.3,'  < ',1pe10.3,' @ i =',i5 &
          ,'j = ',i5 ,' k = ',i5)
     write(6,*) 'x :: ',x(imin),'y :: ',y(jmin),'z :: ',z(kmin)
     write(6,*) 'ro :: ',ro(imin,jmin,kmin),'pr :: ',pr(imin,jmin,kmin)
     write(6,*) 'vx :: ',vx(imin,jmin,kmin),'vy :: ',vy(imin,jmin,kmin)
     write(6,*) 'vz :: ',vz(imin,jmin,kmin)
     write(6,*) 'bxc:: ',bx(imin,jmin,kmin),'byc:: ',by(imin,jmin,kmin)
     write(6,*) 'bzc:: ',bz(imin,jmin,kmin)
     write(6,*) 'eta:: ',eta(imin,jmin,kmin)
  endif

  end subroutine getNewdt__glm


  subroutine getNewdt__glmcyl(margin,safety,dtmin,ix,jx,kx,gm,ro,pr &
                             ,vx,vy,vz,bx,by,bz,x,dx,y,dy,z,dz,eta &
                             ,dt,ch,min_dx)
 
  integer,intent(in) :: ix,jx,kx,margin
  real(8),intent(in) :: safety ! CFL number
  real(8),intent(in) :: dtmin ! minimum dt
  real(8),intent(in) :: min_dx
  real(8),intent(in) :: gm ! specific heat retio
  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: y,dy
  real(8),dimension(kx),intent(in) :: z,dz  
  real(8),dimension(ix,jx,kx),intent(in) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(in) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(in) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(in) :: eta

  real(8),intent(inout) :: dt
! GLM-MHD's wave velocity
! this wave velocity shuld be maximum velocity in system.
  real(8),intent(out) :: ch 
  
  integer :: imin,jmin,kmin
  integer :: i,j,k
  real(8) :: roinverse,vsq,cssq,bsq
  real(8) :: cfxsq,cfysq,cfzsq
  real(8) :: b1,b2,b3,dtmaxi
  real(8) :: v1,v2,v3
  real(8) :: temp, temp2
  real(8) :: temp_sum,temp_dif
  real(8) :: beforedt,dtg

! calculate time step dt in resistive glmmhd @ cylindrical coordinate
  dtmaxi = 0.0d0
  beforedt = dt

  do k=1,kx
     do j=1,jx
        do i=1,ix
           roinverse = 1.0d0/ro(i,j,k)
           vsq = vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2
           
           v1 = vx(i,j,k)
           v2 = vy(i,j,k)
           v3 = vz(i,j,k)
           
           b1 = dabs(bx(i,j,k))
           b2 = dabs(by(i,j,k))
           b3 = dabs(bz(i,j,k))
           bsq = b1**2 + b2**2 + b3**2
           
           cssq = gm*pr(i,j,k)*roinverse
           
           temp_sum = bsq*roinverse + cssq
           temp_dif = bsq*roinverse - cssq

! calculate wave velocity in each direction
           cfxsq = 0.5d0*(temp_sum + sqrt(temp_dif*temp_dif &
                + 4.0d0*cssq*(b2*b2 + b3*b3)*roinverse))
           cfysq = 0.5d0*(temp_sum + sqrt(temp_dif*temp_dif &
                + 4.0d0*cssq*(b3*b3 + b1*b1)*roinverse))
           cfzsq = 0.5d0*(temp_sum + sqrt(temp_dif*temp_dif &
                + 4.0d0*cssq*(b1*b1 + b2*b2)*roinverse))
           
! absolute wave velocity
           temp = max((dabs(v1)+sqrt(cfxsq))/dx(i) &
                ,(dabs(v2)+sqrt(cfysq))/(x(i)*dy(j)) &
                ,(dabs(v3)+sqrt(cfzsq))/dz(k))
           
! diffusion velocity
           temp2 = max(eta(i,j,k)/dx(i)**2 &
                ,eta(i,j,k)/(x(i)*dy(j))**2 &
                ,eta(i,j,k)/dz(k)**2)

           temp = max(temp,temp2)
           if ( temp > dtmaxi) then
              dtmaxi = temp
              imin = i
              jmin = j
              kmin = k
           endif
        end do
     end do
  end do

  dtg = min(2.0d0*beforedt,safety/dtmaxi)
  call mpi_allreduce(dtg,dt,1,mdp,mmin,mcomw,merr)
  ch = safety*min_dx/dt

! Exception print

  if (dt < dtmin) then
     write(6,*) '  ### stop due to small dt, less than dtmin ###'
     write(6,*) dt,dtmin,imin,jmin,kmin
620  format('   dt = ',1pe10.3,'  < ',1pe10.3,' @ i =',i5 &
          ,'j = ',i5 ,' k = ',i5)
     write(6,*) 'x :: ',x(imin),'y :: ',y(jmin),'z :: ',z(kmin)
     write(6,*) 'ro :: ',ro(imin,jmin,kmin),'pr :: ',pr(imin,jmin,kmin)
     write(6,*) 'vx :: ',vx(imin,jmin,kmin),'vy :: ',vy(imin,jmin,kmin)
     write(6,*) 'vz :: ',vz(imin,jmin,kmin)
     write(6,*) 'bxc:: ',bx(imin,jmin,kmin),'byc:: ',by(imin,jmin,kmin)
     write(6,*) 'bzc:: ',bz(imin,jmin,kmin)
     write(6,*) 'eta:: ',eta(imin,jmin,kmin)
  endif

  end subroutine getNewdt__glmcyl
 

end module getNewdt
