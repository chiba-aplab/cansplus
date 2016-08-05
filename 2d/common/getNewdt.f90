module getNewdt

  use mpi_setup
  implicit none
  private

  public :: getNewdt__glm

! for GLM-MHD
  real(8), parameter :: cr=0.18d0


contains


  subroutine getNewdt__glm(margin,safety,dtmin,ix,jx,gm,ro,pr &
                          ,vx,vy,vz,bx,by,bz,x,dx,y,dy,eta &
                          ,dt,ch,cp,min_dx)

  integer,intent(in) :: ix,jx,margin
  real(8),intent(in) :: safety
  real(8),intent(in) :: dtmin
  real(8),intent(in) :: min_dx
  real(8),intent(in) :: gm
  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: y,dy
  real(8),dimension(ix,jx),intent(in) :: ro,pr
  real(8),dimension(ix,jx),intent(in) :: vx,vy,vz
  real(8),dimension(ix,jx),intent(in) :: bx,by,bz
  real(8),dimension(ix,jx),intent(in) :: eta
  real(8), intent(inout) :: dt
! GLM-MHD's wave velocity
! this wave velocity shuld be maximum velocity in system.
  real(8), intent(out) :: ch, cp
  
  integer :: imin,jmin
  integer :: i,j
  real(8) :: roinverse,vsq,cssq,bsq
  real(8) :: cfxsq,cfysq
  real(8) :: b1,b2,b3,dtmaxi
  real(8) :: v1,v2,v3
  real(8) :: temp, temp2
  real(8) :: temp_sum,temp_dif
  real(8) :: beforedt,dtg

! calculate time step dt in resistive glmmhd @ cartesian grid
  dtmaxi = 0.0d0
  beforedt = dt

  do j=1,jx
     do i=1,ix
        roinverse = 1.0d0/ro(i,j)
        vsq = vx(i,j)**2+vy(i,j)**2+vz(i,j)**2
           
        v1 = vx(i,j)
        v2 = vy(i,j)
        v3 = vz(i,j)
           
        b1 = dabs(bx(i,j))
        b2 = dabs(by(i,j))
        b3 = dabs(bz(i,j))
        bsq = b1**2 + b2**2 + b3**2
           
        cssq = gm*pr(i,j)*roinverse
           
        temp_sum = bsq*roinverse + cssq
        temp_dif = bsq*roinverse - cssq

! calculate wave velocity in each direction
        cfxsq = 0.5d0*(temp_sum + sqrt(temp_dif*temp_dif &
             + 4.0d0*cssq*(b2*b2 + b3*b3)*roinverse))
        cfysq = 0.5d0*(temp_sum + sqrt(temp_dif*temp_dif &
             + 4.0d0*cssq*(b3*b3 + b1*b1)*roinverse))
           
! absolute wave velocity
        temp = max( (dabs(v1)+sqrt(cfxsq))/dx(i) &
                   ,(dabs(v2)+sqrt(cfysq))/dy(j) )

! diffusion velocity
        temp2 = max( eta(i,j)/dx(i)**2 &
                    ,eta(i,j)/dy(j)**2)

        temp = max(temp,temp2)
        if ( temp > dtmaxi) then
           dtmaxi = temp
           imin = i
           jmin = j
        endif
     end do
  end do

  dtg = min(2.0d0*beforedt,safety/dtmaxi)
  call mpi_allreduce(dtg,dt,1,mdp,mmin,mcomw,merr)
  ch = safety*min_dx/dt
  cp = sqrt(ch*cr)

! Exception print

  if (dt < dtmin) then
     write(6,*) '  ### stop due to small dt, less than dtmin ###'
     write(6,*) dt,dtmin,imin,jmin
620  format('   dt = ',1pe10.3,'  < ',1pe10.3,' @ i =',i5 &
          ,'j = ',i5)
     write(6,*) 'x :: ',x(imin),'y :: ',y(jmin)
     write(6,*) 'ro :: ',ro(imin,jmin),'pr :: ',pr(imin,jmin)
     write(6,*) 'vx :: ',vx(imin,jmin),'vy :: ',vy(imin,jmin)
     write(6,*) 'vz :: ',vz(imin,jmin)
     write(6,*) 'bxc:: ',bx(imin,jmin),'byc:: ',by(imin,jmin)
     write(6,*) 'bzc:: ',bz(imin,jmin)
     write(6,*) 'eta:: ',eta(imin,jmin)
     call MPI_ABORT(mcomw, 9, merr)
     call MPI_FINALIZE(merr)
  endif

  end subroutine getNewdt__glm
 

end module getNewdt
