!======================================================================|
subroutine perturb_vector &
     & (perturb_type, i_random_seed &
     & , ix, jx, kx, vx, vy, vz, x, dx, dy, dz, v0)
!======================================================================!
! NAME : perturtb_v
!
! PURPOSE: Add perturbation to variables
!
! INPUTS
!    perturb_type   [integer] : flag of perturbation type
!    i_random_seed  [integer] : seed of random number
!    ix, jx, kx     [integer] : grid numbers
!    x(ix)          [real(8)] : x-coord
!    dx(ix)         [real(8)] : cell center interval in x-coord
!    dy(jx)         [real(8)] : cell center interval in y-coord
!    dz(kx)         [real(8)] : cell center interval in z-coord
!    v0             [real(8)] :
!
! INOUTPUTS
!    vx(ix, jx, kx) [real(8)] : x-component of vector to be perturbed
!    vy(ix, jx, kx) [real(8)] : y-component of vector to be perturbed
!    vz(ix, jx, kx) [real(8)] : z-component of vector to be perturbed
!    
! OUTPUTS
!
! INTERNAL VARIABLES
!
! UPDATE LOG:
! 2012/12/30 initial coding (by H. ODA)
!======================================================================!

! perturbation
!  use mpi_domain_xz
  use cans_type

  implicit none

  integer,intent(in) :: ix,jx,kx
  integer,intent(in) :: perturb_type, i_random_seed

  !type(mpidomain) :: mpid

  real(8),intent(in) :: v0

  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: dy
  real(8),dimension(kx),intent(in) :: dz

  real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz

  real(8),dimension(ix,jx,kx) :: dax,day,daz
  real(8),dimension(ix,jx,kx) :: dvx,dvy,dvz
  real(8),dimension(ix,jx,kx) :: dvxc,dvyc,dvzc

  integer :: n
  integer,dimension(:),allocatable :: seed

  real(8) :: amp
  real(8) :: dtody,dtodz,dtodx
  real(8) :: line1,line2
  real(8) :: temp
  real(8) :: dvxmax,dvymax,dvzmax
  real(8) :: vmax

  integer :: i,j,k
!----------------------------------------------------------------------|

  amp = 0.01d0*v0

  call random_seed(size = n)
  allocate(seed(n))
  !seed = 123+mpid%mpirank
  seed = 123 + i_random_seed
  call random_seed(PUT = seed)
  call random_number(dax)
  call random_number(day)
  call random_number(daz)

  if(perturb_type .eq. 1)then
     do k=1,kx
        do j=1,jx
           do i=1,ix
              dvx(i,j,k) = 0.0d0
              dvy(i,j,k) = 0.0d0
              dvz(i,j,k) = 0.0d0
           enddo
        enddo
     enddo


     dvxmax = 0.0d0
     do k=2,kx-1
        do j=2,jx-1
           do i=1,ix
              if(vy(i,j,k) .ne. 0.0d0)then
                 dtodz = 1.0d0/dz(k)
                 dtody = 1.0d0/(x(i)*dy(j))
                 dvx(i,j,k) = dvx(i,j,k)-dtodz*(day(i,j,k-1)-day(i,j,k))&
                      -dtody*(daz(i,j,k)-daz(i,j-1,k))
                 temp = dabs(dvx(i,j,k))
                 if(temp > dvxmax)then
                    dvxmax = temp
                 endif
              endif
           enddo
        enddo
     enddo

     dvymax = 0.0d0
     do k=2,kx-1
        do j=1,jx
           do i=2,ix-1
              if(vy(i,j,k) .ne. 0.0d0)then
                 dtodz = 1.0d0/dz(k)
                 dtodx = 1.0d0/dx(i)
                 dvy(i,j,k) = dvy(i,j,k)-(dtodz)*(dax(i,j,k)-dax(i,j,k-1)) &
                      -(dtodx)*(-daz(i,j,k)+daz(i-1,j,k))
                 temp = dabs(dvy(i,j,k))
                 if(temp > dvymax)then
                    dvymax = temp
                 endif
              endif
           enddo
        enddo
     enddo

     dvzmax = 0.0d0
     do k=1,kx
        do j=1,jx
           do i=2,ix-1
              if(vy(i,j,k) .ne. 0.0d0)then
                 line1 = dabs(x(i)+0.5d0*dx(i))
                 line2 = dabs(x(i)-0.5d0*dx(i))
                 dtodx = 1.0d0/(x(i)*dx(i))
                 dtody = 1.0d0/(x(i)*dy(j))
                 dvz(i,j,k) = dvz(i,j,k)-(dtodx) &
                      *(line1*day(i,j,k)-line2*day(i-1,j,k)) & 
                      -dtody*(dax(i,j-1,k)-dax(i,j,k))
                 temp = dabs(dvz(i,j,k))
                 if(temp > dvzmax)then
                    dvzmax = temp
                 endif
              endif
           enddo
        enddo
     enddo

     do k=1,kx
        do j=1,jx
           do i=1,ix
              if(dvxmax .ne. 0.0d0)then
                 dvx(i,j,k) = amp*(dvx(i,j,k)/dvxmax)
              else
                 dvx(i,j,k) = 0.0d0
              endif
              if(dvymax .ne. 0.0d0)then
                 dvy(i,j,k) = amp*(dvy(i,j,k)/dvymax)
              else
                 dvy(i,j,k) = 0.0d0
              endif
              if(dvzmax .ne. 0.0d0)then
                 dvz(i,j,k) = amp*(dvz(i,j,k)/dvzmax)
              else
                 dvz(i,j,k) = 0.0d0
              endif
           enddo
        enddo
     enddo

     do k=1,kx-1
        do j=1,jx-1
           do i=1,ix-1
              dvxc(i,j,k) = 0.5d0*(dvx(i,j,k)+dvx(i-1,j,k))
              dvyc(i,j,k) = 0.5d0*(dvy(i,j,k)+dvy(i,j-1,k))
              dvzc(i,j,k) = 0.5d0*(dvz(i,j,k)+dvz(i,j,k-1))
           enddo
        enddo
     enddo

     do k=1,kx-1
        do j=1,jx-1
           do i=1,ix-1
              vx(i,j,k) = vx(i,j,k) + dvxc(i,j,k)
              vy(i,j,k) = vy(i,j,k) + dvyc(i,j,k)
              vz(i,j,k) = vz(i,j,k) + dvzc(i,j,k)
           enddo
        enddo
     enddo
  else if(perturb_type .eq. 2)then
     do k=1,kx
        do j=1,jx
           do i=1,ix
              if(dabs(vy(i,j,k)) > 0.0d0)then
                 vx(i,j,k) = vx(i,j,k) + 2.0d0*amp*(dax(i,j,k)-0.5d0)
                 vy(i,j,k) = vy(i,j,k) + 2.0d0*amp*(day(i,j,k)-0.5d0)
                 vz(i,j,k) = vz(i,j,k) + 2.0d0*amp*(daz(i,j,k)-0.5d0)
              endif
           enddo
        enddo
     enddo
  else
     do k=1,kx
        do j=1,jx
           do i=1,ix
              dvx(i,j,k) = 0.0d0
              dvy(i,j,k) = 0.0d0
              dvz(i,j,k) = 0.0d0
           enddo
        enddo
     enddo
  endif

  return
end subroutine perturb_vector
