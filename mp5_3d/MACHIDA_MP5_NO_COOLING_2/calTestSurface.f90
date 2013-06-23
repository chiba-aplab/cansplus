subroutine calTestSurface(mdir,mpid,margin,ix,jx,kx,x,dx,y,dy,z,dz,qq,qqval)
  use mpi_domain_xz
  implicit none

  type(mpidomain),intent(in) :: mpid

  integer,intent(in) :: mdir
  integer,intent(in) :: margin,ix,jx,kx

  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: y,dy
  real(8),dimension(kx),intent(in) :: z,dz

  real(8),dimension(ix,jx,kx),intent(in) :: qq

  real(8),intent(inout) :: qqval

  real(8) :: ds,rr
 
  integer :: i,j,k

  qqval = 0.0d0
  if(mdir == 0)then
     if(mpid%mpirank_2d(1) .eq. (mpid%mpisize_2d(1)-1))then
        i = ix-margin
        rr = sqrt(x(i)**2)
        do k=margin+1,kx-margin
           do j=margin+1,jx-margin
              ds = rr*dy(j)*dz(k)
              qqval = qqval + qq(i,j,k)*ds
           end do
        end do
     endif
  else
     if(mpid%mpirank_2d(2) .eq. (mpid%mpisize_2d(2)-1)) then  
        do j=margin+1,jx-margin
           do i=margin+1, ix-margin
              rr = sqrt(x(i)**2)
              ds = rr*dx(i)*dy(j)
              
              k= margin+1
              qqval = qqval + qq(i,j,k)*ds
              
              k= kx-margin
              qqval = qqval + qq(i,j,k)*ds
           end do
        end do
     end if
  end if
  return 
end subroutine calTestSurface
