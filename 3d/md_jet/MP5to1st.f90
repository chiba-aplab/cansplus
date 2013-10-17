subroutine MP5to1st(mdir,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi &
     ,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw)
  implicit none

  integer,intent(in) :: ix,jx,kx,mdir
  real(8),dimension(ix,jx,kx),intent(in) :: ro,pr,vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(in) :: bx,by,bz,phi
  real(8),dimension(ix,jx,kx,2),intent(inout) :: row,prw,vxw,vyw,vzw
  real(8),dimension(ix,jx,kx,2),intent(inout) :: bxw,byw,bzw,phiw

  integer :: i,j,k
  real(8) :: minvalue

  if (mdir == 1)then
     do k=3,kx-2
        do j=3,jx-2
           do i=3,ix-2
              minvalue = min(row(i-1,j,k,2),row(i,j,k,1),prw(i-1,j,k,2),prw(i,j,k,1))
              if(minvalue <= 0d0)then
                 row(i-1,j,k,2) = ro(i,j,k)
                 row(i  ,j,k,1) = ro(i,j,k)
                 prw(i-1,j,k,2) = pr(i,j,k)
                 prw(i  ,j,k,1) = pr(i,j,k)
                 vxw(i-1,j,k,2) = vx(i,j,k)
                 vxw(i  ,j,k,1) = vx(i,j,k)
                 vyw(i-1,j,k,2) = vy(i,j,k)
                 vyw(i  ,j,k,1) = vy(i,j,k)
                 vzw(i-1,j,k,2) = vz(i,j,k)
                 vzw(i  ,j,k,1) = vz(i,j,k)
                 bxw(i-1,j,k,2) = bx(i,j,k)
                 bxw(i  ,j,k,1) = bx(i,j,k)
                 byw(i-1,j,k,2) = by(i,j,k)
                 byw(i  ,j,k,1) = by(i,j,k)
                 bzw(i-1,j,k,2) = bz(i,j,k)
                 bzw(i  ,j,k,1) = bz(i,j,k)
                 phiw(i-1,j,k,2) = phi(i,j,k)
                 phiw(i  ,j,k,1) = phi(i,j,k)
              endif
           enddo
        enddo
     enddo
  else if(mdir == 2)then
     do k=3,kx-2
        do j=3,jx-2
           do i=3,ix-2
              minvalue = min(row(i,j-1,k,2),row(i,j,k,1),prw(i,j-1,k,2),prw(i,j,k,1))
              if(minvalue <= 0d0)then
                 row(i,j-1,k,2) = ro(i,j,k)
                 row(i,j  ,k,1) = ro(i,j,k)
                 prw(i,j-1,k,2) = pr(i,j,k)
                 prw(i,j  ,k,1) = pr(i,j,k)
                 vxw(i,j-1,k,2) = vx(i,j,k)
                 vxw(i,j  ,k,1) = vx(i,j,k)
                 vyw(i,j-1,k,2) = vy(i,j,k)
                 vyw(i,j  ,k,1) = vy(i,j,k)
                 vzw(i,j-1,k,2) = vz(i,j,k)
                 vzw(i,j  ,k,1) = vz(i,j,k)
                 bxw(i,j-1,k,2) = bx(i,j,k)
                 bxw(i,j  ,k,1) = bx(i,j,k)
                 byw(i,j-1,k,2) = by(i,j,k)
                 byw(i,j  ,k,1) = by(i,j,k)
                 bzw(i,j-1,k,2) = bz(i,j,k)
                 bzw(i,j  ,k,1) = bz(i,j,k)
                 phiw(i,j-1,k,2) = phi(i,j,k)
                 phiw(i,j  ,k,1) = phi(i,j,k)
              endif
           enddo
        enddo
     enddo
  else
     do k=3,kx-2
        do j=3,jx-2
           do i=3,ix-2
              minvalue = min(row(i,j,k-1,2),row(i,j,k,1),prw(i,j,k-1,2),prw(i,j,k,1))
              if(minvalue <= 0d0)then
                 row(i,j,k-1,2) = ro(i,j,k)
                 row(i,j,k  ,1) = ro(i,j,k)
                 prw(i,j,k-1,2) = pr(i,j,k)
                 prw(i,j,k  ,1) = pr(i,j,k)
                 vxw(i,j,k-1,2) = vx(i,j,k)
                 vxw(i,j,k  ,1) = vx(i,j,k)
                 vyw(i,j,k-1,2) = vy(i,j,k)
                 vyw(i,j,k  ,1) = vy(i,j,k)
                 vzw(i,j,k-1,2) = vz(i,j,k)
                 vzw(i,j,k  ,1) = vz(i,j,k)
                 bxw(i,j,k-1,2) = bx(i,j,k)
                 bxw(i,j,k  ,1) = bx(i,j,k)
                 byw(i,j,k-1,2) = by(i,j,k)
                 byw(i,j,k  ,1) = by(i,j,k)
                 bzw(i,j,k-1,2) = bz(i,j,k)
                 bzw(i,j,k  ,1) = bz(i,j,k)
                 phiw(i,j,k-1,2) = phi(i,j,k)
                 phiw(i,j,k  ,1) = phi(i,j,k)
              endif
           enddo
        enddo
     enddo
  endif

  return
end subroutine MP5to1st
