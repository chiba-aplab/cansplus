subroutine MP5to1st(mdir,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi &
     ,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw)
  implicit none

  integer,intent(in) :: ix,jx,kx,mdir
  real(8),dimension(ix,jx,kx),intent(in) :: ro,pr,vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(in) :: bx,by,bz,phi
  real(8),dimension(ix,jx,kx,2),intent(inout) :: row,prw,vxw,vyw,vzw
  real(8),dimension(ix,jx,kx,2),intent(inout) :: bxw,byw,bzw,phiw

  integer :: i,j,k
  real(8) :: minvalue,smv,msmv,psmv

  if (mdir == 1)then
     do k=3,kx-2
        do j=3,jx-2
           do i=3,ix-2
              minvalue = min(row(i-1,j,k,2),row(i,j,k,1),prw(i-1,j,k,2),prw(i,j,k,1))
              smv = sign(1d0,minvalue)
              psmv = max(0d0,smv)
              msmv = max(0d0,-smv)

              row(i-1,j,k,2) = row(i-1,j,k,2)*psmv + ro(i,j,k)*msmv
              row(i  ,j,k,1) = row(i  ,j,k,1)*psmv + ro(i,j,k)*msmv
              prw(i-1,j,k,2) = prw(i-1,j,k,2)*psmv + pr(i,j,k)*msmv
              prw(i  ,j,k,1) = prw(i  ,j,k,1)*psmv + pr(i,j,k)*msmv
              vxw(i-1,j,k,2) = vxw(i-1,j,k,2)*psmv + vx(i,j,k)*msmv
              vxw(i  ,j,k,1) = vxw(i  ,j,k,1)*psmv + vx(i,j,k)*msmv
              vyw(i-1,j,k,2) = vyw(i-1,j,k,2)*psmv + vy(i,j,k)*msmv
              vyw(i  ,j,k,1) = vyw(i  ,j,k,1)*psmv + vy(i,j,k)*msmv
              vzw(i-1,j,k,2) = vzw(i-1,j,k,2)*psmv + vz(i,j,k)*msmv
              vzw(i  ,j,k,1) = vzw(i  ,j,k,1)*psmv + vz(i,j,k)*msmv
              bxw(i-1,j,k,2) = bxw(i-1,j,k,2)*psmv + bx(i,j,k)*msmv
              bxw(i  ,j,k,1) = bxw(i  ,j,k,1)*psmv + bx(i,j,k)*msmv
              byw(i-1,j,k,2) = byw(i-1,j,k,2)*psmv + by(i,j,k)*msmv
              byw(i  ,j,k,1) = byw(i  ,j,k,1)*psmv + by(i,j,k)*msmv
              bzw(i-1,j,k,2) = bzw(i-1,j,k,2)*psmv + bz(i,j,k)*msmv
              bzw(i  ,j,k,1) = bzw(i  ,j,k,1)*psmv + bz(i,j,k)*msmv
              phiw(i-1,j,k,2) = phiw(i-1,j,k,2)*psmv + phi(i,j,k)*msmv
              phiw(i  ,j,k,1) = phiw(i  ,j,k,1)*psmv + phi(i,j,k)*msmv
           enddo
        enddo
     enddo
  else if(mdir == 2)then
     do k=3,kx-2
        do j=3,jx-2
           do i=3,ix-2
              minvalue = min(row(i,j-1,k,2),row(i,j,k,1),prw(i,j-1,k,2),prw(i,j,k,1))
              smv = sign(1d0,minvalue)
              psmv = max(0d0,smv)
              msmv = max(0d0,-smv)

              row(i,j-1,k,2) = row(i,j-1,k,2)*psmv + ro(i,j,k)*msmv
              row(i,j  ,k,1) = row(i,j  ,k,1)*psmv + ro(i,j,k)*msmv
              prw(i,j-1,k,2) = prw(i,j-1,k,2)*psmv + pr(i,j,k)*msmv
              prw(i,j  ,k,1) = prw(i,j  ,k,1)*psmv + pr(i,j,k)*msmv
              vxw(i,j-1,k,2) = vxw(i,j-1,k,2)*psmv + vx(i,j,k)*msmv
              vxw(i,j  ,k,1) = vxw(i,j  ,k,1)*psmv + vx(i,j,k)*msmv
              vyw(i,j-1,k,2) = vyw(i,j-1,k,2)*psmv + vy(i,j,k)*msmv
              vyw(i,j  ,k,1) = vyw(i,j  ,k,1)*psmv + vy(i,j,k)*msmv
              vzw(i,j-1,k,2) = vzw(i,j-1,k,2)*psmv + vz(i,j,k)*msmv
              vzw(i,j  ,k,1) = vzw(i,j  ,k,1)*psmv + vz(i,j,k)*msmv
              bxw(i,j-1,k,2) = bxw(i,j-1,k,2)*psmv + bx(i,j,k)*msmv
              bxw(i,j  ,k,1) = bxw(i,j  ,k,1)*psmv + bx(i,j,k)*msmv
              byw(i,j-1,k,2) = byw(i,j-1,k,2)*psmv + by(i,j,k)*msmv
              byw(i,j  ,k,1) = byw(i,j  ,k,1)*psmv + by(i,j,k)*msmv
              bzw(i,j-1,k,2) = bzw(i,j-1,k,2)*psmv + bz(i,j,k)*msmv
              bzw(i,j  ,k,1) = bzw(i,j  ,k,1)*psmv + bz(i,j,k)*msmv
              phiw(i,j-1,k,2) = phiw(i,j-1,k,2)*psmv + phi(i,j,k)*msmv
              phiw(i,j  ,k,1) = phiw(i,j  ,k,1)*psmv + phi(i,j,k)*msmv
           enddo
        enddo
     enddo
  else
     do k=3,kx-2
        do j=3,jx-2
           do i=3,ix-2
              minvalue = min(row(i,j,k-1,2),row(i,j,k,1),prw(i,j,k-1,2),prw(i,j,k,1))
              smv = sign(1d0,minvalue)
              psmv = max(0d0,smv)
              msmv = max(0d0,-smv)

              row(i,j,k-1,2) = row(i,j,k-1,2)*psmv + ro(i,j,k)*msmv
              row(i,j,k  ,1) = row(i,j,k  ,1)*psmv + ro(i,j,k)*msmv
              prw(i,j,k-1,2) = prw(i,j,k-1,2)*psmv + pr(i,j,k)*msmv
              prw(i,j,k  ,1) = prw(i,j,k  ,1)*psmv + pr(i,j,k)*msmv
              vxw(i,j,k-1,2) = vxw(i,j,k-1,2)*psmv + vx(i,j,k)*msmv
              vxw(i,j,k  ,1) = vxw(i,j,k  ,1)*psmv + vx(i,j,k)*msmv
              vyw(i,j,k-1,2) = vyw(i,j,k-1,2)*psmv + vy(i,j,k)*msmv
              vyw(i,j,k  ,1) = vyw(i,j,k  ,1)*psmv + vy(i,j,k)*msmv
              vzw(i,j,k-1,2) = vzw(i,j,k-1,2)*psmv + vz(i,j,k)*msmv
              vzw(i,j,k  ,1) = vzw(i,j,k  ,1)*psmv + vz(i,j,k)*msmv
              bxw(i,j,k-1,2) = bxw(i,j,k-1,2)*psmv + bx(i,j,k)*msmv
              bxw(i,j,k  ,1) = bxw(i,j,k  ,1)*psmv + bx(i,j,k)*msmv
              byw(i,j,k-1,2) = byw(i,j,k-1,2)*psmv + by(i,j,k)*msmv
              byw(i,j,k  ,1) = byw(i,j,k  ,1)*psmv + by(i,j,k)*msmv
              bzw(i,j,k-1,2) = bzw(i,j,k-1,2)*psmv + bz(i,j,k)*msmv
              bzw(i,j,k  ,1) = bzw(i,j,k  ,1)*psmv + bz(i,j,k)*msmv
              phiw(i,j,k-1,2) = phiw(i,j,k-1,2)*psmv + phi(i,j,k)*msmv
              phiw(i,j,k  ,1) = phiw(i,j,k  ,1)*psmv + phi(i,j,k)*msmv
           enddo
        enddo
     enddo
  endif

  return
end subroutine MP5to1st
