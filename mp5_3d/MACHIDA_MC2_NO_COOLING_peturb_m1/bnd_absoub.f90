subroutine bnd_absoub(ix,jx,kx,x,z,xin,roi,pri,vxi,vyi,vzi,bxi,byi,bzi &
     ,ro,pr,vx,vy,vz,bx,by,bz)

  implicit none

  integer,intent(in) :: ix,jx,kx

  real(8),intent(in) :: xin

  real(8),dimension(ix),intent(in) :: x
  real(8),dimension(kx),intent(in) :: z

  real(8),dimension(ix,jx,kx),intent(in) :: roi,pri
  real(8),dimension(ix,jx,kx),intent(in) :: vxi,vyi,vzi
  real(8),dimension(ix,jx,kx),intent(in) :: bxi,byi,bzi

  real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz

!  call outgoing_vv(ix,jx,kx,x,z,vx,vz,xin)

  call absorb(ix,jx,kx,x,z,xin,ro,roi)
  call absorb(ix,jx,kx,x,z,xin,pr,pri)

  call absorb(ix,jx,kx,x,z,xin,vx,vxi)
  call absorb(ix,jx,kx,x,z,xin,vy,vyi)
  call absorb(ix,jx,kx,x,z,xin,vz,vzi)

  call absorb(ix,jx,kx,x,z,xin,bx,bxi)
  call absorb(ix,jx,kx,x,z,xin,by,byi)
  call absorb(ix,jx,kx,x,z,xin,bz,bzi)

end subroutine bnd_absoub
