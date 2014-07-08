module model

  use const
  use mpi_setup, only : mpid

  implicit none
  private

  public :: model_setup


contains

  
  subroutine model_setup(ro,pr,vx,vy,vz,bx,by,bz,phi        &
                      ,roi,pri,vxi,vyi,vzi,bxi,byi,bzi,phii &
                      ,x,dx,xm,y,dy,ym,z,dz,zm &
                      ,gx,gz,eta,min_dx)

!---Input & Output
  real(8),dimension(ix)      ,intent(out) :: x,dx
  real(8),dimension(0:ix)    ,intent(out) :: xm  
  real(8),dimension(jx)      ,intent(out) :: y,dy
  real(8),dimension(0:jx)    ,intent(out) :: ym  
  real(8),dimension(kx)      ,intent(out) :: z,dz
  real(8),dimension(0:kx)    ,intent(out) :: zm  
  real(8),dimension(ix,jx,kx),intent(out) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(out) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(out) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(out) :: eta,phi
  real(8),dimension(ix,jx,kx),intent(out) :: gx,gz
  real(8),dimension(ix,jx,kx),intent(out) :: roi,pri
  real(8),dimension(ix,jx,kx),intent(out) :: vxi,vyi,vzi
  real(8),dimension(ix,jx,kx),intent(out) :: bxi,byi,bzi,phii
  real(8),                    intent(out) :: min_dx

  integer :: i,j,k
  integer :: ig,jg,kg
  integer :: izero,jzero,kzero
  real(8) :: af1,af2
  real(8) :: dxmax,dymax,dzmax
  real(8),dimension(ix,jx,kx) :: curx,cury,curz
  real(8),dimension(igx) :: xg,dxg
  real(8),dimension(0:igx) :: xmg
  real(8),dimension(jgx) :: yg,dyg
  real(8),dimension(0:jgx) :: ymg
  real(8),dimension(kgx) :: zg,dzg
  real(8),dimension(0:kgx) :: zmg
  real(8),dimension(igx,jgx,kgx) :: gxg,gzg

  real(8),dimension(kgx) :: gzm0,tem0,rbeta,rbetacr,den0,pre0,bmx0   

!---Step 1a.-------------------------------------------------------------|
! set global x-grid 
!

  do i=1,igx
    dxg(i)=dxg0
  enddo

! X origin
  izero=igx/2-1
  xmg(izero)= -dxg0*0.5d0
  do i=izero,igx-1
    xmg(i+1) = xmg(i)+dxg(i+1)
  enddo
  do i=izero-1,0,-1
    xmg(i) = xmg(i+1)-dxg(i+1)
  enddo
  do i=1,igx
     xg(i) = 0.5d0*(xmg(i)+xmg(i-1))
  enddo
  xg(1) = 0.5d0*(xmg(1) + (xmg(1)-dxg0))

!---Step 1b.-------------------------------------------------------------|
! set global y-grid
  do j=1,jgx
     dyg(j) = dyg0
  enddo

! Y origin
  jzero=jgx/2-1
  ymg(jzero)=-dyg0*0.5d0
  do j=jzero,jgx-1
     ymg(j+1) = ymg(j)+dyg(j+1)
  enddo

  do j=jzero-1,0,-1
     ymg(j) = ymg(j+1)-dyg(j+1)
  enddo

  do j=1,jgx
     yg(j) = 0.5d0*(ymg(j)+ymg(j-1))
  enddo
  yg(1) = 0.5d0*(ymg(1) + (ymg(1)-dyg0))

!---Step 1c.-------------------------------------------------------------|
! set global z-grid
  do k=1,kgx
     dzg(k)=dzg0
  enddo

! Z origin
  kzero=kgx/2    
  zmg(kzero)=0.d0
  do k=kzero,kgx-1
    zmg(k+1) = zmg(k)+dzg(k+1)
  enddo
  do k=kzero-1,0,-1
    zmg(k) = zmg(k+1)-dzg(k+1)
  enddo

  do k=2,kgx
    zg(k) = 0.5d0*(zmg(k)+zmg(k-1))
  enddo

  zg(1) = 0.5d0*(zmg(1) + (zmg(1)-dzg0))

!---Step 2a.-------------------------------------------------------------|
! set individual x-grid 
  do i=1,ix
     ig = mpid%mpirank_3d(1)*(ix-2*margin)+i
     x(i)=xg(ig)
     dx(i) = dxg(ig)
  enddo
  do i=0,ix
     ig = mpid%mpirank_3d(1)*(ix-2*margin)+i
     xm(i) = xmg(ig)
  end do

!---Step 2b.-------------------------------------------------------------|
! set individual y-grid 
  do j=1,jx
     jg = mpid%mpirank_3d(2)*(jx-2*margin)+j
     y(j) = yg(jg)
     dy(j) = dyg(jg)
  enddo
  do j=0,jx
     jg = mpid%mpirank_3d(2)*(jx-2*margin)+j
     ym(j) = ymg(jg)
  end do

!---Step 2c.-------------------------------------------------------------|
! set individual z-grid 
  do k=1,kx
     kg = mpid%mpirank_3d(3)*(kx-2*margin)+k
     z(k) = zg(kg)
     dz(k) = dzg(kg)
  enddo
  do k=0,kx
     kg = mpid%mpirank_3d(3)*(kx-2*margin)+k
     zm(k) = zmg(kg)
  end do

! calculate min_dx
  min_dx = min(minval(dxg),minval(dyg),minval(dzg))

!---Step 3a.----------------------------------------------------------|
! set gravitation
  do k=1,kgx
     do j=1,jgx
        do i=1,igx
        gxg(i,j,k)=0.0d0
        gzg(i,j,k)= - g0 * tanh(zg(k)/wg0)
        end do
     end do
        gzm0(k)= - g0 * tanh((zg(k)+0.5d0*dzg(k))/wg0)
  end do

!--Step 3b.----------------------------------------------------------|
! set individual gravitation
  do k=1,kx
    do j=1,jx
        do i=1,ix
           ig = mpid%mpirank_3d(1)*(ix-2*margin)+i
           jg = mpid%mpirank_3d(2)*(jx-2*margin)+j
           kg = mpid%mpirank_3d(3)*(kx-2*margin)+k
           gx(i,j,k) = gxg(ig,jg,kg)
           gz(i,j,k) = gzg(ig,jg,kg)
        enddo
     enddo
  enddo
 
!-------------------------
!   temperature
  do k=1,kgx
    tem0(k)=1.d0 + &
           0.5d0*(tcor-1.d0)*(tanh((abs(zg(k))-ztr)/wtr)+1.d0)
  enddo

!-------------------------
!   pressure ratio
  do k=1,kgx
          af1 = (  tanh((zg(k)-zf1)/wf0) + 1.d0)*0.5d0
          af2 = ( -tanh((zg(k)-zf2)/wf0) + 1.d0)*0.5d0
    rbeta(k) = abs(alpha*af1*af2)
  enddo
 
!----------------------------------------------------------------------|
! set initial model  
  den0(kzero) = ro0
  pre0(kzero) = pr0 * tem0(kzero)

  ! set density pressure
  do k=kzero+1,kgx
    den0(k) = den0(k-1)   &
   *((1+rbeta(k-1))*tem0(k-1)+0.5d0*gm*gzm0(k-1)*dzg(k-1)) &
   /((1+rbeta(k)  )*tem0(k)  -0.5d0*gm*gzm0(k-1)*dzg(k-1))
    pre0(k) = pre0(kzero) &
   *(den0(k)/den0(kzero))*(tem0(k)/ tem0(kzero))
  enddo
  do k=kzero-1,1,-1
    den0(k) = den0(k+1)  &
   *((1+rbeta(k+1))*tem0(k+1)-0.5d0*gm*gzm0(k)*dzg(k)) &
   /((1+rbeta(k)  )*tem0(k)  +0.5d0*gm*gzm0(k)*dzg(k))
    pre0(k) = pre0(kzero) &
   *(den0(k)/den0(kzero))*(tem0(k)/ tem0(kzero))
  enddo
  ! set magnetic field 
  do k=1,kgx
    bmx0(k)  = sqrt(2*pre0(k)*rbeta(k))
  enddo

  do k=1,kx
     do j=1,jx
        do i=1,ix
     kg=mpid%mpirank_3d(3)*(kx-2*margin)+k
     jg=mpid%mpirank_3d(2)*(jx-2*margin)+j
     ig=mpid%mpirank_3d(1)*(ix-2*margin)+i
         ro(i,j,k)=den0(kg)
         pr(i,j,k)=pre0(kg)
         vx(i,j,k)=0.0d0
         vy(i,j,k)=0.0d0 
         vz(i,j,k)=0.0d0
         bx(i,j,k)=0.0d0
         by(i,j,k)=bmx0(kg)
         bz(i,j,k)=0.0d0
         eta(i,j,k)=0.0d0
         phi(i,j,k)=0.0d0

         roi(i,j,k)=den0(kg)
         pri(i,j,k)=pre0(kg)
         vxi(i,j,k)=0.0d0
         vyi(i,j,k)=0.0d0
         vzi(i,j,k)=0.0d0
         bxi(i,j,k)=0.0d0
         byi(i,j,k)=bmx0(kg)
         bzi(i,j,k)=0.0d0
         phii(i,j,k)=0.0d0
        enddo
     enddo
  enddo

  call pertub(vx,vy,vz,x,y,z)

  end subroutine model_setup


  subroutine pertub(vx,vy,vz,x,y,z)
  real(8),dimension(ix)      ,intent(in) :: x
  real(8),dimension(jx)      ,intent(in) :: y
  real(8),dimension(kx)      ,intent(in) :: z
  real(8),dimension(ix,jx,kx),intent(out) :: vx,vy,vz
  integer :: i,j,k

  do k=1,kx
     do j=1,jx
        do i=1,ix
         vx(i,j,k)=0.0d0
         vy(i,j,k)= amp*sin(2*pi*y(j)/ylamd)                      &
            *((tanh((z(k)-yptb1)/wptb1)-tanh((z(k)-yptb2)/wptb2)) &
            - (tanh((z(k)-yptb3)/wptb2)-tanh((z(k)-yptb4)/wptb1)))
         vz(i,j,k)=0.0d0
        enddo
     enddo
  enddo

  end subroutine pertub


end module model
