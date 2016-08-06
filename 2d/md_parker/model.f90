module model

  use const
  use mpi_setup, only : mpid

  implicit none
  private

  public :: model_setup


contains

  
  subroutine model_setup(ro,pr,vx,vy,vz,bx,by,bz,phi        &
                      ,roi,pri,vxi,vyi,vzi,bxi,byi,bzi,phii &
                      ,x,dx,xm,y,dy,ym,gy,eta,min_dx)

!---Input & Output
  real(8),dimension(ix)   ,intent(out) :: x,dx
  real(8),dimension(0:ix) ,intent(out) :: xm  
  real(8),dimension(jx)   ,intent(out) :: y,dy
  real(8),dimension(0:jx) ,intent(out) :: ym  
  real(8),dimension(ix,jx),intent(out) :: ro,pr
  real(8),dimension(ix,jx),intent(out) :: vx,vy,vz
  real(8),dimension(ix,jx),intent(out) :: bx,by,bz
  real(8),dimension(ix,jx),intent(out) :: eta,phi
  real(8),dimension(ix,jx),intent(out) :: gy
  real(8),dimension(ix,jx),intent(out) :: roi,pri
  real(8),dimension(ix,jx),intent(out) :: vxi,vyi,vzi
  real(8),dimension(ix,jx),intent(out) :: bxi,byi,bzi,phii
  real(8),                    intent(out) :: min_dx

  integer :: i,j,k
  integer :: ig,jg
  integer :: izero,jzero
  real(8) :: af1,af2
  real(8) :: dxmax,dymax
  real(8),dimension(ix,jx) :: curx,cury,curz
  real(8),dimension(igx) :: xg,dxg
  real(8),dimension(0:igx) :: xmg
  real(8),dimension(jgx) :: yg,dyg
  real(8),dimension(0:jgx) :: ymg
  real(8),dimension(igx,jgx) :: gyg
  real(8),dimension(jgx) :: gym0,tem0,rbeta,den0,pre0,bmx0   


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

!---Step 2a.-------------------------------------------------------------|
! set individual x-grid 
  do i=1,ix
     ig = mpid%mpirank_2d(1)*(ix-2*margin)+i
     x(i)=xg(ig)
     dx(i) = dxg(ig)
  enddo
  do i=0,ix
     ig = mpid%mpirank_2d(1)*(ix-2*margin)+i
     xm(i) = xmg(ig)
  end do

!---Step 2b.-------------------------------------------------------------|
! set individual y-grid 
  do j=1,jx
     jg = mpid%mpirank_2d(2)*(jx-2*margin)+j
     y(j) = yg(jg)
     dy(j) = dyg(jg)
  enddo
  do j=0,jx
     jg = mpid%mpirank_2d(2)*(jx-2*margin)+j
     ym(j) = ymg(jg)
  end do

! calculate min_dx
  min_dx = min(minval(dxg),minval(dyg))

!---Step 3a.----------------------------------------------------------|
! set gravitation
  do j=1,jgx
     do i=1,igx
     gyg(i,j)= - g0 * tanh(yg(j)/wg)
     end do
     gym0(j)= - g0 * tanh((yg(j)+0.5d0*dyg(j))/wg)
  end do

!--Step 3b.----------------------------------------------------------|
! set individual gravitation
 do j=1,jx
     do i=1,ix
        ig = mpid%mpirank_2d(1)*(ix-2*margin)+i
        jg = mpid%mpirank_2d(2)*(jx-2*margin)+j
        gy(i,j) = gyg(ig,jg)
     enddo
  enddo
 
!-------------------------
!   temperature
  do j=1,jgx
    tem0(j)=1.d0 + &
           0.5d0*(tcor-1.d0)*(tanh((abs(yg(j))-ytr)/wtr)+1.d0)
  enddo

!-------------------------
!   pressure ratio
  do j=1,jgx
          af1 = (  tanh((yg(j)-yf1)/wf0) + 1.d0)*0.5d0
          af2 = ( -tanh((yg(j)-yf2)/wf0) + 1.d0)*0.5d0
    rbeta(j) = abs(alpha*af1*af2)
  enddo
 
!----------------------------------------------------------------------|
! set initial model  
  den0(jzero) = ro0
  pre0(jzero) = pr0 * tem0(jzero)

  ! set density pressure
  do j=jzero+1,jgx
    den0(j) = den0(j-1)   &
   *((1+rbeta(j-1))*tem0(j-1)+0.5d0*gm*gym0(j-1)*dyg(j-1)) &
   /((1+rbeta(j)  )*tem0(j)  -0.5d0*gm*gym0(j-1)*dyg(j-1))
    pre0(j) = pre0(jzero) &
   *(den0(j)/den0(jzero))*(tem0(j)/ tem0(jzero))
  enddo
  do j=jzero-1,1,-1
    den0(j) = den0(j+1)  &
   *((1+rbeta(j+1))*tem0(j+1)-0.5d0*gm*gym0(j)*dyg(j)) &
   /((1+rbeta(j)  )*tem0(j)  +0.5d0*gm*gym0(j)*dyg(j))
    pre0(j) = pre0(jzero) &
   *(den0(j)/den0(jzero))*(tem0(j)/ tem0(jzero))
  enddo
  ! set magnetic field 
  do j=1,jgx
    bmx0(j)  = sqrt(2*pre0(j)*rbeta(j))
  enddo

  !$OMP PARALLEL DO PRIVATE(i,j,ig,jg)
  do j=1,jx
     do i=1,ix
  jg=mpid%mpirank_2d(2)*(jx-2*margin)+j
  ig=mpid%mpirank_2d(1)*(ix-2*margin)+i
      ro(i,j)=den0(jg)
      pr(i,j)=pre0(jg)
      vx(i,j)=0.0d0
      vy(i,j)=0.0d0 
      vz(i,j)=0.0d0
      bx(i,j)=bmx0(jg)
      by(i,j)=0.0d0
      bz(i,j)=0.0d0
      eta(i,j)=0.0d0
      phi(i,j)=0.0d0

      roi(i,j)=den0(jg)
      pri(i,j)=pre0(jg)
      vxi(i,j)=0.0d0
      vyi(i,j)=0.0d0
      vzi(i,j)=0.0d0
      bxi(i,j)=bmx0(jg)
      byi(i,j)=0.0d0
      bzi(i,j)=0.0d0
      phii(i,j)=0.0d0
     enddo
  enddo
  !$OMP END PARALLEL DO

  call pertub(vx,vy,vz,x,y)

  end subroutine model_setup


  subroutine pertub(vx,vy,vz,x,y)
  real(8),dimension(ix)      ,intent(in) :: x
  real(8),dimension(jx)      ,intent(in) :: y
  real(8),dimension(ix,jx),intent(out) :: vx,vy,vz
  integer :: i,j,k

  !$OMP PARALLEL DO
  do j=1,jx
     do i=1,ix
      vx(i,j)= amp*sin(2*pi*x(i)/xlamd)                      &
         *((tanh((y(j)-yptb1)/wptb1)-tanh((y(j)-yptb2)/wptb2)) &
         - (tanh((y(j)-yptb3)/wptb2)-tanh((y(j)-yptb4)/wptb1)))
      vy(i,j)=0.0d0
      vz(i,j)=0.0d0
     enddo
  enddo
  !$OMP END PARALLEL DO

  end subroutine pertub


end module model
