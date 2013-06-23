!======================================================================|
  subroutine bd_sppy(mbnd,margin,qq,ix,jx,kx)
!======================================================================|
!
! NAME  bdsppy
!
! PURPOSE
!    apply symmetric boundary condition. 
!    The symmetric point is between the grid point.
!    The values out of the boundary have same sign.
!
! INPUTS & OUTPUTS
!    qq(ix,jx): [double] variable
!
! OUTPUTS
!    None
!
! INPUTS
!    ix,jx: [integer] dimension size
!    margin: [integer] margin, i.e. # of grid points outside the boundary
!    mbnd: [integer] If mbnd=0, smaller 'i' side. 
!                    If mbnd=1, larger  'i' side.
!
! HISTORY
!    written 2002-3-1 T. Yokoyama
!
!----------------------------------------------------------------------|
!----------------------------------------------------------------------|

      implicit none
      integer,intent(in) :: margin,ix,jx,kx
      integer,intent(in) :: mbnd

      real(8),dimension(ix,jx,kx) :: qq

      integer :: i,j,k,kbnd


      if (mbnd.eq.0) then
        kbnd=1+margin
        do k=1,margin
		do j=1,jx
        do i=1,ix
          qq(i,j,kbnd-k) = qq(i,j,kbnd-1+k)
        enddo
        enddo
		enddo
      else
        kbnd=kx-margin
        do k=1,margin
		do j=1,jx
        do i=1,ix
          qq(i,j,kbnd+k) = qq(i,j,kbnd+1-k)
        enddo
        enddo
		enddo
      endif

      return
      end subroutine bd_sppy
