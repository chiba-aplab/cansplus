subroutine chkNaN(qq,nanflag)
  implicit none

  real(8) :: qq

  integer :: nanflag

  nanflag = 0
  if(dabs(qq) < 0.0d0 .or. dabs(qq) > huge(qq))then
     nanflag = 1
  end if

  return
end subroutine chkNaN

