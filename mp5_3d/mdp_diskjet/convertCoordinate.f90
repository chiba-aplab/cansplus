!----test1012
module convertCoordinate
  implicit none

contains

  subroutine CylinderToCartesian(x,y,z,x1,x2,x3)
    implicit none
    
    real(8),intent(in) :: x,y,z

    real(8),intent(inout) :: x1,x2,x3

    x1 = x*cos(y)
    x2 = x*sin(y)
    x3 = z

    return
  end subroutine CylinderToCartesian

  subroutine CartesianToCylinder(x1,x2,x3,x,y,z)

    real(8),intent(in) :: x1,x2,x3

    real(8),intent(inout) :: x,y,z

    x = sqrt(x1**2 + x2**2)
    y = acos(x1/x)
    z = x3
    
    return
  end subroutine CartesianToCylinder
end module convertCoordinate
