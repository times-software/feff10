module rotation_mod
contains
  subroutine cross(a,b,c)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Calculate c = a X b for 3-vectors a and b
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real*8, dimension(3), intent(in) :: a, b
    real*8, dimension(3), intent(out) :: c

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end subroutine

  subroutine rotation_axis_angle(a,b, axis, theta)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Find rotation axis and angle that transforms a into b
    !
    ! a,b - real 3 vectors
    !
    ! Returns:
    !   axis - rotation axis, real 3 vector
    !   theta - rotation angle in radians
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real*8, dimension(3), intent(in) :: a, b
    real*8, dimension(3), intent(out) :: axis
    real*8, intent(out) :: theta

    call cross(a,b,axis)
    theta = asin(sqrt(sum(axis**2) / (sum(a**2)*sum(b**2))))
  end subroutine

  subroutine rotation_matrix(axis, theta, R)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Build a rotation matrix from a rotion axis and angle
    !
    ! axis  - real 3 vector
    ! theta - angle in radians
    !
    ! Returns:
    !   R - 3x3 matrix
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real*8, dimension(3), intent(in) :: axis
    real*8, intent(in) :: theta
    real*8, dimension(3,3), intent(out) :: R

    real*8, dimension(3) :: u
    real*8 :: c, d, s
    
    u = axis / sqrt(sum(axis**2))

    c = cos(theta)
    d = 1 - c
    s = sin(theta)

    R(1,1) = c + u(1)**2 * d 
    R(1,2) = u(1) * u(2) * d - u(3) * s
    R(1,3) = u(1) * u(3) * d + u(2) * s
    R(2,1) = u(2) * u(1) * d + u(3) * s
    R(2,2) = c + u(2)**2 * d
    R(2,3) = u(2) * u(3) * d - u(1) * s
    R(3,1) = u(3) * u(1) * d - u(2) * s
    R(3,2) = u(3) * u(2) * d + u(1) * s
    R(3,3) = c + u(3)**2 * d

  end subroutine

  subroutine rotate(R, v, Rv)
    real*8, dimension(3,3), intent(in) :: R
    real*8, dimension(3), intent(in) :: v
    real*8, dimension(3), intent(out) :: Rv

    integer i, j

    do i=1,3
      Rv(i) = 0.0
      do j=1,3
        Rv(i) = Rv(i) + R(i,j) * v(j)
      end do
    end do
  end subroutine

  subroutine rotate_in_place(R, v)
    real*8, dimension(3,3), intent(in) :: R
    real*8, dimension(3), intent(inout) :: v
    
    real*8, dimension(3) :: Rv

    call rotate(R, v, Rv)
    v(:) = Rv(:)
  end subroutine

end module
