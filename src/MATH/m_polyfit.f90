module polyfit_mod

contains
  subroutine polyfit(x0, y0, order, coeffs)
    ! fit a polynomial of order "order" to data x0,y0 using linear least-squares
    double precision, intent(in) :: x0(:)
    double complex, intent(in) :: y0(:)
    integer, intent(in) :: order
    double complex, intent(out) :: coeffs(0:order)

    double complex, allocatable :: F(:,:)
    integer i, ier

    allocate(F(size(x0),0:order))
    do i=0,order
      F(:,i) = x0**i
    end do

    call leastsq(F, y0, coeffs, ier)

    deallocate(F)
  end subroutine

  subroutine polyval(coeffs, x, y)
    double complex, intent(in) :: coeffs(:)
    double precision, intent(in) :: x(:)
    double complex, intent(out) :: y(:)

    double complex, allocatable :: F(:,:)
    integer i, order

    order = size(coeffs) - 1

    allocate(F(size(x),0:order))
    do i=0,order
      F(:,i) = x**i
    end do

    y = matmul(F, coeffs)
    deallocate(F)
  end subroutine

  subroutine leastsq(F, y, coeffs, ier)
    ! solve linear least squares
    ! y(i) ~= sum_j F(i,j) * coeffs(j)
    double complex, intent(in) :: F(:,:)
    double complex, intent(in) :: y(:)
    double complex, intent(out) :: coeffs(:)

    double complex, allocatable :: Ft(:,:), FtF(:,:), Fty(:), tmp(:,:)
    integer, allocatable :: ipiv(:)
    integer :: info, ier

    integer N, M
    integer i, j

    N = size(F,1)
    M = size(F,2)

    coeffs(:) = 0

    if (size(y) .ne. N .or. size(coeffs) .ne. M) then
      ier = -3
      return
    end if

    allocate(Ft(M,N), FtF(M,M), Fty(M), ipiv(M), tmp(M,M))

    Ft = transpose(F)
    FtF = matmul(Ft, F)
    Fty = matmul(Ft, y)

    ! use LU decomposition to solve for 1 (FtF)^-1 Fty
    ier = 0
    info = 0
    call zgetrf(M, M, FtF, M, ipiv, info)
    if (info < 0) then
      ier = -1
      goto 999
    end if

    call zgetrs('N', M, 1, FtF, M, ipiv, Fty, M, info)
    if (info < 0) then
      ier = -2
      goto 999
    end if

    coeffs(:) = Fty(:)

  999 continue
    deallocate(tmp, FtF, Fty, ipiv)
  end subroutine  

end module
