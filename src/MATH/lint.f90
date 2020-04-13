
!  Linear interpolation subroutine. Similiar to one in Numerical Recipes.
      subroutine lint(xa,ya,n,flag,khi,klo,x,y)
      implicit none
      logical flag
      integer n
!  try double precision here
       real*8 x, y, xa(n), ya(n)
!      double precision x, y, xa(n), ya(n)
      integer k, khi, klo
       real*8 a, b, h
!      double precision a, b, h

! Arguments:
!      xa = array of x-values of function to be interpolated
!      ya =   "   "  y   "    "      "     "  "      "
!      n = diminsion of xa and ya
!      flag -> false if lint was recently called with a smaller
!              value of x
!      khi = ?
!      klo = ?
!      hi micah

!     if x is outside the domain of the function, return 0.
!     if(x.le.xa(1).or.x.ge.xa(n)) then 
      if(x.le.xa(1)) then 
          y = 0.
      else
         

!     if we can not save work by locating array points on either side
!     of the x-value we want the function interpolated at, find such
!     points. Otherwise we expect klo and khi sent as arguments.
          if(flag) then 
              klo = 1
              khi = n
!     if klo and khi do not yet index adjacent array entries, move them 
!     closer together.
    1         if(khi-klo.gt.1) then
                  k = (khi+klo)/2
                  if(xa(k).gt.x) then
                      khi = k
                  else
                      klo = k
                  end if
!      loop until khi and klo differ by one.
              goto 1
              end if
              flag = .false.
          else
!       if the value we want to interpolate the function at is not in the range
!       between x(khi) and x(klo), change klo and khi to make this the case.
    6         if(xa(khi).lt.x) then
                  klo = khi
                  khi = khi + 1
              goto 6
              end if
          end if
          h = xa(khi) - xa(klo)
! Modified by FDVL
! Avoid warnings in gfortran.
! We should not be using pause, all problems should cause stops
!         if(h.eq.0.) pause 'bad xa input in lint'
          if (h .eq. 0.) then
            call wlog('bad xa input in lint')
            stop
          end if

!         find the interpolated value.
          a = (xa(khi) - x)/h
          b = (x - xa(klo))/h
          y = a*ya(klo) + b*ya(khi)
      end if
      return
      end

