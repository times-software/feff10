!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: xx.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function xx (j)
      implicit double precision (a-h, o-z)
!     x grid point at index j, x = log(r), r=exp(x)
      parameter (delta = 0.050000000000000)
      parameter (c88   = 8.800000000000000)
!     xx = -8.8 + (j-1)*0.05
      xx = -c88 + (j-1)*delta
      return
      end

      double precision function rr(j)
      implicit double precision (a-h, o-z)
!     r grid point at index j
      rr = exp (xx(j))
      return
      end

      function ii(r)
      implicit double precision (a-h, o-z)
!     index of grid point immediately below postion r
      parameter (delta = 0.050000000000000)
      parameter (c88   = 8.800000000000000)
!     ii = (log(r) + 8.8) / 0.05 + 1
      ii = (log(r) + c88) / delta + 1
      return
      end
