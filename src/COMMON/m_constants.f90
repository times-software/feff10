!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_constants.f90,v $:
! $Revision: 1.3 $
! $Author: jorissen $
! $Date: 2010/12/16 18:30:30 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*************************************************************************
      module constants

      implicit none

      ! Symbolic names for kind types of 4-, 2-, and 1-byte integers:
      integer, parameter :: I4 = selected_int_kind(9)
      integer, parameter :: I2 = selected_int_kind(4)
      integer, parameter :: I1 = selected_int_kind(2)
      ! Symbolic names for kind types of single- and double-precision reals:
      integer, parameter :: SP  = kind(1.0)
      integer, parameter :: DP  = kind(1.0D0)
      ! Symbolic names for kind types of single- and double-precision complex:
      integer, parameter :: SZ  = kind((1.0,1.0))
      integer, parameter :: DZ  = kind((1.0D0,1.0D0))
      ! Symbolic name for kind type of default logical:
      integer, parameter :: LGT = kind(.true.)
      private i1,i2,i4,sp,dp,sz,dz,lgt

      real(dp), parameter :: pi2 = 6.283185307179586476925286766559_dp
      real(dp), parameter :: pi  = 3.1415926535897932384626433832795_dp
      real(dp), parameter :: one = 1.0_dp
      real(dp), parameter :: zero  = 0.0_dp
      real(dp), parameter :: third = 1.0_dp/3.0_dp
      real(dp), parameter :: raddeg = 180.0_dp/pi
!     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      real(dp), parameter :: fa = 1.919158292677512811_dp
      complex(dz), parameter :: coni = (0.0_dp,1.0_dp)
      real(dp), parameter :: bohr = 0.529177249_dp
      real(dp), parameter :: ryd  = 13.605698_dp
      real(dp), parameter :: hart = 2.0_dp*ryd
      real(dp), parameter :: alpinv = 137.03598956_dp
      real(dp), parameter :: alphfs = 1.0_dp/alpinv

      ! from moduleseels.f 
      ! conversion from eV to Ry :
      real(dp), parameter :: ev2Ry = 1.0_dp/13.6058_dp
      !  h/2pi c in units eV a.u.
      real(dp), parameter :: hbarc_eV = 1973.2708_dp/0.529177_dp
	  real(dp), parameter :: hbarc_atomic = 137.04188_dp  ! i.e. in Ha, hence 27.2 times smaller than above
      ! electron rest mass times c^2 in au (ie, 1 * alfa * alfa), times eV/Ha (27.2)
      real(dp), parameter :: MeC2 =  511004.0_dp
      REAL(dp), parameter :: HOnSqrtTwoMe = 23.1761_dp
      ! Me c / hbar = 2.5896 10^12 m^(-1) = 137.04188 a.u.^(-1)
      real(dp), parameter :: MeCOnHbar = 137.04188_dp

      end module constants
