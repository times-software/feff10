!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fxc.f90,v $:
! $Revision: 1.7 $
! $Author: jorissen $
! $Date: 2011/12/10 23:17:11 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================================================================
!     LDA EXCHANGE
!=======================================================================
subroutine ldafxc(ilast, ri, edens, exchange, fxc)
  
  use DimsMod, only: nrptx
  use constants
  implicit none

  integer ilast, exchange
  double precision ri(nrptx), edens(251)
  double precision fxc(251)
  double precision rs
  integer i

  do i= 1, ilast
     !       also calculate the local exchange term fxc = vxc*r_s**3/r**2
     if (edens(i) .le. 0.d0) then
        rs = 100.d0
        fxc(i) = 0.d0
     else
        rs = (edens(i)/3)**(-third)
        fxc(i) = rs**3 / ri(i)**2 / 6 * (-1.222d0/rs -0.75924d0/(11.4d0+rs))

        if (exchange .eq. 2) then
           fxc(i) = rs**3 / ri(i)**2 / 6 * (-1.222/rs)
        end if
     endif
     !     vvbh from Von Barth Hedin paper, 1971
     !     see eq.60-61 of Gross&Kohn, Adv. Quant. Chem. 21, p.255(1990).
     !     eq.60 fxc = d V_xc / d rho * (2\ell + 1) /(4*pi*r**2)
     !     where second factor comes from \delta(r-r')
     !     for dipole field \ell=1 
     !     fxc(i) = rs**3 / ri(i)**2 / 6 * (-1.22177412/rs -1.512/(30+rs))
     !     c below are coefficients in Zangwill/Soven paper
  end do

  return
end subroutine ldafxc
