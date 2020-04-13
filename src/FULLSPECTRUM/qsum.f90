
!     givnen eps2, numden computes n_eff via eps_2 sumrule.
      subroutine qsum (numden,eps2,omega,iepts,neff)
      use constants
      implicit none
      include 'HEADERS/params.h'

      character*512 slog
      real*8 emin, emax
      real*8 B, enrgy, m, neff, n,                                        &
     &     c, na, e, rho, am,                                           &
     &     mu, refl, eloss, eps1, e2curr, k,e2prev,                     &
     &     omega(fullpts), eprev,                                       &
     &     factor,numden,eps2(fullpts)
      integer iepts, ii, ios,i
      character comment*100
      
        B=1/(2*pi**2*numden) !for opcons.dat
        
        neff=0.
        do i=1,iepts-1
          neff=neff                                                     &
     &         +(omega(i+1)*eps2(i+1)+omega(i)*eps2(i))/2               &
     &         *(omega(i+1)-omega(i))
        enddo
        neff=neff*B

      end
