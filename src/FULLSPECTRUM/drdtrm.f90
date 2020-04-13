      subroutine drdtrm (omega, eps, iepts,tau,numden)
      !calculates epsilon for a free electron gas of density numden
      !and lifetime tau on energy grid omega
      use constants
      implicit none
      include 'HEADERS/params.h'

      real*8 omega(fullpts) , tau, numden,wp2
      complex eps(fullpts)
      real*8 trnbeg, trnend, gam, eps2, eps1, denom, hbar
      parameter(hbar = 6.58E-16) !eV*seconds, for conversions
      integer i, iepts
      character*512 slog


!     gam=hc/(2.0*pi**tau*c) !drude width in eV
      gam=hbar/tau/hart
!     gam=1/(10/hart) !10 eV for a test
      wp2=4*pi*numden
!     write(slog,fmt="(' wp2: ',e20.10)") wp2
!     call  wlog (slog)
!     write(slog,fmt="(' gam: ',e20.10)") gam
!     call  wlog (slog)


      do i=1,iepts
        denom=omega(i)**2+gam**2
        eps2=wp2*gam/omega(i)/denom
        eps1=-1.0*wp2/denom
        eps(i)=eps1+coni*eps2
      enddo

      open (unit=77,file='drude.dat')
      write(unit=77,fmt="('# gam (eV):',e20.10)") gam*hart
      write(unit=77,fmt="('# wp (eV):',e20.10)") sqrt(wp2)*hart
      do i=1,iepts
        write(unit=77,fmt="(3e20.10)") omega(i), eps(i)
      enddo
      close (unit=77)



      end
