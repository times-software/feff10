      subroutine egrid_lin(npts,emin,ehi,omega)
      !Generates energy grid for full spectrum.
      !Currently we follow the prescriptions of Rivas
      !(i.e. a 1000 point exponential grid for 0-100 eV,
      !followed by a linear grid above that).
      !Now we have switched to a simple linear grid
      use constants
      implicit none
      include 'HEADERS/params.h'
      real*8 trans,emin,eps
      parameter ( trans = 100.0/hart) !energy to switch from log -> linear grid
!     parameter ( trans = 10.0/hart) !energy to switch from log -> linear grid
!     parameter ( emin = 0.01/hart) !smallest energy in grid
      parameter ( eps = 0.01/hart) !smallest energy we will use
      integer npts,i
      real*8 ehi, omega(fullpts), de,dle
      character*512 slog

      goto 100 !skip old grid and make a new, simple, linear one

!     sanity test for energy grid
      if (npts.gt.fullpts) then
        call wlog ('  ...too many points in energy grid.')
        call wlog ('  Check EGRID card.')
        call wlog ('  If you need more points, increase param fullpts.')
        write(slog,fmt="('      npts = ',i10)") npts
        call wlog (slog)
        write(slog,fmt="('   fullpts = ',i10)") fullpts
        call wlog (slog)
        stop
      endif
      if (nlog.gt.npts) then
        call wlog ('  ...more points in logarithmic section ')
        call wlog ('     of the energy grid than in the whole')
        call wlog ('     grid (nlog>npts).')
        call wlog ('  Check EGRID card.')
        call wlog ('  If you need fewer points in the logarithmic')
        call wlog ('  section, decrease param nlog.')
        write(slog,fmt="('   nlog = ',i10)") nlog 
        call wlog (slog)
        write(slog,fmt="('   npts = ',i10)") npts
        call wlog (slog)
        stop
      endif

!     initialized grid
      do i=1,fullpts
        omega(i)=0.0
      enddo

!     energy steps: dle is the linear step in log(omega) for the low energy
!     portion of the spectrum; de is the linear step in omega for the
!     high energy part of the spectrum
      dle=(log(trans)-log(emin))/nlog
      de=(ehi-trans)/(npts-nlog-1) 

      write(slog,fmt="('   emin (eV): ',e20.10)") emin*hart
      call wlog (slog)

!     set grid
      omega(1)=emin 
      do i=2,npts
        if (i.le.nlog) then
!         omega(i)=exp(dle)*omega(i-1)
          omega(i)=exp((i-1)*dle+log(emin))
        else
          omega(i)=omega(i-1)+de
        endif 
      enddo 
!     write(slog,fmt="('   emin (eV): ',e20.10)") (omega(1))*hart
!     call wlog (slog)

100   continue 

      !back to a linear grid ...
      if (emin.le.0.0) emin=eps !make sure no non-positive energies
      de=(ehi-emin)/(npts-1)
      omega(1)=emin 
      do i=2,npts
        omega(i)=omega(i-1)+de
      enddo 
!     open(unit=20,file='egrid.dat')
!     do i=1,npts
!       write(20,fmt="(i6,2e20.10)") i,omega(i),omega(i)*hart
!     enddo 
!     close(unit=20)
      
      
      return
      end
