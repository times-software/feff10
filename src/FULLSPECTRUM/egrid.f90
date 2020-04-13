      subroutine egrid(emin,ehi,nz,ncomps,omega,npts)
      !Generate energy grid spanning the interval (emin, ehi). Grid is
      !regular in k for the lowest energy electron that contributes to absorption at each energy. 

      !inputs:
      !   ehi: highest energy in grid
      !   emin: lowest energy in grid
      !   nz: atomic numbers of all elements in commpound
      !   ncomps: number of different elements in compound

      !outputs:
      !   omega(maxpts): energy grid
      !   npts: actual number of points filled in the grid.

      !Grid is such that omega(1)=emin and omega(npts)=ehi.
      !The coarsness of the grid is controlled by the parameter xkstep in HEADERS/params.h.
      use constants
      use dimsmod, only: nheadx,nex,nphx=>nphu

      implicit none

      include 'HEADERS/params.h'

      integer npts,i,nz(nphx),ncomps
      real*8 ehi, emin, omega(fullpts)
      !ckcurr -- momentum of current energy point
      !cknext -- momentum of next edge for the electron of the current edge
      real*8 ckcurr, cknext
      !ecurr, enext -- energies of the previous and current edges
      real*8 ecurr, enext
      character*512 slog

      !avoid point at omega=0
      if (emin.lt.0.001/hart) emin=0.001/hart

      omega(1)=emin
      call nexted(emin,enext,nz, ncomps)
!     write(slog,fmt="('picked enext before loop. enext: ',e20.10)") enext
!     call wlog(slog)
      call preved(emin,ecurr,nz, ncomps)
!     write(slog,fmt="('picked ecurr before loop. ecurr: ',e20.10)") ecurr
!     call wlog(slog)
      ckcurr=sqrt(2*(emin-ecurr))
      cknext=sqrt(2*(enext-ecurr))

      do i=2,fullpts
        ckcurr=ckcurr+xkstep
!       write(slog,fmt="('i: ',i20)") i
!       call wlog (slog)
!       write(slog,fmt="('xkstep: ',e20.10)") xkstep
!       call wlog (slog)
!       write(slog,fmt="('ckcurr: ',e20.10)") ckcurr
!       call wlog (slog)
!       write(slog,fmt="('cknext: ',e20.10)") cknext
!       call wlog (slog)
        if(ckcurr.ge.cknext) then
          ckcurr=0.0
          ecurr=enext
          omega(i)=ecurr
          call nexted(ecurr,enext,nz,ncomps)
          cknext=sqrt(2*(enext-ecurr))
!         call wlog ('switched edges')
!         write(slog,fmt="('ecurr, enext, cknext:', 3e20.10)") ecurr, enext,cknext
!         call wlog(slog)
        else
          omega(i)=ecurr+ckcurr**2/2
        end if
!       write(slog,fmt="('omega: ',e20.10)") omega(i)
!       call wlog (slog)
        if (omega(i).gt.ehi) then
          omega(i)=ehi
          exit
        end if
      enddo

      npts=i
      if (npts.eq.fullpts) then
        call wlog ('Energy grid clipped at upper end.')
        call wlog ('Adjust parameters fullpts (to get more points in the grid) or xkstep (to get a coarser grid).')
      end if
      
      return
      end
