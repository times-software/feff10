!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: getedg.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine getedg (ihole, iz, emu)
        use constants, only: bohr, ryd, hart
        use Elam_Table_of_Energies, only:exmu
        implicit none
        integer, intent(in)  :: ihole, iz
        real(8), intent(out) :: emu

      ! Insert corrected edges from G.P. Williams' and W.T. Elam's tables.
      ! J Kas - if iz > 100, don't use tables
      IF(iz > 100) THEN
         PRINT*, 'Elam energies only exist up to Z=100: ' //  &
     &           'ignoring SETEDGE' 
      ELSE
         if(exmu(iz,ihole).gt.0.) emu = exmu(iz,ihole) / hart
      END IF

      return
      end


      subroutine preved (ecurr, eprev, nz, ncomps)
      !returns the energy of the edge in the spectrum of a material
      !composed of elements of atomic numbers nz that is less than
      !ecurr by the smallest amount (i.e. this subroutine finds the 'previous' edge).

        use dimsmod, only: nheadx,nex,nphx=>nphu
        use Elam_Table_of_Energies, only:exmu
        use constants
        implicit none
      real*8 ecurr, eprev
      double precision emu
      integer nz(nphx), ncomps, i, ihole, iz

      eprev=0.0 
      do i=1,ncomps
        iz=nz(i)
        do ihole=1,29
          emu=max(exmu(iz,ihole)/hart,0d0)
          if (emu.lt.ecurr.and.emu.gt.eprev) eprev=emu
        enddo
      enddo

      return
      end


      subroutine nexted (ecurr, enext, nz, ncomps)
      !returns the energy of the edge in the spectrum of a material
      !composed of elements of atomic numbers nz that is greater than
      !ecurr by the smallest amount (i.e. this subroutine finds the 'next' edge).
        use dimsmod, only: nheadx,nex,nphx=>nphu
        use Elam_Table_of_Energies, only:exmu
        use constants
        implicit none
      real*8 ecurr, enext, emu
      integer nz(nphx), ncomps, i, ihole, iz
      character*512 slog

      enext=1.0e8 !energy higher than any physical edge energy
      do i=1,ncomps
        iz=nz(i)
        do ihole=1,29
          emu=dble(max(exmu(iz,ihole)/hart,0d0))
          if (emu.gt.ecurr.and.emu.lt.enext) enext=dble(emu)
        enddo
      enddo

      return
      end




