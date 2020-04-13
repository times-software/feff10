!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: feff.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program feff
      integer iabs, nabs, nss !KJ added nss 7-06
      logical ceels !KJ added ceels 5-06

      call rdinp(nabs,nss, ceels)
      call ffmod1
      call ffmod7
      call ffmod8
      call ffmod2
      do 900 iabs = 1, nabs
        if (nabs.gt.1) call ffsort (iabs,nss,ceels) !KJ 5-6 added second argument
        call ffmod3
        call ffmod4
        call ffmod5
        call ffmod6 (iabs)
        call ffmod9
        call eelsmod  !KJ added 3-06 - Josh, changed name to eelsmod 5-08
 900  continue
      stop
      end
