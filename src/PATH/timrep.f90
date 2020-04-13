!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: timrep.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine timrep (npat, ipat, rx, ry, rz, dhash,                 &
     &            ipol, ispin, evec, xivec,ica)   !KJ added ica 5/06

!     subroutine timrev(...) is modified for polarization case 
!     Time-orders path and returns path in standard order,
!     standard order defined below.
!     Input:  npat, ipat
!     Output: ipat in standard order (time reversed if necessary)
!             rx, ry, rz   contain x,y,z coordinates of the path atoms,
!             where z-axis is along polarization vector or first leg, if
!               running usual feff,
!             x-axis is chosen so that first atom, which does not lie on
!               z-axis, lies in xz-plane,
!               for elliptically polarized light, x-axis is along the
!               incidence direction
!             y-axis is cross product of two previos unit vectors
!             Standarrd order is defined so that first nonzero x,y and z
!             coords are positive.(Otherwise we use the inversion of
!             the corresponding unit vector)
!             dhash double precision hash key for path in standard
!                order

      use dimsmod, only: natx, npatx
	  use eels_inp,only: eels
	  use global_inp,only: do_nrixs
      double precision evec(3), xivec(3)
      common /atoms/ rat(3,0:natx), ipot(0:natx), ilb(0:natx)
      dimension ipat(npatx+1), rx(npatx), ry(npatx), rz(npatx)
      dimension ipat0(npatx+1), rx0(npatx), ry0(npatx), rz0(npatx)
      integer icase,ica  !KJ added 5/06
      double precision dhash, dhash0

!     Time reverses path if time reversing it will put it
!     in standard order.  Standard order is defined by min hash
!     number, using path hash algorithm developed for the path
!     degeneracy checker.  See subroutine phash for details.
!     Symmetrical paths are, of course, always standard ordered.
!     Also returns hash number for standard ordered path.

!     Use suffix 0 for (') in variable names

!KJ next block : prepare new calling argument for mpprmp
      if(ica.gt.0.and.ica.lt.8) then
         icase=ica   ! passed as input option from FEFF.INP
      elseif(eels.eq.1 .or. do_nrixs.eq.1) then !KJ 7-09 added NRIXS - Aleksi used to have sth like this in mpprmp
         icase=7  ! force no symmetry for eels
      else
         icase=-1 ! let mpprmp figure out what symmetry to use - default
      endif
!KJ end new block. 5/06


!     If no time-reversal standard ordering needed, make hash number
!     and return.  No timrev needed if 2 leg path (symmetrical).
      nleg = npat + 1
      ipat(nleg) = 0
      do 10 i = 1, npatx
         rx(i)   = 0
         ry(i)   = 0
         rz(i)   = 0
         rx0(i)   = 0
         ry0(i)   = 0
         rz0(i)   = 0
   10 continue
      call mpprmp(npat, ipat, rx, ry, rz, ipol, ispin, evec, xivec,icase)   !KJ added icase 5/06
      call phash (npat, ipat, rx, ry, rz, dhash)

      if (npat .le. 1) return

!     Make time reversed path

      ipat0(nleg) = ipat(nleg)
      do i = 1, npat
         ipat0(i) = ipat(nleg-i)
      enddo
      call mpprmp(npat, ipat0, rx0, ry0, rz0, ipol, ispin, evec, xivec,icase)  !KJ added icase 5/06
      call phash (npat, ipat0, rx0, ry0, rz0, dhash0)

!     turn off path reversal in special cases (make dhash0>dhash)
      if (ispin.ne.0 .and. ipol.ne.0) dhash0 = dhash+1

!     Do the comparison using hash numbers
!     Want representation with smallest hash number
      if (dhash0 .lt. dhash)  then
!        time reversed representation is smaller, so return
!        that version of the path
         dhash = dhash0
         do 300  i = 1, npat
            ipat(i) = ipat0(i)
            rx(i)   = rx0(i)
            ry(i)   = ry0(i)
            rz(i)   = rz0(i)
  300    continue
      endif

      return
      end
