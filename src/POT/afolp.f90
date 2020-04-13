!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: afolp.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine afolp ( nph, nat, iphat, rat, iatph, xnatph,           &
     &                novr, iphovr, nnovr, rovr, folp, folpx, iafolp,   &
     &                edens, edenvl,                                    &
     &                dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,  &
     &                ixc, rhoint, vint, rs, xf, xmu, xmunew,           &
     &                rnrmav, qtotel, inters, totvol)

!     find folp(iph) automatically and recalculates
!     interstitial parameters, rmt, vint, etc.
!     written by ala 11.97
      use constants
      use DimsMod, only: natx, nphx=>nphu, novrx
        
      implicit double precision (a-h, o-z)

      dimension iphat(natx)
      dimension rat(3,natx)
      dimension iatph(0:nphx)
      dimension xnatph(0:nphx)
      dimension novr(0:nphx)
      dimension iphovr(novrx,0:nphx)
      dimension nnovr(novrx,0:nphx)
      dimension rovr(novrx,0:nphx)
      dimension folp(0:nphx), folpx(0:nphx)
      dimension edens(251,0:nphx), edenvl(251,0:nphx)
      dimension dmag(251,0:nphx+1)
      dimension vclap(251,0:nphx)
      dimension vtot (251,0:nphx), vvalgs (251,0:nphx)
      dimension imt(0:nphx)
      dimension inrm(0:nphx)
      dimension rmt(0:nphx), rmtx(0:nphx)
      dimension rnrm(0:nphx)
      character*512 slog

! Debug: FDV
!       write(6,*) ' Entering afolp'
!       do iph=0,nph
!         write(6,fmt='(a,i4,3f16.10)'), &
!         'rmt, folpx, folp: ', iph, rmt(iph)*bohr, folpx(iph), folp(iph)
!       end do

      do 5 iph=0,nph
         rmtx(iph) = rmt(iph) / folp(iph)
   5  continue

      !call wlog(' iph, rnrm(iph)*bohr, rmt(iph)*bohr, folp(iph)')
      call wlog('type, norman radius, muffin tin, overlap factor')
      if (iafolp.ge.0) then
         do 400  iph = 0, nph
!          old algorithm for automatic overlap
!          folp(iph) = 1 + 0.7*(rnrm(iph)/rmt(iph) - 1)
           folp(iph) = folpx(iph)
           rmt(iph) = folp(iph) * rmtx(iph)
! Debug: FDV
!          write(6,fmt='(i4,3f12.8)') iph, rmt(iph)*bohr, folp(iph), rmtx(iph)

  398      format(i5, 1p, 3e13.5)
           write(slog,398) iph, rnrm(iph)*bohr, rmt(iph)*bohr, folp(iph)
           call wlog(slog)
  400    continue

! Debug: FDV
!       do iph=0,nph
!         write(6,fmt='(a,i4,f16.10)'), 'rmt: ', iph, rmt(iph)*bohr
!       end do

         idmag = 0
!        write(6,fmt='(a)') 'Before istprm'
         call istprm (nph, nat, iphat, rat, iatph, xnatph,              &
     &               novr, iphovr, nnovr, rovr, folp, folpx, iafolp,    &
     &               edens, edenvl, idmag,                              &
     &               dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,   &
     &               ixc, rhoint,vint, rs, xf, xmu, xmunew,             &
     &               rnrmav, qtotel, inters, totvol)
!        write(6,fmt='(a)') 'After istprm'
! Debug: FDV
!        stop

      endif

      return
      end
