!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: scmt.f90,v $:
! $Revision: 1.17 $
! $Author: hebhop $
! $Date: 2012/11/29 23:21:14 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine scmt (iscmt, ecv, nph, nat, vclap, edens, edenvl, vtot, vvalgs, rmt, rnrm,qnrm,      &
     &                ixc, rhoint, vint, xmu, jumprm, xnferm, xnvmu, xnval,      &
     &                x0, ri, dx, xnatph, xion, iunf, iz, adgc, adpc, dgc,dpc, ihole,      &
     &                rat,iatph,iphat, lmaxsc, rhoval, xnmues, ok, rgrd, nohole, nscmt, icoul, ca1, rfms1, lfms1 & !)
                      ,edos,scfdos)

!     Finds new Fermi level (xmu), electron counts (qnrm) and new valence densities (rhoval).
      use constants
      use DimsMod, only: nphx=>nphu, nrptx, lx, natx
      implicit double precision (a-h, o-z)

!     input
      dimension dmagx(nrptx), dmag0(251)
      dimension vclap(251,0:nphx)
      dimension vtot (251,0:nphx), vvalgs (251,0:nphx)
      dimension xnval (30,0:nphx)
      dimension qnrm(0:nphx), dq(0:nphx)
      dimension ri(nrptx), ri05(251), nr05(0:nphx)
      dimension xnatph(0:nphx), iz(0:nphx), xion(0:nphx)
      dimension rat(3,natx),iatph(0:nphx),iphat(natx), lmaxsc(0:nphx)
      real  rfms1
      real*8 :: rmt(0:nphx),rnrm(0:nphx)

      real*8, intent(inout) :: xnvmu(0:lx,0:nphx+1)
      real*8, intent(inout) :: xnmues(0:lx, 0:nphx)

!     input and output
      dimension edens(251,0:nphx), edenvl(251,0:nphx)
      dimension rhoval(251,0:nphx+1)

!     work space
      dimension dum(nrptx), vtotph(nrptx),vvalph(nrptx)
      dimension dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1)
      dimension adgc(10,30,0:nphx+1), adpc(10,30,0:nphx+1)
      dimension dgcn(nrptx,30), dpcn(nrptx,30)
      complex*16 yrhoce(251,0:nphx), yrhocp(251,0:nphx)

      integer iph
!     complex energy grid emg is decomposed into em and eref
      parameter (negx = 80)
      complex*16 emg(negx), em, eref, ee, ep, fl, fr, fxa
!     nflrx should be odd and defines the max of Im energy for the countour 
      parameter (nflrx = 17)
      dimension step(nflrx)
!     stuff from feff.f for rdinp, pathfinder and genfmt
      logical wnstar, upok, ok
!     Following passed to pathfinder, which is single precision.
      character*512 slog

      complex,    allocatable :: gtr(:,:)
      complex*16, allocatable, dimension(:,:):: xrhoce, xrhocp, xrhole, ph
      complex*16, allocatable :: yrhole(:,:,:)

      real*8 scfdos(negx,0:lx,0:nphx),edos(negx)

      integer ient
      data ient /0/

!     save staff from rdinp, so no need to call it again
      save   ri05, ient

      ! Allocate local variables
      allocate(gtr(0:lx, 0:nphx),xrhoce(0:lx,0:nphx), xrhocp(0:lx,0:nphx))
      allocate(xrhole(0:lx,0:nphx),yrhole(251,0:lx,0:nphx),ph(lx+1, 0:nphx))


      ient = ient + 1
      if (ient.eq.1) then
         do 15 i= 1,251
  15     ri05(i) = exp (-8.8+0.05*(i-1))
      endif

!      write (slog,10) iscmt, nscmt
!  10  format('              SCF ITERATION NUMBER',i3,'  OUT OF ',i3)
      write (slog,10) iscmt
  10  format('SCF ITERATION NUMBER',i3)
      call wlog(slog)
      write(29,*) trim(slog)

      !call wlog (' Calculating energy and space dependent l-DOS ....')

!     initialize new valence density
      rhoval(:,0:nphx) = 0.d0

!     polarization average in scmt and ldos
      call grids (ecv, xmu, negx, neg, emg, step, nflrx)

!     ie - is number of energy points calculated
      ie = 0
      ee = emg(1)
      ep = dble(ee)
	  xrhoce(:,:)=dcmplx(0)
	  xnmues(:,:)=0.d0
	  yrhoce(:,:)=dcmplx(0)
      iflr = nflrx
      iflrp = nflrx

!     Start the cycle over energy points (ie)
  25  continue
      ie = ie + 1

      xrhocp(:,0:nph)=xrhoce(:,0:nph)
	  yrhocp(:,0:nph)=yrhoce(:,0:nph)

      if (ie.eq.1 .or. mod(ie,20).eq.0) then
         write(slog,30) ie, dble(ee)*hart
   30    format('     point # ', i3, '  energy = ', f7.3)
         call wlog(slog)
      endif

      do iph = 0, nph
         dmag0(:) = 0.d0
!c       use spin-unpolarized case to get SCF. set dmagx to zero
!c       may want to replace dmag0 with dmag(1,iph) for spin-dependent extension of SCF procedure.
         call fixvar (rmt(iph),edens(1,iph),vtot(1,iph),dmag0, vint, rhoint, dx, rgrd, jumprm, vjump, ri, vtotph, dum, dmagx) !makes vtotph,dum,dmagx out of vtot,edens,dmag0
         if (mod(ixc,10) .ge.5) then
            if (jumprm .gt. 0) jumprm = 2
            call fixvar (rmt(iph), edenvl(1,iph), vvalgs(1,iph), dmag0, &
              vint, rhoint, dx, rgrd , jumprm, vjump, ri, vvalph, dum, dmagx)
            if (jumprm .gt. 0) jumprm = 1
         endif

         call fixdsx (iph, dx, rgrd , dgc, dpc, dgcn, dpcn)  !makes dgcn, dpcn out of dgc, dpc
         jri = (log(rmt(iph)) + x0) / rgrd + 2
         jri1 = jri+1
         eref = vtotph(jri1)
         do 40 i = 1, jri1
  40     vtotph(i) = vtotph(i) - eref
         if (ixc.ge.5) then
            do 50 i = 1, jri1
  50        vvalph(i) = vvalph(i) - eref
         else
            do 60 i = 1, jri1
  60        vvalph(i) = vtotph(i)
         endif

         itmp = 0
         if (iph.eq.0 .and. nohole.lt.0) itmp = ihole
         call rholie( ri05, nr05(iph), rgrd, x0, ri, ee, ixc, rmt(iph), rnrm(iph),    &     ! vtotph,vvalph,(a)dg/pcn,xnval? => xrhole,xrhoce,yrhole,yrhoce,ph
                vtotph, vvalph, xnval(1,iph), dgcn, dpcn, eref, adgc(1,1,iph), adpc(1,1,iph), xrhole(0,iph),           &
                xrhoce(0,iph),yrhole(1,0,iph),yrhoce(1,iph),ph(1,iph), iz(iph), xion(iph), iunf, itmp,lmaxsc(iph), iph) !KJ iph
		 !solves Dirac equation to obtain orbitals and phase shifts for FMS below
      enddo ! iph

!     transform neg,emg to em,ne,eref first
      em= dble(ee)
      eref=dble(eref)-coni*dimag(ee)


!c    call fms for a cluster around central atom
      gtr(:,0:nph)=cmplx(0)
      if (rfms1 .gt. 0) then
        if (lfms1 .ne. 0) then
          call fmsie( 0, nph, lmaxsc, ie,  em, eref, ph, iz, rfms1, lfms1, nat, iphat, rat, gtr,iscmt.eq.1)
        else
          do iph0 = 0, nph 
              call fmsie( iph0, nph, lmaxsc, ie, em, eref, ph, iz, rfms1, lfms1, nat, iphat, rat, gtr,iscmt.eq.1)
          enddo
        endif
      endif


      xntot = 0
      fl = 0
      fr = 0
      do iph = 0,nph
!       calculate density and integrated number of electrons in each
!       channel for each type of atoms density, etc., find xntot. 
        call ff2g (gtr(0,iph), iph, ie, nr05(iph), xrhoce, xrhole(0,iph), xrhocp, ee, ep,                               &
           yrhole(1,0,iph), yrhoce(1,iph), yrhocp(1,iph), rhoval(1,iph),&
           xnmues(0,iph), xnatph(iph), xntot, iflr, iflrp, fl, fr, iunf)
		   ! gtr, xrhoce etc. => rhoval, yrhoce, xrhoce, xntot, xnmues
           if (ie.le.negx) then
              edos(ie)=dble(ee)
              scfdos(ie,0:lx,iph)=dimag(xrhoce(0:lx,iph))
           endif
      enddo
      if (ie.ne.1) xndifp = xndif
      xndif = xntot - xnferm

!     decide on next energy point; there are nflrx floors, defined
!     by the magnitude of Im part. Each floor has its height and
!     horizontal step to search for Fermi level associated with it.
!     The driver below will decide whether to go left or right on
!     the current floor, go one floor up or down.

      if ((ie.lt.neg .and. ient.gt.1) .or. (ient.eq.1.and.ie.lt.nflrx)) then
         ep = ee
         ee = emg(ie+1)
         if (ie.eq.neg-1) then
!          reset iflr variables
           iflrp = 2
           iflr  = 1
         endif
         goto 25
      elseif (ient.eq.1 .and. ie.eq.nflrx) then
         upok = .false.
         idir = 1
         if (xntot.gt. xnferm) idir = -1
         ep = ee
         ee = ee + idir * step(iflr)
         goto 25
      elseif (ient.gt.1 .and. ie.eq.neg) then
         upok = .true.
         iflrp = 1
         iflr  = 1
         idir = -1
         if (xntot.lt. xnferm) idir = 1
         ep = ee
         ee = ee + idir * step(iflr)
         goto 25
      else
!       check if the fermi level is found
        if (iflrp.eq.1 .and. iflr.eq.1 .and. xndifp*xndif .le. 0.e0) then
!          Fermi level is found ; do not goto 25
           if (xndif.eq.0) then
              xmunew = dble(ee)
              a=0
           else
              a = xndif/(xndif-xndifp)
              do i = 1,4
                fxa = a*fl + (1-a)*fr
                bb = dimag((ep-ee)*(fr+fxa)/2 + coni*dimag(ee)*(fr-fl))
                xndif1 = xndif + a * bb
                a = a - xndif1 / bb
              enddo
              xmunew = dble((1-a)*ee+a*ep)
           endif

!          add end cap corrections to the configuration and density
!          factor 2 for spin degeneracy
           do iph = 0,nph
              do il = 0,lx
               if (il.le.2 .or. iunf.ne.0) then
                fl = xrhocp(il,iph) * 2
                fr = xrhoce(il,iph) * 2
                fxa = a*fl + (1-a)*fr
                bb = dimag((ep-ee)*(fr+fxa)/2 + coni*dimag(ee)*(fr-fl))
                xnmues(il,iph) = xnmues(il,iph) + a * bb
               endif
              enddo
              do ir = 1,nr05(iph)
                fl = yrhocp(ir,iph) * 2
                fr = yrhoce(ir,iph) * 2
                fxa = a*fl + (1-a)*fr
                bb = dimag((ep-ee)*(fr+fxa)/2 + coni*dimag(ee)*(fr-fl))
                rhoval(ir,iph) = rhoval(ir,iph) + a * bb
              enddo
             enddo
        else
!          continue search ; goto 25 eventually
           if (iflr.eq.iflrp) then
!            previous step was horizontal
             if (xndifp*xndif.le.0) then
!               need to step down
                upok =.false.
                iflrp = iflr
                iflr = iflr - 1
                ep = ee
                ee = dble(ee) + coni*4*step(iflr)
             elseif (abs(xndif).gt.10.d0*abs(xndif-xndifp)              &
     &          .and. upok) then
!               need to go up one floor since too far from fermi level
                iflrp = iflr
                if (iflr.lt.nflrx) then
                  iflr = iflr+1
                  ep = ee
                  ee = dble(ee) +  coni*4*step(iflr)
                else
                  ep = ee
                  ee = ee + idir* step(iflr)
                endif
             else
!               keep the same floor and direction
                ep = ee
                ee = ee + idir* step(iflr)
             endif
           else
!            previous step was up or down (vertical)
!            check the direction of search
             idir = -1
             if (xndif.lt.0) idir = 1
             iflrp = iflr
             ep = ee
             ee = ee + idir* step(iflr)
           endif
           goto 25
        endif
      endif
!     END of the loop over energy in complex plane.
!     new fermi level and densities are calculated.

!     report configuration; repeat iteration if found bad counts.
      ok = .true.
      !call wlog('  Electronic configuration')
      !call wlog('  type     l     N_el')
      write(29,*) '  Electronic configuration'
      write(29,*) '  type     l     N_el'
 310  format (2i6, f9.3)
      do 320 ip= 0,nph
      do 320 il = 0,lx
         write (slog,310) ip,il,xnmues(il,ip)
         !call wlog(slog)
         write(29,*) trim(slog)
!        check that occupation numbers are consistent with those set in getorb.f
         diff = abs(xnmues(il,ip) - xnvmu(il,ip))
         if (diff.gt.13.1 .or. (il.eq.2 .and. diff.gt. 9.1) .or.        &
     &   (il.eq.1 .and. diff.gt.5.1) .or. (il.eq.0 .and. diff.gt.1.95)) then
            call wlog (' Found bad counts.')
            write (slog,311) xnvmu(il,ip)
  311       format('  Occupation number in getorb is ', f9.3)
            call wlog(slog)
            call wlog ('  Will repeat this iteration. ')
            if (ient.gt.1) ok = .false.
         endif
 320  continue

!     if (.not. ok) then will restart SCF loop 
      if (ok) then
         xmu = xmunew
!        find rhoval via Broyden algorithm
         call broydn( iscmt, ca1, nph, xnvmu, nr05 , xnatph, rnrm, qnrm, edenvl, rhoval, dq) ! xnatph,qrnm,edenvl,rhoval => rhoval,dq

!        calculate new vclap - overlap coulomb potential
         call coulom (icoul, nph, nr05 , rhoval, edenvl, edens, nat, rat, iatph, iphat, rnrm, dq, iz, vclap) ! everything => vclap

!       update array edens
        do 350 ip=0,nph
           do 330 ir=1,nr05 (ip)
             edens(ir,ip)=edens(ir,ip)-edenvl(ir,ip)+rhoval(ir,ip)
  330      continue
           do 340 ir=nr05 (ip)+1,251
             edens(ir,ip)=0.0d0
             edenvl(ir,ip)=0.0d0
  340      continue
  350   continue
      endif

      ! Deallocate local variables
      deallocate(gtr, xrhoce, xrhocp)
      deallocate(xrhole, yrhole, ph)


      return
      end
