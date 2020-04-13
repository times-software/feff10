!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: szlz.f90,v $:
! $Revision: 1.7 $
! $Author: jorissen $
! $Date: 2012/02/17 07:39:12 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine szlz (ispin, ecv, nph, nat, rgrd, nohole, rfms2, lfms2,&
     &           lmaxph, edens, edenvl, dmag, vtot, vvalgs, rmt, rnrm,  &
     &           ixc, rhoint, vint, xmu, jumprm,                        &
     &           xnval, iorb, x0, dx, xion, iunf, iz,                   &
     &           adgc, adpc, dgc, dpc, ihole, rat, iphat, corr)


!     Finds new Fermi level (xmu), electron counts 

      use constants
      use DimsMod, only: nphx=>nphu, natx, nrptx, lx
      implicit double precision (a-h, o-z)

      integer ispin

!     input
      dimension dmagx(nrptx), dmag(251,0:nphx) !KJ 12-2011 removed "+1" from dmag(,0:nphx+1) - extra field not used and conflicts with declaration in other routines
      dimension vtot (251,0:nphx), vvalgs (251,0:nphx)
      dimension rmt(0:nphx),rnrm(0:nphx)
      dimension xnval (30,0:nphx), iorb (-4:3,0:nphx)
      dimension ri(nrptx)
      dimension iz(0:nphx), xion(0:nphx), lmaxph(0:nphx)
      dimension rat(3,natx),iphat(natx)
      real  rfms, rfms2
!     input and output
      dimension edens(251,0:nphx), edenvl(251,0:nphx)

!     work space

      dimension dum(nrptx), vtotph(nrptx),vvalph(nrptx)
      dimension dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1)
      dimension adgc(10,30,0:nphx+1), adpc(10,30,0:nphx+1)
      dimension dgcn(nrptx,30), dpcn(nrptx,30)
      complex*16 xrhoce(-4:3, -4:3, 0:nphx), xrhole(-4:3, -4:3, 0:nphx)

      integer iph
!     complex energy grid emg is decomposed into em and eref
      parameter (negx = 80)
      complex*16 emg(negx), em, eref, ee, ep, cchi, de
!     nflrx should be odd and defines the max of Im energy for
!     the countour 
      parameter (nflrx = 17)
      dimension step(nflrx)
      character*512 slog

      real,       allocatable :: amat(:,:,:,:,:), gctr(:,:,:,:,:)
      real*8,     allocatable :: xnmues(:,:,:)
      complex,    allocatable :: gtr(:,:,:,:,:)
      complex*16, allocatable :: fl(:,:,:), fr(:,:,:), ph(:,:)

      ! Allocate locale variables
      allocate(xnmues(3,0:lx,0:nphx))
      allocate(fl(3,0:lx,0:nphx), fr(3,0:lx,0:nphx))
      allocate(gtr(2,2, 3,0:lx, 0:nphx))
      allocate(amat(-lx:lx,2,2, 3,0:lx), gctr(2,2,3,0:lx,0:nphx))
      allocate(ph(lx+1, 0:nphx))

      call setkap(ihole, kinit, linit)

      if (ispin.eq.0) then
        write (slog,8)
   8    format('              N_l, N_j- and N_j+ calculation')
        write (slog,9)
   9    format('              ONLY central atom contribution! ')
      elseif (abs(ispin).le.1) then
        write (slog,10)
  10    format('              S_z, L_z and t_z calculation')
      else 
        write (slog,11)
  11    format('              S_z, N_l and N_j calculation')
      endif
      call wlog(slog)

      call wlog (' Calculating energy and space dependent l-DOS.')
      call wlog (' It takes time ...')

!     calculate energy independent matrix of angular coefficients
      call acoef(ispin, amat)

      call grids (ecv, xmu, negx, neg, emg, step, nflrx)

!     ie - is number of energy points calculated
      ie = 0
      ee = emg(1)
      ep = dble(ee)
      do 22 iph=0,nphx
      do 22 il=0,lx
      do 22 i=1,3
        xnmues(i, il,iph) = 0
  22  continue

!     Start the cycle over energy points (ie)
  25  continue
      ie = ie + 1

      if (ie.eq.1 .or. mod(ie,20).eq.0) then
         write(slog,30) ie, dble(ee)*hart
   30    format('     point # ', i3, '  energy = ', f7.3)
         call wlog(slog)
      endif

      do 100  iph = 0, nph

         call fixvar (rmt(iph),edens(1,iph),vtot(1,iph),dmag(1,iph),    &
     &                vint, rhoint, dx, rgrd, jumprm,                   &
     &                vjump, ri, vtotph, dum, dmagx)
         if (mod(ixc,10) .ge.5) then
            if (jumprm .gt. 0) jumprm = 2
            call fixvar (rmt(iph), edenvl(1,iph), vvalgs(1,iph),        &
     &                dmag(1,iph), vint, rhoint, dx, rgrd , jumprm,     &
     &                vjump, ri, vvalph, dum, dmagx)
            if (jumprm .gt. 0) jumprm = 1
         endif

         call fixdsx (iph, dx, rgrd , dgc, dpc, dgcn, dpcn)
        jri = (log(rmt(iph)) + x0) / rgrd + 2
        jri1 = jri+1
        eref = vtotph(jri1)
        do 40 i = 1, jri1
  40    vtotph(i) = vtotph(i) - eref
        if (ixc.ge.5) then
           do 50 i = 1, jri1
  50       vvalph(i) = vvalph(i) - eref
        else
           do 60 i = 1, jri1
  60       vvalph(i) = vtotph(i)
        endif

         itmp = 0
!        icount=1 for Renormalized atom counts
!        icount=2 for Mulliken counts
         icount = 0
         if (iph.eq.0 .and. nohole.lt.0) itmp = ihole
         if (icount.gt.0) then
            call rholat( icount, rgrd, x0, ri, ee,                      &
     &           ixc, rmt(iph), rnrm(iph),                              &
     &           vtotph, vvalph, xnval(1,iph), iorb(-4,iph),            &
     &           dgcn, dpcn, eref,                                      &
     &           adgc(1,1,iph), adpc(1,1,iph), xrhole(-4,-4,iph),       &
     &           xrhoce(-4,-4,iph), ph(1,iph),                          &
     &           iz(iph), xion(iph), iunf, itmp,3, iph) !KJ iph
         else
            call rholsz( rgrd, x0, ri, ee,                              &
     &           ixc, rmt(iph), rnrm(iph),                              &
     &           vtotph, vvalph, xnval(1,iph), dgcn, dpcn, eref,        &
     &           adgc(1,1,iph), adpc(1,1,iph), xrhole(-4,-4,iph),       &
     &           xrhoce(-4,-4,iph), ph(1,iph),                          &
     &           iz(iph), xion(iph), iunf, itmp,3, iph) !KJ iph
         endif
  100 continue

!     Write out phases for fmssz
!     transform neg,emg to em,ne,eref first
      em= dble(ee)
      eref=dble(eref)-coni*dimag(ee)

!c    call fms for a cluster around central atom
      do 195 iph0 = 0,nph
      do 195 il = 0, lx
      do 195 i = 1, 3
      do 195 i2= 1, 2
      do 195 i1= 1, 2
         gtr( i1,i2, i, il, iph0) = 0
         gctr(i1,i2, i, il, iph0) = 0
  195 continue

      rfms = 0
!     only central atom contribution for ispin = 0
!temp if (ispin.ne.0)  rfms = rfms2
      rfms = rfms2

      if (lfms2 .ne. 0) then
        iph0 = 0
        call fmssz( iph0, ie,  em, eref, ph, iz, nph,                   &
     &        rfms, lfms2, nat, iphat, rat, amat, lmaxph, gctr, gtr)
      else
        do 190 iph0 = 0, nph 
  190   call fmssz( iph0,  ie, em, eref, ph, iz, nph,                   &
     &        rfms, lfms2, nat, iphat, rat, amat, lmaxph, gctr, gtr)
      endif

      de = ee-ep
      do 300 iph = 0,nph
      do 300 lpp = 0,lx
      do 300 iop = 1,3
!       calculate density and integrated number of electrons in each
!       channel for each type of atoms density, etc.
        if (ie.gt.1) fl(iop,lpp,iph) = fr( iop,lpp,iph)
        fr( iop,lpp,iph) = 0
        call kfromi (1, lpp, j1, kk1)
        call kfromi (2, lpp, j1, kk2)
        do 200 i1=1,2
        do 200 i2=1,2
          call kfromi (i1, lpp, j1, k1)
          call kfromi (i2, lpp, j1, k2)
          if (k1.eq.0 .or. k2.eq.0) goto 200

          cchi =  dble( real( gtr(i1,i2, iop,lpp,iph) )) +              &
     &           coni* dble(aimag( gtr(i1,i2, iop,lpp,iph) ))
!         fr( iop,lpp,iph) = fr( iop,lpp,iph) + cchi * xrhole(k1,k2,iph)
!         use above kk1,kk1 for j- value, kk2,kk2 for j+ value
          if (ispin.ne.0 .or. iop.eq.1) then
            fr( iop,lpp,iph) = fr( iop,lpp,iph) + cchi*xrhole(k1,k2,iph)
          elseif(iop.eq.2) then
            fr( iop,lpp,iph) = fr( iop,lpp,iph)+cchi*xrhole(kk1,kk1,iph)
          elseif(iop.eq.3) then
            fr( iop,lpp,iph) = fr( iop,lpp,iph)+cchi*xrhole(kk2,kk2,iph)
          endif

!         add central atom part
          cchi =  dble(  gctr(i1,i2, iop,lpp,iph) ) 
          if (ispin.ne.0 .or. iop.eq.1) then
            fr( iop,lpp,iph) = fr( iop,lpp,iph) + cchi*xrhoce(k1,k2,iph)
!           use above k1,k1 for j- value, k2,k2 for j+ value
          elseif(iop.eq.2) then
            fr( iop,lpp,iph) = fr( iop,lpp,iph)+cchi*xrhoce(kk1,kk1,iph)
          elseif(iop.eq.3) then
            fr( iop,lpp,iph) = fr( iop,lpp,iph)+cchi*xrhoce(kk2,kk2,iph)
          endif
 200    continue

!       do integral over energy with trapezoidal rule
        if (ie.eq.1)  fl( iop,lpp,iph) = fr( iop,lpp,iph)
        xnmues(iop,lpp,iph) =  xnmues(iop,lpp,iph) +                    &
     &  dimag((fl(iop,lpp,iph) + fr(iop,lpp,iph)) * de /2)
        if (ie.eq.neg) then
!          end point correction
           xnmues(iop,lpp,iph) =  xnmues(iop,lpp,iph) +                 &
     &     dimag( fr(iop,lpp,iph) * (dble(ee)-ee) )
        endif

  300 continue

!     next energy point
      if (ie.lt.neg) then
	   ! write(*,*) 'did point ',ie,', energy ',emg(ie)
         ep = ee
         ee = emg(ie+1)
         goto 25
      endif

!     report configuration; repeat iteration if found bad counts.
      call wlog('  Electronic configuration')
      call wlog('  Electronic configuration:Mulliken counts')
      if (ispin.eq.0) then
         call wlog('   iph    il      N_l   N_j-  N_j+')
      elseif (abs(ispin).eq.1) then
         call wlog('   iph    il      S_z   L_z   T_z')
      else
         call wlog('   iph    il      S_z   N_l   N_j')
      endif
 310  format (2i6, 3f9.4)
      do 320 ip= 0,nph
      do 320 il = 0,lx
         write (slog,310) ip,il,(xnmues(i,il,ip), i=1,3)
         call wlog(slog)
 320  continue
      corr = 1.d0
      if (ispin.eq.0 .and. kinit.ne.-1) then
!       calculation  changes in counts due to spin-orbit interaction
        ip = 2
        if (kinit.lt.0) ip = 3
        il = linit + 1
        if (linit.eq.3) il = linit - 1
        corr = xnmues(1,il,0) /xnmues (ip,il, 0)
      endif

      ! Deallocate locale variables
      deallocate(xnmues,fl,fr,gtr,amat,gctr,ph)

      return
      end
