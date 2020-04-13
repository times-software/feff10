!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: scmtmp.f90,v $:
! $Revision: 1.13 $
! $Author: jorissen $
! $Date: 2012/10/23 20:08:40 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!KJ I have added a discussion at the bottom of this source file.  1-2012
      subroutine scmtmp (npr, iscmt, ecv, nph, nat, vclap,              &
     &                edens, edenvl, vtot, vvalgs, rmt, rnrm,qnrm,      &
     &                ixc, rhoint, vint, xmu, jumprm,                   &
     &                xnferm, xnvmu, xnval,                             &
     &                x0, ri, dx, xnatph, xion, iunf, iz,               &
     &                adgc, adpc, dgc,dpc, ihole,                       &
     &                rat,iatph,iphat, lmaxsc, rhoval, xnmues, ok,      &
     &                rgrd, nohole, nscmt, icoul, ca1, rfms1, lfms1,    &
     &                gtr, xrhole, xrhoce, yrhole, yrhoce )

!     Finds new Fermi level (xmu), electron counts (qnrm) 
!     and new valence densities (rhoval).
!     KJ 10-2012: added "if master" to most stdout commands to get rid of mess ...

      use DimsMod, only: nrptx, nphx=> nphu, natx, lx
	  use par
	  use constants
      implicit double precision (a-h, o-z)
      real*8 wall_commend, wall_commst

!     input
      dimension dmagx(nrptx), dmag0(251)
      dimension vclap(251,0:nphx)
      dimension vtot (251,0:nphx), vvalgs (251,0:nphx)
      dimension rmt(0:nphx),rnrm(0:nphx)
      dimension xnval (30,0:nphx)
      dimension qnrm(0:nphx), dq(0:nphx)
      dimension ri(nrptx), ri05(251), nr05(0:nphx)
      dimension xnatph(0:nphx), iz(0:nphx), xion(0:nphx)
      dimension rat(3,natx),iatph(0:nphx),iphat(natx), lmaxsc(0:nphx)
      real  rfms1
      
      real*8, intent(inout) :: xnvmu(0:lx,0:nphx+1)
      real*8, intent(inout) :: xnmues(0:lx, 0:nphx)
      complex*16, intent(inout) :: xrhoce(0:lx,0:nph,npr)
      complex*16, intent(inout) :: xrhole(0:lx,0:nph,npr)
      complex*16, intent(inout) :: yrhole(251,0:lx,0:nph,npr)
      complex, intent(inout) :: gtr(0:lx, 0:nph, npr)

!     input and output
      dimension edens(251,0:nphx), edenvl(251,0:nphx)
      dimension rhoval(251,0:nphx+1)

!     work space

      dimension dum(nrptx), vtotph(nrptx),vvalph(nrptx)
      dimension dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1)
      dimension adgc(10,30,0:nphx+1), adpc(10,30,0:nphx+1)
      dimension dgcn(nrptx,30), dpcn(nrptx,30)

      complex*16 yrhocp(251,0:nphx)

!     special dimension for MPI
!     Maxprocs = max number of processors for parallel execution
!--   This is set in parallel.h
!     npr - actual number of processors is passed to this subroutine

      complex*16 yrhoce(251,0:nph,npr)

      integer iph
!     complex energy grid emg is decomposed into em and eref
      parameter (negx = 80)
      complex*16 emg(negx), em, eref, ee, ep, fl, fr, fxa
!     nflrx should be odd and defines the max of Im energy for
!     the countour 
      parameter (nflrx = 17)
      dimension step(nflrx)
!     stuff from feff.f for rdinp, pathfinder and genfmt
      logical wnstar, ok
!     Following passed to pathfinder, which is single precision.
      character*512 slog

      complex*16, allocatable :: ph(:,:)
      complex*16, allocatable :: xrhocp(:,:)

      integer ient, ixl(1), ixlc(1), ixly(1)
      data ient /0/

!     save stuff from rdinp, so no need to call it again
      save   ri05, ient

      ! Allocate local variables
      allocate(ph(lx+1, 0:nphx), xrhocp(0:lx,0:nphx))

      ient = ient + 1  ! Counts the number of times scmtmp is called, starting at 1
      if (ient.eq.1) then
         ! The initial guess at xmu is typically positive while the actualy value should be negative. In this case, this MP algorithm is quite slow at finding the correct Fermi level. So, here we pick an arbitrary negative value that is typically closer to the actual Fermi level.
         xmu = -0.25d0
         do 15 i= 1,251
  15     ri05(i) = exp (-8.8+0.05*(i-1))
      endif

      if (master) write (slog,10) iscmt
!  10  format('              SCF ITERATION NUMBER',i3)
  10  format('SCF ITERATION NUMBER',i3)
      call wlog(slog)
      write(29,*) trim(slog)

      !call wlog (' Calculating energy and space dependent l-DOS ...')

!     initialize new valence density
      rhoval(:,:) = 0

      call grids (ecv, xmu, negx, neg, emg, step, nflrx)

!     ie - is number of energy points calculated
      ietot0 = 1
      ee = emg(1)
      ep = dble(ee)

      xrhoce(:,:,:) = 0
      xnmues(:,:) = 0
      yrhoce(:,:,:) = 0
      iflr = nflrx
      iflrp = nflrx

      nproc = npr
      n1 = 1
      n2 = min(neg, nproc)

!     Start the cycle over energy points (ie)
  25  continue  ! This is really a loop over sets of energy points (i.e. every time you come here, a batch of energy points goes to the processor grid)
                ! The total energy grid consists of an initial grid of neg points, to which points are added in grids of nproc points until the Fermi energy is found. ("current grids")
				! Every time we pass this point here, a batch of min(nproc,points left in current grid) is sent to the nproc processors.
				! The arrays only hold the quantities calculated for the current grid of energy points.
				! However, the integrated quantities "remember" the result of all the previous sets and add/subtract the contribution from the current set.
				! Once the first grid of neg points is all done, we are allowed to check for the presence of the Fermi level.  If not found, we add a new grid of nproc points.

!     slow loop for MPI execution
      ie = this_process + n1    ! every processor calculates one energy point
      ietot = ietot0 + this_process
      if (ie .gt. n2) then     ! meaning there's no energy point available for the current processor (so don't calculate anything)
							   ! careful : n2 can change from set to set.
           !KJ 1-2012 If we don't set the grid properly, things will go wrong in the call to ff2g way below (bugfix for nproc > neg)
           do iph=0,nph
              nr05(iph)=(log(rnrm(iph)) + x0) / 0.05d0 + 5
              if (nr05(iph) .gt. 251) nr05(iph)=251
           enddo
           go to 200
      endif

        ipr = 1 + ie - n1  ! "ie" is an index in the current energy grid (>= 1 sets) whereas ipr is an index in the current set
        if (worker) par_type = 3

        if (ietot.eq.1 .or. mod(ietot,20).eq.0) then
           write(slog,30) ietot, dble(emg(ie))*hart
   30      format('     point # ', i3, '  energy = ', f7.3)
           call wlog(slog) !here, exceptionally, the "worker" is allowed to talk to the user
        endif

        do iph = 0, nph
          dmag0(:) = 0.d0
!c        use spin-unpolarized case to get SCF. set dmagx to zero
!c        may want to replace dmag0 with dmag(1,iph) for spin-dependent
!c        extension of SCF procedure.
          call fixvar (rmt(iph),edens(1,iph),vtot(1,iph),dmag0,vint, rhoint, dx, rgrd, jumprm, vjump, ri, vtotph, dum, dmagx)
          if (mod(ixc,10) .ge.5) then
            if (jumprm .gt. 0) jumprm = 2
            call fixvar (rmt(iph),edenvl(1,iph),vvalgs(1,iph),dmag0, vint, rhoint, dx, rgrd , jumprm, vjump, ri, vvalph, dum, dmagx)
            if (jumprm .gt. 0) jumprm = 1
          endif

          call fixdsx (iph, dx, rgrd , dgc, dpc, dgcn, dpcn)
          jri = (log(rmt(iph)) + x0) / rgrd + 2
          jri1 = jri+1
          eref = vtotph(jri1)
          do 40 i = 1, jri1
  40      vtotph(i) = vtotph(i) - eref
          if (ixc.ge.5) then
            do 50 i = 1, jri1
  50        vvalph(i) = vvalph(i) - eref
          else
            do 60 i = 1, jri1
  60        vvalph(i) = vtotph(i)
          endif

           itmp = 0
           if (iph.eq.0 .and. nohole.lt.0) itmp = ihole
           call rholie( ri05, nr05(iph), rgrd, x0, ri, emg(ie), ixc,    &  ! from potential, solve and calculate densities and phases
                 rmt(iph), rnrm(iph), vtotph, vvalph, xnval(1,iph),     &
                 dgcn, dpcn, eref, adgc(1,1,iph), adpc(1,1,iph),        &
                 xrhole(0,iph,ipr), xrhoce(0,iph,ipr),                  &
                 yrhole(1,0,iph,ipr), yrhoce(1,iph,ipr),                &
                 ph(1,iph), iz(iph), xion(iph), iunf, itmp,lmaxsc(iph), iph) !KJ iph
        enddo   ! iph=0,nph

!       Write out phases for fmsie
!       transform neg,emg to em,ne,eref first
        em= dble(emg(ie))
        eref=dble(eref)-coni*dimag(emg(ie))

!c      call fms for a cluster around central atom
        gtr(:,:,ipr) = 0
        if (rfms1 .gt. 0) then
          if (lfms1 .ne. 0) then
            iph0 = 0
!           set logic to call yprep on every processor
            lfms = lfms1
            if (ietot0.eq.1) lfms = 2
            call fmsie( iph0, nph, lmaxsc, ietot, em, eref, ph, iz, rfms1, lfms, nat, iphat, rat, gtr(0,0,ipr),iscmt.eq.1)  ! calculate Green's function from phases
          else
            do 190 iph0 = 0, nph 
  190       call fmsie( iph0, nph, lmaxsc, ietot, em, eref, ph, iz, rfms1, lfms1, nat, iphat, rat, gtr(0,0,ipr),iscmt.eq.1)
          endif
        endif
  200 continue
!     end of slow loop for MPI execution

      ietot0 = ietot0 + n2 - n1 + 1  ! count total number of energy points.
      if (worker) par_type = 2

      ixl(1) = (lx + 1) * (nph + 1)
      ixly(1) = ixl(1) * 251
      ixlc(1) = (nph + 1) * 251
      if (nproc .gt. 1) then  ! broadcast arrays to/from all nodes now (including the ones that didn't calculate anything)
        call seconds(wall_commst)
        if (worker .and. (ie .le. n2)) then
!-- Send pointers for gtr buffer to master
          call par_send_int(ixl,1,0,this_process)
          call par_send_int(ixly,1,0,this_process)
          call par_send_int(ixlc,1,0,this_process)
!-- Send buffer
          if (ixl(1) .ne. 0) then
            call par_send_cmplx(gtr(0,0,ipr),ixl(1),0,this_process)
            call par_send_dc(xrhoce(0,0,ipr),ixl(1), 0, this_process)
            call par_send_dc(xrhole(0,0,ipr),ixl(1), 0, this_process)
          endif
          if (ixly(1) .ne. 0)  call par_send_dc(yrhole(1,0,0,ipr),ixly(1), 0, this_process)
          if (ixlc(1) .ne. 0)  call par_send_dc(yrhoce(1,0,ipr),ixlc(1), 0, this_process)
        else if (master) then
          do i = 1,n2-n1
!-- Receive pointers for gtr buffer from i
            call par_recv_int(ixl,1,i,i)
            call par_recv_int(ixly,1,i,i)
            call par_recv_int(ixlc,1,i,i)
!-- Receive buffer from i
            if (ixl(1) .ne. 0) then
              call par_recv_cmplx(gtr(0,0,i+1),ixl(1),i,i)
              call par_recv_dc(xrhoce(0,0,i+1),ixl(1),i,i)
              call par_recv_dc(xrhole(0,0,i+1),ixl(1),i,i)
            endif
            if (ixly(1) .ne. 0)  call par_recv_dc(yrhole(1,0,0,i+1),ixly(1),i,i)
            if (ixlc(1) .ne. 0)  call par_recv_dc(yrhoce(1,0,i+1),ixlc(1),i,i)
          enddo
        endif
!-- Broadcast gtr
        ilen = ixl(1) * (n2 - n1 + 1)
        ileny = ilen * 251
        ilenc = (nph + 1) * (n2 - n1 + 1) * 251
        call par_bcast_cmplx(gtr(0,0,1),ilen,0)
        call par_bcast_dc(xrhoce(0,0,1),ilen,0)
        call par_bcast_dc(xrhole(0,0,1),ilen,0)
        call par_bcast_dc(yrhole(1,0,0,1),ileny,0)
        call par_bcast_dc(yrhoce(1,0,1),ilenc,0)
        call seconds(wall_commend)
        wall_comm = wall_comm + wall_commend - wall_commst
      endif
      
	  ! Now all nodes have the complete gtr, xrhoce, xrhole, yrhole, yrhoce.  I think only min(nproc,neg) points have been calculated the first time we get here. 
		 
!     fast loop (does not need parallel execution)
!     uses results of the above loop to find Fermi level and to decide on next set of energy points
!     every node does the same work here
      do 300 ie = n1, n2  !Sum over energy points in current set
        ipr = 1+ ie -n1   !matching index in the processor grid
        ee = emg(ie)

        if (ie.eq.1 .and. iflrp.ne.1) then
!         the absolutely first point on energy grid - get the integral started
          xrhocp(0:lx,0:nph) = xrhoce(0:lx,0:nph, ipr)
          yrhocp(1:251,0:nph) = yrhoce(1:251,0:nph, ipr)
        endif

        xntot = 0
        if (ie.eq.neg .and. iflrp.gt.1) iflr = 1  !Means we're done calculating the first grid of neg points.
        fl = 0
        fr = 0
        do iph = 0,nph
!         Calculate density and integrated number of electrons in each channel for each type of atoms density, etc.  Find total charge xntot. 
          call ff2g (gtr(0,iph,ipr), iph,ie, nr05(iph), xrhoce(0,0,ipr), xrhole(0,iph,ipr), xrhocp, ee, ep, yrhole(1,0,iph,ipr),     &
            yrhoce(1,iph,ipr),yrhocp(1,iph),rhoval(1,iph), xnmues(0,iph), xnatph(iph), xntot, iflr, iflrp, fl, fr,iunf)
        enddo

!       check whether Fermi level is found between points n1 and n2, and decide on next set of energy points;
        if (ie.ne.1 .or. iflrp.eq.1) xndifp = xndif
        xndif = xntot - xnferm

!       check if the fermi level is found
        if ( iflr.eq.1) then ! Can't do this test until the first grid of neg points is calculated completely.
          if (xndifp*xndif .le. 0.e0) then
!         Fermi level is found ; exit from energy loop
             if (xndif.eq.0) then
               xmunew = dble(emg(ie))
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

!            add end cap corrections to the configuration and density
!            factor 2 for spin degeneracy
             do iph = 0,nph
               do il = 0,lx
                if (il.le.2 .or. iunf.ne.0) then
                 fl = xrhocp(il,iph) * 2
                 fr = xrhoce(il,iph,ipr) * 2
                 fxa = a*fl + (1-a)*fr
                 bb = dimag((ep-ee)*(fr+fxa)/2 + coni*dimag(ee)*(fr-fl))
                 xnmues(il,iph) = xnmues(il,iph) + a * bb
                endif
               enddo
               do ir = 1,nr05(iph)
                 fl = yrhocp(ir,iph) * 2
                 fr = yrhoce(ir,iph,ipr) * 2
                 fxa = a*fl + (1-a)*fr
                 bb = dimag((ep-ee)*(fr+fxa)/2 + coni*dimag(ee)*(fr-fl))
                 rhoval(ir,iph) = rhoval(ir,iph) + a * bb
               enddo
             enddo !iph

!            exit from the energy loop since Fermi level is found
             goto 305
          endif
        endif
        ep = emg(ie)
        xrhocp(0:lx,0:nph) = xrhoce(0:lx,0:nph, ipr)
        yrhocp(1:251,0:nph) = yrhoce(1:251,0:nph, ipr)

 300  continue

      ! Loop over energy points finished ; Fermi level not found.
	  ! Prepare the next set of energy points.
      if (n2.lt.neg .and. iflrp.gt.1) then
	    ! We're in the first grid and   nproc < neg : let's continue calculating the first grid.
        n1 = n2+1
        n2 = min(neg, n2+nproc)
      else
	    ! Exhausted all points in the current grid ; set up a new grid for the next set of energy points
		! In this new grid, make sure to use all processors, i.e. neg = nproc.  Might as well.
		! The new grid starts from the last point of the previous grid.  It searches on nproc points step(1) apart either above or below last energy point.  
        iflr = 1 !Means from now on, we can evaluate for the Fermi level in every set of energy points
        iflrp = 1 !Means we're in at least the second grid of energy points
        idir = -1  !       set direction of search
        if (xndif.lt.0) idir = 1
        n1 = 1
        n2 = min(nproc, negx)
        do ie = n1, n2
           emg(ie) = ep+ idir*step(iflr) * ie
        enddo
      endif
      goto 25  ! Start the whole process over again with a new set of energy points.

!     END of the loop over energy in complex plane.  Fermi level found and densities are calculated.
 305  continue

!     report configuration; repeat iteration if found bad counts.
      ok = .true.
      if(master) write(29,*) '  Electronic configuration'
      if(master) write(29,*) '  type     l     N_el'
 310  format (2i6, f9.3)
      do 320 ip= 0,nph
      do 320 il = 0,lx
         write (slog,310) ip,il,xnmues(il,ip)
         write(29,*) trim(slog)
!        check that occupation numbers are consistent with those
!        set in getorb.f
         diff = abs(xnmues(il,ip) - xnvmu(il,ip))
         if (diff.gt.13.1 .or. (il.eq.2 .and. diff.gt. 9.1) .or.        &
               (il.eq.1 .and. diff.gt.5.1) .or.                               &
               (il.eq.0 .and. diff.gt.1.95)) then
            call wlog (' Found bad counts.')
            write (slog,311) xnvmu(il,ip)
  311       format('  Occupation number in getorb is ', f9.3)
            call wlog(slog)
            call wlog ('  Will repeat this iteration ')
            if (ient.gt.1) ok = .false.
         endif
 320  continue

!     if (.not. ok) then will restart SCF loop in potsub
      if (ok) then
         xmu = xmunew
!        find rhoval via Broyden algorithm (valence density)
         call broydn( iscmt, ca1, nph, xnvmu, nr05 , xnatph, rnrm, qnrm, edenvl, rhoval, dq)

!        calculate new vclap - overlap coulomb potential
         call coulom (icoul, nph, nr05 , rhoval, edenvl, edens, nat, rat, iatph, iphat, rnrm, dq, iz, vclap)

!       update array edens
        do ip=0,nph
           do ir=1,nr05 (ip)
             edens(ir,ip)=edens(ir,ip)-edenvl(ir,ip)+rhoval(ir,ip)
           enddo
           do ir=nr05 (ip)+1,251
             edens(ir,ip)=0.0d0
             edenvl(ir,ip)=0.0d0
           enddo
        enddo
      endif

      ! Deallocate local variables
      deallocate(ph, xrhocp)

      return
      end
	  
	  
! NOTES ON THE MPI SCF LOOP
! Kevin Jorissen, Jan 2012, feff 9.5.1
! Tested on feff90/examples/DANES/GeCl_4
! mpirun -np 58 -host n01,n02,n03,….  ~/feff90/bin/MPI/pot
! Using output in .scfconvergence-feff (and, while debugging, various debugging outputs not present in the release version)
! 
! 
! * I have fixed the bug that caused calculations with many processors (nproc > neg) to hang indefinitely.  (Problem : radial grid in scmtmp not initialized -> integrals zero -> surplus nodes didn't realize the calculation had converged and launched into new iterations by themselves!)
! 
! * I have added loads of comments to scmtmp.f90, so that I can now understand the flow of that routine almost perfectly.
! 
! * There's this "nflr / step(:)" thing which is somewhat misleading : only 2 "ladders" are actually used ; the second one is used as many times as needed.
! 
! * Results are dependent on number of processors used.  E.g. for GeCl4, the SCF loop gives Fermi level 
! -5.354 eV  (nproc = 1, 4, 14, 58)
! -5.356 eV  (nproc = 59, 64)
! -5.405 eV  (nproc = 80,128, 256)
! where it is interesting to note that neg=58, and negx=80 for all calculations.
! 
! The convergence criteria are 0.001 charge and 0.003 Ha on the Fermi level, which equals 0.08 eV.  So some of the differences above are within that criterium.  In all cases, it is the charge convergence which is the stricter criterium (takes longer to converge).
! 
! * The parameter negx sets a limit on the number of concurrent processes in scmtmp.f90.  Since it is currently by default set to negx=80, any number of processors above 80 will NOT be used unless the user recompiles.  We may want to revisit that number?
! Note that negx also occurs in grids, where it can have some effect on the energy grid emg.
! 
! * The influence of negx for nproc > negx is illustrated here.  When changing negx=300, the above results change to  
! -5.354 eV  (nproc = 1, 4, 14, 58) (negx=80)
! -5.356 eV  (nproc = 59, 64) (negx=80)
! -5.405 eV  (nproc = 80) (negx=80 or negx=300)
! -5.352 eV  (nproc = 128, 256) (negx=300)
! which brings the results closer and still within tolerance (0.08 eV).  I believe the variation of results with negx is not a bug (since for nproc=80 there's no difference when varying negx) but simply due to the variation in the energy grid.  In this case, the results with higher negx seem "better", but it's hard to substantiate that.
! 
! * While nproc < negx, the energy grids are the same for all nproc (except nproc=1, which is processed in the scmt.f90 routine).  The grid is simply processed in different-sized chunks (4, 14, 58, 59).  I think this holds for nproc > negx also (since the calculation simply limits nproc anyway!), but I haven't tested.
! 
! 
! * I feel complete with this investigation, except for the following question : 
! 
! - In what circumstances would it be useful to parallellize feff over a larger number of processors?
!   --- k-space : implement loop over k-points
!   --- fms module : use a processor for each energy point, if possible
	  
