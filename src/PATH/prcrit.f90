!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: prcrit.f90,v $:
! $Revision: 1.6 $
! $Author: jorissen $
! $Date: 2012/12/11 23:20:30 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine prcrit (neout, nncrit, ik0out, cksp, fbeta, ckspc,     &
     &                   fbetac, potlb0, xlam, xlamc)
!KJ from ffmod4.f90 :         call prcrit (ne, nncrit, ik0, cksp, fbeta, ckspc,              &
!KJ     &                fbetac, potlbl, xlam, xlamc)

!     Prepare fbeta arrays, etc., for pathfinder criteria
!
!     Note that path finder is single precision, so be sure that
!     things are correct precision in calls and declarations!
!     See declarations below for details.
!     
!     Inputs:  Reads phase.bin
!     Output:  neout   'ne', number of energy grid points
!              ik0out  index of energy grid with k=0
!              cksp    |p| at each energy grid point in single precision
!              fbeta   |f(beta)| for each angle, npot, energy point, sp
!              ckspc   |p| at each necrit point in single precision
!              fbetac  |f(beta)| for each angle, npot, nncrit point, sp
!              potlb0  unique potential labels
!              xlam    mean free path for each energy point in Ang, sp
!              xlamc   mean free path for each nncrit point in Ang, sp

      use dimsmod, only: ltot, nex, nphx=>nphu, nspx=>nspu
	  use constants
      use global_inp,only: do_nrixs
	  use nrixs_inp,only: jinit,jmax,kfinmax !KJ 7-09 for NRIXS
      implicit none
      character*6  potlbl(0:nphx)

!     staff originally kept in common blocks of pdata.h
      double precision rnrmav, xmu, edge
      complex*16 ph( nex, ltot+1, 0:nphx), eref(nex), em(nex)
      integer lmax(nex,0:nphx), iz(0:nphx)

!     Output variables SINGLE PRECISION for use with path finder.
!     BE CAREFUL!!
      integer, parameter :: necrit=9, nbeta=40
      real fbetac(-nbeta:nbeta,0:nphx,necrit), ckspc(necrit)
      real fbeta(-nbeta:nbeta,0:nphx,nex), cksp(nex)
      real xlamc(necrit)
      real xlam(nex)
      character*6  potlb0(0:nphx)

!     Local variables
      complex*16 cfbeta, tl, cktmp
      real*8 dcosb(-nbeta:nbeta)
      real*8 pl(ltot+1)
      integer iecrit(necrit)
      real*8, parameter :: eps = 1.0e-16
      complex*16 eref2(nex,nspx)
      complex*16 ph4(nex, -ltot:ltot, nspx, 0:nphx)
	  real*8 krange(9) !KJ 7-09 NRIXS
	  complex*16, allocatable :: rkk(:,:,:) !KJ 7-09
	  integer ie,il,iii,ibeta,ik0out,neout,iph,lmaxp1,indmax,icrit,nncrit,ne,ne1,ne3,npot,ihole,ik0

!     Need stuff from phase.bin
!     Read phase calculation input

      if (do_nrixs .eq. 1) then  !NRIXS
	     allocate(rkk(nex,kfinmax,nspx))
         call rdxsphjas (ne, ne1, ne3, npot, ihole, &
           rnrmav, xmu, edge, ik0, em, eref2, iz, potlbl, ph4, rkk, lmax, lmaxp1)
	  else
		 allocate(rkk(nex,8,nspx))
         call rdxsph (ne, ne1, ne3, npot, ihole,    &
           rnrmav, xmu, edge, ik0, em, eref2, iz, potlbl, ph4, rkk, lmax, lmaxp1)
	  endif

      eref(1:ne)=eref2(1:ne,1)
      do iph = 0, npot
      do ie = 1, ne
      do il = 0, lmax(ie, iph)
      ph(ie,il+1, iph) = ph4(ie, -il, 1, iph)
      enddo
	  enddo
	  enddo

      neout = ne1
      ik0out = ik0
      potlb0(0:nphx)=potlbl(0:nphx)

!     |p| at each energy point (path finder uses invA, convert here)
!     Also make mfp (xlam) in Ang
      do ie = 1, ne
         cktmp = sqrt (2*(em(ie) - eref(ie)))
         cksp(ie) = dble (cktmp) / bohr
!        xlam code lifted from genfmt
         xlam(ie) = 1.0e10
         if (abs(dimag(cktmp)) .gt. eps) xlam(ie) = 1/dimag(cktmp)
         xlam(ie) = xlam(ie) * bohr
      enddo

!     Make the cos(beta)'s
!     Grid is from -40 to 40, 81 points from -1 to 1, spaced .025
      do ibeta = -nbeta, nbeta
         dcosb(ibeta) = 0.025 * ibeta
      enddo
!     watch out for round-off error
      dcosb(-nbeta) = -1
      dcosb(nbeta)  =  1

!     make fbeta (f(beta) for all energy points
      do ibeta = -nbeta, nbeta
         call cpl0 (dcosb(ibeta), pl, lmaxp1)
         do iii = 0, npot
            do ie = 1, ne
               cfbeta = 0
               do il = 1, lmax(ie,iii)+1
                  tl = (exp (2*coni*ph(ie,il,iii)) - 1) / (2*coni)
                  cfbeta = cfbeta + tl*pl(il)*(2*il-1)
               enddo
               fbeta(ibeta,iii,ie) = abs(cfbeta)
            enddo
         enddo
      enddo

!     Make similar arrays for only the icrit points

!     Use 9 points at k=0,1,2,3,4,6,8,10,12 invA
!     See phmesh for energy gid definition.  These seem to work fine, 
!     and results aren't too sensitive to choices of k.  As few as 4
!     points work well (used 0,3,6,9), but time penalty for 9 points
!     is small and increased safety seems to be worth it.
      iecrit(1) = ik0
      iecrit(2) = ik0 + 5
      iecrit(3) = ik0 + 10
      iecrit(4) = ik0 + 15
      iecrit(5) = ik0 + 20
      iecrit(6) = ik0 + 30
      iecrit(7) = ik0 + 34
      iecrit(8) = ik0 + 38
      iecrit(9) = ik0 + 40

!KJ 7-09 Aleksi had code here so that one could in principle use different iecrit and different energy grid.
!        However, he never implemented this seriously and the code in phmes***.f90 simply writes the same defaults
!        as specified above.  So I've just taken it out.
!        After all, Josh has implemented various methods to control the energy grid.  I'd rather stick with those
!        and not make the code any messier than it already is.  Sorry, Aleksi.
!      if (do_nrixs .eq. 1) then
!!		JAS changed so that we do not all have to use the same grid
!!		trying to set up iecrit so that it is done in phmeshjas.f  
!		open (unit=44, file='emesh.dat', status='unknown')
!		read(44,*)
!		read(44,*) (iecrit(ie),ie=1,9)
!		close(44)
!	  endif
!KJ

!     make sure that we have enough energy grid points to use all
!     9 iecrits
      nncrit = 0
      do ie = 1, necrit
         if (iecrit(ie) .gt. ne)  goto 295
         nncrit = ie
      enddo
  295 continue
      if (nncrit .eq. 0) call par_stop('bad nncrit in prcrit')
            
      do icrit = 1, nncrit
         ie = iecrit(icrit)
         ckspc(icrit) = cksp(ie)
         xlamc(icrit) = xlam(ie)
         do ibeta = -nbeta, nbeta
            do iii = 0, npot
               fbetac(ibeta,iii,icrit) = fbeta(ibeta,iii,ie)
            enddo
         enddo
      enddo

      deallocate(rkk)

      return
      end
