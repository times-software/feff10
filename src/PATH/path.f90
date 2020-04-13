!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: path.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2012/05/15 21:29:59 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     sub-program exchange
      program  ffmod4
!     subroutine ffmod4

!     makes paths list using cluster geometry and phase shifts
!     written by a.ankudinov 2000 using earlier subroutines
!     written by s.zabinsky
!     modified by a.ankudinov 2001 for new I/O structure

!     INPUT FILES
!       global.dat, geom.dat - global infomation file is read here 
!       mod4.inp - specific information for present module
!       phase.bin - output of XSPH module is read using subroutine 
!                  'rdxsph' inside subroutine 'prcrit'.
!                   needed  data: (list of variables)
!                  (ne, ne1, npot, ik0, em, eref2, potlbl, ph4)
!     OUTPUT FILE
!       paths.dat - list of filtered paths

      use dimsmod, only: nphx=>nphu, nex, nspx=>nspu, init_dimensions
	  use constants
	  use par
	  use eels_inp
	  use paths_inp
	  use global_inp
	  use atoms_inp
	  use errorfile

	  implicit none !double precision (a-h, o-z)
      character*30 fname
	  integer ios,ne,ik0,iph,ie,ibeta
	  real*8 cosb,angle
	  character*6, allocatable ::  potlbl(:)
!     Following passed to pathfinder, which is single precision.
!     Be careful to always declare these!
      integer, parameter :: necrit=9, nbeta=40
      real, allocatable :: fbetac(:,:,:),  fbeta(:,:,:), cksp(:), xlam(:)
      real xlamc(necrit), ckspc(necrit)

      call par_begin
      if (worker) go to 400
	  call OpenErrorfileAtLaunch('path')
      call init_dimensions

      
!     open the log file, unit 11.  See subroutine wlog.
      open (unit=11, file='log4.dat', status='unknown', iostat=ios)
      call chopen (ios, 'log4.dat', 'feff')

!     INPUT: read geom.dat, global.dat, mod4.inp
      call repath

	  allocate(potlbl(0:nphx), fbetac(-nbeta:nbeta,0:nphx,necrit),  &
           fbeta(-nbeta:nbeta,0:nphx,nex), cksp(nex), xlam(nex))


      if (nspx.gt.1) ispin = abs(ispin)

      if (ms.eq.1  .and.  mpath.eq.1)  then
         call wlog('Pathfinder: finding scattering paths...')
         call wlog('Preparing plane wave scattering amplitudes')
         call prcrit (ne, nncrit, ik0, cksp, fbeta, ckspc, fbetac, potlbl, xlam, xlamc)

!        Dump out fbetac for central atom and first pot
         if (ipr4 .ge. 3 .and. ipr4.ne.5)  then
            do iph = 0, 1
               do ie = 1, nncrit
                  write(fname,200)  ie, iph
  200             format ('fbeta', i1, 'p', i1, '.dat')
                  open (unit=1, file=fname, status='unknown')
                  write(1,210)  iph, ie, ckspc(ie)
  210             format ('# iph, ie, ckspc(ie) ', 2i5, 1pe20.6, /      &
     &                    '#  angle(degrees), fbeta/|p|,  fbeta')
                  do ibeta = -nbeta, nbeta
                     cosb = .025 * ibeta
                     if (cosb .gt.  1)  cosb =  1
                     if (cosb .lt. -1)  cosb = -1
                     angle = acos (cosb)
                     write(1,230)  angle*raddeg, fbetac(ibeta,iph,ie)/ckspc(ie), fbetac(ibeta,iph,ie)
  230                format (f10.4, 1p, 2e15.6)
                  enddo
                  close (unit=1)
               enddo
            enddo
         endif

         call wlog('Searching for paths')
         call paths (ckspc, fbetac, xlamc, pcritk, pcrith, critpw, nncrit, rmax, nlegxx, rfms2, nat, rat, iphat, ibounc) 
         call wlog('Eliminating path degeneracies')
         call pathsd (ckspc, fbetac, xlamc, ne, ik0, cksp, fbeta, xlam, critpw, ipr4, nncrit, potlbl, ipol, ispin, evec, xivec)
         call wlog('Done with module: pathfinder.'//char(13)//char(10))
      endif
      close (unit=11)

  400 call par_barrier
      call par_end
      if(master)call WipeErrorfileAtFinish
!     sub-program exchange
      stop
!     return

      end
