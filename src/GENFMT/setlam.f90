!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: setlam.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setlam (icalc, ie)

      use dimsmod, only: lamtot, mtot, ntot
	  use constants
      use lambda
      use pdata
      implicit double precision (a-h, o-z)

!     Set lambda array based on icalc and ie
!     icalc  what to do
!      0     i0, ss exact
!      1     i1, ss exact
!      2     i2, ss exact
!     10     cute algorithm
!     <0     do exactly as told, decode as:
!               icalc = -(nmax + 100*mmax + 10 000*(iord+1))
!               Note that iord=0 <=> nmax=mmax=0, so use
!                  icalc = -10 000 for this case.
!               iord = 2*nmax + mmax, so if you want iord to control,
!               set nmax and mmax large enough-- if you want nmax and
!               mmax to control, set iord = 2*nmax + mmax...

!     inputs: ie used for cute algorithm
!             nsc used from /pdata/ to recognize ss paths
!     output: variables in /lambda/ set

      dimension mlam0(lamtot), nlam0(lamtot)

!     one degree in radians
      parameter (onedeg = .01745329252)
      character*512 slog

!     Set iord, nmax and mmax based on icalc
      if (icalc .lt. 0)  then
!        decode it and do what user wants
         icode = -icalc
         nmax = mod(icode,100)
         mmax = mod(icode,10000)/100
         iord = icode/10000 -1
      elseif (nsc .eq. 1)  then
         mmax = ilinit
         nmax = ilinit
         iord = 2*nmax + mmax
      elseif (icalc .lt. 10)  then
         iord = icalc
         mmax = iord
         nmax = iord/2
      elseif (icalc .eq. 10)  then
!        do cute algorithm
!        set mmax = L0 if straight line path, otherwise set mmax = 3
         mmax = ilinit
         do 10  ileg = 1, nleg
            mag1 = abs(beta(ileg))
            mag2 = abs(mag1 - pi)
!           if beta is not 0 or pi, path is non-linear
            if (mag1.gt.onedeg .and. mag2.gt.onedeg) mmax = 3
   10    continue
!        Set nmax based on ie and l0.
!        k <= 12 invA (ie=41)  nmax = L0
!        k >= 13 invA (ie=42)  nmax =  9
         nmax = ilinit
         if (ie .ge. 42)  nmax = 9
         iord = 2*nmax + mmax
      else
         write(slog,'(a,i8)') ' undefined icalc ', icalc
         call wlog(slog)
         call par_stop('setlam')
      endif

!-----construct index lambda (lam), (mu, nu) = mlam(lam), nlam(lam)
!     lamtot, ntot, mtot are maximum lambda, mu and nu to consider
!     Use ...0 for making indices, then sort into arrays with no
!     trailing 0 so laml0x is minimimized. (note: this is a crude
!     n**2 sort -- can 'improve' to nlog_2(n) if necessary)
      lam = 0
      do 20 in = 1, nmax+1
         n = in - 1
         do 20  im = 1, mmax+1
            m = im-1
            jord = 2*n+m
            if (jord .gt. iord)  goto 20
            if (lam .ge. lamtot)  then
               call wlog(' Lambda array filled, some order lost')
               goto 21
            endif
            lam = lam+1
            mlam0(lam) = -m
            nlam0(lam) = n
            if (m .eq. 0)  goto 20
            if (lam .ge. lamtot)  then
               call wlog(' Lambda array filled, some order lost')
               goto 21
            endif
            lam = lam+1
            mlam0(lam) = m
            nlam0(lam) = n
   20 continue
   21 continue
      lamx=lam
!     lamx must be less than lamtot
      if (lamx .gt. lamtot) call par_stop('SETLAM lamx > lamtot')

!     laml0x is biggest lam for non-zero fmatrix, also set mmax and nmax
!     Sort mlam0 and nlam0 to use min possible laml0x
      lam = 0
      do 30  lam0 = 1, lamx
         if ((nlam0(lam0).le.ilinit) .and.                              &
     &       (iabs(mlam0(lam0)).le.ilinit)) then
            lam = lam+1
            nlam(lam) = nlam0(lam0)
            mlam(lam) = mlam0(lam0)
            nlam0(lam0) = -1
         endif
   30 continue
      laml0x = lam
      do 40  lam0 = 1, lamx
         if (nlam0(lam0) .ge. 0)  then
            lam = lam+1
            nlam(lam) = nlam0(lam0)
            mlam(lam) = mlam0(lam0)
         endif
   40 continue

      mmaxp1 = 0
      nmax = 0
      do 50  lam = 1, lamx
         if (mlam(lam)+1 .gt. mmaxp1)  mmaxp1 = mlam(lam)+1
         if (nlam(lam) .gt. nmax)  nmax = nlam(lam)
   50 continue

      if (nmax.gt.ntot .or. mmaxp1.gt.mtot+1)  then
   52    format (a, 4i8)
         write(slog,52) ' mmaxp1, nmax, mtot, ntot ',                   &
     &                    mmaxp1, nmax, mtot, ntot
         call wlog(slog)
         write(slog,52) ' icalc ', icalc
         call wlog(slog)
         call par_stop('setlam')
      endif

      return
      end
