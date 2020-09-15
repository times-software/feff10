!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: getmat.f90,v $:
! $Revision: 1.3 $
! $Author: jorissen $
! $Date: 2012/06/29 01:05:24 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getmat(ihole, lin, nlp1, nlm1, jinit, minit, kinit,    &
     &     jfin, mfin, kfin, ncore, nph, matsize, kappa, xnval, ibasis)
!   input:
!     lin - orbital momentum of initial state
!     nlp1 - number of (lin+1) orbitals in basis set
!     nlm1 - number of (lin-1) orbitals in basis set
!    output:
!     matsize - matrix size, and for each matrix index
!     jinit, jfin - initial and final j (*2)
!     minit, mfin - initial and final m (*2) mfin - minit = 1 only 
!     kinit, kfin - initial and final kappa
!     nph - index of the orbital in the basis set (dgcn and dgcnp)
!     array nph is initialized for the use with basis set, but
!     it can be overwritten later (see example in subroutine getchi0)

      implicit double precision (a-h, o-z)
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
      integer,parameter :: maxx = 78
      integer,intent(in) :: lin,kappa
      integer,intent(inout) :: nlp1,nlm1
      integer,intent(out) :: jinit,jfin,minit,mfin,kinit,kfin,nph,matsize
      dimension jinit(maxx), minit(maxx), jfin(maxx), mfin(maxx)
      dimension kinit(maxx), kfin(maxx), nph(maxx), ncore(maxx)
      dimension kappa(41)
      real*8 xnval(41)

      np = 3 * (2*lin + 1)
      nm = 3 * (2*lin - 1)
      if (nm.lt.0) nm = 0
!     there should be at least one lin+1 orbital
      if (nlp1.lt.1) nlp1=1
      if (nlm1.lt.0) nlm1 = 0
      matsize = nlp1*np + nlm1*nm
!     check whether dimension is big enough
      if (matsize.gt.maxx) stop 'Increase dimension maxx in getmat'
!     initialize to zero all arrays
       jinit(:) = 0
       minit(:) = 0
       kinit(:) = 0
       jfin(:)  = 0
       mfin(:)  = 0
       kfin(:)  = 0
       nph(:)  = 0

!     i is current matrix index
      i = 0
!     cycle over  lin+1 basis set orbitals
      do 50 j = 0, nlp1-1
!       cycle over 3 dipole allowed lin-->lin+1 transitions
        do 40 idip = 1,3
          if (idip.eq.1) then
            kin = lin
            kfi = lin + 1
          elseif (idip.eq.2) then
            kin = -lin-1
            kfi = lin+1
          else
            kin = -lin-1
            kfi = -lin-2
          endif
          if (kin.eq.0) goto 40

          jin = 2*abs(kin) - 1
          jfi = 2*abs(kfi) - 1
          nbs = -2*j-1
          if (kfi.gt.0) nbs = nbs-1

!         cycle over minit
          do 30 mini = -jin, jin, 2
            mfi = mini + 2
            if (mfi.lt.-jfi .or. mfi.gt.jfi) goto 30

!           set variables for matrix index i
            i = i + 1
            jinit(i) = jin
            minit(i) = mini
            kinit(i) = kin
            jfin(i)  = jfi
            mfin(i)  = mfi
            kfin(i)  = kfi
            nph(i)   = nbs
  30      continue
  40    continue
  50  continue

!     add lin-->lin-1 transitions
!     cycle over lin-1 basis set orbitals
      do 90 j = 0, nlm1-1
!       cycle over 3 dipole allowed lin-->lin-1 transitions
        do 80 idip = 1,3
          if (idip.eq.1) then
            kin = lin
            kfi = lin - 1
          elseif (idip.eq.2) then
            kin = lin
            kfi = -lin
          else
            kin = -lin-1
            kfi = -lin
          endif
          if (kin.eq.0 .or.kfi.eq.0) goto 80

          jin = 2*abs(kin) - 1
          jfi = 2*abs(kfi) - 1
          nbs = -2*(j+nlp1)-1
          if (kfi.gt.0) nbs = nbs-1

!         cycle over minit
          do 70 mini = -jin, jin, 2
            mfi = mini + 2
            if (mfi.lt.-jfi .or. mfi.gt.jfi) goto 70

!           set variables for matrix index i
            i = i + 1
            jinit(i) = jin
            minit(i) = mini
            kinit(i) = kin
            jfin(i)  = jfi
            mfin(i)  = mfi
            kfin(i)  = kfi
            nph(i)   = nbs
  70      continue
  80    continue
  90  continue
!     last sanity check for matrix dimension
      if (i.ne. matsize) then
        print*, i, matsize
        stop 'FAILED matrix size check in subroutine getmat'
      endif

!c    Manual Input array nph
!c    automate it later
!c    L3 edge- inorb = 4, L2 edge - inorb = 3
      do im = 1, matsize
        if (kinit(im) .lt. 0) then
          ncore(im) = ihole
        else
          ncore(im) = ihole - 1
        endif
        if (ibasis.eq.0) then
         lfin = kfin(im)
         if (lfin.lt.0) lfin = abs(lfin) -1
         do iorb = 1, 41 
           lorb = kappa(iorb)
           if (lorb.lt.0) lorb = abs(lorb) - 1
           if (xnval(iorb).gt.0 .and. lfin.eq.lorb) then
             if (kfin(im).eq.kappa(iorb) .or. kfin(im).lt.0) nph(im)= iorb
!            above logic relies on the order or orbitals from subroutine
!            getorb: for the same lorb first appears the orbit
!            with j=lorb-1/2 and second with j=lorb+1/2
           endif
         enddo
        endif
!c     manual input
!c     for diamond and other 2p elements
!       if (kinit(im) .eq. -1) ncore(im) = 1
!       if (kfin(im) .eq. 1) nph(im) = 3
!       if (kfin(im) .eq. -2) nph(im) = 4
!      for Mg and other  3p
!       if (kinit(im) .eq. -1) ncore(im) = 1
!       if (kfin(im) .eq. 1 .and. im.le.np) nph(im) = 6
!       if (kfin(im) .eq. -2 .and. im.le.np) nph(im) = 6
!       if (kfin(im).eq. 1 .and. im.gt.np .and. im.le.2*np) nph(im) = 6
!       if (kfin(im).eq.-2 .and. im.gt.np .and. im.le.2*np) nph(im) = 6
!c     for 3d transition metal series L2,3 edges
!       if (kinit(im) .eq. 1) ncore(im) = 3
!       if (kinit(im) .eq. -2) ncore(im) = 4
!       if (kfin(im) .eq. 2) nph(im) = 8
!       if (kfin(im) .eq. -3) nph(im) = 9
!c     for 4d transition metal series N2,3 edges
!       if (kinit(im) .eq. 1) ncore(im) = 6
!       if (kinit(im) .eq. -2) ncore(im) = 7
!       if (kfin(im) .eq. 2) nph(im) = 13
!       if (kfin(im) .eq. -3) nph(im) = 14
!c     for Xe and 4f series N4,5 edges
!       if (kinit(im) .eq. 2) ncore(im) = 13
!       if (kinit(im) .eq. -3) ncore(im) = 14
!c     for W (tungsten) M4,5
!       if (kinit(im) .eq. 2) ncore(im) = 8
!       if (kinit(im) .eq. -3) ncore(im) = 9
!       nph is already set in getmat.f
      enddo

      return
      end

      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) stop 'singular matrix in ludcmp'    !pause to stop  KJ 6-2012
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software .


      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software .
