!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ovp2mt.f90,v $:
! $Revision: 1.6 $
! $Author: jorissen $
! $Date: 2011/02/11 04:14:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ovp2mt( nph, vtot, lrewr, qtot,ri,xnatph,lnear,        &
     &             inrm, imt, rnrm, rmt, cmovp, ipiv, vint, inters)
!  INPUT: nph - number of diferent potentials
!   vtot(i,iph) - potential OR density at point i for potential iph
!   lrewr       - if lrewr .gt. 0 potential will be overwritten
!                   density is never overwritten (lcoul.lt.0)
!                  lrewr=0 density calculation
!                  lrewr=1 potential calculation, vint estimated
!                  lrewr=2 potential calculation, vint is fixed
!   lcoul       -  .gt.0  (potential only) calculate charge for each iph
!                  .eq.0  (potential only) flat interstitial potential 
!                  .lt.0  (density only) calc charge inside MT spheres
!   qtot       -  for density only, total electron charge of cluster
!   ri         -  loucks radial grid
!   xnatph     -  number of atoms of type iph in the cluster
!   cmovp      -  LU decomposed overlapped matrix from movrlp.f
!   ipiv       -  pivoting indices for matrix cmovp
!  OUTPUT
!    vtot    if lrewr.gt.0  decomposed overlapped potential
!            if lrewr.le.0  old prescription for potential inside MT
!              spheres or don't want to overwrite densities
!    vint    mt zero level for potentials; charge outside mt spheres for
!            densities

      use constants
      use DimsMod, only: nphx=>nphu, novp
      implicit double precision (a-h, o-z)

      dimension vtot(251,0:nphx), xnatph(0:nphx)
      dimension inrm(0:nphx), imt(0:nphx), rmt(0:nphx), rnrm(0:nphx)
      dimension vtotav(0:nphx)
!     work space for linear algebra
!      parameter (novp=40)
      complex cmovp(novp*(nphx+1)+1,novp*(nphx+1)+1)
      complex cvovp(novp*(nphx+1)+1)
      integer ipiv(novp*(nphx+1)+1)
      dimension  ri(251)
      character*13 trans
      dimension  crho(251)
      logical lnear
      dimension lnear(0:nphx)
!pot      character*30 fname

!      get ipot and irav from inters
      ipot = mod(inters,2)
      irav = (inters-ipot) / 2
!     prepare cvovp and bvec from vtot
      ncp=0
      do 25 ip1=0,nph
      do 25 i=1,novp
        ncp = ncp + 1
        ix1 = imt(ip1)-novp + i
        cvovp(ncp)= real( vtot(ix1,ip1) )
       if (lrewr.eq.2) cvovp(ncp) = cvovp(ncp) - vint
  25  continue
      do 27 ip1=0,nph
         if (irav .eq. 1) then
           rav = (rmt(ip1) + rnrm(ip1)) / 2
         elseif(irav.eq.0) then
           rav =  rnrm(ip1)
         else
           rav = ri(imt(ip1)+1)
         endif
         if (lnear(ip1)) rav = ri(imt(ip1)+1)
         call terp(ri,vtot(1,ip1),inrm(ip1)+2,3,rav,vtotav(ip1))
  27  continue
      istx=novp*(nphx+1)+1
      trans = 'NotTransposed'
      nrhs = 1

!     find parameters for interstitial potential
      if (lrewr.gt.0) then
!        dealing with potentials
         if (lrewr.eq.1) then
!           additional equation to find vint
            ncp = ncp + 1
            cvovp(ncp) = 0
            bsum = 0
!           switch from average equation for vint to the local one
            nphlst = 0
            if (ipot .eq. 0) nphlst = nph
            do 430 iph=0,nphlst
               cvovp(ncp) = cvovp(ncp) + vtotav(iph)*xnatph(iph)
               bsum = bsum + xnatph(iph)
  430       continue
            cvovp(ncp) = cvovp(ncp) / bsum
         endif

         call cgetrs(trans, ncp, nrhs, cmovp, istx, ipiv, cvovp, istx, info)
         if (info.lt.0) then
             call par_stop('    *** Error in cgetrf')
         endif

         if (lrewr.eq.1) vint = dble(real(cvovp(ncp))) /100.0

!        rewrite vtot
         do 550 iph=0,nph
 
!pot  to write out ovp tot pot and its mt approxim, comment out cpot
!pot         write(fname,172)  iph
!pot  172    format('potp', i2.2, '.dat')
!pot         open (unit=1, file=fname, status='unknown', iostat=ios)
!pot         call chopen (ios, fname, 'wpot')

            do 500 i=1,novp
              index1=imt(iph)-novp + i
              index2=i+novp*iph

!pot            write(1,176) i, ri(index1), 
!pot     1             vtot(index1,iph),  dble(real(cvovp(index2)))+vint
!pot  176       format (1x, i4, 1p, 3e12.4)

              vtot(index1,iph) = dble(real(cvovp(index2)))+vint
  500       continue

!pot         close (unit=1)

!           use second order extrapolation
            j=imt(iph)+1
            call terp (ri,vtot(1,iph),imt(iph),2,ri(j),vtot(j,iph))
            do 505 j=imt(iph)+2, 251
  505       vtot(j,iph) = vint
  550    continue
      else
!        dealing with  density calculations. vint is the total charge inside mt spheres.
!        Divided by interstitial volume in istprm

         call cgetrs(trans, ncp, nrhs, cmovp, istx, ipiv, cvovp, istx, info)
         if (info.lt.0) then
             call par_stop('    *** Error in cgetrf')
!            stop
         endif

         vint = 0
         do 450 iph=0,nph
            do 440 i=1,imt(iph)+2
               if (i.lt.imt(iph)-novp+1) then
                 crho(i) =  vtot(i,iph)*ri(i)**2
               elseif (i.le. imt(iph)) then
                 ix1 = novp*iph +i-imt(iph)+novp
                 crho(i) = real(cvovp(ix1)) * ri(i)**2
!                crho(i) =  vtot(i,iph)*ri(i)**2
               else
                 call terp(ri,crho,imt(iph),2,ri(i), crho(i) )
               endif
  440       continue
            np = imt(iph) + 2
            cdum = 0
            dpas = 0.05d0
            call somm2 (ri,crho,dpas,cdum,rmt(iph),0,np)
            vint = vint + xnatph(iph) * cdum
  450    continue
         vint=qtot-vint
      endif

      return
      end
