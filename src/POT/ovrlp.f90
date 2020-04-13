!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ovrlp.f90,v $:
! $Revision: 1.5 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ovrlp (iph, iphat, rat, iatph, novr, iphovr,           &
     &                nnovr, rovr, iz, nat, rho, dmag,                  &
     &                rhoval, vcoul, edens, edenvl, vclap, rnrm)

!     Overlaps coulomb potentials and electron densities for current
!     unique potential
      use constants
      use DimsMod, only: natx, nphx=> nphu, novrx, nrptx

      implicit double precision (a-h, o-z)


      dimension iphat(natx)
      dimension rat(3,natx)
      dimension iatph(0:nphx)
      dimension novr(0:nphx)
      dimension iphovr(novrx,0:nphx)
      dimension nnovr(novrx,0:nphx)
      dimension rovr(novrx,0:nphx)
      dimension iz(0:nphx)
      dimension rho(251,0:nphx+1), dmag(251,0:nphx+1)
      dimension vcoul(251,0:nphx+1), rhoval(251,0:nphx+1)
      dimension edens(251,0:nphx), edenvl(251,0:nphx)
      dimension vclap(251,0:nphx)
      dimension rnrm(0:nphx)
!#mn
       external dist

!     start with free atom values for current atom
      do 100  i = 1, 251
         vclap(i,iph) = vcoul(i,iph)
         edens(i,iph) = rho  (i,iph)
         
!c       investigate effect of central atom spin only
!        if (iph.ge.1) dmag(i,iph) = 0.0

         edenvl(i,iph) = rhoval  (i,iph)
  100 continue

      if (novr(iph) .gt. 0)  then
         do 104  iovr = 1, novr(iph)
            rnn  = rovr(iovr,iph)
            ann  = nnovr(iovr,iph)
            infr = iphovr(iovr,iph)
            call sumax (rnn, ann, vcoul(1,infr), vclap(1,iph))
            call sumax (rnn, ann, rho  (1,infr), edens(1,iph))
            call sumax (rnn, ann, rho  (1,infr), edenvl(1,iph))
  104    continue
      else
!        Do overlapping from geometry with model atom iat
         iat = iatph(iph)

!        overlap with all atoms within r overlap max (rlapx)
!        12 au = 6.35 ang  This number pulled out of a hat...
         rlapx = 12
!        inat is Index of Neighboring ATom
         do 110  inat = 1, nat
!           don't overlap atom with itself
            if (inat .eq. iat)  goto 110

!           if neighbor is too far away, don't overlap it
            rnn = dist (rat(1,inat), rat(1,iat))
            if (rnn .gt. rlapx)  goto 110

            infr = iphat(inat)
            call sumax (rnn, one, vcoul(1,infr), vclap(1,iph))
            call sumax (rnn, one, rho  (1,infr), edens(1,iph))
            call sumax (rnn, one, rho  (1,infr), edenvl(1,iph))
!ala        call sumax (rnn, one, rhoval(1,infr), edenvl(1,iph))
  110       continue
      endif

!     set norman radius
!     set norman radius
      IF(iz(iph).eq.0) THEN
         call frnrm (edens(1,iph), 1, rnrm(iph))
         PRINT '(A,I2,A,F20.10)', 'Norman radius for empty cell. iph = ', iph, ': ', rnrm(iph)
      ELSE
         call frnrm (edens(1,iph), iz(iph), rnrm(iph))
      END IF

!     remember ratio dmag/edens , not dmag itself
      do 200 i = 1,251
        if (edens(i,iph) .gt. 0.d0) then
          dmag(i,iph) = dmag(i,iph) / edens(i,iph)
        else
          dmag(i,iph) = 0.d0
        endif
 200  continue

      return
      end
