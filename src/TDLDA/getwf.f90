!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: getwf.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2011/02/11 04:14:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getwf(ibasis, jj, nlp,nlm, nlpoc,nlmoc, rmt, rint, jri,&
     &           jint, em, eref, dx, x0, ri, v, vval, pat, qat,         &
     &           dgcn, dpcn, adgc, adpc, dgcnp, dpcnp, xnval,           &
     &           iz, ihole, xion, iunf, kinitm, kfinm, nph, matsize, iph) !KJ iph

      use dimsmod, only: nrptx
	  use constants
      implicit double precision (a-h, o-z)
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
!  information specific for finding basis set orbitals
!     jj - index for the basis set orbital (starts from 0)
!     nlp -number of (L+1) orbitals in basis (L-initial state orb.mom.)
!     nlm -number of (L-1) orbitals in basis
!     nlpoc - number of completely occcupied (L+1) orbitals
!     nlmoc - number of completely occcupied (L-1) orbitals
!  information for dfovrg call, to get radial solution for some energy
!     rmt, rint - muffin-tin and Norman radius
!     jri, jint - indeces for muffin-tin and Norman radii on radial grid
!     em, eref - the energy of photoelectron
!     dx, x0, ri - radial grid information
!     v, vval - potential
!     pat, qat - orbitals
!     dgcn, dpcn, adgc, adpc - fully and partially occupied orbitals
!     dgcnp, dpcnp -  unoccupied orbitals (to be calculated and stored)
!     xnval - atomic occupations for valence orbitals
!     iz - nuclei charge
!     ihole - index of core-hole orbital
!     xion - ionicity (zero by default)
!     iunf - flag to freeze/unfreeze f-electrons
!  information about matrix indeces for K and Chi_0
!     kinitm - initial state kappa for each matrix index
!     kfinm - final state kappa for each matrix index
!     nph - index for the basis set orbital for each matrix index
!     matsize - matrix dimension

      complex*16 ptz
      dimension ptz(-1:1, -1:1)
      complex*16 xrcold(nrptx) , xncold(nrptx)

      dimension ri(nrptx), vtot(nrptx), edens(nrptx),dmag(nrptx)
      dimension vvalgs(nrptx), edenvl(nrptx)
      dimension dgcn(nrptx,41), dpcn(nrptx,41)
      dimension adgc(10,41), adpc(10,41), xnval(41)
      dimension dgcnp(nrptx,41), dpcnp(nrptx,41)
   
      dimension xp(nrptx), xq(nrptx)

!     work space for xcpot
      dimension vxcrmu(nrptx), vxcimu(nrptx), gsrel(nrptx)
      dimension vvxcrm(nrptx), vvxcim(nrptx)

!     work space for fovrg
      complex*16 p(nrptx), q(nrptx), pn(nrptx), qn(nrptx)

      complex*16  p2, ck, xkmt, xkmtp
      complex*16  pu, qu, dum1, factor
      complex*16  xfnorm, xirf
      complex*16  temp, aa, bb, cc, phold
      complex*16  phx(8), ph0
      complex*16  xm1, xm2, xm3, xm4
      complex*16 jl,jlp1,nl,nlp1
      complex*16  v(nrptx), vval(nrptx)
      complex*16  xrc(nrptx), xnc(nrptx)
      character*512 slog

!     nesvi:  
      parameter (maxsize = 78)
      dimension kinitm(maxsize), kfinm(maxsize)
      dimension nph(maxsize)
      complex*16 eref
      dimension pat(nrptx), qat(nrptx)
!     orbital input file imension
      parameter (norbx=5000)
      dimension phi(norbx), rip(norbx)
      

!     set initial orb. mom. and dimension in matrix per each final orbital
      ll = kinitm(1)
      if (ll.le.0) ll = -ll -1
      n1 = 3*(2*ll+1)
      n2 = 3*(2*ll-1)
!     set ilast larger than jri for better interpolation for pu
      ilast = jint + 1
      if (ilast.le.jri) ilast = jri+1

      if (jj.lt.nlm+nlp) then
!       set initial and final indeces for this orbital in matrix
        if (jj.lt.nlp) then
          imi = 1+n1*jj
          imf = n1*(jj+1)
        else
          imi = 1 + n1*nlp + n2*(jj-nlp)
          imf = imi + n2 - 1
        endif
        
        do im = imi, imf
!         notice that actually don't need to cycle over all im's 
!         but rather only over positive and negative kappa (ik)
!         iforb = 2*jj + ik + 1
!         however neglect it for now since it is not a bottleneck

!         check whether really need to calculate the orbital
!         iforb = - nph(im)
          iforb = 2*jj+1
          if (kfinm(im).gt.0) iforb = iforb+1
          kfin = kfinm (im)
          lfin = kfin
          if (kfin.le.0) lfin = abs(kfin) - 1
          do i=1,nrptx
            pat(i) = 0
            qat(i) = 0
          enddo

          if (ibasis.eq.0)  then
!           copy orbital from array dgcn to dgcnp
            do i=1, nrptx
              pat(i) = dgcn(i,nph(im))
            enddo
            do i=1, nrptx
              qat(i) = dpcn(i,nph(im))
            enddo
            nph(im) = - iforb
          elseif (ibasis.eq.1) then
!           read orbital from file, neglecting SO in final state
            nph(im) = - iforb
            if (jj.eq.2) then
              open (unit=3, file='Vila/Orbs/mg.4p.dat', status='old')
            elseif (jj.eq.1) then
              open (unit=3, file='Vila/Orbs/mg.4p.dat', status='old')
            else 
              open (unit=3, file='Vila/Orbs/mg.3p.dat', status='old')
            endif
            n=0
  10        n = n+1
              read(3,*, end=20) rip(n), phi(n)
!             need \psi(r) * r and distance in bohrs
              phi(n) = phi(n)*rip(n)
              rip(n) = rip(n) /bohr
              goto 10
  20        continue
            n = n-1
            close (unit=3)

!           interpolate on our radial grid
            do i = 1, ilast
              call terp (rip, phi, n, 1, ri(i), pat(i))
            enddo
            do i = 1, nrptx
              qat(i) = 0
              if (i.gt.ilast) pat(i) = 0
            enddo
          else
!           ibasis=2
!           find orbital for required number of nodes and zero at R_int

!           note that we neglect SO interaction for projection operator
!           thus need to calculate orbital only for one of j=l +/- 1/2
!           p2 is (complex momentum)**2 referenced to energy dep xc
!           if the initial state is p1/2(L2) edge, then subtract spin-orbit
!           splitting, because em
!           is linked to p3/2 energy origin
            p2 = em - (dble(eref))      
            ck = sqrt (2*p2 + (p2*alphfs)**2)
            xkmt = rmt * ck

            irr = -1
            ncycle = 0
            ic3 = 1
            call dfovrg (ncycle, kfin, rmt, ilast, jri, p2, dx,         &
     &               ri, v,vval, dgcn, dpcn, adgc, adpc,                &
     &               xnval, pu, qu, p, q,                               &
     &               iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph

            ilp = lfin - 1
            if (kfin .lt. 0) ilp = lfin + 1
            call exjlnl (xkmt, lfin, jl, nl)
            call exjlnl (xkmt, ilp, jlp1, nlp1)
            call phamp(rmt,pu,qu, ck, jl,nl,jlp1,nlp1, kfin, ph0,temp)

            sign = -1.0
            if (kfin.gt.0) sign = 1.0
            factor = ck*alphfs 
            factor = sign * factor/(1+sqrt(1+factor**2))
            dum1 = 1/ sqrt(1+factor**2)
            xfnorm = 1 / temp *dum1
!           normalization factor
!           xfnorm = dum1*rmt*(jl*cos(delta) - nl*sin(delta))/ Rl(rmt)
!           dum1 is relativistic correction to normalization
!           normalize regular solution
            do i = 1,ilast
              p(i)=p(i)*xfnorm
              q(i)=q(i)*xfnorm
            enddo

!           cut solutions beyond ii-th zero
            inul = ilast
!           set ii to the number of nodes to be found (including one at
!           norman radius) In case of Xe simple connection with jj index.
            ii = jj+ 1 + nlpoc
            if (jj.ge.nlp) ii = (jj-nlp) +1 + nlmoc
            do  i = 5, ilast
              if (dble(p(i-1))* dble(p(i)) .le. 0) then
                inul = i -1
                if (dble(p(i)).ne.0) ii = ii - 1
                if (ii.eq.0) goto 30
              endif
            enddo
   30       continue
            !print*, 'getwf: second index should be 0; '
            !print*, '   and distance close to norman radius'
            !print*,  jj,ii,kfin,  'r(inul)' , ri(inul)*bohr
            !print*,  jj, em*hart

            do i = 1,inul 
              pat(i)=dble(p(i))
              qat(i)=dble(q(i))
            enddo
          endif
!         orthogonalize pat to previous f-orbitals
          if(jj.le.nlp-1) then
            jin = 0
          else
            jin = nlp
          endif
          do jp = jin, jj-1
            ifp = iforb-2*((jp-jin)+1)
            do i = 1, ilast
              xp(i) = pat(i)*dgcnp(i,ifp) + qat(i)*dpcnp(i,ifp)
              xq(i) = 0
            enddo
            xinorm = 2*lfin + 2
            i0 = jint + 1
            call somm2 (ri, xp, dx, xinorm, rint, 0, i0)
            do i=1,nrptx
              pat(i) = pat(i) - xinorm*dgcnp(i,ifp)
              qat(i) = qat(i) - xinorm*dpcnp(i,ifp)
            enddo
          enddo

!         normalize pat and qat in the Norman radius sphere: <n|n>=1,
!         (renormalized atomic sphere method)
          do i = 1, ilast
            xp(i) = pat(i)**2 + qat(i)**2
            xq(i) = 0
          enddo
!         nb, xinorm is used for exponent on input to somm 
          xinorm = 2*lfin + 2
          i0 = jint + 1
          call somm2 (ri, xp, dx, xinorm, rint, 0, i0)
        
          xinorm = sqrt(xinorm)
          do i=1,nrptx
            pat(i) = pat(i) / xinorm
            qat(i) = qat(i) / xinorm
          enddo
          do  i=1, nrptx
            dgcnp(i,iforb) = pat(i)
            dpcnp(i,iforb) = qat(i)
          enddo
        enddo
      endif

      return 
      end 
