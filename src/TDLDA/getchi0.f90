module getchi0_workspace
      use dimsmod,only: nex
!     to pass energy levels and projected DOS
      integer neg(30)
      real*8 rhoj(nex,30), eng(nex,30)

      integer norbp
end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: getchi0.f90,v $:
! $Revision: 1.7 $
! $Author: jorissen $
! $Date: 2011/03/29 01:59:33 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getchi0(ie, em, eref, edge, emu, refsh, omega,         &
     &           ipmbse, ifxc, sfun,                                    &
     &           rmt, rint, jri, jint, dx, x0, ri, edens, v, vval, vch, &
     &           eorb, kappa, dgcn, dpcn, adgc, adpc,                   &
     &           dgcnp, dpcnp, xnval, iz, ihole, xion, iunf, nfo, npo,  &
     &           jinit,minit,kinitm,jfin, mfin, kfinm, ncore, nph, matsize,&
     &           dipmatl, dipmat, chi0im,                               &
     &           wmat, xkmat, xkmatp, phf, iph) !KJ iph
     
     
      use dimsmod, only: nrptx
	  use constants
      use getchi0_workspace
      implicit double precision (a-h, o-z)

      complex*16 xrcold(nrptx) , xncold(nrptx)

      real*8 ri(nrptx), edens(nrptx),fxcim(nrptx)
      real*8 fxc(nrptx), fxc0(nrptx)
      real*8 eorb(30)
      integer kappa(30)
      real*8 dgcn(nrptx,30), dpcn(nrptx,30)
      real*8 dgcnp(nrptx,30), dpcnp(nrptx,30)
      real*8 adgc(10,30), adpc(10,30), xnval(30)
      integer iorb(-4:3)
   
      real*8 xp(nrptx)

!     work space for xcpot
      real*8 vxcrmu(nrptx), vxcimu(nrptx), gsrel(nrptx)
      real*8 vvxcrm(nrptx), vvxcim(nrptx)

!     work space for fovrg
      complex*16 p(nrptx), q(nrptx), pn(nrptx), qn(nrptx), fscf(nrptx)
      complex*16 pp(nrptx), qp(nrptx), pnp(nrptx), qnp(nrptx)
!!     to pass energy levels and projected DOS
!      integer neg(30), eng(nex,30)
!      real*8 rhoj(nex,30)

      complex*16  p2, ck, xkmt, xkmtp
      complex*16  pu, qu, dum1, factor
      complex*16  xfnorm, xirf, xirf1
      complex*16  temp, aa, bb, cc, rkk1, rkk0, phold
      complex*16  phx(8), ph0
      complex*16  xm1, xm2, xm3, xm4
      complex*16 jl,jlp1,nl,nlp1
      complex*16  v(nrptx), vval(nrptx)
      real*8 vch(nrptx),  vcx(nrptx)
      complex*16  xrc(nrptx), xnc(nrptx)
      character*512 slog

!     nesvi:  
      complex*16 apsm(10,matsize), aqsm(10,matsize),p2m(matsize)
      parameter (maxsize = 78)
      real*8 chi0im(maxsize, maxsize)
      integer jinit(maxsize), minit(maxsize)
      integer jfin(maxsize), mfin(maxsize)
      integer kinitm(maxsize), kfinm(maxsize)
      real*8 dipmat(maxsize), dipmatl(maxsize)
      real*8 phf(maxsize), pref(maxsize)
      real*8 bf(0:2, nrptx)
      real*8 refsh(maxsize)
      complex*16 eref
      integer ncore(maxsize)
      real*8 dgc0(nrptx), dpc0(nrptx)
      complex*16 pf(nrptx,matsize), qf(nrptx,maxsize)
      complex*16 ptot(nrptx,matsize), qtot(nrptx,maxsize)
      real*8 pc(nrptx,maxsize), qc(nrptx,maxsize)
      complex*16 wmat(maxsize,maxsize), yvec(nrptx, maxsize)
      complex*16 var(nrptx), xkmat(maxsize,maxsize)
      complex*16 xkmatp(maxsize,maxsize), rabcd, rabcdp
      complex*16 ykgr(nrptx), ykgrex(nrptx) 
      real*8 rnorm1(maxsize)
      integer nph(maxsize)
      real*8 pat(nrptx), qat(nrptx)
      real*8 ovrl(maxsize)

      real*8, external :: ellpi
      

      xkmat(:,:) = 0
      xkmatp(:,:) = 0
!write(*,*) 'getchi0 at energy ',ie

!     remember the bessel functions for multipole matrix elements
      xk0 = omega * alphfs
      ilast = jri+6
      if (ilast.le.jint) ilast = jint + 1
      if (ilast.gt.nrptx) ilast = nrptx
      do 50 i = 1, ilast
        temp = xk0 * ri(i)
        if (abs(temp).lt.1.d0) then
!         use series expansion
          do ll = 0,2
            call bjnser(temp,ll, xirf, dum1,1)
            bf(ll,i) = dble(xirf)
          enddo
        else
!         use formula
          x = dble(temp)
          sinx = sin(x)
          cosx = cos(x)
          bf(0,i) = sinx/x
          bf(1,i) = sinx/x**2 - cosx/x
          bf(2,i) = sinx*(3/x**3-1/x) - 3*cosx/x**2
        endif

!       also calculate the local exchange term fxc = vxc*r_s**3/r**2
        if (edens(i).le.0) then
          if(mod(i,10).eq.0) then
            write(slog, 149) 'negative dens ', i, ' - usually due to harmless precision errors, but check DOS to make sure'
  149       format (a, 2i3)
            call wlog(slog)
          endif
          rs = 100
        else
          rs = (4*pi*edens(i)/3)**(-third)
        endif
!       vvbh from Von Barth Hedin paper, 1971
!       see eq.60-61 of Gross&Kohn, Adv. Quant. Chem. 21, p.255(1990).
!       eq.60 fxc = d V_xc / d rho * (2\ell + 1) /(4*pi*r**2)
!       where second factor comes from \delta(r-r')
!       for dipole field \ell=1 
        fxc0(i) = rs**3 / ri(i)**2 / 6 * (-1.22177412/rs -1.512/(30+rs))
!c       below are coefficients in Zangwill/Soven paper
!       fxc0(i) = rs**3 / ri(i)**2 / 6 * (-1.222/rs -0.75924/(11.4+rs))
!       from eq.61 fxc = f_inf * (2\ell + 1) /(4*pi*r**2)
        rsx = rs / 30
        fxc(i) = rs**3 / ri(i)**2 / 6 * (-1.22177412*0.6/rs -1.008* (1.0/3.0 - rsx/2 + rsx**2- (rsx**3+0.22)*log(1.0+1/rsx) ))

!       choose ifxc=0 for RPA, 1 for Zangwill Soven, 2 for Gross kohn
!       3 for kernel in our paper 4 ZS for diag, GK for off diagonal
        if (ifxc .eq.0) then
          fxc(i) = 0
          fxcim(i) = 0
          fxc0(i) = 0
        elseif (ifxc.eq.1) then
          fxc(i) = fxc0(i)
          fxcim(i) = 0
        elseif (ifxc.eq.2 .or. ifxc.eq.4) then
          cc = 23*pi/15 / ( 4./3.*pi*ri(i)**2)
          gam = 3.6256099082**2/ (4*sqrt(2*pi))
          bxc = gam/cc*(fxc(i)-fxc0(i))
          axc = -cc * bxc**(5./3.)
          bxc = bxc**(4./3.)
          vso = - refsh(1)
          if (ifxc.eq.2) vso = omega
          fxcim(i) = axc*vso / (1+bxc*vso**2)**(5./4.)
          ss = sqrt(1+bxc*vso**2)
          s1 = (1-ss)/2
          s2 = (1+ss)/2
          fxc(i) = fxc(i) + axc/pi/ss**2 *sqrt(8./bxc)* (2*1.350644 - s2*ellpi(s1) - s1*ellpi(s2))
!         to compare with Gross Kohn paper use line below
!         fxc(i) = fxc(i) * ( 4./3.*pi*ri(i)**2)/rs**3
        elseif (ifxc.eq.3 .or. ifxc.eq.5) then
          fxc(i) = 0
          fxcim(i) = 0
        endif
 50   continue


      chi0im(:,:) = 0
      dipmatl(:) = 0
      dipmat(:) = 0
      ovrl(:) = 0

! -------------------- separate into K^zs and K^pm

!------------------------
!c              calculate K^zs(r,r')*\phi_i(r') phi_f(r')
      do 29 imj = 1, matsize
!       calculate contribution to K due to core-hole potential
        yvec (:,imj) = 0
!       v = 1
        nu = 1
        nu2 = 2*nu

        call yzktd (ncore(imj),nu,flps, pf(1,imj), qf(1,imj), apsm(1,imj), aqsm(1,imj), ykgr, 0)

!       ykgr is Y(r) = r * int dr' U(r,r') Psi_k'(r') R_l'(r') 
!       where R_l(r) is pf (qf)
!       Psi_k(r) is dgc0(ncore) (dpc0)
!       Y(B,D,r), where B = k', D = l', see Grant's paper

        do i = 1, ilast
!c          need Yk(r) /r
            yvec(i,imj) =  dble(ykgr(i)) / ri(i) 
!           add xc term here
            if (ifxc.ne.2) then
              yvec(i,imj) = yvec(i,imj) + fxc0(i) * dble(pc(i,imj)*pf(i,imj) + qc(i,imj)*qf(i,imj))
            else
              yvec(i,imj) = yvec(i,imj) + (fxc(i)+coni*fxcim(i)) * dble(pc(i,imj)*pf(i,imj) + qc(i,imj)*qf(i,imj))
            endif
!           multiply by separation function
            yvec(i,imj) = yvec(i,imj) * sfun
        enddo
              
!         A = k, C = l
!         d^v prefactor (angular part of integral)

          jb2 = jinit(imj)
          mb2 = minit(imj)
          jd2 = jfin(imj)
          md2 = mfin(imj)

          ind2 = 2
!         Grant's condition (7.13)
!         ja+v+jc is odd if Aa = Ac, even otherwise
!         jd+v+jb is odd if Ad = Ab, even otherwise
!         in our case it turns out j+v+j' should be even => v is odd
!         v =1 is the leading term

          dnudb= cwig3j(jd2,nu2,jb2,1,0,2)
          dnudb= dnudb * cwig3j(jd2,nu2,jb2,-md2,(md2-mb2),ind2)
          x2 = (jd2+1) * (jb2+1)
          pref(imj) = (-1)**((md2 + 1)/2.0) * dnudb * sqrt(x2)
  29  continue
!---------------------

!     dipole only
      ks = 1  
!c    calculate screened dipole field
      ww = dble(emu+p2-edge)
      if (ks.eq.1) then
        p2 = em - dble(eref) + coni*1.d-8
        p2f = edge - dble(eref)
        if (ie.eq.1) call correorb(iz, ihole, rmt, jri, dx,ri,          &
     &               p2f,edge, v, dgcn, dpcn, adgc, adpc,               &
     &               eorb, neg, eng, rhoj, kappa, norbp, iph) !KJ
!       print*,'ie=', ie
! Debug: Fer
!write(*,*) 'norbp',norbp
!write(*,*) 'neg',neg
!write(*,*) 'kappa',kappa
!write(*,*) 'iph',iph

        call phiscf (ifxc, rmt, ilast, jri, p2, p2f, emu, dx,           &
     &              ri, v, edens, dgcn, dpcn, adgc, adpc,               &
     &              iz, ihole, neg, eng, rhoj,kappa, norbp, fscf,       &
     &              yvec, maxsize, matsize, sfun, iph) !KJ

        do imi = 1, matsize
           do imj = 1, matsize
             do i = 1,ilast
               var(i) = yvec(i,imj)* dble(pc(i,imi)*pf(i,imi) + qc(i,imi)*qf(i,imi))
             enddo
!            integrate over r
             rabcd = 0
             do i = 2, ilast
               rabcd=rabcd + (var(i)+var(i-1))*(ri(i)-ri(i-1)) / 2
             enddo
             wmat(imi,imj) = rabcd * pref(imi) * pref(imj)
           enddo
        enddo
        wse = dble(p2-eng(1,ihole))
!     endif
      else
        fscf(:) = 1.d0
        wmat(:,:) = 0
        wse = ww
      endif
      ww = sqrt(wse/ww)
      idim = matsize
!c    change idim here to remove some channels
!     idim = 60
      do 20 im = 1, idim
        kinit = kinitm(im)
        kfin = kfinm(im)
        kdif = kfin - kinit

        lfin = kfin
        if (kfin.le.0) lfin = abs(kfin) - 1


!       p2 is (complex momentum)**2 referenced to energy dep xc

!       if the initial state is p1/2(L2) edge, then subtract spin-orbit
!       splitting, because em
!       is linked to p3/2 energy origin
        p2 = em - (dble(eref) - refsh(im))      
        ck = sqrt (2*p2 + (p2*alphfs)**2)
        xkmt = rmt * ck
        p2m(im) = p2 

        if (dble(p2).le.0.d0) goto 20
        if (dble(p2).le.dble(v(jri+1))) goto 20

      
!       check that orbital momentum does not exceed max allowed
!       if (lfin .gt. lx) then
!c        set final j and l to unphysical values
!         lfin = -1 
!       endif
       
!       if (ltolm1.eq.0 .and. ((kinit.lt.0 .and. ind.ge.3) .or.
!     1          (kinit.gt.0 .and. ind.ne.3)) ) goto 300

        ikap = kfin

        irr = -1
        ncycle = 0
        ic3 = 0
!       set ilast larger than jri for better interpolation for pu
!       also need 5 points after jri for irregular solution

        p2 = p2 + coni*0.0001/hart
        call dfovrg (ncycle, ikap, rmt, ilast, jri, p2, dx,             &
     &               ri, v,vval, dgcn, dpcn, adgc, adpc,                &
     &               xnval, pu, qu, p, q,                               &
     &               iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph

        ilp = lfin - 1
        if (ikap .lt. 0) ilp = lfin + 1
        call exjlnl (xkmt, lfin, jl, nl)
        call exjlnl (xkmt, ilp, jlp1, nlp1)
        call phamp(rmt,pu,qu, ck, jl,nl,jlp1,nlp1, ikap, ph0,temp)
        phf(im) = dble(ph0)

        sign = -1.0
        if (ikap.gt.0) sign = 1.0
        factor = ck*alphfs 
        factor = sign * factor/(1+sqrt(1+factor**2))
        dum1 = 1/ sqrt(1+factor**2)
        xfnorm = 1 / temp *dum1
!       normalization factor
!       xfnorm = dum1*rmt*(jl*cos(delta) - nl*sin(delta))/ Rl(rmt)
!       dum1 is relativistic correction to normalization
!       normalize regular solution
        do i = 1,ilast
          p(i)=p(i)*xfnorm
          q(i)=q(i)*xfnorm
        enddo

        do ii = 1, nrptx
          dgc0(ii) = dgcn(ii,ncore(im))
          dpc0(ii) = dpcn(ii,ncore(im))
        enddo
      
!       prepare stuff for K matrix
        do i = 1, ilast
          pc(i,im) = dgc0(i)
          qc(i,im) = dpc0(i)
        enddo

        apsm(:,im) = 0
        aqsm(:,im) = 0
        flps = 1.d0
        apsm(1,im) = pc(1,im) / ri(1)
        aqsm(1,im) = qc(1,im) / ri(1)
  
!       pat, qat - atomic functions that we make projection on.
        do i=1,nrptx
         if (nph(im).gt.0) then
          pat(i) = dgcn(i,nph(im))
          qat(i) = dpcn(i,nph(im))
         else
          pat(i) = dgcnp(i,-nph(im))
          qat(i) = dpcnp(i,-nph(im))
         endif
        enddo

!       normalize pat and qat in the Norman radius sphere: <n|n>=1,
!       (renormalized atomic sphere method)
     
        do i = 1, ilast
          xp(i) = pat(i)**2 + qat(i)**2
        enddo
!       nb, xinorm is used for exponent on input to somm 
        xinorm = 2*lfin + 2
        i0 = jint + 1
        call somm2 (ri, xp, dx, xinorm, rint, 0, i0)
      
        xinorm = sqrt(xinorm)
        do i=1,ilast
          pat(i) = pat(i) / xinorm
          qat(i) = qat(i) / xinorm
        enddo
  
!       calculate overlap integral between f and atomic function 
        do i = 1, ilast
          xp(i) = pat(i)* p(i) + qat(i)*q(i)
        enddo
!       nb, xinorm is used for exponent on input to somm 
        xinorm = 2*lfin + 2
        i0 = jint + 1
        call somm2 (ri, xp, dx, xinorm, rint, 0, i0)
      
        ovrl(im) = xinorm
!       ploc, qloc - localized part of the functions.
        do i=1,ilast
          pf(i,im) = pat(i)
          qf(i,im) = qat(i)
          ptot(i,im) = p(i)
          qtot(i,im) = q(i)
        enddo
        do i=1+ilast, nrptx
          pf(i,im) = 0
          qf(i,im) = 0
          ptot(i,im) = 0
          qtot(i,im) = 0
        enddo
      
        mult = 0
        iold = 0
!       call radint(1, mult, bf, kinit, dgc0, dpc0, ikap, p, q,
!    1     pn, qn,ri, dx, ilast, iold, xrc, xnc, xrcold, xncold, xirf)
!c      calculate xirf including fscf - TDLDA result
        do id = 1, 2
          if (id.eq.1) then
            do j = 1,ilast
              pp(j)  = p(j)*dble(fscf(j))
              qp(j)  = q(j)*dble(fscf(j))
            enddo
          else
            do j = 1,ilast
              pp(j)  = p(j)*dimag(fscf(j))
              qp(j)  = q(j)*dimag(fscf(j))
            enddo
          endif
          ifl = -1
!         ifl = 1
          xirf1 = 0
          call radint(ifl, mult, bf, kinit, dgc0,dpc0, ikap, pp,qp, pn,qn,ri,dx, ilast,iold, xrc,xnc, xrcold,xncold, xirf1)
          if (ifl.lt.0) xirf1 = xirf1 * xk0 * ww
          if (id.eq.1) then
            xirf = xirf1
          else
           if (abs(xirf1) .lt. abs(xirf)) then
              dum = abs(xirf1) / abs(xirf)
              xirf = xirf * sqrt(1.d0 + dum**2)
            else
              dum = abs(xirf) / abs(xirf1)
              xirf = xirf1 * sqrt(1.d0 + dum**2)
            endif
          endif
        enddo

        dum = dimag(xirf)
!       dipmat(im) = dimag(xirf)
!       note that for real potential  xirf is real or reduced matrix
!       element for dipole transition is pure imaginary.
            
!       add (-1)^(j-m) (j L j') factor 
!                         -m pol m'                                 
!       this factor averages to 1/3 after summation over m,m'
!       note, that j,m already have factor 2
!       L = 1 in dipole approximation 
        l2 = 2
!       polarization -1, which corresponds to m'-m = 1  
!       ipol = minit(im) - mfin(im)  
        ipol = mfin(im) - minit(im)
        ind2 = 2
!       angpart= cwig3j(jinit(im),l2,jfin(im),-minit(im),ipol,ind2)
        angpart= cwig3j(jfin(im),l2,jinit(im),-mfin(im),ipol,ind2)
 
!       full dipole m.e.
!       dipmat(im)=dipmat(im) * angpart*(-1)**((jinit(im)-minit(im))/2)
!       dipmat(im)=dipmat(im) * angpart*(-1)**((jfin(im)-mfin(im))/2)
        dum=dum * angpart*(-1)**((jfin(im)-mfin(im))/2)
!       manual changes required: for Xe p-->d  final m.e. is in
!       1--15 positions while p-->s in 16--24 positions
!       logic is also fine for 3d, 4d, and 2p elements calculations
        lin = kinitm(1)
        if (lin.lt.0) lin = abs(lin) - 1
        matszp = 3*(2*lin+1)
        matszm = 3*(2*lin-1)
        if (lin.eq.0) matszm = 0

        if (nfo.gt.0) then
          if (im.le.matszp) dipmat(im) = dum
          imp = im - nfo*matszp
          if (imp.gt.0 .and. imp.le.matszm)  dipmat(imp+matszp) = dum
        else
          if (im.le.matszm)  dipmat(im) = dum
        endif


!c      localized part only
!       call radint(1, mult, bf, kinit,dgc0,dpc0,ikap,pf(1,im),qf(1,im),
!    1       pn, qn, ri, dx, ilast, iold, xrc, xnc, xrcold,xncold, xirf)
        do id = 1, 2
          if (id.eq.1) then
            do j = 1,ilast
              pp(j)  = pf(j,im)*dble(fscf(j))
              qp(j)  = qf(j,im)*dble(fscf(j))
            enddo
          else
            do j = 1,ilast
              pp(j)  = pf(j,im)*dimag(fscf(j))
              qp(j)  = qf(j,im)*dimag(fscf(j))
            enddo
          endif
          ifl = -1
!         ifl = 1
          xirf1 = 0
          call radint(ifl, mult, bf, kinit, dgc0,dpc0, ikap, pp,qp,pn,qn,ri,dx, ilast,iold, xrc,xnc, xrcold,xncold, xirf1)
          if (ifl.lt.0) xirf1 = xirf1 * xk0 * ww
          if (id.eq.1) then
            xirf = xirf1
          else
           if (abs(xirf1) .lt. abs(xirf)) then
              dum = abs(xirf1) / abs(xirf)
              xirf = xirf * sqrt(1.d0 + dum**2)
            else
              dum = abs(xirf) / abs(xirf1)
              xirf = xirf1 * sqrt(1.d0 + dum**2)
            endif
          endif
        enddo

        dipmatl(im) = dimag(xirf)
!       dipmatl(im)=dipmatl(im)* angpart*(-1)**((jinit(im)-minit(im))/2)
        dipmatl(im)=dipmatl(im)* angpart*(-1)**((jfin(im)-mfin(im))/2)

!       if overlap integral < 0, then no localized part
!temp   if (ovrl(im) .le. 0) then
!temp     dipmatl(im) = 0.0
!temp     ovrl(im) = 0.0
!temp   endif

!       selection rules for dipole matrix elements and chi0
!       if ( abs(kdif) .ne. ks .and. kfin .ne. - kinit ) goto 20   
     
        chi0im(im,im) = - 2*dble(ck) * (ovrl(im))**2
        if (nfo.gt.1) then
          do i=1,nfo-1
            imp = im-i*matszp
            if (imp.gt.0 .and. em .ge.(edge - refsh(im))) then
              chi0im(im,imp) = - 2*dble(ck) * ovrl(im)* ovrl(imp)
              chi0im(imp,im) = chi0im(im,imp) 
            endif
          enddo
        endif
        if (npo.gt.1) then
          do i=1,npo-1
            imp = im-i*matszm
            if (imp.gt.nfo*matszp.and. em .ge.(edge - refsh(im))) then
              chi0im(im,imp) = - 2*dble(ck) * ovrl(im)* ovrl(imp)
              chi0im(imp,im) = chi0im(im,imp) 
            endif
          enddo
        endif

!       only unoccupied part 
        if (em .lt. (edge - refsh(im))) then
          chi0im(im,im) = 0
          dipmatl(im) = 0
          dipmat(im) = 0
        endif
         
   20 continue

!     write out projected Im chi0 for testing: e.g.
!     write (44, *) em, chi0im(1,1), chi0im(4,4), chi0im(7,7)

!------------------------
!c              calculate K

!     first fix multipliers for direct and xc terms
      sfx = 1-sfun
      sxc = 1
      if (ipmbse.ge.3) sxc = -sfun
! test logic
      do i = 1, ilast
        vcx(i) = vch(i)*sfx
      enddo
! end test
!     scale  xc term according to to ipmbse
!     scale direct term later 
      do i = 1, ilast
        fxc0(i) = sxc * fxc0(i)
        fxc(i)  = sxc * fxc(i)
        fxcim(i)= sxc * fxcim(i)
      enddo

      do 30 imj = 1, idim
!     do 30 imj = 1, matsize
        if (dble(p2m(imj)).le.0.d0) goto 30

!       calculate contribution to K due to core-hole potential
        do 241 i = 1, jint
  241   var(i) = vcx(i)* ( pf(i,imj)**2 + qf(i,imj)**2)
        rabcd = 0
        do 242 i = 2, jint
  242   rabcd = rabcd + (var(i)+var(i-1))*(ri(i)-ri(i-1)) / 2
        do 243 i = 1, jint
  243   var(i) = vcx(i)*( pf(i,imj)*ptot(i,imj)+qf(i,imj)*qtot(i,imj) )
        rabcdp = 0
        do 244 i = 2, jint
  244   rabcdp = rabcdp + (var(i)+var(i-1))*(ri(i)-ri(i-1)) / 2
        xkmat(imj, imj) = rabcd
        if (imj.le.nfo*matszp) then
          imp = mod(imj,matszp)
          if (imp.eq.0) imp = matszp
        else
          imp = matszp + mod((imj-nfo*matszp), matszm)
          if (imp.eq.matszp) imp = matszp + matszm
        endif
        xkmatp(imp, imj) = rabcdp

        if (nfo.gt.1) then
          do j=1, nfo-1
            imp = imj-j*matszp
            if (imp.gt.0) then
!c           calculate off-diagonal matrix elements of V_0
             do  i = 1, jint
              var(i)=vcx(i)* (pf(i,imj)*pf(i,imp) +qf(i,imj)*qf(i,imp))
             enddo
             rabcd = 0
             do  i = 2, jint
              rabcd = rabcd + (var(i)+var(i-1))*(ri(i)-ri(i-1)) / 2
             enddo
             xkmat(imj, imp) = rabcd
             xkmat(imp, imj) = rabcd
            endif
          enddo
        endif
        if (npo.gt.1) then
          do j=1,npo-1
            imp = im-j*matszm
            if (imp.gt.nfo*matszp.and. em .ge.(edge - refsh(im))) then
!             notice: same as above
!c           calculate off-diagonal matrix elements of V_0
             do  i = 1, jint
              var(i)=vcx(i)* (pf(i,imj)*pf(i,imp) +qf(i,imj)*qf(i,imp))
             enddo
             rabcd = 0
             do  i = 2, jint
              rabcd = rabcd + (var(i)+var(i-1))*(ri(i)-ri(i-1)) / 2
             enddo
             xkmat(imj, imp) = rabcd
             xkmat(imp, imj) = rabcd
            endif
          enddo
        endif

!       v = 1
        nu = 1
        nu2 = 2*nu

        call yzktd (ncore(imj),nu,flps, pf(1,imj),qf(1,imj), apsm(1,imj), aqsm(1,imj),ykgr, 0)

!       ykgr is Y(r) = r * int dr' U(r,r') Psi_k'(r') R_l'(r') 
!       where R_l(r) is pf (qf)
!       Psi_k(r) is dgc0(ncore) (dpc0)
!       Y(B,D,r), where B = k', D = l', see Grant's paper

        do 35 imi = 1, idim
!       do 35 imi = 1, matsize
          if (dble(p2m(imi)).le.0.d0) goto 35

          do i=1,ilast
!c          calculate int Yk(r) /r * (Psi_k(r)*R_l(r)+qc*qf) in Grant
            var(i)=dble(pc(i,imi)*pf(i,imi) + qc(i,imi)*qf(i,imi)) 
            var(i) = var(i) / ri(i) * dble(ykgr(i)) * sfx
!           add xc term here
            if (kinitm(imi).eq.kinitm(imj) .and. ifxc.ne.2) then
              var(i) = var(i) + fxc0(i) *                               &
     &         dble(pc(i,imi)*pf(i,imi) + qc(i,imi)*qf(i,imi))*         &
     &         dble(pc(i,imj)*pf(i,imj) + qc(i,imj)*qf(i,imj))
            elseif (kinitm(imi).gt.0 .or. ifxc.eq.2) then
              var(i) = var(i) + (fxc(i)+coni*fxcim(i)) *                &
     &         dble(pc(i,imi)*pf(i,imi) + qc(i,imi)*qf(i,imi))*         &
     &         dble(pc(i,imj)*pf(i,imj) + qc(i,imj)*qf(i,imj))
            else
              var(i) = var(i) + (fxc(i)-coni*fxcim(i)) *                &
     &         dble(pc(i,imi)*pf(i,imi) + qc(i,imi)*qf(i,imi))*         &
     &         dble(pc(i,imj)*pf(i,imj) + qc(i,imj)*qf(i,imj))
            endif
          enddo
              
!c        integration by trapezoid method
!         this gives R_ABCD (see Grant's review)

!         A = k, C = l
          rabcd = 0
          do i=2,ilast
            rabcd=rabcd + (var(i)+var(i-1))*(ri(i)-ri(i-1)) / 2
          enddo

          do i=1,ilast
!c          calculate int ( Yk(r) /r * (Psi_k(r)*R_l(r)+qc*qf) in Grant 
            var(i)=dble(pc(i,imi)*ptot(i,imi) + qc(i,imi)*qtot(i,imi)) 
            var(i) = var(i) / ri(i) * dble(ykgr(i)) * sfx
!           add xc term here
            if (kinitm(imi).eq.kinitm(imj) .and. ifxc.ne.2) then
              var(i) = var(i) + fxc0(i) *                               &
     &         dble(pc(i,imi)*ptot(i,imi) + qc(i,imi)*qtot(i,imi)) *    &
     &         dble(pc(i,imj)*pf(i,imj) + qc(i,imj)*qf(i,imj))
            elseif (kinitm(imi).gt.0 .or. ifxc.eq.2) then
              var(i) = var(i) + (fxc(i)+coni*fxcim(i)) *                &
     &         dble(pc(i,imi)*ptot(i,imi) + qc(i,imi)*qtot(i,imi)) *    &
     &         dble(pc(i,imj)*pf(i,imj) + qc(i,imj)*qf(i,imj))
            else
              var(i) = var(i) + (fxc(i)-coni*fxcim(i)) *                &
     &         dble(pc(i,imi)*ptot(i,imi) + qc(i,imi)*qtot(i,imi)) *    &
     &         dble(pc(i,imj)*pf(i,imj) + qc(i,imj)*qf(i,imj))
            endif
          enddo
              
!c        integration by trapezoid method
!         this gives R_ABCD (see Grant's review)

!         A = k, C = l
          rabcdp = 0
          do i=2,ilast
            rabcdp=rabcdp + (var(i)+var(i-1))*(ri(i)-ri(i-1)) / 2
          enddo

!         d^v prefactors

!         ja2 = jinit(imi)
!         ma2 = minit(imi)
!         jc2 = jfin(imi)
!         mc2 = mfin(imi)
!         jd2 = jinit(imj)
!         md2 = minit(imj)
!         jb2 = jfin(imj)
!         mb2 = mfin(imj)
          jc2 = jinit(imi)
          mc2 = minit(imi)
          ja2 = jfin(imi)
          ma2 = mfin(imi)
          jb2 = jinit(imj)
          mb2 = minit(imj)
          jd2 = jfin(imj)
          md2 = mfin(imj)

          if( (ma2+mb2) .ne. (mc2+md2) ) goto 35 

          ind2 = 2

!         Grant's condition (7.13)
!         ja+v+jc is odd if Aa = Ac, even otherwise
!         jd+v+jb is odd if Ad = Ab, even otherwise
!         in our case it turns out j+v+j' should be even => v is odd
!         v =1 is the leading term

          dnuca= cwig3j(ja2,nu2,jc2,1,0,2)
          dnuca= dnuca * cwig3j(ja2,nu2,jc2,-ma2,(ma2-mc2),ind2)

          dnudb= cwig3j(jd2,nu2,jb2,1,0,2)
          dnudb= dnudb * cwig3j(jd2,nu2,jb2,-md2,(md2-mb2),ind2)

          prefac = (-1)**((ma2 + 1)/2.0) * (-1)**((md2 + 1)/2.0)
          prefac = prefac * dnuca * dnudb
          x1 = (ja2+1) * (jd2+1)
          x2 = (jc2+1) * (jb2+1)
          prefac = prefac * sqrt(x1 * x2) 
!         multiply by separation function
          xkmat(imi,imj) = xkmat(imi,imj) + prefac * rabcd 
          xkmatp(imi,imj) = xkmatp(imi,imj) + prefac * rabcdp

          if (ifxc.eq.5 .and. kinitm(imi).ne.kinitm(imj)) then
!           add dominant nonlocal exchange term (nu=2)
            nu = 2
            call yzktd (ncore(imj),nu,flps, pf(1,imj),qf(1,imj), apsm(1,imj), aqsm(1,imj),ykgrex, ncore(imi))

            do i=1,ilast
!c            calculate int Yk(r) /r * (Psi_k(r)*R_l(r)+qc*qf) in Grant
              var(i)=dble(pf(i,imj)*pf(i,imi) + qf(i,imj)*qf(i,imi)) 
              var(i) = var(i) / ri(i) * dble(ykgrex(i)) * sfx
            enddo
              
!c          integration by trapezoid method
!           this gives R_ABCD (see Grant's review)
!           A = k, C = l
            rabcd = 0
            do i=2,ilast
              rabcd=rabcd + (var(i)+var(i-1))*(ri(i)-ri(i-1)) / 2
            enddo

            do i=1,ilast
!c            calculate int ( Yk(r) /r * (Psi_k(r)*R_l(r)+qc*qf) in Grant 
              var(i)=dble(pf(i,imj)*ptot(i,imi) + qf(i,imj)*qtot(i,imi)) 
              var(i) = var(i) / ri(i) * dble(ykgrex(i)) * sfx
            enddo
              
!c          integration by trapezoid method
!           this gives R_ABCD (see Grant's review)
!           A = k, C = l
            rabcdp = 0
            do i=2,ilast
              rabcdp=rabcdp + (var(i)+var(i-1))*(ri(i)-ri(i-1)) / 2
            enddo

!           d^v prefactors

            jd2 = jinit(imi)
            md2 = minit(imi)
            ja2 = jfin(imi)
            ma2 = mfin(imi)
            jb2 = jinit(imj)
            mb2 = minit(imj)
            jc2 = jfin(imj)
            mc2 = mfin(imj)
            if( (ma2+mb2) .ne. (mc2+md2) ) goto 35 

            ind2 = 2

!           Grant's condition (7.13)
!           ja+v+jc is odd if Aa = Ac, even otherwise
!           jd+v+jb is odd if Ad = Ab, even otherwise
!           in our case it turns out j+v+j' should be even => v is odd
!           v =2 is the leading term

            nux = 2*nu
            dnuca= cwig3j(ja2,nux,jc2,1,0,2)
            dnuca= dnuca * cwig3j(ja2,nux,jc2,-ma2,(ma2-mc2),ind2)
  
            dnudb= cwig3j(jd2,nux,jb2,1,0,2)
            dnudb= dnudb * cwig3j(jd2,nux,jb2,-md2,(md2-mb2),ind2)

            prefac = (-1)**((ma2 + 1)/2.0) * (-1)**((md2 + 1)/2.0)
            prefac = prefac * dnuca * dnudb
            x1 = (ja2+1) * (jd2+1)
            x2 = (jc2+1) * (jb2+1)
            prefac = prefac * sqrt(x1 * x2) 
            xkmat(imi,imj) = xkmat(imi, imj) - prefac * rabcd 
            xkmatp(imi,imj) = xkmatp(imi, imj) - prefac * rabcdp
          endif
!         fix later: need generalization
!         if (imi.gt.69)  xkmatp(imi-54,imj) = xkmatp(imi,imj)
!         if (imi.gt.15)  xkmatp(imi,imj) = 0
!c        first 15 for d->f transitions and next 9 for d->p
          if (nfo.gt.0) then
           if (imi.gt.matszp .and. imi.le.nfo*matszp) xkmatp(imi,imj)=0
           imp = imi - nfo*matszp
           if (imp.le.matszm.and.imp.gt.0)                              &
     &           xkmatp(imp+matszp,imj) = xkmatp(imi,imj)
           if (imp.gt.0) xkmatp(imi,imj) = 0
          else
           if (imi.gt.matszm)  xkmatp(imi,imj) = 0
          endif
            
!         test importance of non-diagonal elements
!         if (imi .ne. imj) xkmat(imi,imj) = 0

!         test mixing between two edges - no mixing
!         if (imi .le. 4 .and. imj .gt. 4) xkmat(imi,imj) = 0
!         if (imj .le. 4 .and. imi .gt. 4) xkmat(imi,imj) = 0

!         only mixing between L3-L2 (no L2-L2 or L3-L3)
!         if (imi .gt. 2 .and. imj .gt. 2) xkmat(imi,imj) = 0
!         if (imj .le. 2 .and. imi .le. 2) xkmat(imi,imj) = 0

!         xkmat(imi,imj) = 0.0


   35   continue
   30 continue   
!---------------------
                
      return 
      end 
                     

