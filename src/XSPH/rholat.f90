!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: rholat.f90,v $:
! $Revision: 1.6 $
! $Author: bmattern $
! $Date: 2012/07/01 16:29:35 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rholat ( icount, dx, x0, ri, em,                       &
     &                  ixc, rmt, rnrm,                                 &
     &                  vtot, vvalgs, xnval, iorb, dgcn, dpcn, eref,    &
     &                  adgc, adpc, xrhole, xrhoce, ph,                 &
     &                  iz, xion, iunf, ihole, lmaxsc, iph) !KJ iph

      use constants
      use DimsMod, only: nrptx, lx
      implicit double precision (a-h, o-z)

!     INPUT
!     dx, x0, ri(nr)
!                  Loucks r-grid, ri=exp((i-1)*dx-x0)
!     ne, em(ne)   number of energy points,  complex energy grid
!     ixc          0  Hedin-Lunqist + const real & imag part
!                  1  Dirac-Hara + const real & imag part
!                  2  ground state + const real & imag part
!                  3  Dirac-Hara + HL imag part + const real & imag part
!                  5  Dirac-Fock exchange with core electrons +
!                     ixc=0 for valence electron density
!     rmt          r muffin tin
!     rnrm         r norman
!     vtot(nr)     total potential, including gsxc, final state
!     dgcn(dpcn)   large (small) dirac components for central atom
!     adgc(adpc)   their development coefficients
!
!     OUTPUT
!     xrhole(0:lx)  integral over r of density function
!     xrhoce(0:lx)  the same integral for embedded atom only

!     max number allowed in xsect r-grid
      parameter (nrx = nrptx)

!     output
      complex*16  xrhole(-4:3,-4:3)
      complex*16  xrhoce(-4:3, -4:3)


      complex*16, intent(inout) ::  ph(lx+1)

      dimension ri(nrptx), ri05(251)
      dimension  vtot(nrptx), vvalgs(nrptx)
      complex*16 vtotc(nrptx), vvalc(nrptx)
      dimension xnval(30), iorb(-4:3), dgcn(nrptx,30), dpcn(nrptx,30)
      dimension adgc(10,30), adpc(10,30)

!     energy grid in complex e-plane
      complex*16 em, eref

!     work space for dfovrg: regular and irregular solutions
      complex*16 pr(nrx,2,2), qr(nrx,2,2), pn(nrx,2,2), qn(nrx,2,2)

      complex*16  p2, xkmt, ck, xck
      complex*16  pu, qu
      complex*16  xfnorm, xirf, xmp, xpm
      complex*16  temp,  phx, phm(2,2), factor

      complex*16 jl,jlp1,nl,nlp1
      complex*16  xpc(nrx)

!     nesvi
      dimension pat(nrx,2,2),qat(nrx,2,2)
      complex*16 intr(nrx,2,2),var(nrx) 
      dimension xq(nrptx),xp(nrptx)
      complex*16 xkmi

!     initialize
      lmax=lmaxsc
      if (lmax.gt.lx) lmax = lx
      do 20 i = 1, nrptx
         vtotc(i)=vtot(i)
         vvalc(i)= vvalgs(i)
  20  continue
!     set imt and jri (use general Loucks grid)
!     rmt is between imt and jri (see function ii(r) in file xx.f)
      imt  = (log(rmt) + x0) / dx  +  1
      jri  = imt+1
      if (jri .gt. nrptx)  call par_stop('jri .gt. nrptx in phase')

      inrm = (log(rnrm) + x0) / dx  +  1
      jnrm = inrm+1

!     define ilast1,rlast
      rlast=rnrm
      if (icount.eq.2) rlast=10*rnrm
      jlast1=(log(rlast) + x0)/ dx + 2
      ilast1=jlast1 + 6

!c    nesvi
!c    dgcn and dpcn should be normalized <n|n>=1, check this here
     
      do 440 j = -4, 3
        jj = iorb(j)
        if (jj.le.0) goto 440

        do 420  i = 1, jlast1
         xp(i) = dpcn(i,jj)**2 + dgcn(i,jj)**2
         xq(i) = 0
  420   continue
!c      nb, xinorm is used for exponent on input to somm
        lfin = j
        if (j.lt.0) lfin = -j - 1
        xinorm = 2*lfin + 2
        i0 = jnrm + 1
        call somm2 (ri, xp, dx, xinorm, rnrm, 0, i0)
        if (xinorm.lt.0.99 .and. icount.eq.2) then
           call wlog                                                    &
     &     ('  WARNING: small overlap integral for Mulliken count')
        endif
    
        xinorm = 1.d0 / sqrt(xinorm)
        do 430 i=1,nrptx
          dpcn(i,jj)=dpcn(i,jj) * xinorm
          dgcn(i,jj)=dgcn(i,jj) * xinorm
  430   continue
  440 continue

!     set limits for tabulations
      nr05= (log(rnrm) + x0) / 0.05d0 + 5
      if (nr05.gt.251) nr05 = 251
!     ilast is the last integration point
!     it is larger than jnrm for better interpolations
      ilast = nint( (nr05-1) *0.05d0 / dx ) + 1
      if (ilast.gt.nrptx) ilast=nrptx

      if (ilast1.gt.nrptx) ilast1=nrptx

      do 10 i = -4,3
      do 10 j = -4,3
         xrhole(i,j) = 0
         xrhoce(i,j) = 0
  10  continue
      do 15 i=1,lx+1
  15  ph(i) = 0

!     p2 is 0.5*(complex momentum)**2 referenced to energy dep xc
!     need hartree units for dfovrg
      p2 = em - eref
      if (mod(ixc,10) .lt. 5) then
        ncycle = 0
      else
        ncycle = 3
      endif
      ck = sqrt(2*p2 + (p2*alphfs)**2)
      xkmt = rmt * ck

      do 200 lll=0,lmax
        do 199 jd = 0,1
          ikap = (lll+jd)* (-1)**jd
          if (ikap.eq.0) goto 199

          ilp = lll + 1
          if (ikap.gt.0) ilp = lll - 1
          im = 1+ jd

          do 150 j = 1, 2
            ic3 = j-1
            if (lll.eq.0 .and. ic3.eq.1) goto 150

            irr = -1
            call dfovrg ( ncycle, ikap, rmt, ilast1, jri, p2, dx,       &
     &                ri, vtotc, vvalc, dgcn, dpcn, adgc, adpc,         &
     &                xnval, pu, qu, pn(1,im,j), qn(1,im,j),            &
     &                iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph
            
            call exjlnl (xkmt, lll, jl, nl)
            call exjlnl (xkmt, ilp, jlp1, nlp1)
            call phamp (rmt, pu, qu, ck,  jl, nl, jlp1, nlp1, ikap,     &
     &                  phx, temp)
            if (lll.eq.0)  ph(1)=phx
            phm(im,j) = phx

!           Normalize final state  at rmt to
!           rmt*(jl*cos(delta) - nl*sin(delta))
            xfnorm = 1 / temp
!           normalize regular solution
            do 133  i = 1,ilast1
              pr(i,im,j)=pn(i,im,j)*xfnorm
              qr(i,im,j)=qn(i,im,j)*xfnorm
  133       continue

!-----------------------
!           nesvi            

!c           add solution beyond Rmt:
!             do 1010 i=jri+1, ilast1
!                xkmi=ri(i)*ck
!                call exjlnl(xkmi,lll,jl,nl)
!                pr(i,im,j)=(jl*cos(phx)-nl*sin(phx))*ri(i)
!                qr(i,im,j)=0.0d0
!1010         continue

!             chose atomic function for making projection.
!             Project on corresponding atomic states. 
              jj = iorb(ikap)

!             make corresponding atomic functions
              if (jj.eq.0) then
                do 397 i=1,nrptx
                  pat(i,im,j)=0
                  qat(i,im,j)=0
  397           continue
              else
                do 398 i=1,nrptx
                  pat(i,im,j)=dgcn(i,jj)
                  qat(i,im,j)=dpcn(i,jj)    
  398           continue
              endif

            open(unit=3,file='wfat.dat',status='unknown')    

!         only central atom contribution needs irregular solution
            do 194  i = 1, ilast1           
                write(3,1019) ri(i)/rnrm, dgcn(i,6),dgcn(i,8),          &
     &          dgcn(i,10),dgcn(i,12)
 1019           format(f10.5,1x,e10.4,1x,e10.4,1x,e10.4,1x,e10.4)

 194         continue

            close(3)




!             calculate overlap integral between f and atomic function
!             (integral Rl(r)*Psi_at(r)dr from 0 till r') 
!             intr(i) is that overlap integral. Later it
!             will be multiplied by pr(i)*Psi_at(r') and integrated till
!             r=infinity (ideal case), but actually till rlast.

              do 400 i=1,ilast1
                var(i)=pat(i,im,j)*pr(i,im,j)+qat(i,im,j)*qr(i,im,j)
!             factor of 2 -integration r< r>  -->2 r r'
  400         continue

!             integration by trapezoid method
              
              intr(1,im,j)=var(1)*ri(1)
   
              do 410 i=2,ilast1
                intr(i,im,j)=intr(i-1,im,j)+                            &
     &                       (var(i)+var(i-1))*(ri(i)-ri(i-1))
  410         continue 

!c              old way, no double integration 
!              do 415 i=1,ilast1
!                 intr(i,im,j)=intr(ilast1,im,j)/2.0                
!  415         continue    

              
!----------------


!          find irregular solution
            irr = 1
            pu = ck*alphfs
            factor = pu/(1+sqrt(1+pu**2))
            if (ikap.lt.0) factor = -factor
!           set pu, qu - initial condition for irregular solution at ilast
!           qu=(nlp1*cos(phx)+jlp1*sin(phx))*pu *rmt
!           pu = (nl*cos(phx)+jl*sin(phx)) *rmt
            qu=(nlp1*cos(phx)+jlp1*sin(phx))* factor *rmt 
            pu = (nl*cos(phx)+jl*sin(phx)) *rmt 

            call dfovrg (ncycle, ikap, rmt, ilast1, jri, p2, dx,        &
     &              ri, vtotc,vvalc, dgcn, dpcn, adgc, adpc,            &
     &              xnval, pu, qu, pn(1,im,j), qn(1,im,j),              &
     &              iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph
!c            set N- irregular solution , which is outside
!c            N=(nlp1*cos(ph0)+jlp1*sin(ph0))*factor *rmt * dum1
!c            N = i*R - H*exp(i*ph0)
              temp = exp(coni*phx)
              do i = 1, ilast
                pn(i,im,j) = coni * pr(i,im,j) - temp * pn(i,im,j)
                qn(i,im,j) = coni * qr(i,im,j) - temp * qn(i,im,j)
              enddo

 150      continue

!         combine all constant factors to temp
!         add relativistic correction to normaliz. and factor 2*lll+1
          temp = 2*ck / (1+factor**2) / pi

!         nesvi add irregular solution beyond Rmt          
!           do 1020 i=(jri+1), ilast1
!                xkmi=ri(i)*ck
!                call exjlnl(xkmi,lll,jl,nl)
!                pn(i,im,j)=(nl*cos(phx)+jl*sin(phx))*ri(i)
!                qn(i,im,j)=0.0d0
!1020        continue
     
!          open(unit=2,file='wfunc1.dat',status='unknown')    
!         ic3 = 0, j= ic3+1
          j = 1
!         calculate diagonal radial integrals R(k1,k1) - xrhoce and xrhole
            do 190  i = 1, ilast1
              xpc(i) = pr(i,im,j)*pat(i,im,j)*intr(i,im,j)+             &
     &              qr(i,im,j)*qat(i,im,j)*intr(i,im,j)

!            if (ikap .eq. -3 .and. (dble(em) +12.0/hart)
!     1          .le. 1.0/hart) then
!                write(2,1015) ri(i)/rnrm, dble(pr(i,im,j)),
!     1          pat(i,im,j), dble(intr(i,im,j)),dble(xpc(i))
! 1015           format(f10.6,1x,e10.4,1x,e10.4,1x,e10.4,1x,e10.4)
!             endif

 190        continue
            xirf = lll*2 + 2
            i0=jlast1+1
            call csomm2 (ri, xpc, dx, xirf, rlast, i0)
            xrhole(ikap,ikap) =xirf*temp*exp(coni*(phm(im,j)+phm(im,j)))

!            close(2)
            open(unit=2,file='wfunc.dat',status='unknown')    

!         only central atom contribution needs irregular solution
            do 195  i = 1, ilast1
              xpc(i) = pn(i,im,j)*pat(i,im,j)*intr(i,im,j)+             &
     &              qn(i,im,j)*qat(i,im,j)*intr(i,im,j)
              xpc(i) = xpc(i) -                                         &
     &              coni*(pr(i,im,j)*pat(i,im,j)*intr(i,im,j) +         &
     &              qr(i,im,j)*qat(i,im,j)*intr(i,im,j))

!         for test purposes
 
!           do 195  i = 1, ilast1
!              xpc(i) = pn(i,im,j)*pat(i,im,j)*intr(i,im,j)
!            xpc(i) = -1.0*coni*(pr(i,im,j)*pat(i,im,j)*intr(i,im,j))
           
             if (ikap .eq. 1 .and. (dble(em) +12.0/hart)                &
     &          .lt. 1.0/hart) then
                write(2,1016) ri(i)/rnrm, dble(pr(i,im,j)),             &
     &          pat(i,im,j), dble(intr(i,im,j)),-dimag(xpc(i))
 1016           format(f10.4,1x,e10.4,1x,e10.4,1x,e10.4,1x,e10.4)
             endif

 195        continue

            close(2)

            xirf =  1
            call csomm2 (ri, xpc, dx, xirf, rlast, i0)
            xrhoce(ikap,ikap) = - xirf * temp

!         calculate cross terms
          if (ikap.lt.-1) then
            k1 = ikap + 2*lll + 1
            do 290  i = 1, ilast1
              xpc(i) = pr(i,1,j)*pat(i,1,j)*intr(i,2,j) +               &
     &                 qr(i,1,j)*qat(i,1,j)*intr(i,2,j) 
 290        continue
            xirf = lll*2 + 2
!           i0 should be less or equal to  ilast
            i0=jlast1+1
            call csomm2 (ri, xpc, dx, xirf, rlast, i0)
!            xrhole (ikap, k1) = xirf*temp* exp(coni*(phm(1,j)+phm(2,j)))
!             xrhoce(ikap,k1)=0.0d0           
             xrhole (k1, ikap) = xrhole (ikap, k1)
!            nesvi: checked that cross-terms are not important for N_h 
            
!           ic3 = 1, j= ic3+1
            j = 2
            xpm =  exp(coni*(phm(1,j)-phm(2,j))) / 2
            xmp =  exp(coni*(phm(2,j)-phm(1,j))) / 2
            do 295  i = 1, ilast1
              xpc(i) = (pn(i,1,j)*pat(i,1,j)*intr(i,2,j)+               &
     &                  qn(i,1,j)*qat(i,1,j)*intr(i,2,j)) * xmp +       &
     &                 (pn(i,2,j)*pat(i,2,j)*intr(i,1,j)+               &
     &                  qn(i,2,j)*qat(i,2,j)*intr(i,1,j)) * xpm
              xpc(i) = xpc(i) - coni*(xpm+xmp) *                        &
     &                 (pr(i,1,j)*pat(i,1,j)*intr(i,2,j) +              &
     &                  qr(i,1,j)*qat(i,1,j)*intr(i,2,j))
 295        continue
            xirf =  1
            call csomm2 (ri, xpc, dx, xirf, rlast, i0)
            xrhoce(ikap,k1) = - xirf * temp
!        cross term not important for N_h
!            xrhoce(ikap,k1)=0.0d0
            xrhoce(k1,ikap) =  xrhoce(ikap,k1)
          endif
 199    continue 
 200  continue 


           
          if ((dble(em) +12.0/hart) .lt. 1.0/hart) then
          open(unit=4,file='xrhocet.dat',status='unknown')  
          do 1195  i=-4,3
              do 1195 j=-4,3
                write(4,1018) i,j,dimag(xrhoce(i,j))
 1018           format(i3,1x,i3,1x,f10.4)
 1195        continue
          close(4)
          endif
         



!     calculate phase shift in old way (ic3=1) test new one
!     which is commented out above later
      do 300 lll = 1,lmax
          im = 1
          ikap = -lll-1
          irr = -1
          ic3 = 1
          call dfovrg ( ncycle, ikap, rmt, ilast1, jri, p2, dx,         &
     &                ri, vtotc, vvalc, dgcn, dpcn, adgc, adpc,         &
     &                xnval, pu, qu, pr(1,im,1), qr(1,im,1),            &
     &                iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph
            
          call exjlnl (xkmt, lll, jl, nl)
          call exjlnl (xkmt, lll+1, jlp1, nlp1)
          call phamp (rmt, pu, qu, ck,  jl, nl, jlp1, nlp1, ikap,       &
     &                  phx, temp)
          ph(1+lll)=phx
 300  continue

      return
      end
