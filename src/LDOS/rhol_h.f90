      subroutine rhol_h (iph, dx, x0, ri, ne, em,                       &
                       ixc, rmt, rnrm,                                 &
                       vtot, vvalgs, dgcn, dpcn, eref,                 &
                       adgc, adpc, xrhole, xmrhole, xrhoce,              &
                       xmrhoce, ph,aph,                                &
                       iz, xion, iunf, ihole, lmaxsc,          &
                       xnval,Vnlm,i_opt,iss,xmu)

      use constants
      use DimsMod, only: nrptx, nex, lx, nspx=>nspu
      use hubbard_inp
      implicit none
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
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
!     xrhole(0:lx,nex)  integral over r of density function
!     xrhoce(nex)  integral over r of density function for embedded atom

!     max number allowed in xsect r-grid
      integer, parameter :: nrx = nrptx

!     output
      real*8, intent(out) ::   xrhoce(0:lx, nex)
      real*8, intent(out) ::   xmrhoce(0:lx, nspx*(lx+1)**2,nex)
      complex*16, intent(out) ::  ph(nex,lx+1)
      complex*16, intent(out) ::  aph(nex,lx+1,(lx+1)**2)
!     input
      integer, intent(in) :: ixc, ihole, iz, lmaxsc, iunf, i_opt, iph, ne, iss
      real*8, intent(in) :: dx, x0, rmt, rnrm, xmu, xion
      !energy grid in complex e-plane
      complex*16, intent(in) ::  em(nex), eref
      real*8, intent(in) ::   ri(nrptx)
      real*8, intent(in) ::   vtot(nrptx), vvalgs(nrptx)
      real*8, intent(in) ::   Vnlm(0:lx,(lx+1)**2,2)
      real*8, intent(in) ::   xnval(41), dgcn(nrptx,41), dpcn(nrptx,41)
      real*8, intent(in) ::   adgc(10,41), adpc(10,41)
      real*8 ri05(251)
      complex*16  xrhole(0:lx, nex)
      complex*16  xmrhole(0:lx,nspx*(lx+1)**2, nex)
      complex*16  vtotc(nrptx), vvalc(nrptx)
      complex*16  vtotc_tmp(nrptx), vvalc_tmp(nrptx)
!     work space for dfovrg: regular and irregular solutions
      complex*16 pr(nrx), qr(nrx), pn(nrx), qn(nrx)

      complex*16  p2, xkmt, ck, xck
      complex*16  pu, qu
      complex*16  xfnorm, xirf
      complex*16  temp,  phx
      complex*16 jl,jlp1,nl,nlp1
      complex*16  xpc(nrx)
      character*30 f1name
      integer jri, imt, ie, ncycle, jnrm, inrm, i, lmax, ilast, lll, imm, ikap, irr, ic3, ios, &
         mmm, i0, ir00
      real*8 xmuO, dx05

      xmuO=xmu-fermi_shift/hart
!     make  radial grid with 0.05 step
      dx05=0.05d0
      do i=1,251
         ri05(i) = exp(-8.8+dx05*(i-1))
      enddo
      lmax=lmaxsc
      if (lmax.gt.lx) lmax = lx
!KJ    Why the next 2 statements?  Commenting out.
!      if (iz.le.4) lmax=2
!      if (iz.le.2) lmax=1
      vtotc=vtot
      vvalc=vvalgs
!     set imt and jri (use general Loucks grid)
!     rmt is between imt and jri (see function ii(r) in file xx.f)
      imt  = (log(rmt) + x0) / dx  +  1
      jri  = imt+1
      if (jri .gt. nrptx)  call par_stop('jri .gt. nrptx in phase')
      inrm = (log(rnrm) + x0) / dx  +  1
      jnrm = inrm+1
!     ilast is the last integration point
!     it is larger than jnrm for better interpolations
      ilast = jnrm + 6
      if (ilast.gt.nrptx) ilast = nrptx

      if (lmax.lt.lx) then
          xmrhoce(lmax+1:lx,:,:)=0
          xmrhole(lmax+1:lx,:,:)=0
          ph(:,lmax+1:lx)=0
          xrhoce(lmax+1:lx,:) = 0
          xrhole(lmax+1:lx,:) = 0
      endif
      if ((lmax+1).lt.lx) then
          aph(:,lmax+2:lx,:) = 0
      endif

      do ie = 1, ne
!       p2 is (complex momentum)**2 referenced to energy dep xc
        p2 = em(ie) - eref
        ! if (mod(ixc,10) .lt. 5) then
        if ((mod(ixc,10).lt.5).OR.(ixc.EQ.6).OR.(ixc.EQ.7)) then ! TTS
          ncycle = 0
        else
          ncycle = 3
        endif
        ck = sqrt(2*p2+ (p2*alphfs)**2)
        xkmt = rmt * ck
        do lll=0,lmax
           !KJ changed the loop upper index to "lmax"
           ! For higher values of lll, used to just set to 0 and skip to end of loop.
           ! I now do that before we start the loop.
           ! Cleaner and more efficient.

            ! Begin(Towfiq): filling ph(ie,l,m,s) 05/04/2010
            if(i_opt.eq.2) then
                 do imm=(lll**2+1),(lll+1)**2
                      ikap = -1-lll
                      irr = -1
                      ic3 = 1
                      if (lll.eq.0) ic3 = 0
                      do ir00=1,nrptx
                          vtotc_tmp(ir00)=vtotc(ir00)+1.0*Vnlm(lll,imm,iss)
                          vvalc_tmp(ir00)=vvalc(ir00)+1.0*Vnlm(lll,imm,iss)
                      enddo
                      call dfovrg ( ncycle, ikap, rmt, ilast, jri, p2, dx, ri, vtotc_tmp, vvalc_tmp, dgcn, dpcn, adgc, adpc,    &
                                xnval, pu, qu, pn, qn, iz, ihole, xion, iunf, irr, ic3,iph)
                      call exjlnl (xkmt, lll, jl, nl)
                      call exjlnl (xkmt, lll+1, jlp1, nlp1)
                      call phamp (rmt, pu, qu, ck,  jl, nl, jlp1, nlp1, ikap, phx, temp)
                      aph(ie,lll+1,imm)=phx

                      xfnorm = 1 / temp
                      do i = 1,ilast
                        pr(i)=pn(i)*xfnorm
                        qr(i)=qn(i)*xfnorm
                      enddo
                      irr = 1
                      pu = ck*alphfs
                      pu = - pu/(1+sqrt(1+pu**2))
                      qu=(nlp1*cos(phx)+jlp1*sin(phx))*pu *rmt 
                      pu = (nl*cos(phx)+jl*sin(phx)) *rmt 
                      call dfovrg (ncycle, ikap, rmt, ilast, jri, p2, dx, ri, vtotc_tmp,vvalc_tmp, dgcn, dpcn, adgc, adpc,    &
                               xnval, pu, qu, pn, qn, iz, ihole, xion, iunf, irr, ic3,iph)
                      temp = exp(coni*phx)
                      qu = 2 * alpinv * temp * ( pn(jri)*qr(jri) - pr(jri)*qn(jri) )
                      qu = 1 /qu / ck
                      do i = 1, ilast
                         pn(i) = coni * pr(i) - temp * pn(i)*qu
                         qn(i) = coni * qr(i) - temp * qn(i)*qu
                      enddo
                      pu = ck*alphfs
                      pu = - pu/(1+sqrt(1+pu**2))
                      do i=jri,ilast
                         xck = ck * ri(i)
                         temp = xck
                         call exjlnl (temp, lll, jl, nl)
                         call exjlnl (temp, lll+1, jlp1, nlp1)
                         pr(i)= (jl*cos(phx)-nl*sin(phx)) *ri(i)
                         qr(i)=(jlp1*cos(phx)-nlp1*sin(phx))*pu *ri(i)
                         pn(i)= (nl*cos(phx)+jl*sin(phx)) *ri(i)
                         qn(i)=(nlp1*cos(phx)+jlp1*sin(phx))*pu *ri(i)
                      enddo
                      pu = ck*alphfs
                      pu = - pu/(1+sqrt(1+pu**2))
                      temp = (2*lll+1.0d0)/(1+pu**2) /pi *ck*4 / hart
                      do i = 1, ilast
                        xpc(i) = pr(i) * pr(i) + qr(i) * qr(i) 
                      enddo
                      xirf = lll*2 + 2
                      i0=jnrm+1
                      
                      call csomm2 (ri, xpc, dx, xirf, rnrm, i0)
                      xmrhole(lll,imm,ie) = xirf*temp
                      xrhole(lll,ie) = xirf*temp
                      do i = 1, ilast
                        xpc(i) = pn(i)*pr(i)-coni*pr(i)*pr(i) + qn(i)*qr(i)-coni*qr(i)*qr(i)
                      enddo
                      xirf =  1
                      call csomm2 (ri, xpc, dx, xirf, rnrm, i0)
                      xmrhoce(lll,imm,ie) =  - dimag(xirf* temp)
                      xrhoce(lll,ie) =  - dimag(xirf* temp)
                 enddo !imm
           endif ! i_opt=2
           ! End(Towfiq)

           if(i_opt.eq.1) then

              ikap = -1-lll
              irr = -1
              ic3 = 1
              if (lll.eq.0) ic3 = 0

              call dfovrg ( ncycle, ikap, rmt, ilast, jri, p2, dx,          &
                         ri, vtotc, vvalc, dgcn, dpcn, adgc, adpc,         &
                         xnval, pu, qu, pn, qn,                            &
                         iz, ihole, xion, iunf, irr, ic3,iph)
                
              call exjlnl (xkmt, lll, jl, nl)
              call exjlnl (xkmt, lll+1, jlp1, nlp1)
              call phamp (rmt, pu, qu, ck,  jl, nl, jlp1, nlp1, ikap, phx, temp)
              ph(ie,lll+1)=phx
          
              !     Begin-0
              !c    Normalize final state  at rmt to
              !c    rmt*(jl*cos(delta) - nl*sin(delta))
              xfnorm = 1 / temp
              !c        xfnorm = 1
              !c    normalize regular solution
              do i = 1,ilast
                pr(i)=pn(i)*xfnorm
                qr(i)=qn(i)*xfnorm
              enddo
              !c     find irregular solution
              irr = 1
              pu = ck*alphfs
              pu = - pu/(1+sqrt(1+pu**2))
              !set pu, qu - initial condition for irregular solution at ilast
              !qu=(nlp1*cos(phx)+jlp1*sin(phx))*pu *rmt
              !pu = (nl*cos(phx)+jl*sin(phx)) *rmt
              qu=(nlp1*cos(phx)+jlp1*sin(phx))*pu *rmt 
              pu = (nl*cos(phx)+jl*sin(phx)) *rmt 
              call dfovrg (ncycle, ikap, rmt, ilast, jri, p2, dx, ri, vtotc,vvalc, dgcn, dpcn, adgc, adpc,            &
                       xnval, pu, qu, pn, qn, iz, ihole, xion, iunf, irr, ic3,iph)

              ! set N- irregular solution , which is outside
              ! N=(nlp1*cos(ph0)+jlp1*sin(ph0))*factor *rmt * dum1
              ! N = i*R - H*exp(i*ph0)
              temp = exp(coni*phx)
              ! calculate wronskian (qu)
              qu = 2*alpinv*temp*( pn(jri)*qr(jri) - pr(jri)*qn(jri) )
              qu = 1 /qu / ck
              !  qu should be close to 1
              do i = 1, ilast
                  !pr(i) = pr(i)*qu
                  !qr(i) = qr(i)*qu
                  pn(i) = coni * pr(i) - temp * pn(i)*qu
                  qn(i) = coni * qr(i) - temp * qn(i)*qu
              enddo

              !c    ATOM,  dgc0 is large component, ground state hole orbital
              !c    .      dpc0 is small component, ground state hole orbital
              !c    FOVRG, p    is large component, final state photo electron
              !c    .      q    is small component, final state photo electron
                
              !   Use exact solution to continue solutions beyond rmt
              pu = ck*alphfs
              pu = - pu/(1+sqrt(1+pu**2))
              do i=jri,ilast
                 xck = ck * ri(i)
                 temp = xck
                 call exjlnl (temp, lll, jl, nl)
                 call exjlnl (temp, lll+1, jlp1, nlp1)
                 pr(i)= (jl*cos(phx)-nl*sin(phx)) *ri(i)
                 qr(i)=(jlp1*cos(phx)-nlp1*sin(phx))*pu *ri(i)
                 pn(i)= (nl*cos(phx)+jl*sin(phx)) *ri(i)
                 qn(i)=(nlp1*cos(phx)+jlp1*sin(phx))*pu *ri(i)
              enddo

              ! combine all constant factors to temp
              ! add relativistic correction to normalization, factor 2*lll+1,
              ! 2*ck for G.F., factor 2 for spin, and hart to transform to eV
              pu = ck*alphfs
              pu = - pu/(1+sqrt(1+pu**2))
              temp = (2*lll+1.0d0)/(1+pu**2) /pi *ck*4 / hart
              do  i = 1, ilast
                xpc(i) = pr(i) * pr(i) + qr(i) * qr(i) 
              enddo
              
              xirf = lll*2 + 2
              ! i0 should be less or equal to  ilast
              i0=jnrm+1
              call csomm2 (ri, xpc, dx, xirf, rnrm, i0)
              xrhole(lll,ie) = xirf*temp

              ! only central atom contribution needs irregular solution
              do i = 1, ilast
                xpc(i) = pn(i)*pr(i)-coni*pr(i)*pr(i) + qn(i)*qr(i)-coni*qr(i)*qr(i)
              enddo

              xirf =  1
              call csomm2 (ri, xpc, dx, xirf, rnrm, i0)
              xrhoce(lll, ie) =  - dimag(xirf* temp)
            endif ! i_opt=1

        enddo ! lll loop over angular momentum index
      enddo ! ie big loop over energy


      if (i_opt.eq.2) then
         !KJ added the "if" :
         ! this file is not read again until the xsph module
         ! so, while it is harmless to write it in the first pass of rhol_h,
         ! it is also completely useless.
 2011    format('aphase_dos', i2.2, '.dat')
         write(f1name,2011) iph
         open (unit=27, file=f1name, status='unknown', iostat=ios)   
         call chopen (ios,'aphase.dat','phase')
         do ie=1,ne
            do lll=0,lx
               do mmm=(lll**2+1),(lll+1)**2
                  write(27,*) dble(em(ie))*hart,lll,mmm,aph(ie,lll+1,mmm)
               enddo
            enddo
         enddo
         close(27)
      endif

      return
      end

