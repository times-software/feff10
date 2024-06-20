!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: screensub.f90,v $:
! $Revision: 1.6 $
! $Author: jorissen $
! $Date: 2011/02/11 04:14:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================================================================
!     get_pe_dens - calculate photoelectron density as an integrol over
!     energy.
!=======================================================================
! inputs: geom.dat => 1st line, ldos.inp => 2nd, scrn.inp => 3rd
!         pot.bin => 4,5th, grids => 6th, atomic and potentials => 6,7th
! output: totden_pe
subroutine get_pe_dens( nat, nph, iphat, rat,                          &
     &                   toler1, toler2, lmaxph, rdirec, rfms2, lfms2,  &
     &                   maxl, irrh, iend, lfxc, nrptx0,                &
     &                   ihole, xmu, adgc, adpc, iz,                    &
     &                   xion, iunf, xnval,edens, dgc0, dpc0,           &
     &                   ri, dx, x0, rmt, rnrm, em, ne,                 &
     &                   ixc0, iph, vtot, vvalgs, eref, dgcn, dpcn, ph, &
     &                   ilast, vch, totden_pe, rcutin, l_pe)
  !     linear response calculation
  !     calculate static screening

  use DimsMod, only: nrptx, lx, natx, nphx=>nphu, ltot, nex, nspx=>nspu
  use constants
  implicit none
!Changed the dimensions to 41 to account for superheavy elements. Pavlo Baranov 07/2016

  external dgetrf, dgetrs
  !     function getiat
  integer getiat
  !     pot.bin --------------------------------------------------------
  integer ihole, iz, iunf, l_pe
  double precision xmu, adgc(10,41), adpc(10,41)
  double precision xion, xnval(41), edens(251)
  double precision dgc0(251), dpc0(251)
  !     -----------------------------------------------------------------
  double precision vtot(nrptx), vvalgs(nrptx)
  double precision dgcn(nrptx,41), dpcn(nrptx,41)
  complex*16 eref
  !

  integer ixc0, ispin
  !     -----------------------------------------------------------------
  !     radial grid
  double precision rmt, rnrm, rcut, rcutin, rcut0, rnew, rpnew, x0, dx, ri(nrptx)
  integer kinit, linit, jri, jri1, jnrm, ilast, iph, ic3,irr, ikap
  integer i, j, k, m, n, ir1, ir2, ll, ll2, lfin, ilp, ie, ncycle, ios
  !     work space for xcpot
  complex*16 v(nrptx), vval(nrptx)
  !     work space for fovrg
  complex*16 pr(nrptx,0:lx), qr(nrptx,0:lx), pn(nrptx,0:lx), qn(nrptx,0:lx)
  complex*16 p2, ck, xkmt, xkmtp, xck
  complex*16 pu, qu, dum1, factor, de
  complex*16 xfnorm, xirf, temp, ph0
  complex*16 bessj(8), bessh(8)
  complex*16 jl, jlp1, nl, nlp1

  !     fms routine (variables are in single precision!)
  integer nat
  real rfms2, rdirec, toler1, toler2, rpart, aipart, rat(3,natx)
  integer il, il2, ix, ix2, im, im2, nsp, minv, iverb, nph, ipp, ill, lfms2

  complex cks(nspx)
  complex conis
  parameter ( conis = (0.0,1.0) )

  integer lmaxph(0:nphx)
  integer inclus, iphat(natx)
  !     phase shift      
  complex*16 ph(nex, ltot+1, 0:nphx)

  !     Linear Response
  !     scrn.inp -------------------------------------------------------
  integer maxl, irrh, iend, lfxc, nrptx0
  !     -----------------------------------------------------------------
  integer nrx
  parameter ( nrx = 251 )
  double precision percent
  !     energy mesh
  complex*16 em(nex)
  integer ne
  !     matrix inversion
  character*1 trans
  parameter ( trans = 'N' )
  integer info, nrhs, lda, ldb
  integer ipiv(nrx), ll_pe
  !     work space
  !     core hole potential
  double precision vch(nrx), U_Hub, n_occ, UTot, vbare(nrx), U_Bare
  !     Coulomb Exchange      
  double precision Kmat(nrx, nrx)
  !     LDA fxc      
  double precision fxc(nrx)
  !     chi_0(r,r')
  complex*16 chi0r(nrx, nrx)
  !     chi_0(r,r',e)
  complex*16 chi0re(nrx, nrx), pratom, pnatom, gcProj
  !     (1 - K Chi0)
  double precision Amat(nrx, nrx)
  !     W, the result of calculation
  double precision wscrn(nrx), den_pe(nrx,nex), totden_pe(nrx), emin, emax, atnorm, xnel

  complex*16, allocatable :: gtrl(:,:,:)
  complex, allocatable :: gg(:,:,:), gtr(:,:,:,:), xphase(:,:,:)
  logical, allocatable :: lcalc(:)
  logical UseProjection
  character*300 message
  ! Allocate local variables
  allocate(gtrl(0:lx, 0:lx, nex))
  allocate(gg(nspx*(lx+1)**2,nspx*(lx+1)**2,0:nphx))
  allocate(gtr(0:lx, 0:lx, 0:nphx, nex))
  allocate(xphase(nspx, -lx:lx, 0:nphx))
  allocate(lcalc(0:lx))
  UseProjection = .TRUE.
  rcut = rnrm*rcutin !wscrn(1)*rnrm
  rcut0 = rnrm
  PRINT*, 'rcut/rcut0', rcut/rcut0
  emin = -1.d30
  emax = -1.d30
  UTot = 0.d0
  xnel = 0.d0
  PRINT*, 1
  den_pe = 0.d0
  PRINT*, 2
  totden_pe = 0.d0
  PRINT*, 3
  gtrl(:,:,:) = 0.d0 ! JK Need to initialize for ifort. 04/2010
  PRINT*, 4
  !     =================================================================      
  call setkap(ihole, kinit, linit)

  !     set imt and jri (use general Loucks grid)
  !     rmt is between imt and jri (see function ii(r) in file xx.f)
  jri   = getiat(x0, dx, rmt) + 1
  jri1  = jri+1
  if ( jri1 .gt. nrptx )  call par_stop('jri .gt. nrptx in phase')
  !     nesvi: define jnrm
  jnrm  = getiat(x0, dx, rnrm) + 1
  !      ilast is the last integration point
  !     it is larger than jnrm for better interpolations
  ilast = jnrm + 6 + iend
  if (ilast .gt. nrx) ilast = nrx
  
  !     initialize values to zero
  do i = 1, nrx
     do j = 1, nrx
        chi0r (i,j) = 0.0d0
     end do
  end do
  percent  = 0

  call yprep(iph, nat, inclus, nph, iphat, rfms2, rat)

  if (inclus .gt. 1) then
     write(message,'(a,i4,a,i3)') ' Doing FMS for a cluster of ', &
          &         inclus, ' atoms around iph = ', iph
     call wlog(message)
  end if

  !===============================================================
  !     start energy cycle (big loop!)
  !===============================================================
  do i = 1, nrptx
     v(i)    = vtot(i)
     vval(i) = vvalgs(i)
  end do
  PRINT*, 2
  do ie = 1, ne
     chi0re = 0.0d0
     !       set the method to calculate atomic cross section
     !       p2 is (complex momentum)**2 referenced to energy dep xc
     p2   = em(ie) - eref
     ck   = sqrt(2*p2 + (p2*alphfs)**2)
     rpart  = real( dble(ck))
     aipart = real(dimag(ck))
     cks(1) = cmplx(rpart, aipart)
     xkmt = rmt * ck
     !write(99,*) dreal(em(ie))*hart, dimag(em(ie))*hart, ie

     if ( mod(ixc0,10) .lt. 5 ) then
        ncycle = 0
     else
        !         fix later . may be ncycle can be less
        ncycle = 3
     endif
     
     !=============================================================
     !         FMS
     !=============================================================
     if (inclus .gt. 1) then
        do ipp = 0, nph
           do ill = -lmaxph(ipp), lmaxph(ipp)
              rpart  = dble( ph(ie, abs(ill)+1, ipp))
              aipart = dimag(ph(ie, abs(ill)+1, ipp))
              xphase(1, ill, ipp) = cmplx(rpart, aipart)
           end do
        end do
        
        iverb = 0
        nsp = 1
        ispin = 0
        do i = 0, lx
           lcalc(i) = .true.
        end do
        minv = 0

        gtr(:,:,:,ie) = 0.d0
        
        call fms(lfms2, nsp, ispin, inclus, nph, cks,lmaxph, &
             & xphase,ie,iverb,minv,rdirec,toler1,toler2,lcalc, gg)
        
        do il = 0, lmaxph(iph)
           do il2 = il, il
              ix = il**2
              ix2 = il2**2
              do im = 1, 2*il+1
                 do im2 = im, im
                    gtr(il,il2,iph,ie) = gtr(il,il2,iph,ie) + gg(ix+im,ix2+im2,iph)
                 end do
              end do
              gtrl(il,il2,ie) = ( dble(real( gtr(il,il2,iph,ie))) &
                & + coni*dble(aimag(gtr(il,il2,iph,ie))) )&
                & * exp(coni*(ph(ie,il+1,iph)+ph(ie,il2+1,iph)))/SQRT((2*il+1.0d0)*(2*il2+1.d0))
              IF(il.ne.il2) gtrl(il,il2,ie) = 0.d0
           end do
        end do
     end if

     !       the maximum angular momentum        
     ll_pe = l_pe ! choose state
     do k = ll_pe + 1, ll_pe + 1 ! JK - only need pe wavefunction.

        ll = k - 1
        ikap = -k

        ic3 = 1
        if (ll .eq. 0) ic3 = 0

        !=========================================================
        !         get ''regular'' solution
        !=========================================================
        irr = -1
         
        call dfovrg( ncycle, ikap, rmt, ilast, jri, p2, dx,           &
             &                  ri , v, vval, dgcn, dpcn,  adgc, adpc, xnval,   &
             &                  pu, qu, pr(:,ll), qr(:,ll), iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph

        call exjlnl(xkmt, ll, jl, nl)
        call exjlnl(xkmt, ll+1, jlp1, nlp1)
        call phamp(rmt,pu,qu,ck,jl,nl,jlp1,nlp1,ikap,ph0,temp)

        factor = ck*alphfs
        factor = - factor/(1+sqrt(1+factor**2))
        dum1   = 1 / sqrt(1+factor**2)
        IF(temp.ne.0.d0) THEN
           xfnorm = 1.d0 / temp * dum1
        ELSE
           xfnorm = 0.d0
        END IF

        !         normalization factor
        !         dum1 is relativistic correction to normalization
        !         normalize regular solution
        do i = 1, ilast
           pr(i,ll) = pr(i,ll)*xfnorm
           qr(i,ll) = qr(i,ll)*xfnorm
        end do
        !===========================================================
        !         get ''irregular'' solution
        !===========================================================
        irr = 1
        !         pu, qu are upper and lower components at muffin tin
        !         set pu, qu - initial condition for irregular solution

        pu =   (nl*cos(ph0) +   jl*sin(ph0)) * rmt * dum1
        qu = (nlp1*cos(ph0) + jlp1*sin(ph0)) * rmt * dum1 * factor

        if (irrh .eq. 1) then
           call besjh(xkmt, maxl+1, bessj, bessh)
           pu = bessh(ll+1)* exp(coni*ph0) *rmt * dum1
           qu = bessh(ll+2)* exp(coni*ph0) *rmt * dum1* factor
        end if

        call dfovrg( ncycle, ikap, rmt, ilast, jri, p2, dx,           &
             &                  ri, v, vval, dgcn, dpcn, adgc, adpc, xnval,     &
             &                  pu, qu, pn(:,ll), qn(:,ll), iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph

        !         set N- irregular solution , which is outside
        !         N=(nlp1*cos(ph0)+jlp1*sin(ph0))*factor *rmt * dum1
        !         N = i*R - H*exp(i*ph0)
        temp = exp(coni*ph0)
        !         calculate wronskian (qu)
        qu = 2*alpinv *temp*(pn(jri,ll)*qr(jri,ll)-pr(jri,ll)*qn(jri,ll))
        IF(qu.ne.0.d0) THEN
           qu = 1 /qu / ck
        ELSE
           qu = 0.d0
        END IF
        !         qu should be close to 1
        do i = 1, ilast
           pn(i,ll) = temp * pn(i,ll)*qu
           qn(i,ll) = temp * qn(i,ll)*qu
        end do

        !         Use exact solution to continue solutions beyond rmt
        do j = jri, ilast
           xck  = ck * ri(j)
           call exjlnl(xck, ll, jl, nl)
           call exjlnl(xck, ll+1, jlp1, nlp1)
           call besjh(xck, maxl+1, bessj, bessh)

           pr(j,ll) =  (jl*cos(ph0) -   nl*sin(ph0))*ri(j)*dum1
           qr(j,ll) =(jlp1*cos(ph0) - nlp1*sin(ph0))*ri(j)*dum1*factor
           pn(j,ll) = bessh(ll+1)*exp(coni*ph0) *ri(j)*dum1
           qn(j,ll) = bessh(ll+2)*exp(coni*ph0) *ri(j)*dum1*factor
        end do
     end do

     IF(xnel.GT.10.d0) THEN
        EXIT
     END IF
     PRINT*, xmu, DBLE(em(ie)), xnel
     
     if (ie .eq. 1) then
        de = (em(ie+1) - em( ie ))/2.0d0
     else if (ie .eq. ne) then
        de = (em( ie ) - em(ie-1))/2.0d0
     else
        de = (em(ie+1) - em(ie-1))/2.0d0
     end if
     do m = 1, ilast
        ! Get the spectral density of angular momentum ll_pe
        den_pe(m,ie) = den_pe(m,ie) + DIMAG( (pr(m,ll_pe)*pn(m,ll_pe) &
       &   + pr(m,ll_pe)**2*gtrl(ll_pe,ll_pe,ie))*ck*4)*(2*ll_pe+1.0d0)/pi
     end do

     do ir1 = 1, ilast
        do i = ir1, ilast
           chi0r( ir1,i) = chi0r(ir1,i) + chi0re(ir1,i)*de
           chi0re(ir1,i) = (0.0d0, 0.0d0)
        end do
     end do

     !       print the current progress to the wscrn
     if ( 100.0*ie/ne .gt. percent) then

        write(message,'(a,f5.1,a)') '  ', percent, &
            &     '% of energy integration'
        call wlog(message)
        percent = percent + 10.0
     end if
     
     IF(ie.GT.1) THEN
        do m = 1, ilast
           totden_pe(m) = totden_pe(m) + den_pe(m,ie)*de
           IF(m.le.jnrm) xnel = xnel + den_pe(m,ie)*de*ri(m)*dx
        end do
     END IF
  end do
  call wlog('100.% -- end of energy integration --')
  !=================================================================
  !     end of energy cycle
  !=================================================================

  ! Normalize  total pe-density. Using den_pe(1,1) to store normalization.
  den_pe(1,1) = 0.d0
  DO i = 1, ilast
     IF(UseProjection) THEN
        rnew = MIN(MAX(rcut0,ri(i)),rcut)
        rnew = rnew - rcut0
        WRITE(77,*) ri(i), rnew
        totden_pe(i) = totden_pe(i)*COS(rnew/(rcut-rcut0)*pi/2.d0)**4
        den_pe(1,1) = den_pe(1,1) + totden_pe(i)*ri(i)*dx
     ELSE
        den_pe(1,1) = den_pe(1,1) + totden_pe(i)*ri(i)*dx
     END IF
  END DO
  totden_pe(:) = totden_pe(:)/den_pe(1,1)

  deallocate(gtrl,gg,gtr,xphase,lcalc)

  return
  !      end subroutine screen
end subroutine get_pe_dens
!=======================================================================
!     END wscrn
!=======================================================================
