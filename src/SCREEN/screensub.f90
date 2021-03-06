!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: screensub.f90,v $:
! $Revision: 1.8 $
! $Author: jorissen $
! $Date: 2011/12/10 23:17:11 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================================================================
!     wscrn
!=======================================================================
! inputs: geom.dat => 1st line, ldos.inp => 2nd, scrn.inp => 3rd
!         pot.bin => 4,5th, grids => 6th, atomic and potentials => 6,7th
! output: wscrn
subroutine screen( nat, nph, iphat, rat,                          &
     &                   toler1, toler2, lmaxph, rdirec, rfms2, lfms2,  &
     &                   maxl, irrh, iend, lfxc, nrptx0,                &
     &                   ihole, xmu, adgc, adpc, iz,                    &
     &                   xion, iunf, xnval,edens, dgc0, dpc0,           &
     &                   ri, dx, x0, rmt, rnrm, em, ne,                 &
     &                   ixc0, iph, vtot, vvalgs, eref, dgcn, dpcn, ph, &
     &                   ilast, vch, wscrn)
  !     linear response calculation
  !     calculate static screening

  use DimsMod, only: lx, nrptx, natx, nphx=>nphu, nspx=>nspu, ltot, nex
  use constants
  implicit none
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
  external dgetrf, dgetrs
  !     function getiat
  integer getiat
  !     pot.bin --------------------------------------------------------
  integer ihole, iz, iunf
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
  double precision rmt, rnrm, x0, dx, ri(nrptx)
  integer kinit, linit, jri, jri1, jnrm, ilast, iph, ic3,irr, ikap
  integer i, j, k, m, n, ir1, ir2, ll, lfin, ilp, ie, ncycle, ios
  !     work space for xcpot
  complex*16 v(nrptx), vval(nrptx)
  !     work space for fovrg
  complex*16 pr(nrptx), qr(nrptx), pn(nrptx), qn(nrptx)
  complex*16 p2, ck, xkmt, xkmtp, xck
  complex*16 pu, qu, dum1, factor, de
  complex*16 xfnorm, xirf, temp, ph0
  complex*16 bessj(8), bessh(8)
  complex*16 jl, jlp1, nl, nlp1

  !     fms routine (variables are in single precision!)
  integer nat
  real rfms2, rdirec, toler1, toler2, rpart, aipart, rat(3,natx)
  integer il, ix, im, nsp, minv, iverb, nph, ipp, ill, lfms2

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
  integer ipiv(nrx)
  !     work space
  !     core hole potential
  double precision vch(nrx)
  !     Coulomb Exchange      
  double precision Kmat(nrx, nrx)
  !     LDA fxc      
  double precision fxc(nrx)
  !     chi_0(r,r')
  complex*16 chi0r(nrx, nrx)
  !     chi_0(r,r',e)
  complex*16 chi0re(nrx, nrx)
  !     (1 - K Chi0)
  double precision Amat(nrx, nrx)
  !     W, the result of calculation
  double precision wscrn(nrx)

  complex*16, allocatable :: gtrl(:,:)
  complex, allocatable :: gg(:,:,:), gtr(:,:,:), xphase(:,:,:)
  logical, allocatable :: lcalc(:)
  ! String for log messages:
  character*512 slog
  ! Allocate local variables
  allocate(gtrl(0:lx, nex))
  allocate(gg(nspx*(lx+1)**2,nspx*(lx+1)**2,0:nphx))
  allocate(gtr(0:lx, 0:nphx, nex))
  allocate(xphase(nspx, -lx:lx, 0:nphx))
  allocate(lcalc(0:lx))


  gtrl(:,:) = 0.d0 ! JK Need to initialize for ifort. 04/2010
  !     =================================================================      
  call setkap(ihole, kinit, linit)

  !     set imt and jri (use general Loucks grid)
  !     rmt is between imt and jri (see function ii(r) in file xx.f)
  rmt   = rmt
  rnrm  = rnrm
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
        chi0re(i,j) = 0.0d0
     end do
  end do
  percent  = 0


  call yprep(iph, nat, inclus, nph, iphat, rfms2, rat)

  if (inclus .gt. 1) then
     write(slog,'(a,i4,a,i3)') 'FMS for a cluster of ', inclus, ' atoms around iph = ', iph
     call wlog(slog)
     !print *, 'FMS for a cluster of ', inclus, ' atoms around iph = ', iph
  end if

  !===============================================================
  !     start energy cycle (big loop!)
  !===============================================================
  do i = 1, nrptx
     v(i)    = vtot(i)
     vval(i) = vvalgs(i)
  end do

  do ie = 1, ne
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

        do i = 0, lx
           do j = 0, nphx
              gtr(i,j,ie) = 0.0d0
           end do
        end do
        
        call fms(lfms2, nsp, ispin, inclus, nph, cks,lmaxph, &
             & xphase,ie,iverb,minv,rdirec,toler1,toler2,lcalc, gg)
        
        do il = 0, lmaxph(iph)
           ix = il**2
           do im = 1, 2*il+1
              gtr(il,iph,ie) = gtr(il,iph,ie) + gg(ix+im,ix+im,iph)
           end do
           
           gtrl(il,ie) = ( dble(real( gtr(il,iph,ie))) &
                & + coni*dble(aimag(gtr(il,iph,ie))) )&
                & * exp(2*coni*ph(ie,il+1,iph))/(2*il+1.0d0)
        end do
     end if

     !       the maximum angular momentum        
     do k = 1, maxl

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
             &                  pu, qu, pr, qr, iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph

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
           pr(i) = pr(i)*xfnorm
           qr(i) = qr(i)*xfnorm
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
             &                  pu, qu, pn, qn, iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph

        !         set N- irregular solution , which is outside
        !         N=(nlp1*cos(ph0)+jlp1*sin(ph0))*factor *rmt * dum1
        !         N = i*R - H*exp(i*ph0)
        temp = exp(coni*ph0)
        !         calculate wronskian (qu)
        qu = 2*alpinv *temp*(pn(jri)*qr(jri)-pr(jri)*qn(jri))
        IF(qu.ne.0.d0) THEN
           qu = 1 /qu / ck
        ELSE
           qu = 0.d0
        END IF
        !         qu should be close to 1
        do i = 1, ilast
           pn(i) = temp * pn(i)*qu
           qn(i) = temp * qn(i)*qu
        end do

        !         Use exact solution to continue solutions beyond rmt
        do j = jri, ilast
           xck  = ck * ri(j)
           call exjlnl(xck, ll, jl, nl)
           call exjlnl(xck, ll+1, jlp1, nlp1)
           call besjh(xck, maxl+1, bessj, bessh)

           pr(j) =  (jl*cos(ph0) -   nl*sin(ph0))*ri(j)*dum1
           qr(j) =(jlp1*cos(ph0) - nlp1*sin(ph0))*ri(j)*dum1*factor
           pn(j) = bessh(ll+1)*exp(coni*ph0) *ri(j)*dum1
           qn(j) = bessh(ll+2)*exp(coni*ph0) *ri(j)*dum1*factor
        end do
        !=============================================================
        !         setup prefactors
        !=============================================================
        ll = k - 1

        !         the factor includes averaged over spins
        factor = -1.0d0/(2*(pi**2))*(2*ll+1.0d0)*(2.0d0*ck)**2*dx**2
!         !=============================================================
!         !         FMS
!         !=============================================================
!         if (inclus .gt. 1) then
!            do ipp = 0, nph
!               do ill = -lmaxph(ipp), lmaxph(ipp)
!                  rpart  = dble( ph(ie, abs(ill)+1, ipp))
!                  aipart = dimag(ph(ie, abs(ill)+1, ipp))
!                  xphase(1, ill, ipp) = cmplx(rpart, aipart)
!               end do
!            end do

!            iverb = 0
!            nsp = 1
!            ispin = 0
!            do i = 0, lx
!               lcalc(i) = .true.
!            end do
!            minv = 0

!            do i = 0, lx
!               do j = 0, nphx
!                  gtr(i,j,ie) = 0.0d0
!               end do
!            end do

!            call fms(lfms2, nsp, ispin, inclus, nph, cks,lmaxph,        &
!                 &              xphase,ie,iverb,minv,rdirec,toler1,toler2,lcalc, gg)

!            do il = 0, lmaxph(iph)
!               ix = il**2
!               do im = 1, 2*il+1
!                  gtr(il,iph,ie) = gtr(il,iph,ie) + gg(ix+im,ix+im,iph)
!               end do

!               gtrl(il,ie) = (        dble(real( gtr(il,iph,ie)))      &
!                    &                          + coni*dble(aimag(gtr(il,iph,ie)))     )&
!                    &                        * exp(2*coni*ph(ie,il+1,iph))/(2*il+1.0d0)
!            end do
!         end if
        !=============================================================
        !         construct response function
        !=============================================================         
        !         Construct G_l(r,r')G_l(r',r) for W0 calculation, G(r,r') = R(r<)*H(r>)
        do m = 1, ilast
           do n = m, ilast
              chi0re(m,n) = chi0re(m,n) + factor * ri(m)*ri(n)          &
                   &                             * pr(m)*pr(m) * pn(n)*pn(n)
           end do

        end do

        if (inclus .gt. 1) then
!		write(*,*) 'size ri, size pr, size pn',size(ri),size(pr),size(pn)
!		write(*,*) 'size chi0re',size(chi0re,1),size(chi0re,2)
!		write(*,*) 'size gtrl',size(gtrl,1),size(gtrl,2)
!		write(*,*) 'll,ie,jnrm',ll,ie,jnrm
!		write(*,*) 'lx,nex,nrx',lx,nex,nrx
           do m = 1, jnrm
              do n = m, jnrm
                 chi0re(m,n) = chi0re(m,n)+factor * ri(m)*ri(n) * (      &
                      &              2 * gtrl(ll,ie)* pr(m)*pr(m) * pr(n)*pn(n)          &
                      &            + gtrl(ll,ie)**2 * pr(m)*pr(m) * pr(n)*pr(n)   )
              end do

           end do
        end if

     end do

     !===============================================================
     !       integration over energy
     !===============================================================
     if (ie .eq. 1) then
        de = (em(ie+1) - em( ie ))/2.0d0
     else if (ie .eq. ne) then
        de = (em( ie ) - em(ie-1))/2.0d0
     else
        de = (em(ie+1) - em(ie-1))/2.0d0
     end if

     do ir1 = 1, ilast
        do i = ir1, ilast
           chi0r( ir1,i) = chi0r(ir1,i) + chi0re(ir1,i)*de
           chi0re(ir1,i) = (0.0d0, 0.0d0)
        end do
     end do

     !       print the current progress to the wscrn
     if ( 100.0*ie/ne .gt. percent) then
        if(percent.ge.1.d0) then
           write(slog,'(2x,i3,a)') nint(percent), '%'
        else
           write(slog,'(2x,i3,a)') nint(percent), '% of energy integration'
        endif
        call wlog(slog)
        percent = percent + 20.0
     end if
  end do

  call wlog('  100%')
  !=================================================================
  !     end of energy cycle
  !=================================================================

  !=================================================================
  !     setupt K-matrix used for calculating screening effect
  !=================================================================
  do m = 1, ilast
     do n = m, ilast
        Kmat(m,n) =  4.0d0*pi * 1 * 1/ri(n)
     end do

  end do
  !=================================================================
  !     LDA Exchange Correlation (actually this probably doesn't work)
  !=================================================================  
  if (lfxc .gt. 0) then
     call ldafxc(ilast, ri, edens, lfxc, fxc)
     call wlog('Using TDLDA kernel.')
     do i = 1, ilast
        Kmat(i,i) = Kmat(i,i) + 4.0d0*pi * fxc(i)
     end do
  end if
  !=================================================================
  !     get response function
  !=================================================================
  call wlog('Preparing response function.')

  do ir1 = 2, ilast
     do i = 1, ir1-1
        chi0r(ir1,i) = chi0r(i,ir1)
        Kmat(ir1,i) = Kmat(i,ir1)
     end do
  end do

  !     here performs v_ch(r) = \int (dgc0'^2 + dpc0'^2)*1/r> dr'
  do i = 1, ilast
     vch(i)  = (dpc0(i)**2 + dgc0(i)**2) * dx * ri(i)
  end do
  !     print *, 'Check normalization', sum(vch)

  do i = 1, ilast
     wscrn(i) = 0.0d0
     do j = 1, i
        wscrn(i) = wscrn(i) + vch(j)
     end do
     wscrn(i) = wscrn(i)/ri(i)
     do j = i + 1, ilast
        wscrn(i) = wscrn(i) + vch(j)/ri(j)
     end do
  end do

  do i = 1, ilast
     vch(i) = wscrn(i)
  end do

  !=================================================================
  !     start matrix inversion
  !=================================================================
  !     the matrix to be inverted ( .I. N = N x N Identity Matrix )
  !call wlog('Start matrix inversion')
  !     Amat = Amat - matmul(Kmat, Chi0)
  do i = 1, ilast
     do j = 1, ilast
        if (i .eq. j) then
           Amat(i,j) = 1.0d0
        else
           Amat(i,j) = 0.0d0
        end if
        do k = 1, nrx
           Amat(i,j) = Amat(i,j)-Kmat(i,k)*dimag(chi0r(k,j))
        end do
     end do
  end do

  call wlog('Computing (1 - K Chi0)^-1 v_ch')

  nrhs = 1
  !     size(Amat,1)
  lda  = nrx
  !     size(wscrn)
  ldb  = nrx
  !     actual size of nxn matrix
  n    = ilast

  call dgetrf(n, n, Amat, lda, ipiv, info)

  if (info .eq. 0) then
     call dgetrs(trans, n, nrhs, Amat, lda, ipiv, wscrn, ldb, info)

     !if (info .ne. 0) print *, info
  else
     call wlog('The factor U is singular')
  end if
  !=================================================================
  !     end of matrix inversion
  !=================================================================

  ! Deallocate local variables
  deallocate(gtrl,gg,gtr,xphase,lcalc)

  return
  !      end subroutine screen
end subroutine screen
!=======================================================================
!     END wscrn
!=======================================================================
