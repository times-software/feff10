!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: screensub.f90,v $:
! $Revision: 1.6 $
! $Author: jorissen $
! $Date: 2011/02/11 04:14:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================================================================
!     wscrn
!=======================================================================
! inputs: geom.dat => 1st line, ldos.inp => 2nd, scrn.inp => 3rd
!         pot.bin => 4,5th, grids => 6th, atomic and potentials => 6,7th
! output: wscrn
subroutine chi_crpa( nat, nph, iphat, rat,                          &
     &                   toler1, toler2, lmaxph, rdirec, rfms2, lfms2,  &
     &                   maxl, irrh, iend, lfxc, nrptx0,                &
     &                   ihole, xmu, adgc, adpc, iz,                    &
     &                   xion, iunf, xnval,edens, dgc0, dpc0,           &
     &                   ri, dx, x0, rmt, rnrm, em, ne,                 &
     &                   ixc0, iph, vtot, vvalgs, eref, dgcn, dpcn, ph, &
     &                   ilast, vch, wscrn, rcutin, l_crpa)
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
  integer ihole, iz, iunf, l_crpa
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
  parameter ( nrx = 1251 )
  double precision percent
  !     energy mesh
  complex*16 em(nex)
  integer ne
  !     matrix inversion
  character*1 trans
  parameter ( trans = 'N' )
  integer info, nrhs, lda, ldb
  integer ipiv(nrx), ll_CRPA
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
  double precision wscrn(nrx), den_CRPA(nrx,nex), totden_CRPA(nrx), emin, emax, atnorm, xnel

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
  den_CRPA = 0.d0
  totden_CRPA = 0.d0
  gtrl(:,:,:) = 0.d0 ! JK Need to initialize for ifort. 04/2010
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
     ll_CRPA = l_crpa ! choose state
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
     IF(DBLE(em(ie)).GT.xmu) GOTO 105
        !=============================================================
        !         setup prefactors
        !=============================================================
     do k = 1, maxl
        ll = k - 1
       

        !         the factor includes averaged over spins
        !factor = -1.0d0/(2.d0*(pi**2))*(2*ll+1.0d0)*(2.0d0*ck)**2*dx**2
        factor = -1.0d0/(2.d0*(pi**2))*(2.0d0*ck)**2*dx**2
        pratom = 0.d0
        gcProj = 0.d0 
        IF(ll_CRPA.EQ.ll) THEN
           ! find projection onto normalized atomic state.
           ! First normalize inside the norman radius.
           atnorm = 0.d0
           DO m = 1, jnrm
              atnorm = atnorm + (dgc0(m)**2+dpc0(m)**2)*ri(m)*dx
           END DO
           dgc0(:) = dgc0(:)/SQRT(atnorm)
           dpc0(:) = dpc0(:)/SQRT(atnorm)
           dgc0(jnrm+1:) = 0.d0
           dpc0(jnrm+1:) = 0.d0
           DO m = 1, jnrm
              pnatom = 0.d0
              IF(m.LE.jnrm) THEN
                 pratom = pratom + (dgc0(m)*pr(m,ll)+dpc0(m)*qr(m,ll))*ri(m)*dx
              END IF
              DO n = 1, m
                 pnatom = pnatom + (dgc0(n)*pr(n,ll)+dpc0(n)*qr(n,ll))*ri(n)*dx
              END DO
              pnatom = pnatom*(dgc0(m)*pn(m,ll) + dpc0(m)*qn(m,ll))
              DO n = m+1, jnrm
                 pnatom = pnatom + (pr(m,ll)*dgc0(m)+qr(m,ll)*dpc0(m))*(dgc0(n)*pn(n,ll)+dpc0(n)*qn(n,ll))*ri(n)*dx
              END DO
              ! pnatom now holds \int G(r,r') phi(r) r**2 dr
              gcProj = gcProj + pnatom*ri(m)*dx
           END DO
        END IF
        !pratom = 0.d0
        !gcproj = 0.d0

        do m = 1, ilast
           do n = m, ilast
              !do ll2 = 0, maxl-1
              do ll2 = ll, ll
                 IF((ll_CRPA.EQ.ll2).AND.(ll_CRPA.EQ.ll)) THEN
                    IF(UseProjection) THEN
                       ! subtract projection onto atomic state dgc0.
                       rnew = MIN(MAX(rcut0,ri(m)),rcut)
                       rnew = rnew - rcut0
                       rpnew = MIN(MAX(rcut0,ri(n)),rcut)
                       rpnew = rpnew - rcut0
                       chi0re(m,n) = chi0re(m,n) + factor*sqrt((2.d0*ll+1)*(2.d0*ll2+1))*ri(m)*ri(n)* &
                            &  gtrl(ll,ll2,ie)*gtrl(ll2,ll,ie) *& 
                            &      pr(m,ll)*pr(m,ll) * pr(n,ll2)*pr(n,ll2)  &
                            &     * SIN(rnew/(rcut-rcut0)*pi/2.d0)**4 &
                            &     * SIN(rpnew/(rcut-rcut0)*pi/2.d0)**4
! Modified by FDV
! Added a comment at beggining of line or Solaris Studio gives error
!                           &     !( pr(m,ll)*pr(m,ll) * pr(n,ll2)*pr(n,ll2) - pratom**4*dgc0(m)**2*dgc0(n)**2 )   
                    END IF
                 ELSE
                    chi0re(m,n) = chi0re(m,n) + factor*sqrt((2.d0*ll+1)*(2.d0*ll2+1))*ri(m)*ri(n)* &
                       &  gtrl(ll,ll2,ie)*gtrl(ll2,ll,ie) * pr(m,ll)*pr(m,ll) * pr(n,ll2)*pr(n,ll2)   
                 END IF
              end do
           end do
        end do
        !IF(ll_CRPA.EQ.ll) CYCLE
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
              IF(ll_CRPA.EQ.ll) THEN
                 IF(UseProjection) THEN
                    rnew = MIN(MAX(rcut0,ri(m)),rcut)
                    rnew = rnew - rcut0
                    rpnew = MIN(MAX(rcut0,ri(n)),rcut)
                    rpnew = rpnew - rcut0
                    chi0re(m,n) = chi0re(m,n) + factor*(2.d0*ll+1.d0)* ri(m)*ri(n)          &
                         &   * ( pr(m,ll)*pr(m,ll) * pn(n,ll)*pn(n,ll)) &
                         &   * SIN(rpnew/(rcut-rcut0)*pi/2.d0)**4  &
                         &   * SIN(rnew/(rcut-rcut0)*pi/2.d0)**4 
                    
                         !&             * ( pr(m,ll)*pr(m,ll) * pn(n,ll)*pn(n,ll) - & 
                         !&                      gcProj**2*dgc0(m)**2*dgc0(n)**2 )
                 END IF
              ELSE
                 chi0re(m,n) = chi0re(m,n) + factor*(2.d0*ll+1.d0)* ri(m)*ri(n)          &
                   &                             *pr(m,ll)*pr(m,ll) * pn(n,ll)*pn(n,ll)
              END IF
           end do
        end do

        if (inclus .gt. 1) then
!		write(*,*) 'size ri, size pr, size pn',size(ri),size(pr),size(pn)
!		write(*,*) 'size chi0re',size(chi0re,1),size(chi0re,2)
!		write(*,*) 'size gtrl',size(gtrl,1),size(gtrl,2)
!		write(*,*) 'll,ie,jnrm',ll,ie,jnrm
!		write(*,*) 'lx,nex,nrx',lx,nex,nrx
           do m = 1, ilast !jri !-20 !jnrm
              do n = m, ilast !jri !-20 !jnrm
                 IF(ll_CRPA.EQ.ll) THEN
                    IF(UseProjection) THEN
                       rnew = MIN(MAX(rcut0,ri(m)),rcut)
                       rnew = rnew - rcut0
                       rpnew = MIN(MAX(rcut0,ri(n)),rcut)
                       rpnew = rpnew - rcut0
                       chi0re(m,n) = chi0re(m,n)+factor *(2.d0*ll+1) * ri(m)*ri(n) *       &
                            &           2 *gtrl(ll,ll,ie)* & 
                            & (pr(m,ll)*pr(m,ll)*pr(n,ll)*pn(n,ll)) &
                            &   * SIN(rnew/(rcut-rcut0)*pi/2.d0)**4  &
                            &   * SIN(rpnew/(rcut-rcut0)*pi/2.d0)**4 
                            !& (pr(m,ll)*pr(m,ll)*pr(n,ll)*pn(n,ll) - &
                            !& gcProj*pratom**2*dgc0(m)**2*dgc0(n)**2 )          
                    END IF
                 ELSE
                    chi0re(m,n) = chi0re(m,n)+factor *(2.d0*ll+1) * ri(m)*ri(n) * (      &
                         &              2 * gtrl(ll,ll,ie)* pr(m,ll)*pr(m,ll) * pr(n,ll)*pn(n,ll))          
                 END IF
              end do
           end do
        end if

     end do
     GOTO 103   
! Note by FDV
! According to Solaris Studio, the code below can not be reached.
        !========================================================
        !   Get central atom density for ll_CRPA
        !========================================================
        ikap = -ll_CRPA - 1

        ic3 = 1
        if (ll_CRPA .eq. 0) ic3 = 0

        !=========================================================
        !         get ''regular'' solution 
        !=========================================================
        irr = -1

        call dfovrg( ncycle, ikap, rmt, ilast, jri, p2, dx,           &
             &                  ri , v, vval, dgcn, dpcn,  adgc, adpc, xnval,   &
             &                  pu, qu, pr(:,ll_CRPA), qr(:,ll_CRPA), iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph

        call exjlnl(xkmt, ll_CRPA, jl, nl)
        call exjlnl(xkmt, ll_CRPA+1, jlp1, nlp1)
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
           pr(i,ll_CRPA) = pr(i,ll_CRPA)*xfnorm
           qr(i,ll_CRPA) = qr(i,ll_CRPA)*xfnorm
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
           pu = bessh(ll_CRPA+1)* exp(coni*ph0) *rmt * dum1
           qu = bessh(ll_CRPA+2)* exp(coni*ph0) *rmt * dum1* factor
        end if

        call dfovrg( ncycle, ikap, rmt, ilast, jri, p2, dx,           &
             &                  ri, v, vval, dgcn, dpcn, adgc, adpc, xnval,     &
             &                  pu, qu, pn(:,ll_CRPA), qn(:,ll_CRPA), iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph                                                                                    

        !         set N- irregular solution , which is outside
        !         N=(nlp1*cos(ph0)+jlp1*sin(ph0))*factor *rmt * dum1
        !         N = i*R - H*exp(i*ph0)                             
        temp = exp(coni*ph0)
        !         calculate wronskian (qu)                                                                                                                                       
        qu = 2*alpinv *temp*(pn(jri,ll_CRPA)*qr(jri,ll_CRPA)-pr(jri,ll_CRPA)*qn(jri,ll_CRPA))
        IF(qu.ne.0.d0) THEN
           qu = 1 /qu / ck
        ELSE
           qu = 0.d0
        END IF
        !         qu should be close to 1                                                                                                                                        
        do i = 1, ilast
           pn(i,ll_CRPA) = temp * pn(i,ll_CRPA)*qu
           qn(i,ll_CRPA) = temp * qn(i,ll_CRPA)*qu
        end do

        !         Use exact solution to continue solutions beyond rmt
        do j = jri, ilast
           xck  = ck * ri(j)
           call exjlnl(xck, ll_CRPA, jl, nl)
           call exjlnl(xck, ll_CRPA+1, jlp1, nlp1)
           call besjh(xck, maxl+1, bessj, bessh)

           pr(j,ll_CRPA) =  (jl*cos(ph0) -   nl*sin(ph0))*ri(j)*dum1
           qr(j,ll_CRPA) =(jlp1*cos(ph0) - nlp1*sin(ph0))*ri(j)*dum1*factor
           pn(j,ll_CRPA) = bessh(ll_CRPA+1)*exp(coni*ph0) *ri(j)*dum1
           qn(j,ll_CRPA) = bessh(ll_CRPA+2)*exp(coni*ph0) *ri(j)*dum1*factor
        end do
103 CONTINUE
     !===============================================================
     !       integration over energy
     !===============================================================
105 CONTINUE
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
        ! Get the density of states of angular momentum ll_CRPA
        den_CRPA(m,ie) = den_CRPA(m,ie) + DIMAG( (pr(m,ll_CRPA)*pn(m,ll_CRPA) &
       &   + pr(m,ll_CRPA)**2*gtrl(ll_CRPA,ll_CRPA,ie))*ck*4)*(2*ll_CRPA+1.0d0)/pi
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
           totden_CRPA(m) = totden_CRPA(m) + den_CRPA(m,ie)*de
           IF(m.le.jnrm) xnel = xnel + den_CRPA(m,ie)*de*ri(m)*dx
        end do
     END IF
  end do
  call wlog('100.% -- end of energy integration --')
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
     write(message,'(a)') 'Using TDLDA kernel.'
     call wlog(message)
     do i = 1, ilast
        Kmat(i,i) = Kmat(i,i) + 4.0d0*pi * fxc(i)
     end do
  end if
  !=================================================================
  !     get response function
  !=================================================================
  call wlog('Prepare response function')

  do ir1 = 2, ilast
     do i = 1, ir1-1
        chi0r(ir1,i) = chi0r(i,ir1)
        Kmat(ir1,i) = Kmat(i,ir1)
     end do
  end do

  !     here performs v_ch(r) = \int (dgc0'^2 + dpc0'^2)*1/r> dr'
  !den_CRPA = 0.d0
  !totden_CRPA = 0.d0
  !read(89,*) 
  ! Normalize the total d-density
  den_CRPA(1,1) = 0.d0
  DO i = 1, ilast
     IF(UseProjection) THEN
        rnew = MIN(MAX(rcut0,ri(i)),rcut)
        rnew = rnew - rcut0
        WRITE(77,*) ri(i), rnew
        totden_CRPA(i) = totden_CRPA(i)*COS(rnew/(rcut-rcut0)*pi/2.d0)**4
        den_CRPA(1,1) = den_CRPA(1,1) + totden_CRPA(i)*ri(i)*dx
     ELSE
        den_CRPA(1,1) = den_CRPA(1,1) + totden_CRPA(i)*ri(i)*dx
     END IF
  END DO
  totden_CRPA(:) = totden_CRPA(:)/den_CRPA(1,1)
  !totden_CRPA(jnrm+1:) = 0.d0
  do i = 1, ilast
  !   read(89,*,end=20) den_CRPA(i), totden_CRPA, den_CRPA(i) 
  !   PRINT*, totden_CRPA
  !   den_CRPA(i) = den_CRPA(i)/totden_CRPA
     vch(i)  = totden_CRPA(i) * dx * ri(i)
  end do
  vch(jnrm+1:) = 0.d0
20 continue
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
  !print *, 'Start matrix inversion'
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

  !print *, 'Compute (1 - K Chi0)^-1 v_ch'

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

     if (info < 0) call par_stop(' ')
  else
     call wlog('The factor U is singular in dgetrf')
  end if
  !=================================================================
  !     end of matrix inversion
  !=================================================================

  !=================================================================
  !     Integrate w_ch(r)*rho_c(r) to get U
  !=================================================================
  ! wscrn now holds the integral of W(r,r')*rho_c(r')
  U_Hub = 0.d0
  U_Bare = 0.d0
  n_occ = 0.d0
  vbare(:) = vch(:)
  do i = 1, ilast
     vch(i) = wscrn(i)*den_CRPA(i,ie)
     !write(88,*) ri(i), vch(i), den_CRPA(i,ie)
     U_Hub = U_Hub + wscrn(i)*totden_CRPA(i)*dx*ri(i)
     U_Bare = U_Bare + vbare(i)*totden_CRPA(i)*dx*ri(i)
     n_occ = n_occ + totden_CRPA(i)*ri(i)*dx
  end do
  OPEN(UNIT=44,FILE='crpa.dat',STATUS='REPLACE')
  WRITE(44,'(A)') 'U, n, U_Bare'
  WRITE(44,*) U_Hub, n_occ, U_Bare
  PRINT*, 'U, n, U_Bare'
  PRINT*, U_Hub*hart, n_occ, U_Bare*hart, xnel
  CLOSE(44)

  ! Deallocate local variables
  deallocate(gtrl,gg,gtr,xphase,lcalc)

  return
  !      end subroutine screen
end subroutine chi_crpa
!=======================================================================
!     END wscrn
!=======================================================================
