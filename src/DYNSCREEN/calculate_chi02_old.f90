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
subroutine calculate_chi02( inclus,                           &
     &                   maxl, irrh, iend, nrptx0,            &
     &                   ihole, xmu, adgc, adpc, iz,          &
     &                   xion, iunf, xnval,                   &
     &                   ri, dx, x0, rmt, rnrm, em, ne, omega, ixc0, &
     &                   iph, vtot, vvalgs, eref, dgcn, dpcn, &
     &                   ilast0, gtrl1, gtrl2, chi0r)
  !     linear response calculation
  !     calculate static screening

  use DimsMod, only: lx, nrptx, natx, nphx=>nphu, nspx=>nspu, ltot, nex
  use constants
  implicit none
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016

  ! Input
  !     pot.bin --------------------------------------------------------
  integer ihole, iz, iunf
  COMPLEX*16 gtrl1(0:lx, nex), gtrl2(0:lx, nex)
  double precision xmu, adgc(10,41), adpc(10,41)
  double precision xion, xnval(41)
  !     -----------------------------------------------------------------
  double precision vtot(nrptx), vvalgs(nrptx)
  double precision dgcn(nrptx,41), dpcn(nrptx,41)
  complex*16 eref, omega
  complex*16 em(nex)
  integer ne, ixc0
  !

  !     -----------------------------------------------------------------
  !     radial grid
  double precision rmt, rnrm, x0, dx, ri(nrptx) 
  integer ilast, ilast0, iph

  ! FMS on or not?
  integer inclus
  !     scrn.inp -------------------------------------------------------
  integer maxl, irrh, iend, nrptx0
  
  ! Local
  integer jri, jnrm, jri0
  integer ic3, irr, ikap
  integer i, j, k,  ll, ie, ncycle
  !     work space for xcpot
  complex*16 v(nrptx), vval(nrptx)
  !     work space for fovrg
  complex*16 pr(nrptx,2), qr(nrptx,2), pn(nrptx,2), qn(nrptx,2)
  complex*16 p2, ck, ck1, ck2, xkmt, xck
  complex*16 pu, qu, dum1, factor, de
  complex*16 xfnorm, xirf, temp, ph0
  complex*16 bessj(8), bessh(8)
  complex*16 jl, jlp1, nl, nlp1


  !     phase shift      

  !     Linear Response
  !     -----------------------------------------------------------------
  integer nrx
  parameter ( nrx = 2502 )
  double precision percent
  !     chi_0(r,r')
  complex*16 chi0r(nrx, nrx)
  !     chi_0(r,r',e)
  complex*16 chi0re(nrx, nrx)

  ! String for log messages:
  character*512 slog
  integer ipart
  integer getiat
  external getiat

  LOGICAL,PARAMETER :: FreeElectronTest = .FALSE.
  !     function getiat
  if ( mod(ixc0,10) .lt. 5 ) then
     ncycle = 0
  else
     !         fix later . may be ncycle can be less
     ncycle = 3
  endif

  !     set imt and jri (use general Loucks grid)
  !     rmt is between imt and jri (see function ii(r) in file xx.f)
  !PRINT*, 'x0,dx,rmt=', x0,dx,rmt
  jri   = getiat(x0, dx, rmt) + 1
  if ( jri + 1 .gt. nrptx )  call par_stop('jri .gt. nrptx in phase')
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



  !===============================================================
  !     start energy cycle (big loop!)
  !===============================================================
  do i = 1, nrptx
     v(i)    = vtot(i)
     vval(i) = vvalgs(i)
  end do

  do ie = 1, ne
     !       the maximum kappa
     do k = 1, maxl
        do ipart = 1, 2
           if(ipart.EQ.1) then
              p2 = em(ie) - eref
           else
              p2 = em(ie) - eref + omega
           end if
           !       set the method to calculate atomic cross section
           !       p2 is (complex momentum)**2 referenced to energy dep xc
           
           ck   = sqrt(2*p2 + (p2*alphfs)**2)
           if(ipart.EQ.1) then
              ck1 = ck
           else
              ck2 = ck
           end if
           xkmt = rmt * ck

           ll = k - 1
           ikap = -k
           
           ic3 = 1
           if (ll .eq. 0) ic3 = 0
           
           !=========================================================
           !         get ''regular'' solution
           !=========================================================
           irr = -1
            
           !PRINT*, ncycle, ikap, rmt
           !PRINT*, ilast, jri, p2
           !PRINT*, dx, ri(1), v(1)
           !PRINT*, vval(1), dgcn(1), dpcn(1)
           !PRINT*, adgc(1), adpc(1), xnval(1)
           !PRINT*, pu, qu, pr(1)
           !PRINT*, qr(1), iz, ihole
           !PRINT*, xion, iunf, irr
           !PRINT*, ic3, iph 
           call dfovrg( ncycle, ikap, rmt, ilast, jri, p2, dx,           &
                &                  ri , v, vval, dgcn, dpcn,  adgc, adpc, xnval,   &
                &                  pu, qu, pr(:,ipart), qr(:,ipart), iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph
           !PRINT*, 'pr=', pr(2,ipart)
           !STOP
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
              pr(i,ipart) = pr(i,ipart)*xfnorm
              qr(i,ipart) = qr(i,ipart)*xfnorm
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
                &                  pu, qu, pn(:,ipart), qn(:,ipart), iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph
           
           !         set N- irregular solution , which is outside
           !         N=(nlp1*cos(ph0)+jlp1*sin(ph0))*factor *rmt * dum1
           !         N = i*R - H*exp(i*ph0)
           temp = exp(coni*ph0)
           !         calculate wronskian (qu)
           qu = 2*alpinv *temp*(pn(jri,ipart)*qr(jri,ipart)-pr(jri,ipart)*qn(jri,ipart))
           IF(qu.ne.0.d0) THEN
              qu = 1 /qu / ck
           ELSE
              qu = 0.d0
           END IF
           !         qu should be close to 1
           do i = 1, ilast
              pn(i,ipart) = temp * pn(i,ipart)*qu
              qn(i,ipart) = temp * qn(i,ipart)*qu
           end do
        
           !         Use exact solution to continue solutions beyond rmt
           IF(FreeElectronTest) THEN
              jri0 = 1
              ph0 = 0.d0
           ELSE
              jri0 = jri
           END IF
           IF(ilast0.LE.0) ilast0 = ilast 
           do j = jri0, ilast0
              xck  = ck * ri(j)
              call exjlnl(xck, ll, jl, nl)
              call exjlnl(xck, ll+1, jlp1, nlp1)
              call besjh(xck, maxl+1, bessj, bessh)
              
              pr(j,ipart) =  (jl*cos(ph0) -   nl*sin(ph0))*ri(j)*dum1
              qr(j,ipart) =(jlp1*cos(ph0) - nlp1*sin(ph0))*ri(j)*dum1*factor
              pn(j,ipart) = bessh(ll+1)*exp(coni*ph0) *ri(j)*dum1
              qn(j,ipart) = bessh(ll+2)*exp(coni*ph0) *ri(j)*dum1*factor
           end do
        end do ! end of loop over parts G(E) and G(E-w)
        !=============================================================
        !         setup prefactors
        !=============================================================
        ll = k - 1
        
        !         the factor includes averaged over spins
        factor = -1.0d0/(2*(pi**2))*(2*ll+1.0d0)*(2.0d0)*(2.d0*ck2)*dx**2 ! Need to multiply ck1 below to take imag part.
        !PRINT*, 'factor=', factor
        !PRINT*, 'll=', ll
        !PRINT*, 'ck1=', ck1, 'ck2=', ck2
        !PRINT*, dx
        !PRINT*, pi
        
        !=============================================================
        !         construct response function
        !=============================================================         
        !         Construct G_l(r,r')G_l(r',r) for W0 calculation, G(r,r') = R(r<)*H(r>)
        do i = 1, ilast0
           do j = i, ilast0
              chi0re(i,j) = chi0re(i,j) + factor * ri(i)*ri(j)          &
                   &                             * DIMAG(ck1*pr(i,1)*pn(j,1))*pr(i,2)*pn(j,2)
           end do
        end do
        if (inclus .gt. 1 .AND. (.NOT.FreeElectronTest)) then
           !		write(*,*) 'size ri, size pr, size pn',size(ri),size(pr),size(pn)
           !		write(*,*) 'size chi0re',size(chi0re,1),size(chi0re,2)
           !		write(*,*) 'size gtrl',size(gtrl,1),size(gtrl,2)
           !		write(*,*) 'll,ie,jnrm',ll,ie,jnrm
           !		write(*,*) 'lx,nex,nrx',lx,nex,nrx
           do i = 1, jnrm
              do j = i, jnrm
                 chi0re(i,j) = chi0re(i,j)+factor * ri(i)*ri(j) * (      &
                      &            + DIMAG(ck1*gtrl1(ll,ie)*pr(i,1)*pr(j,1))*pr(i,2)*pn(j,2)          &
                      &            + gtrl2(ll,ie)* pr(i,2)*pr(j,2)*DIMAG(ck1*pr(i,1)*pn(j,1))          &
                      &            + DIMAG(ck1*gtrl1(ll,ie)*pr(i,1)*pr(j,1))*gtrl2(ll,ie)*pr(i,2)*pr(j,2)   )
              end do
              
           end do
        !WRITE(80,*) chi0re(1,2)
        !WRITE(80,*) 'em=', em(ie)
        !WRITE(80,*) 'factor=', factor
        !WRITE(80,*) 'ris=', ri(1), ri(2)
        !WRITE(80,*) 'pr=', pr(1,1), pr(1,2)
        !WRITE(80,*) 'pn=', pn(2,1), pn(2,2)
        !WRITE(80,*) 'pr*pn=', pn(2,2)*pr(1,2),pr(1,1)*pn(2,1)
        !WRITE(80,*) 'pr*pn=', pn(2,2)*pr(1,2)*pr(1,1)*pn(2,1)*ri(1)*ri(2)*factor
        !WRITE(80,*) 'gtrl=', gtrl1(ll,ie), gtrl2(ll,ie)
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

     !WRITE(82,*) ie, ne, em(ie), chi0re(1,2), de
     do i = 1, ilast0
        do j = i, ilast0
           chi0r( i,j) = chi0r(i,j) + chi0re(i,j)*de
           chi0re(i,j) = (0.0d0, 0.0d0)
        end do
     end do
     !PRINT*, 'chi0=', chi0r(1,10)
     !STOP

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

  do i = 2, ilast0
     do j = 1, i-1
        chi0r(i,j) = chi0r(j,i)
     end do
  end do
  
  call wlog('  100%')

  return
  !      end subroutine screen
end subroutine calculate_chi02
