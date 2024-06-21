!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: prep.f90,v $:
! $Revision: 1.14 $
! $Author: hebhop $
! $Date: 2012/01/07 00:55:48 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================================================================
!     PREP
!=======================================================================
subroutine prep (vr0, ixc0, nrx, ri, dx, x0, ilast, vch, wscrn, CRPA)
!KJ removed unused ibounc
  use DimsMod, only: natx, nrptx, nphx=>nphu, nex, lx, ltot
  use constants
  use atoms_inp
  use ldos_inp
  use screen_inp
  use crpa_inp
  implicit none
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
  integer, intent(in) :: ixc0, nrx
  double precision, intent(in) :: vr0, dx, x0, ri(nrptx)
  integer, intent(out) :: ilast
  logical :: CRPA
  !     W, the result of calculation
  double precision,intent(out) :: vch(nrx), wscrn(nrx)

  !     =================================================================
  !     function getiat
  integer getiat
  !     -----------------------------------------------------------------
  real srat(3,natx)
  !     pot.bin --------------------------------------------------------
  integer ntitle
  character*80 title(nheadx)
  integer nohole, ihole, inters, iafolp, jumprm, iunf
  double precision rnrmav, xmu, vint, rhoint
  double precision emu, s02, erelax, wp, ecv, rs, xf, qtotel, totvol
  integer imt(0:nphx), inrm(0:nphx)
  double precision rmt(0:nphx), rnrm(0:nphx)
  double precision folp(0:nphx),folpx(0:nphx)
  double precision dgc0(251), dpc0(251)
  double precision dgc(251, 41, 0:nphx+1), dpc(251, 41, 0:nphx+1)
  double precision adgc(10, 41, 0:nphx+1), adpc(10, 41, 0:nphx+1)
  double precision edens(251, 0:nphx), vclap(251, 0:nphx)
  double precision vtot(251, 0:nphx), edenvl(251, 0:nphx)
  double precision vvalgs(251, 0:nphx), dmag(251, 0:nphx)
  double precision xnval(41,0:nphx)
  double precision eorb(41)
  integer kappa(41)
  double precision qnrm(0:nphx)

  real*8, allocatable :: xnmues(:,:)

  integer iorb(-5:4,0:nphx)
  integer iz(0:nphx)
  double precision xion(0:nphx)
  double precision xnatph(0:nphx)
  !     -----------------------------------------------------------------
  double precision vtotph(nrptx,0:nphx), vvalph(nrptx,0:nphx)
  double precision rhoph(nrptx), rhphvl(nrptx), dmagx(nrptx)
  double precision dgcx(nrptx), dpcx(nrptx)
  double precision dgcn(nrptx,41,0:nphx), dpcn(nrptx,41,0:nphx)
  integer iph, jnew
  double precision edge, vjump

  integer jri, jri1, ne, ne1, ne3, i, j, iholetmp
  integer inclus(0:nphx), lmax(0:nphx)
  complex*16 ph(nex, ltot+1, 0:nphx)
  double precision de, enext
  complex*16 eref(nex), em(nex)

  ! Allocate variables
    allocate(xnmues(0:lx,0:nphx))

  !     read pot.bin
  call rdpot ( ntitle, title, rnrmav, xmu, vint, rhoint,            &
       &                  emu, s02, erelax, wp, ecv,rs,xf, qtotel,        &
       &                  imt, rmt, inrm, rnrm, folp, folpx, xnatph,      &
       &                  dgc0, dpc0, dgc, dpc, adgc, adpc,               &
       &                  edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,&
       &                  eorb, kappa, iorb, qnrm, xnmues, nohole, ihole, &
       &                  inters, totvol, iafolp, xion, iunf, iz, jumprm)

!KJ here, it is necessary to be careful ...  If pot ran in k-space, rfms 
  rfms2  = ScreenI%rfms
  rdirec = ScreenI%rfms*2

  edge = xmu - vr0
  emu  = emu - vr0

  !     make energy grid
  ! JK - replacing the energy grid with an exponential grid from EFermi
  !      to EFermi + i*emax

  ! This sets an exponential grid on the imaginary axis.
  ! This grid seems to have trouble for Mg. Revert to earlier grid for now.
  !ne = ner
  !PRINT*, eimax, ne
  !CALL SetEGrid(eimax, ne, em)
  !PRINT*, em
  ! Now shift grid so that real part is at EFermi.
  !em(1:ne) = em(1:ne) + xmu
  call setegi(xmu+ScreenI%emin, xmu+ScreenI%emax, ScreenI%eimax, ScreenI%ermin, ScreenI%ner, ScreenI%nei, em, ne)
  !     compute phase shifts
  do iph = 0, nph
     IF(iph.EQ.0) THEN ! JK Fixed ihole bug 10/2010
        iholetmp = ihole
     ELSE
        iholetmp = 0
     END IF
     !print *, '     get phase shift for potential type ', iph
     lmax(iph) = lx

     call fixvar (rmt(iph),edens(1,iph),vtot(1,iph),dmag(1,iph),vint,&
          &           rhoint,dx,dx,jumprm,vjump,ri,vtotph(1,iph),rhoph,dmagx)
     call fixdsx (iph, dx, dx, dgc, dpc, dgcn(1,1,iph),dpcn(1,1,iph))

     jri     = getiat(x0, dx, rmt(iph)) + 1
     jri1    = jri + 1
     eref(1) = vtotph(jri1,iph)
     do i = 1, jri1
        vtotph(i,iph) = vtotph(i,iph) - eref(1)
     end do

     if (ixc .ge. 5) then
        do i = 1, jri1
           vvalph(i,iph) = vvalph(i,iph) - eref(1)
        end do
     else
        do i = 1, jri1
           vvalph(i,iph) = vtotph(i,iph)
        end do
     endif

     call getph( ri, dx, x0, rmt(iph), rnrm(iph), ne, em, ixc,       &
          &             vtotph(1,iph), vvalph(1,iph),                        &
          &             dgcn(1,1,iph), dpcn(1,1,iph), eref(1),               &
          &             adgc(1,1,iph), adpc(1,1,iph),                        &
          &             iz(iph), xion(iph), iunf,                            &
          &             iholetmp, lx, xnval(1,iph), ph(1,1,iph), iph) !KJ iph
  end do

  !=================================================================
  !     Main part of calculation
  !=================================================================
  iph = 0
  !     transform to single precision
  srat=real(rat)

  !     only the potential type 0 (central atom) is evaluated
  ! geom.dat => 1st, ldos.inp => 2nd, scrn.inp => 3rd, pot.bin => 4,5th
  
  !KJ 12-2011 bugfix.  It's not because only type=0 is evaluated that you should pass arrays here and expect them to come out as scalars on the other end !!!! Ugh.
  IF(CRPA) THEN     
     !wscrn(1) = ermin
     call chi_crpa( nat, nph, iphat, srat,                               &
       &            toler1, toler2, lmaxph, rdirec, rfms2, lfms2,         &
       &            ScreenI%maxl, ScreenI%irrh, ScreenI%iend, ScreenI%lfxc, ScreenI%nrptx0,                       &
       &            ihole, xmu, adgc, adpc, iz(0),                           &
       &            xion(0), iunf, xnval, edens, dgc0, dpc0,                 &
       &            ri, dx, x0, rmt(0), rnrm(0), em, ne,                        &
       &            ixc0, iph, vtotph, vvalph, eref(1), dgcn, dpcn, ph,      &
       &            ilast, vch, wscrn,CRPAI%rcut, CRPAI%l_crpa)
  ELSE

     call screen( nat, nph, iphat, srat,                               &
       &            toler1, toler2, lmaxph, rdirec, rfms2, lfms2,         &
       &            ScreenI%maxl, ScreenI%irrh, ScreenI%iend, ScreenI%lfxc, ScreenI%nrptx0,                       &
       &            ihole, xmu, adgc, adpc, iz(0),                           &
       &            xion(0), iunf, xnval, edens, dgc0, dpc0,                 &  !KJ added (0) for xion, iz, rmt, rnrm
       &            ri, dx, x0, rmt(0), rnrm(0), em, ne,                        &
       &            ixc0, iph, vtotph, vvalph, eref(1), dgcn, dpcn, ph,      &  !KJ added (1) for eref.  Frankly this one is suspicious ...
       &            ilast, vch, wscrn)
  END IF
  ! Deallocate variables
    deallocate(xnmues)


  return
end subroutine prep
!=======================================================================
!     END PREP
!=======================================================================
