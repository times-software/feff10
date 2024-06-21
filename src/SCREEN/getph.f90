!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: getph.f90,v $:
! $Revision: 1.7 $
! $Author: jorissen $
! $Date: 2011/02/11 04:14:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getph ( ri, dx, x0, rmt, rnrm, ne, em, ixc,            &
     &                  vtot, vvalgs,                                   &
     &                  dgcn, dpcn, eref, adgc, adpc, iz, xion, iunf,   &
     &                  ihole, lmaxsc, xnval, ph, iph) !KJ iph

  use DimsMod, only: ltot, nrptx, lx, nex
  use constants
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
  !     ph

  !     function getiat
  integer getiat

  integer ne, ixc, iunf, ihole, lmaxsc, iz, iph !KJ iph
  double precision xion
  !     output
  complex*16 ph(nex, ltot+1)

  double precision ri(nrptx)
  double precision vtot(nrptx), vvalgs(nrptx)
  complex*16 vtotc(nrptx), vvalc(nrptx)
  double precision xnval(41)
  double precision dgcn(nrptx,41), dpcn(nrptx,41)
  double precision adgc(10,41), adpc(10,41)

  !     energy grid in complex e-plane
  complex*16 em(nex), eref

  !     work space for dfovrg: regular and irregular solutions
  complex*16 pn(nrptx), qn(nrptx)
  complex*16 p2, xkmt, ck
  complex*16 pu, qu
  complex*16 temp, phx

  complex*16 jl, jlp1, nl, nlp1
  double precision dx, rmt, rnrm, x0
  integer i, ie, ilast, irr, ikap, lll, ic3
  integer ncycle, jri, jnrm, lmax

  lmax = lmaxsc
  if (lmax .gt. lx) lmax = lx
  if (iz .le. 4)   lmax = 2
  if (iz .le. 2)   lmax = 1

  do i = 1, nrptx
     vtotc(i) = vtot(i)
     vvalc(i) = vvalgs(i)
  end do

  !     set imt and jri (use general Loucks grid)
  !     rmt is between imt and jri (see function ii(r) in file xx.f)
  jri  = getiat(x0, dx, rmt) + 1
  if (jri .gt. nrptx)  call par_stop('jri .gt. nrptx in phase')
  jnrm = getiat(x0, dx, rnrm) + 1
  !     ilast is the last integration point
  !     it is larger than jnrm for better interpolations
  ilast = jnrm + 6
  if (ilast .gt. nrptx) then
          PRINT*, 'WARNING: ilast > nrptx in getph!'
          ilast = nrptx
  end if

  do ie = 1, ne
     !       p2 is (complex momentum)**2 referenced to energy dep xc
     p2 = em(ie) - eref
     IF(DBLE(p2).GT.-40.d0/hart) THEN 
     if (mod(ixc,10) .lt. 5) then
        ncycle = 0
     else
        ncycle = 3
     endif
     ck = sqrt(2*p2+ (p2*alphfs)**2)
     xkmt = rmt * ck

     do lll = 0, lx
        if (lll .gt. lmax) then
           ph(ie, lll+1)  = 0
           goto 10
        endif

        !         may want to use ihole=0 for new screening. 
        !         don't want to use it now
        !         ihole = 0
        ikap = -1-lll
        irr  = -1
        ic3  = 1
        if (lll .eq. 0) ic3 = 0
        call dfovrg( ncycle, ikap, rmt, ilast, jri, p2, dx,           &
             &                  ri, vtotc, vvalc, dgcn, dpcn, adgc, adpc, xnval,&
             &                  pu, qu, pn, qn, iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph
        call exjlnl(xkmt, lll, jl, nl)
        call exjlnl(xkmt, lll+1, jlp1, nlp1)
        call phamp(rmt, pu, qu, ck, jl,nl,jlp1,nlp1, ikap, phx, temp)
        ph(ie, lll+1) = phx

     end do
10   continue
     ELSE
        ph(ie, :) = 0.d0
     END IF
  end do

  return
end subroutine getph
