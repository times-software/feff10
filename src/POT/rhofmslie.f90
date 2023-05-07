
subroutine rhofmslie(edens, edenvl, vtot, vvalgs, rmt, rnrm, rhoint, vint, &
                     xnval, x0, ri, dx, adgc, adpc, dgc, dpc, ie, xrhoce,  &
                     yrhoce, xrhole, yrhole, ph, gtr, ri05, nr05, ee,iscmt)
  !
  ! Calculate density (including scattering term) for all L values at a given energy point
  !
  ! This was extracted from scmt.f90 for use in m_thermal_scf.f90
  !
  ! Inputs:
  !   edens - electron density
  !   edenvl - valence electron density
  !   vtot - total potential
  !   vvalg -
  !   rmt - muffin tin radii
  !   rnrm - norman radii
  !   rhoint - interstitial density
  !   vint - interstitial potential value
  !   xnval - number of valence electrons
  !   x0 - parameter in exponential spatial grid
  !   ri - exponential spatial grid  ( ri(j) = exp(-x0 + (j-1)*dx) )
  !   dx - delta parameter in exponential spatial grid
  !   ie - energy index (some subroutines initialize if ie.eq.1)
  !   ee - energy
  !   ri05 - radial grid
  !
  ! Outputs:
  !   xrhoce - LDOS
  !   yrhoce - spatial density
  !   xrhole - workspace for scattering portion of LDOS
  !   yrhole - workspace for scattering portion of spatial density
  !   ph - phase shifts
  !   gtr - scattering matrix
  !   nr05 - indices of norman radii

  use constants
  use DimsMod, only: nrptx, lx, nphx=>nphu

  use atoms_inp, only: nat, nph, iphat, rat
  use potential_inp, only: ixc, xion, iunf, iz, ihole, lmaxsc, nohole, rfms1, lfms1, jumprm0=>jumprm, rgrd

  implicit None
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016

  ! Input
  integer, intent(in):: ie, iscmt
  double precision, intent(in):: vtot (251,0:nphx), vvalgs (251,0:nphx)
  double precision, intent(in):: xnval (41,0:nphx)
  double precision:: ri(nrptx)
  double precision, intent(in):: vint, rhoint, dx, x0
  real*8, intent(in) :: rmt(0:nphx),rnrm(0:nphx)
  double precision, intent(in):: ri05(251)

!     input and output
  double precision, intent(in):: edens(251,0:nphx), edenvl(251,0:nphx)
  integer, intent(inout):: nr05(0:nphx)

  ! Internal variables
  double precision:: dmagx(nrptx), dmag0(251), vjump
  double precision:: dum(nrptx), vtotph(nrptx), vvalph(nrptx)
  double precision:: dgc(251,41,0:nphx+1), dpc(251,41,0:nphx+1)
  double precision:: adgc(10,41,0:nphx+1), adpc(10,41,0:nphx+1)
  double precision:: dgcn(nrptx,41), dpcn(nrptx,41)
  complex*16 yrhoce(251,0:nphx)

  integer iph, jri, jri1, i, itmp, iph0, il, ir, jumprm
  complex*16 em, eref, ee
  !character*512 slog

  complex, intent(out):: gtr(0:lx, 0:nphx)
  complex*16, intent(out):: xrhoce(0:lx, 0:nphx), xrhole(0:lx,0:nphx), ph(lx+1, 0:nphx), yrhole(251,0:lx,0:nphx)

  jumprm = jumprm0 ! Avoid changing global jumprm (even though it is changed back)
  do iph = 0, nph
      dmag0(:) = 0.d0
!c       use spin-unpolarized case to get SCF. set dmagx to zero
!c       may want to replace dmag0 with dmag(1,iph) for spin-dependent extension of SCF procedure.
      call fixvar (rmt(iph),edens(1,iph),vtot(1,iph),dmag0, vint, rhoint, dx, rgrd, jumprm, vjump, ri, vtotph, dum, dmagx) !makes vtotph,dum,dmagx out of vtot,edens,dmag0
      ! if (mod(ixc,10) .ge.5) then
      IF ((ixc.EQ.5).OR.(ixc.EQ.9).OR.(ixc.EQ.15)) THEN ! Replace mod(ixc,10).ge.5
        if (jumprm .gt. 0) jumprm = 2
        call fixvar (rmt(iph), edenvl(1,iph), vvalgs(1,iph), dmag0, vint, rhoint, dx, rgrd , jumprm, vjump, ri, vvalph, dum, dmagx)
        if (jumprm .gt. 0) jumprm = 1
      endif

      call fixdsx (iph, dx, rgrd , dgc, dpc, dgcn, dpcn)  !makes dgcn, dpcn out of dgc, dpc
      jri = (log(rmt(iph)) + x0) / rgrd + 2
      jri1 = jri+1
      eref = vtotph(jri1)
      do i = 1, jri1
        vtotph(i) = vtotph(i) - eref
      end do
      ! if (ixc.ge.5) then
      IF ((ixc.EQ.5).OR.(ixc.EQ.9).OR.(ixc.EQ.10).OR.(ixc.EQ.13).OR.(ixc.EQ.15)) THEN ! Replace ixc.ge.5
        do i = 1, jri1
          vvalph(i) = vvalph(i) - eref
        end do
      else
        do i = 1, jri1
          vvalph(i) = vtotph(i)
        end do
      endif

      itmp = 0
      if (iph.eq.0 .and. nohole.lt.0) itmp = ihole
      call rholie( ri05, nr05(iph), rgrd, x0, ri, ee, ixc, rmt(iph), rnrm(iph),    &     ! vtotph,vvalph,(a)dg/pcn,xnval? => xrhole,xrhoce,yrhole,yrhoce,ph
            vtotph, vvalph, xnval(1,iph), dgcn, dpcn, eref, adgc(1,1,iph), adpc(1,1,iph), xrhole(0,iph),           &
            xrhoce(0,iph),yrhole(1,0,iph),yrhoce(1,iph),ph(1,iph), iz(iph), xion(iph), iunf, itmp,lmaxsc(iph), iph) !KJ iph
      !solves Dirac equation to obtain orbitals and phase shifts for FMS below
  enddo ! iph

!     transform neg,emg to em,ne,eref first
  em= dble(ee)
  eref=dble(eref)-coni*dimag(ee)


!c    call fms for a cluster around central atom
  gtr(:,0:nph)=cmplx(0)
  if (rfms1 .gt. 0) then
    if (lfms1 .ne. 0) then
      call fmsie( 0, nph, lmaxsc, ie,  em, eref, ph, iz, rfms1, lfms1, nat, iphat, rat, gtr,iscmt.eq.1)
    else
      do iph0 = 0, nph
          call fmsie( iph0, nph, lmaxsc, ie, em, eref, ph, iz, rfms1, lfms1, nat, iphat, rat, gtr,iscmt.eq.1)
      enddo
    endif
  endif

  ! form total densities (including scattering contribution) from *rhoce, *rhole and gtr
  ! the result is stored back into *rhoce
  xrhoce(:,0:nph) = xrhoce(:,0:nph) + gtr(:,0:nph)*xrhole(:,0:nph)

  do il = 0,lx
    if (il.le.2 .or. iunf.ne.0) then
      do iph = 0,nph
        do ir = 1,nr05(iph)
          yrhoce(ir, iph) = yrhoce(ir, iph) + gtr(il,iph)*yrhole(ir,il,iph)
        end do
      end do
    end if
  end do

end subroutine
