module rhorrp_mod
  !
  ! This module calculates the density matrix rho(r,r') at arbitrary points r and r' within cluster
  !
  ! Required input:
  !   pot.bin   (output by pot)
  !   phase.bin (output by xsph)
  !   gg_slice.bin (output by fms when DENSITY or COMPTON card is present)
  !
  ! Public subroutines:
  !    rhorrp_init - load input files, precalculate wavefunctions
  !    rhorrp - calculate density matrix at r,r'
  !    rhoerrp - calculate Green's function at r,r'
  !    rhorrp_deinit - clean up
  !    atomic_density - calculate 1s atomic density at r,r'
  !
  ! The energy grid for integration over occupied states is specified in
  ! XSPH/phmesh2.f90
  !
  ! For debugging, wavefunctions are saved to wf.bin (this should be made optional)
  !
  ! Definitions of values imported from other modules:
  !   ixc - exchange-correlation type
  !   rgrd - step parameter in exponential grid
  !   scf_temperature - temperature in eV
  !   nph - number of unique potentials
  !   rat - coordinates of each atom
  !   iatph - map from atom index to potential index
  !   iphat - map from potential index to atom index
  !   nat - number of atoms
  !   rfms2 - radius of sphere for atoms to include in FMS calculation
  !   lmaxph - maximum angular momentum value for each potential type
  !
  ! FIXME: currently lmaxph must be the same for all potential types. The
  ! indexing of gg_slice needs to be fixed to handle the case where lmax
  ! differs between potentials

  use dimsmod, only: nrptx, nheadx, nphx=>nphu, ltot, nspx=>nspu, nex, lx, init_dimensions
  use constants
  use potential_inp, only: ixc, rgrd, scf_temperature
  use atoms_inp, only: nph, rat, iatph, iphat, nat
  use xsph_inp, only: rfms2
  use fms_inp, only: lmaxph
  use par
  implicit none
  !Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
  private

  !
  ! public interface
  !
  public :: rhorrp_init, rhorrp_deinit, nearest_atom, rhorrp,&
            rhoerrp, atomic_density
  logical, public :: restrict_r_voronoi = .false.

  !
  ! Internal variables
  !

  ! original radial grid parameter
  real*8, parameter :: dxpot=0.05d0, x0=8.8d0

  ! work array to hold radial grid
  real*8 :: ri(nrptx), ripot(nrptx)
  real*8 :: dx
  integer :: nr

  ! values loaded from pot.bin
  ! (some of these could be moved into rhorrp_init as dummy vars)
  integer :: ntitle
  character*80  :: title(nheadx)
  real*8  :: rnrmav
  real*8  :: xmu,vint,rhoint,emu,s02,erelax,wp,ecv,rs,xf,qtotel
  !integer, dimension(0:nphx) :: imt,inrm,iz
  !real*8,  dimension(0:nphx) :: rmt,folp,folpx,xnatph,qnrm,xion
  !real*8,  dimension(0:nphx), public :: rnrm
  integer, allocatable :: imt(:),inrm(:),iz(:)
  real*8, allocatable :: rmt(:),folp(:),folpx(:),xnatph(:),qnrm(:),xion(:),rnrm(:)
  public rnrm
  real*8,  dimension(251)    :: dgc0,dpc0
  !real*8,  dimension(251,41,0:nphx) :: dgc,dpc
  !real*8,  dimension(10,41,0:nphx) :: adgc,adpc
  !real*8,  dimension(251,0:nphx) :: edens,vclap,vtot,edenvl,vvalgs,dmag
  !real*8,  dimension(41,0:nphx) :: xnval
  real*8, allocatable :: dgc(:,:,:),dpc(:,:,:),adgc(:,:,:),adpc(:,:,:)
  real*8, allocatable :: edens(:,:),vclap(:,:),vtot(:,:),edenvl(:,:),vvalgs(:,:),dmag(:,:),xnval(:,:)
  real*8,  dimension(41) :: eorb
  integer, dimension(41) :: kappa
  !integer,  dimension(-4:3,0:nphx) :: iorb
  integer, allocatable :: iorb(:,:)
  real*8, allocatable :: xnmues(:,:)
  integer :: nohole,ihole,inters,iafolp,iunf,jumprm
  real*8  :: totvol

  ! values loaded from xsph.bin
  integer, public :: ne, ne1, ne3
  complex*16, allocatable :: ph(:,:,:,:) ! phase shifts from phase.bin
  complex*16, allocatable :: ph2(:,:,:) ! phase shifts from init_wavefunctions (same as how ldos does it)
  complex*16, public, allocatable :: em(:)  ! energy grid
  complex*16, allocatable :: eref(:,:) ! reference energy

  complex*16 :: eref0
  ! wavefunctions (calculated in init_wavefunctions)
  ! p and q are large and small Dirac components
  ! r and n are regular and irregular solutions
  ! ("el" refers to the energy and angular momentum indices)
  complex*16, allocatable, dimension(:,:,:,:), target :: prel, pnel, qrel, qnel

  ! fms matrix (fms.bin)
  integer ldim
  complex, allocatable, dimension(:,:,:) :: gg_slice
  complex, allocatable, dimension(:,:,:,:) :: gg_diag

  ! number of atoms within rfms2 of each representative atom
  integer, allocatable :: inclus(:)
contains

subroutine rhorrp_init
  ! reads in required information for calculating density matrix
  use potential_inp, only: potential_read
  use atoms_inp, only: atoms_read
  use xsph_inp, only: xsph_read
  use fms_inp, only: fms_read
  use ifuns

  integer i

  ! initialize and read in input files
  call init_dimensions
  call potential_read
  call atoms_read
  call xsph_read
  call fms_read

  ! transform to code units (FIXME: why isn't this handled in *_read/*_write?)
  rfms2 = rfms2 / bohr
  rat=rat / bohr

  dx = rgrd
  nr = 251
  ldim = nspx*(lx+1)**2

  allocate(inclus(0:nph))
  call init_inclus

  allocate(ph(nex, -ltot:ltot, nspx, 0:nphx))
  allocate(ph2(nex, 0:ltot, 0:nphx))
  allocate(em(nex), eref(nex, nspx))
  allocate(xnmues(0:lx,0:nphx))

!KJ more allocations now that nphx, nspx are dynamic:
  allocate( imt(0:nphx),inrm(0:nphx),iz(0:nphx), &
   rmt(0:nphx),folp(0:nphx),folpx(0:nphx),xnatph(0:nphx),qnrm(0:nphx),xion(0:nphx), &
   rnrm(0:nphx), &
   dgc(251,41,0:nphx),dpc(251,41,0:nphx), &
   adgc(10,41,0:nphx),adpc(10,41,0:nphx), &
   edens(251,0:nphx),vclap(251,0:nphx),vtot(251,0:nphx),edenvl(251,0:nphx),vvalgs(251,0:nphx),dmag(251,0:nphx), &
   xnval(41,0:nphx), &
   iorb(-5:4,0:nphx) )
!KJ


  !call wlog("Read pot.dat")
  call rdpot ( ntitle, title, rnrmav, xmu, vint, rhoint,          &
               emu, s02, erelax, wp, ecv,rs,xf, qtotel,           &
               imt, rmt, inrm, rnrm, folp, folpx, xnatph,         &
               dgc0, dpc0, dgc, dpc, adgc, adpc,                  &
               edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,   &
               eorb, kappa, iorb, qnrm, xnmues, nohole, ihole,    &
               inters, totvol, iafolp, xion, iunf, iz, jumprm )

  !call wlog("Read phase.dat")
  call init_rdxsph

  do i=1,nrptx
    ripot(i) = rrr(i, dxpot)
  end do

  !call wlog("Init wavefunctions")
  allocate(prel(ne, 0:lx, nr, 0:nph), pnel(ne, 0:lx, nr, 0:nph), &
           qrel(ne, 0:lx, nr, 0:nph), qnel(ne, 0:lx, nr, 0:nph))
  call init_wavefunctions

  !call wlog("Read FMS matrix")
  call init_fms_matrix

end subroutine

subroutine init_inclus
  !
  ! Count how many atoms are within FMS radius (rfms2) of each representative atom ("in cluster")
  !
  integer iph, i
  real*8 rmax2, rr2

  IPH_LOOP: do iph = 0,nph
    inclus(iph)=0
    rmax2 = rfms2**2
    NAT_LOOP: do i=1,nat
      ! find squared distance between this ith atom and representative atom for potential iph
      rr2 = sum((rat(:,i) - rat(:,iatph(iph)))**2)
      if (rr2.le.rmax2) then
          inclus(iph) = inclus(iph) + 1
      endif
    enddo NAT_LOOP
  enddo IPH_LOOP
end subroutine

subroutine init_rdxsph
  !
  ! Read values in from xsph.bin
  !
  integer :: nphtmp, ihole
  real*8 :: rnrmav, xmutmp, edge
  integer :: iz(0:nphx)

  character*6, allocatable  :: potlbl(:)
  complex*16, allocatable :: rkk(:,:,:),rkk_nrixs(:,:,:,:)
  integer, allocatable    :: lmax(:,:) ! maximum angular momenta
  integer :: lmaxp1
  integer ik0

  allocate(rkk(nex,8,nspx))
  allocate(potlbl(0:nphx))
  allocate(lmax(nex,0:nphx))

  call rdxsph ( ne, ne1, ne3, nphtmp, ihole, rnrmav,xmu,edge,    &
                ik0, em, eref, iz, potlbl, ph, rkk, lmax, lmaxp1)

  deallocate(rkk, potlbl, lmax)
end subroutine

subroutine init_wavefunctions
  !
  ! Initialize solutions to Dirac equation
  !
  ! This code is a bit messy and has been cobbled together from LDOS code.
  ! You'll notice that chunks of it appear all over the place in FEFF. Nasty...
  ! We should really encapsulate this stuff somehow.
  !
  integer iph
  real*8 :: vjump

  real*8, allocatable, dimension(:) :: dum, vtotph, vvalph
  real*8, allocatable, dimension(:,:) :: dgcn, dpcn
  real*8, allocatable :: xrhoce(:,:,:)
  complex*16, allocatable, dimension(:,:,:) :: xrhole, xbruce
  complex*16 :: vtotc(nrptx), vvalc(nrptx)

  complex*16  p2, xkmt, ck, xck
  complex*16  pu, qu
  complex*16  xfnorm, xirf
  complex*16  temp,  phx

  complex*16 jl,jlp1,nl,nlp1
  integer :: i, ie, jri, jri1, ilast, lll, ncycle

  complex*16 pr(nrptx), qr(nrptx), pn(nrptx), qn(nrptx)
  real*8,  dimension(nrptx) :: dmagx

  integer ikap, irr, ic3

  allocate(dum(nrptx), vtotph(nrptx), vvalph(nrptx), dgcn(nrptx,41), dpcn(nrptx,41))

  IPH_LOOP: do iph = 0, nph

    ! If RGRID card is used, the potentials need to be interpolated onto
    ! the grid specified by that card. (Why the hell isn't this done in
    ! the potentials module instead of EVERYWHERE ELSE that needs to calculate
    ! wavefunctions?)

    call fixvar (rmt(iph),edens(1,iph),vtot(1,iph),dmag(1,iph),  &
          &       vint, rhoint, dxpot, dx, jumprm,                 &
          &       vjump, ri, vtotph, dum, dmagx)
    if (mod(ixc,10) .ge.5) then
        if (jumprm .gt. 0) jumprm = 2
        call fixvar (rmt(iph), edenvl(1,iph), vvalgs(1,iph),        &
            &       dmag(1,iph), vint, rhoint, dxpot, dx, jumprm,  &
            &       vjump, ri, vvalph, dum, dmagx)
        if (jumprm .gt. 0) jumprm = 1
    endif
    call fixdsx (iph, dxpot, dx, dgc, dpc, dgcn, dpcn)

    ! Do some magic adjusting of potentials
    jri = (log(rmt(iph)) + x0) / rgrd + 2
    jri1 = jri + 1
    eref0 = vtotph(jri1)

    do i = 1, jri1
        vtotph(i) = vtotph(i) - eref0
    enddo
    if (ixc.ge.5) then
        do i = 1, jri1
          vvalph(i) = vvalph(i) - eref0
        enddo
    else
        do i = 1, jri1
          vvalph(i) = vtotph(i)
        enddo
    endif


    ! complex potentials are needed for dfovrg routine
    do i = 1, nr
      vtotc(i)=vtotph(i)
      vvalc(i)= vvalph(i)
    enddo

    ! ilast is last point to integrate to.
    ! It is apparantly past the norman radius "for better interpolation"
    ! this may not be needed...
    ilast = (log(rnrm(iph)) + x0) / dx  + 8
    if (ilast.gt.nrptx) ilast = nrptx

    ! XXX LDOS sets ihole=0 here. is this correct?
    ihole = 0

    IE_LOOP: do ie = 1, ne
      ! p2 is (complex momentum)**2 referenced to energy dep xc
      p2 = em(ie) - eref0
      if (mod(ixc,10) .lt. 5) then
          ncycle = 0
      else
          ncycle = 3
      endif
      ck = sqrt(2*p2+ (p2*alphfs)**2)
      xkmt = rmt(iph) * ck
      LLL_LOOP: do lll=0,lmaxph(iph)
          !
          ! calculate regular solution
          !
          ikap = -1-lll
          irr = -1
          ic3 = 1
          if (lll.eq.0) ic3 = 0

          call dfovrg ( ncycle, ikap, rmt(iph), ilast, jri, p2, dx,          &
              &        ri, vtotc, vvalc, dgcn, dpcn, adgc, adpc,       &
              &        xnval, pu, qu, pn, qn,                          &
              &        iz(iph), ihole, xion(iph), iunf, irr, ic3, iph) !KJ iph

          call exjlnl (xkmt, lll, jl, nl)
          call exjlnl (xkmt, lll+1, jlp1, nlp1)
          call phamp (rmt(iph), pu, qu, ck,  jl, nl, jlp1, nlp1, ikap, phx, temp)
          ph2(ie, lll, iph) = phx

          ! Normalize solutions at rmt to rmt*(jl*cos(delta) - nl*sin(delta))
          xfnorm = 1 / temp
          do i = 1,ilast
            pr(i)=pn(i)*xfnorm
            qr(i)=qn(i)*xfnorm
          enddo

          !
          ! calculate irregular solution
          !
          irr = 1
          pu = ck*alphfs
          pu = - pu/(1+sqrt(1+pu**2))
          ! set pu, qu - initial condition for irregular solution at ilast
          qu=(nlp1*cos(phx)+jlp1*sin(phx))*pu *rmt(iph)
          pu = (nl*cos(phx)+jl*sin(phx)) *rmt(iph)

          call dfovrg (ncycle, ikap, rmt(iph), ilast, jri, p2, dx,           &
              &       ri, vtotc,vvalc, dgcn, dpcn, adgc, adpc,         &
              &       xnval, pu, qu, pn, qn,                           &
              &       iz(iph), ihole, xion(iph), iunf, irr, ic3, iph) !KJ iph
          ! set N- irregular solution , which is outside
          ! N=(nlp1*cos(ph0)+jlp1*sin(ph0))*factor *rmt * dum1
          ! N = i*R - H*exp(i*ph0)
          temp = exp(coni*phx)
          ! calculate wronskian (qu)
          qu = 2 * alpinv * temp * ( pn(jri)*qr(jri) - pr(jri)*qn(jri) )
          qu = 1/qu/ck
          !       qu should be close to 1
          do i = 1, ilast
            !         pr(i) = pr(i)*qu (should this qu be put in the line below?)
            !         qr(i) = qr(i)*qu
            pn(i) = coni * pr(i) - temp * pn(i)*qu
            qn(i) = coni * qr(i) - temp * qn(i)*qu
          enddo

          !  Use exact solution to continue solutions beyond rmt
          pu = ck*alphfs
          pu = -pu/(1+sqrt(1+pu**2))

          do i=jri,nr
            xck = ck * ri(i)
            temp = xck
            call exjlnl (temp, lll, jl, nl)
            call exjlnl (temp, lll+1, jlp1, nlp1)
            pr(i)= (jl*cos(phx)-nl*sin(phx)) *ri(i)
            qr(i)=(jlp1*cos(phx)-nlp1*sin(phx))*pu *ri(i)
            pn(i)= (nl*cos(phx)+jl*sin(phx)) *ri(i)
            qn(i)=(nlp1*cos(phx)+jlp1*sin(phx))*pu *ri(i)
          enddo

          ! smooth irregular solutions near origin (TODO: fix this correctly!)
          if (lll == 0) then
            call fix_irreg(pn)
            call fix_irreg(qn)
          end if

          ! store solutions
          prel(ie,lll,:,iph) = pr(:nr)
          pnel(ie,lll,:,iph) = pn(:nr)
          qrel(ie,lll,:,iph) = qr(:nr)
          qnel(ie,lll,:,iph) = qn(:nr)
      enddo LLL_LOOP
    enddo IE_LOOP
  enddo IPH_LOOP

  ! DEBUG: save out wavefunctions for testing
  !open(unit=8, file="wf.bin", form="unformatted")
  !write(8) ne, lx+1, nr, nph+1
  !write(8) prel
  !write(8) pnel
  !write(8) qrel
  !write(8) qnel
  !close(8)

end subroutine

subroutine fix_irreg(y0)
  !
  ! The irregular solutions tend to be oscillatory near the origin for l=0.
  ! This is an extremely hacky fix that simply fits a polynomial in this region.
  ! A better solution is needed...
  !
  use polyfit_mod
  complex*16, intent(inout) :: y0(nrptx)

  complex*16 :: y(100), dy(100)
  complex*16 :: fit(4)
  integer :: i

  call polyfit(ri(50:100), y0(50:100), 3, fit)
  call polyval(fit, ri(:100), y0(:100))
end subroutine

subroutine init_fms_matrix
  integer ldim, istate, ne_file, nat_file
  integer ier

  ! Read gg_slice.bin (contains piece of scattering matrix for r at central site, r' anywhere)
  open(unit=8, file='gg_slice.bin', status='old', iostat=ier, form='unformatted')
  if (ier .eq. 0) then
    read(8) ldim, istate, ne_file
    ! XXX check that these match with values elsewhere
    allocate(gg_slice(ldim, istate, ne_file))
    read(8) gg_slice
  end if
  close(8)

  open(unit=8, file='gg_diag.bin', status='old', iostat=ier, form='unformatted')
  if (ier .eq. 0) then
    read(8) ldim, ldim, nat_file, ne_file

    ! XXX check that these match with values elsewhere
    allocate(gg_diag(ldim, ldim, nat_file, ne_file))
    read(8) gg_diag
  end if
  close(8)
end subroutine

subroutine rhorrp_deinit
  deallocate(ph, em, eref, xnmues, prel, pnel, qrel, qnel, inclus)
  if (allocated(gg_slice)) deallocate(gg_slice)
  if (allocated(gg_diag)) deallocate(gg_diag)
end subroutine

subroutine rhoerrp(v, vp, rhoe)

  !------------------------------------------------------
  !
  ! Calculate energy dependent 1-electron density matrix rho(r, rp, E)
  !
  ! Parameters:
  !       v: 3-vector of r point  (real*8, dimension(3))
  !      vp: 3-vector of r' point (real*8, dimension(3))
  !    rhoe: array to hold result (real*8, dimension(ne))
  !
  !------------------------------------------------------
  use DimsMod, only: lx
  use constants

  real*8, intent(in), dimension(3) :: v, vp
  complex*16, intent(out), dimension(ne) :: rhoe

  real*8 r, rp
  real*8 f, fp
  integer i, ip

  integer ie, l, lp, iL, iLp

  integer iat, iatp ! atomic sites closest to v and v'
  integer iph, iphp ! potential type closest to v and v'
  real*8, dimension(3) :: dv, dvp ! displacements from site center

  complex*16, dimension((lx+1)**2) :: y,yp
  real*8, dimension(0:lx) :: pl0

  complex*16, dimension(ne, 0:lx) :: rhoel
  complex*16, dimension(ne) :: Ge

  complex*16, dimension(ne,0:lx,2) :: Rl, Rlp

  integer i1, i2
  real*8 f1, f2

  real*8 ctheta
  real*8 tmp, tmp2
  real*8 de, epsilon

  complex*16 pref, ck, pu, p2
  complex, dimension(ne, ldim,ldim) :: gg

  logical CentralOnly, NoScatter

  ! These flags are for debugging only
  CentralOnly = .false.
  NoScatter = .false.

  ! determine which site is closest
  r = sqrt(sum(v*v))
  if (CentralOnly) then
    iat = 1
    iph = 0
    dv(:) = v(:)
  else
    call nearest_atom(v, iat, iph, dv, .true.)
    r = sqrt(sum(dv*dv))
  end if

  rp = sqrt(sum(vp*vp))
  if (CentralOnly) then
    iatp = 1
    iphp = 0
    dvp(:) = vp(:)
  else
    call nearest_atom(vp, iatp, iphp, dvp, .true.)
    rp = sqrt(sum(dvp*dvp))
  end if

  ! For Compton calculations (at least in an ordered crystal), we restrict
  ! the r coordinate to the Voronoi cell of the central atom, while allowing
  ! r' to extend elsewhere.
  if (restrict_r_voronoi .and. iat .ne. 1) then
    rhoe(:) = 0.0
    return
  end if

  ! The full scattering matrix is a bit unwieldy for large clusters. So,
  ! we save only a few useful slices of it (see FMS/fmstot.f90)
  !   gg_diag: site-diagonal piece (iat == iatp) - needed for density calculations
  !   gg_slice: r near central atom, r' anywhere
  ! If rho(r,r') is needed for r away from the central atom and in a site different from r',
  ! gg_full in FMS/fmstot.f90 needs to be sent back to the master process, saved out
  ! and reloaded in this module.
  if (iat.eq.iatp) then
    do ie=1,ne
      gg(ie,:,:) = gg_diag(:,:,iat,ie)
    end do
  else
    if (iat.ne.1) then ! we don't have this piece
      rhoe(:) = 0.0
      return
    endif
    do ie=1,ne
      gg(ie,:,:) = gg_slice((iat-1)*ldim+1:iat*ldim, (iatp-1)*ldim+1:iatp*ldim, ie)
    end do
  end if

  ! If either v or vp put us right on top of an atom, offset slightly in the z direction
  ! to avoid dividing by zero. (XXX is there a better way to handle this?)
  epsilon = 1e-3
  if (abs(r) < epsilon) then
    r = epsilon
    dv(3) = dv(3) + epsilon
  end if
  if (abs(rp) < epsilon) then
    rp = epsilon
    dvp(3) = dvp(3) + epsilon
  end if

  ! values for interpolating in r (see interp_wf function)
  ! XXX remove the r,rp.eq.0 cases since they can't be reached
  if (r .eq. 0) then
    f = 1
  else
    f = (log(r) + x0) / dx + 1
  end if
  if (f < 1) f = 1
  if (f > nr) f = nr !XXX we should probably just cut it off here
  i = f
  f = f - i

  ! values for interpolating in r'
  if (rp .eq. 0) then
    fp = 1
  else
    fp = (log(rp) + x0) / dx + 1
  end if
  if (fp < 1) fp = 1
  if (fp > nr) fp = nr
  ip = fp
  fp = fp - ip

  rhoel(:,:) = 0
  Ge(:) = 0

  !--------------------
  !
  ! v, v' in same cell
  !
  !--------------------

  if (iat .eq. iatp) then
    ! find cos of angle between v and v' and calculate legendre poly
    ctheta = sum(dv*dvp) / sqrt(sum(dv*dv)*sum(dvp*dvp))
    call cpl0(ctheta, pl0, lx+1)

    ! lesser and greater indices
    i1 = i
    f1 = f
    i2 = ip
    f2 = fp
    if (i1 .gt. i2) then
      i2 = i
      f2 = f
      i1 = ip
      f1 = fp
    end if

    ! build same site part of G
    rhoel(:,:) = rhoel(:,:) - &
    &  interpwf(prel(:,:,:,iph), i1, f1) * &
    &  (interpwf(pnel(:,:,:,iph), i2, f2) - coni * interpwf(prel(:,:,:,iph), i2, f2)) -&
    &  interpwf(qrel(:,:,:,iph), i1, f1) * &
    &  (interpwf(qnel(:,:,:,iph), i2, f2) - coni * interpwf(qrel(:,:,:,iph), i2, f2))

    ! sum over l
    do l = 0,lx
      Ge(:) = Ge(:) + rhoel(:,l) * pl0(l) * (2*l + 1) / (4*pi)
    end do
  endif

  !----------------------
  !
  ! scattering piece
  !
  !----------------------

  if (.not. NoScatter) then
    call ylm(dv,lx,y)
    call ylm(dvp,lx,yp) ! note that this is relative to iatp'th site

    ! radial parts of wave functions (upper and lower components)
    Rl(:,:,1) = interpwf(prel(:,:,:,iph), i, f)
    Rl(:,:,2) = interpwf(qrel(:,:,:,iph), i, f)

    Rlp(:,:,1) = interpwf(prel(:,:,:,iphp), ip, fp)
    Rlp(:,:,2) = interpwf(qrel(:,:,:,iphp), ip, fp)

    do iL = 1,(lx+1)**2
      l = ceiling(sqrt(1.0*iL)-1)
      do iLp = 1,(lx+1)**2
        lp = ceiling(sqrt(1.0*iLp)-1)
        Ge(:) = Ge(:) + ( &
                 (Rl(:,l,1) * Rlp(:,lp,1) + Rl(:,l,2) * Rlp(:,lp,2)) * &
                 y(iL) * conjg(yp(iLp))  * &
                 coni**l * (-coni)**lp * &
                 !exp(coni * ph(:,l,1,iph)) * exp(coni * ph(:,lp,1,iphp)) * &
                 exp(coni * (ph2(:,l,iph) + ph2(:,lp,iphp))) * &
                 gg(:,iL,iLp))
      end do
    end do
  end if

  !-------------------------------------
  !
  ! Finish up calculation of rho(E,r,r')
  !
  !-------------------------------------

  do ie = 1,ne
    p2 = em(ie) - eref0
    ck = sqrt(2*p2+ (p2*alphfs)**2)
    pu = ck*alphfs
    pu = - pu/(1+sqrt(1+pu**2))

    pref = 1/(1+pu**2) /pi *ck*4

    rhoe(ie) = Ge(ie) * pref / r / rp
  end do
end subroutine

subroutine rhorrp(v, vp,  rho)
  !------------------------------------------------
  !
  ! Calculate 1-electron density at point `v`, `vp`
  !
  !------------------------------------------------
  !use DimsMod, only:
  use constants

  real*8, intent(in), dimension(3) :: v, vp
  real*8, intent(out) :: rho

  complex*16, dimension(ne) :: rhoe
  complex*16 :: xrho, rho0, rho0p
  integer ie, ihstart, nsub, isub
  real*8 :: temperature

  complex*16 :: fermi, fermip
  complex*16 :: ee, de
  real*8 :: emreal(ne)

  ! Calculate Green's function
  call rhoerrp(v, vp, rhoe)

  ! XXX phmesh should save the temperature it used somewhere so we can load it back in
  temperature = max(scf_temperature/hart, 0.001)

  ! Now integrate rho(E,r,rp) over fermi distribution to get occupied density

  ! We start with the piece between the real axis and the first point, assuming rho(ecv)=0
  fermi = fermi_dist(em(1), xmu, temperature)
  rho0 = 0.0d0
  xrho = em(1) * rhoe(1) * fermi
  rho0p = rhoe(1)
  fermip = fermi
  ie = 1
  ! Integrate up vertically
  do ie = 2,ne1
    de = (em(ie) - em(ie-1))

    if (dble(de) > 1e-15) then
      ! ihstart is index of corner connecting vertical and horizontal legs
      ihstart = ie-1
      exit
    end if

    fermi = fermi_dist(em(ie), xmu, temperature)
    rho0 = rhoe(ie)

    xrho = xrho + (rho0p * fermip + rho0* fermi) / 2.0 * de

    rho0p = rho0
    fermip = fermi

  end do

  ! Now, although rhoe is smooth on the sampled grid, the Fermi distribution is not
  ! So, we perform the horizontal integration on a finer grid, interpolating rhoe
  ! For now, hard code nsub. We could be able to change it based on temperature, de
  ! and how close energy is to mu, but this is pretty fast, so it would probably not
  ! be a worthwhile gain.
  nsub = 10
  emreal(:) = dble(em)
  do ie = ihstart+1,ne1
    de = (em(ie) - em(ie-1)) / nsub
    do isub=1,nsub
      ee = em(ie-1) + isub*de

      fermi = fermi_dist(ee, xmu, temperature)
      call terpc(emreal(ihstart:ne1), rhoe(ihstart:ne1), ne1-ihstart+1, 2, dble(ee), rho0)

      xrho = xrho + (rho0p * fermip + rho0 * fermi) / 2.0 * de

      rho0p = rho0
      fermip = fermi
    end do
  end do

  ! Add on the Matsubara poles
  do ie = ne1+1, ne
    xrho = xrho + (-2*pi*coni*temperature*rhoe(ie))
  end do

  ! Finally, take the imaginary part to get the density
  rho = dimag(xrho)
end subroutine

subroutine nearest_atom(v, iat, iph, dv, fmsF)

  !-------------------------------------------------------------------
  !
  ! Return the atomic site closest to v and coordinates realtive to it
  !
  ! Input
  !   v   - vector (x,y,z)
  !   fmsF- flag for considering atoms only in fms radius
  ! Output
  !   iat - index of atomic site closest to v
  !   iph - index of potential of iat'th atom
  !   dv  - v - rat(iat) (relative displacement of v from atom)
  !
  !-------------------------------------------------------------------
  use atoms_inp

  real*8, intent(in), dimension(3) :: v
  logical fmsF
  integer, intent(out) :: iat, iph
  real*8, intent(out), dimension(3) :: dv

  integer i, num
  real*8 dr2, mindr2
  real*8, dimension(3) :: ddv

  mindr2 = -1
  dv = v
  if (fmsF) then
    num = inclus(0)
  else
    num = nat 
  endif
  do i=1,num
    ddv = v - rat(:,i)
    dr2 = sum(ddv**2)

    if (mindr2 < 0 .or. dr2 < mindr2) then
      mindr2 = dr2
      iat = i
      iph = iphat(iat)
      dv(:) = ddv(:)
    endif
  enddo
end subroutine

function interpwf(wf, i, f)

  !---------------------------------------------------------
  !
  ! Interpolate a wavefunction wf at a point a fraction f of
  ! way between point i and point i+1
  !
  ! wf - the wavefunction
  !  i - index of closest point below desired point
  !  f - fraction of distance between ith and (i+1)th point
  !
  !----------------------------------------------------------

  complex*16, intent(in), dimension(ne, 0:lx, nr) :: wf
  integer, intent(in) :: i
  real*8, intent(in) :: f

  complex*16, dimension(ne, 0:lx) :: interpwf

  if (i < 0) then
    interpwf(:,:) = 0
  else if (i == 0) then
    interpwf(:,:) = wf(:,:,i+1) * f
  else
    interpwf(:,:) = wf(:,:,i) * (1-f) + wf(:,:,i+1) * f
  endif

end function

function fermi_dist(energy, mu, T)
  use compton_inp, only:  set_chemical_potential, chemical_potential, temperature
  implicit none

  complex*16, intent(in) :: energy ! in Hartree
  real*8, intent(in) :: mu, T ! in Hartree
  complex*16 fermi_dist
  real*8 temp, mu_set

  if (set_chemical_potential) then
      mu_set = chemical_potential / hart
  else
      mu_set = mu
  endif

  if ( abs( temperature / hart - T )  < 2.0d-2 ) then
     temp = T
  else
     temp = temperature / hart
  endif

  if (T < 1e-5) then
    if (dble(energy) < dble(mu_set)) then
      fermi_dist = 1.0d0
    else
      fermi_dist = 0.0d0
    endif
  else
    fermi_dist = 1.0d0 / (exp ( (energy - mu_set) / T ) + 1.0d0 )
  endif
end function

subroutine atomic_density(r, il, rho)
  ! calculate the real-space density for core from atomic wavefunctions
  !
  ! Parameters:
  !   r:  cooordinate to calculate density at
  !   il: angular momentum to calculate for
  !   rho: output - density at r
  !
  ! For il > 0, this is the spherically averaged density
  !
  double precision, intent(in) :: r(3)
  integer, intent(in) :: il
  double precision, intent(out) :: rho

  integer iat
  double precision :: dr2, dr
  double precision :: p,q, drho

  rho = 0.0d0
  do iat=1,nat
    dr2 = sum((rat(:,iat) - r)**2)
    if (dr2 > 4.0) cycle
    dr = sqrt(dr2)
    if (dr < 1e-4) dr = 1e-4

    call terp(ripot, dgc(:,il,iphat(iat)), nr, 2, dr, p)
    call terp(ripot, dpc(:,il,iphat(iat)), nr, 2, dr, q)

    drho = (p**2 + q**2)/4/pi/dr**2
    rho = rho + drho
  end do
end subroutine

end module
