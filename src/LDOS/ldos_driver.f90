
  subroutine ldos_driver
!  The main workhorse of the LDOS module
!  After the main has taken care of generic tasks, such as setting up the 
!  dimensions module, the parallel definitions, and figuring out whether
!  any work needs to be done at all;
!  this routine dispatches the real tasks and creates data structures.
!  - read input data from pot.bin file (potentials calculated by pot module)
!  - control loop calls "ldos" sub for regular LDOS calculation;
!         or "ldos_h" sub for HUBBARD LDOS calculation.

  use DimsMod, only: nphx=>nphu, novrx, nheadx, init_dimensions, istatx, lx, nclusx
  use rotx
  use stkets
  use lnlm
  use xstruc
  use t3j
  use par
  use ldos_inp, only: lfms2, ixc, ispin, minv, rfms2, emin, emax, eimag, rgrd, rdirec, toler1, toler2, lmaxph
  use atoms_inp, only: rat, nat, iphat, iatph
  use hubbard_inp, only: mldos_hubb
  use potential_inp, only: folp, iafolp, nohole, xnatph, iz, xion, nph, &
      jumprm, ihole, iunf, ecv, inters, totvol, ntitle, title, novr, nnovr, &
      iphovr, rovr
  use errorfile
  use xsph_inp, only: spinph, xsph_read
  implicit none

  ! HUBBARD data
  real*8 vint_sp(2)
  real*8,  allocatable :: vtot_sp(:,:,:), vvalgs_sp(:,:,:)
  real*8,  allocatable :: edens_sp(:,:,:),edenvl_sp(:,:,:)
  real*8,  allocatable :: rmt_sp(:,:)
  !  pot.bin data
  integer, allocatable :: imt(:),inrm(:)
  real*8,  allocatable :: folpx(:),rmt(:),rnrm(:),qnrm(:)
  real*8,  allocatable :: dgc0(:),dpc0(:)
  real*8,  allocatable :: dgc(:,:,:),dpc(:,:,:)
  real*8,  allocatable :: adgc(:,:,:),adpc(:,:,:)
  integer, allocatable :: kappa(:)
  real*8,  allocatable :: edens(:,:),edenvl(:,:),vclap(:,:),vvalgs(:,:),vtot(:,:),dmag(:,:)
  real*8,  allocatable :: xnval(:,:)
  real*8,  allocatable :: eorb(:)
  integer,  allocatable :: iorb(:,:)  !KJ 12-2011 changed from real*8 ; it's read as format (i2) in rdpot.  Not really used here anyway ...
  real*8, allocatable :: xnmues(:,:)
! Variables that exist in both pot.inp and pot.bin are declared/allocated in the pot_inp module
! (initialized in reldos)
! end of pot.bin data

  ! Added to satisfy implicit none:
  integer :: ios, i, is, iph !,ntitle
  real*8  :: dx,dxnew,x0,rnrmav
  real*8  :: xmu,vint,rhoint,emu,s02,erelax,wp,rs,xf,qtotel
  real*8  :: xmunew
  integer :: idmag


    ! Allocate dimensions and module variables
  call init_stkets(istatx)
  call init_rotx(lx,nclusx)
  call init_lnlm(lx,nclusx)
  call init_xstruc(nclusx)
  call init_t3j(lx)

  ! Allocate local variables
  allocate(xnmues(0:lx,0:nphx))
  allocate(imt(0:nphx),inrm(0:nphx))
  allocate(folpx(0:nphx),rmt(0:nphx),rnrm(0:nphx),qnrm(0:nphx))
  allocate(dgc0(251),dpc0(251))
  allocate(dgc(251,30,0:nphx),dpc(251,30,0:nphx))
  allocate(adgc(10,30,0:nphx),adpc(10,30,0:nphx))
  allocate(kappa(30))
  allocate(edens(251,0:nphx),edenvl(251,0:nphx),vclap(251,0:nphx),vvalgs(251,0:nphx),vtot(251,0:nphx),dmag(251,0:nphx))
  allocate(xnval(30,0:nphx))
  allocate(eorb(30))
  allocate(iorb(-4:3,0:nphx))


  ! Read the potential (it is spin-averaged)
  ! KJ this is dangerous/distasteful: many of these variables were already read by reldos->pot_read
  ! Now we read the same variables from pot.bin ...
  ! Potential conflicts; ugly; unnecessary.
  ! rdpot also used in SCREEN/prep and XSPH/xsphsub.
  call rdpot ( ntitle, title, rnrmav, xmu, vint, rhoint,          &
             emu, s02, erelax, wp, ecv,rs,xf, qtotel,           &
             imt, rmt, inrm, rnrm, folp, folpx, xnatph,         &
             dgc0, dpc0, dgc, dpc, adgc, adpc,                  &
             edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,   &
             eorb, kappa, iorb, qnrm, xnmues, nohole, ihole,    &
             inters, totvol, iafolp, xion, iunf, iz, jumprm)

  ! Atom r-grid
  dx = 0.05d0
  x0 = 8.8d0
  ! Phase r-grid
  dxnew = rgrd

  if (mldos_hubb.eq.1) then
    ! Regular, non-Hubbard LDOS calculation
    ! make spin dependent potential only if needed
    if (ispin.ne.0) then
       idmag = 1
       if (ispin.lt.0) idmag = -1
       do iph = 0, nph
          do i = 1, 251
             dmag(i,iph) = dmag(i, iph) * spinph(iph)
          enddo
       enddo

       call  istprm (nph, nat, iphat, rat, iatph, xnatph,               &
                 novr, iphovr, nnovr, rovr, folp, folpx, iafolp,    &
                 edens, edenvl, idmag,                              &
                 dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,   &
                 ixc, rhoint,vint, rs, xf, xmu, xmunew,             &
                 rnrmav, qtotel, inters, totvol)
       xmunew = xmu
    endif

    ! calculate LDOS
    call ldos (nph, edens, edenvl, dmag, vtot, vvalgs, rmt, rnrm,           &
          ixc, rhoint, vint, xmu, jumprm, x0, dx, rgrd, xion, iunf, iz, &
          xnval, adgc, adpc, dgc, dpc,                                  &
          ihole, qnrm, xnmues,                                          &
          emin, emax, eimag, rfms2, lfms2, lmaxph, nat, iphat, rat,     &
          minv, rdirec, toler1, toler2)

  elseif (mldos_hubb.eq.2) then
    ! Hubbard LDOS calculation
    call wlog(' Using HUBBARD approximation')

    ! HUBBARD data
    allocate(vtot_sp(251,0:nphx,2), vvalgs_sp(251,0:nphx,2))
    allocate(edens_sp(251,0:nphx,2), edenvl_sp(251,0:nphx,2), rmt_sp(0:nphx,2))

    ! make spin dependent potential, always
    do iph = 0, nph
       do i = 1, 251
          dmag(i,iph) = dmag(i, iph) * spinph(iph)
       enddo
    enddo
    do is=1,2
       if(is.eq.1) idmag=1
       if(is.eq.2) idmag=-1
       call  istprm (nph, nat, iphat, rat, iatph, xnatph,         &
                 novr, iphovr, nnovr, rovr, folp, folpx, iafolp,  &
                 edens, edenvl, idmag,                            &
                 dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm, &
                 ixc, rhoint,vint, rs, xf, xmu, xmunew,           &
                 rnrmav, qtotel, inters, totvol)
       xmunew = xmu
       vint_sp(is)=vint
       rmt_sp(:,is)=rmt(:)
       vtot_sp(:,:,is)=vtot(:,:)
       vvalgs_sp(:,:,is)=vvalgs(:,:)
       edens_sp(:,:,is)=edens(:,:)
       edenvl_sp(:,:,is)=edenvl(:,:)
    enddo
    write(*,*)'in ldos_driver , I am ',master

    ! calculate LDOS in Hubbard approximation
    call ldos_h_unrolled(nph, edens_sp, edenvl_sp, dmag, vtot_sp, vvalgs_sp, &
      rmt, rmt_sp, rnrm, ixc, rhoint, vint, vint_sp, xmu, jumprm, x0, dx, &
      rgrd, xion, iunf, iz, xnval, adgc, adpc, dgc, dpc, ihole, qnrm, xnmues, &
      emin, emax, eimag, rfms2, lfms2, lmaxph, nat, iphat, rat, minv, rdirec, &
      toler1, toler2)
    deallocate(vtot_sp, vvalgs_sp, edens_sp, edenvl_sp, rmt_sp)
  endif  ! Hubbard or no Hubbard


  ! Deallocate local variables
  deallocate(xnmues)
  ! Deallocate modules
  call kill_stkets
  call kill_rotx
  call kill_lnlm
  call kill_xstruc
  call kill_t3j

  return
end subroutine ldos_driver
